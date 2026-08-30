# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The abundances: one matrix held as both `species × location` and `species × Y × X` views of the
# same memory, plus the saved and cached forms a recorded run produces.

using DimensionalData
using DimensionalData.Lookups: NoLookup, Categorical
using DimensionalData: Ti
using Missings
using Random

"""
    GridLandscape

An ecosystem's abundances: how many individuals of each species are in each cell. The same numbers
are held under two names, as two views of one block of memory, because the simulation asks two
different questions of them — the hot loop walks species against flat cells, while inspecting or
seeding a run wants to say *where*.

Immutable, and that is what makes the two views agree: the only way to change shape is to build a
new `GridLandscape` and reassign the whole field holding it, so they cannot drift apart.

# Fields

  - `matrix`: species against flat grid cell, with no `Y`/`X` structure — the per-timestep access
    pattern, and a plain `Matrix{Int64}`.
  - `grid`: species against `Y` and `X`, a plain `Array{Int64, 3}` sharing `matrix`'s memory.
  - `dimmatrix`: `matrix` as a `(Dim{:species}, Dim{:location})` `DimArray`, carrying real species
    names.
  - `dimgrid`: `grid` as a `(Dim{:species}, Y, X)` `DimArray`, with `Y` and `X` the habitat's own
    dimensions — real coordinates as `Unitful` lengths, or degrees on a geographic grid.

`dimmatrix`'s `location` dimension is deliberately a `NoLookup`: a flat cell index has no coordinate
to carry, and `dimgrid` is where a cell's position is answered.

Every landscape built from a habitat has real coordinates, synthetic or not — a `StudyArea` fixes
them whether or not any data was read. Only the size-only constructor below, which has no habitat to
take them from, produces `NoLookup` positions.

**The raw arrays and the labelled views are the same memory, and both are stored on purpose.** The
hot loop reads `matrix`; asking *which* species or *where* reads `dimmatrix`/`dimgrid`.

**The raw fields are what make the hot loop free, and the reason is not obvious.** `Ecosystem`
declares `abundances::GridLandscape` — the bare `UnionAll`, an abstract field. Julia can still infer
a field of that container concretely **provided the field's declared type does not mention the type
parameters**, which is exactly why `matrix::Matrix{Int64}` is written out rather than left as one of
the parameters. Measured: reaching the labelled view through the abstract field costs about 176 bytes
**per cell** and grows with the grid; reading the raw field is **0**. Declaring these two fields in
terms of `Am`/`Ag` would silently undo that.
"""
struct GridLandscape{Am <:
                     DimensionalData.AbstractDimArray{Int64, 2,
                                                      <:Tuple{<:Dim{:species},
                                                              <:Dim{:location}}},
                     Ag <:
                     DimensionalData.AbstractDimArray{Int64, 3,
                                                      <:Tuple{<:Dim{:species},
                                                              <:Y, <:X}}}
    matrix::Matrix{Int64}
    grid::Array{Int64, 3}
    dimmatrix::Am
    dimgrid::Ag

    function GridLandscape(abun::Matrix{Int64}, names::Vector{String},
                           yx::Tuple{<:Y, <:X})
        length(names) == size(abun, 1) ||
            throw(DimensionMismatch("`names` has $(length(names)) entries but `abun` has $(size(abun, 1)) species"))
        prod(length.(yx)) == size(abun, 2) ||
            throw(DimensionMismatch("grid `yx` has $(prod(length.(yx))) cells but `abun` has $(size(abun, 2))"))
        # `reshape` of an `Array` shares its memory, and a `DimArray` wraps rather than copies, so
        # all four fields below are views of one buffer. That is what the immutability above is for:
        # the only way to change shape is to build a new `GridLandscape`, so they cannot drift.
        g = reshape(abun, (length(names), length.(yx)...))
        sp = Dim{:species}(Categorical(names))
        dm = DimArray(abun, (sp, Dim{:location}(NoLookup())))
        dg = DimArray(g, (sp, yx...))
        return new{typeof(dm), typeof(dg)}(abun, g, dm, dg)
    end
end

"""
    SavedLandscape

A [`GridLandscape`](@ref) reduced to what has to survive a round trip through disk: the bare
abundance matrix, without its dimensions, and a snapshot of the per-species random number streams.

The streams are what make a resumed run reproducible — a run draws from one generator per species
rather than from a shared one, so restoring the abundances alone would restart every stream from
wherever the loading process happened to be. The dimensions are not saved because the ecosystem being
restored into already has them.

# Fields

  - `matrix`: the abundances, species by flat grid cell.
  - `rngs`: one generator state per species, copied at the moment of saving.
"""
struct SavedLandscape
    matrix::Matrix{Int64}
    rngs::Vector{Random.Xoshiro}
end

# `matrix` is annotated ABSTRACTLY ON PURPOSE, and must stay that way. The obvious repair —
# parameterising the struct on the array type — does not remove the abstraction, it moves it:
# `CachedEcosystem` declares `abundances::CachedGridLandscape`, which is a concrete annotation only
# while this type is concrete. Parameterising here makes that a `UnionAll` instead, so the same
# imprecision reappears one level up, and removing it there would give an exported type a fourth
# type parameter. Measured, not assumed: `isconcretetype(fieldtype(CachedEcosystem, :abundances))`
# is `true` as written and `false` once this is parameterised.
#
# The cost is bounded because this is the disk-cache landscape: every read is
# `cache.abundances.matrix[Ti(At(tm))]` on the save/load path, where the I/O dominates. Nothing in
# the hot loop touches it.
"""
    CachedGridLandscape

A recorded run's abundances over time: one [`GridLandscape`](@ref) per timepoint, or `missing` where
that point has not been computed yet.

How often checkpoints reach disk is separate from how often the simulation steps, and does not affect
the answer: the run always advances by `timestep`, so `saveinterval` chooses only how much is kept.

# Fields

  - `matrix`: a `DimArray` over a `Ti` axis holding `GridLandscape`s and `missing`s. Indexed by
    **time value**, which needs a selector — `matrix[Ti(At(t))]`, not `matrix[t]`, which is
    positional and an error for a `Unitful.Time`.
  - `outputfolder`: where the JLD2 cache files are written.
  - `saveinterval`: how often a checkpoint reaches disk. A multiple of `timestep`.
  - `timestep`: the simulation step, and the granularity of the time axis.
"""
struct CachedGridLandscape
    matrix::DimensionalData.AbstractDimVector{Union{GridLandscape, Missing}}
    outputfolder::String
    saveinterval::Unitful.Time
    timestep::Unitful.Time
end
# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# Without these the default prints the whole abundance matrix twice over — `matrix` and `grid` are
# two views of one buffer — measured at 291 648 characters on a single line for 12 species over a
# 60 x 100 grid.
#
# The total abundance is a genuine reduction over the data, which every other `show` in the package
# avoids. It is here deliberately: it is the one number a reader wants from a landscape, and it is a
# sum over an `Int64` matrix rather than anything allocating. Nothing else here touches the values.
function Base.show(io::IO, l::GridLandscape)
    nsp, ny, nx = size(l.grid)
    return print(io,
                 "GridLandscape($(nsp) species × $(ny) × $(nx), $(sum(l.matrix)) individuals)")
end

function Base.show(io::IO, ::MIME"text/plain", l::GridLandscape)
    nsp, ny, nx = size(l.grid)
    println(io, "GridLandscape")
    println(io, "  species   ", nsp)
    println(io, "  grid      $(ny) × $(nx) cells")
    print(io, "  total     ", sum(l.matrix), " individuals")
    return nothing
end

# Build one from a plain size, for callers with no `SpeciesList` or `GridHabitat` to hand. Species
# are named by their number, matching `GridHabitat`'s own default-name convention, and `Y`/`X` are
# fresh `NoLookup` dimensions, because a bare size carries no positions to attach. This is the only
# route that produces unpositioned cells: a landscape built from a habitat always has real
# coordinates, since a `StudyArea` decides them whether or not any data was read.
function GridLandscape(abun::Matrix{Int64}, dimension::Tuple)
    names = string.(axes(abun, 1))
    yx = (Y(NoLookup(Base.OneTo(dimension[2]))),
          X(NoLookup(Base.OneTo(dimension[3]))))
    return GridLandscape(abun, names, yx)
end

"""
    GridLandscape(sl::SavedLandscape, names::Vector{String}, yx::Tuple{<:Y, <:X})

Restore a `GridLandscape` from a [`SavedLandscape`](@ref), putting its bare abundance matrix back
onto dimensions.

# Arguments

  - `sl`: the saved abundances.
  - `names`: the species names, from the ecosystem being restored into.
  - `yx`: that ecosystem's `Y` and `X` dimensions.

The saved random number streams are restored separately, by the caller that already holds the
ecosystem to put them in.
"""
function GridLandscape(sl::SavedLandscape, names::Vector{String},
                       yx::Tuple{<:Y, <:X})
    return GridLandscape(sl.matrix, names, yx)
end

"""
    SavedLandscape(gl::GridLandscape, rngs::Vector{Random.Xoshiro})

Reduce a [`GridLandscape`](@ref) to its serialisable form.

# Arguments

  - `gl`: the landscape whose abundances are to be kept.
  - `rngs`: the ecosystem's per-species generators, copied rather than aliased so that continuing the
    run does not alter what was saved.
"""
function SavedLandscape(gl::GridLandscape, rngs::Vector{Random.Xoshiro})
    return SavedLandscape(gl.matrix, copy.(rngs))
end

"""
    CachedGridLandscape(file::String, times::StepRangeLen;
                        saveinterval::Unitful.Time = step(times))

Construct a `CachedGridLandscape` backed by a folder, with every timepoint still `missing`.

# Arguments

  - `file`: the folder the cache files are written to.
  - `times`: every timepoint the run will cover. Its step size is the simulation timestep.
  - `saveinterval`: how often a checkpoint reaches disk. Must be a multiple of the timestep, and
    defaults to saving every step.
"""
function CachedGridLandscape(file::String, times::StepRangeLen;
                             saveinterval::Unitful.Time = step(times))
    timestep = step(times)
    iszero(mod(saveinterval, timestep)) ||
        error("saveinterval ($saveinterval) must be a multiple of the timestep ($timestep)")
    v = Vector{Union{GridLandscape, Missing}}(undef, length(times))
    fill!(v, missing)
    a = DimArray(v, (Ti(times),))
    return CachedGridLandscape(a, file, saveinterval, timestep)
end

"""
    emptygridlandscape(habitat::GridHabitat, spplist::SpeciesList)

Create a [`GridLandscape`](@ref) of the right shape with every abundance zero, taking the species
names from `spplist` and the grid dimensions from `habitat`.
"""
function emptygridlandscape(habitat::GridHabitat, spplist::SpeciesList)
    mat = zeros(Int64, counttypes(spplist, true), countsubcommunities(habitat))
    return GridLandscape(mat, spplist.names, dims(habitat.active, (Y, X)))
end
