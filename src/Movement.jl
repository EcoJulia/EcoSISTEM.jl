# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How offspring disperse: the dispersal kernels, the movement rules that draw from them, and the
# per-species lookup table of dispersal probabilities that the hot loop actually reads.

using Unitful

using DataFrames: DataFrame, nrow

using SpecialFunctions: loggamma

"""
    AbstractKernel

The shape of the probability distribution over where an individual's offspring lands, relative to its
own cell. [`GaussianKernel`](@ref) and [`LongTailKernel`](@ref) are the concrete kernels.

A kernel says only *how far and how likely*. What happens when a draw crosses an edge is the grid's
[`AbstractTopology`](@ref), and what becomes of one aimed at a dead cell is `disperse_safely` on the
movement - three separate questions, answered by three separate things.
"""
abstract type AbstractKernel end

"""
    GaussianKernel(; meandist, sd, thresh = 1.0e-9)

A dispersal kernel whose displacement is Gaussian along each axis independently, so that the
distance travelled is Rayleigh distributed. The density of landing at distance `r` is
`exp(-(r / scale)^2) / (pi * scale^2)` per unit area.

Give exactly one of `meandist` and `sd`. They describe the same kernel: `meandist = sd *
sqrt(pi / 2)`.

# Keywords

  - `meandist`: the mean distance travelled.
  - `sd`: the standard deviation of the displacement along each axis.
  - `thresh`: the probability below which a destination is dropped, which is what bounds the lookup
    table to a finite neighbourhood.

# Fields

  - `scale`: the distance scale in the density above, in kilometres: `sqrt(2) * sd`, the
    root-mean-square distance travelled.
  - `thresh`: as above.
"""
struct GaussianKernel <: AbstractKernel
    scale::typeof(1.0km)
    thresh::Float64

    function GaussianKernel(;
                            meandist::Union{Nothing, Unitful.Length} = nothing,
                            sd::Union{Nothing, Unitful.Length} = nothing,
                            thresh::Real = 1.0e-9)
        _onedistance(GaussianKernel, (meandist = meandist, sd = sd))
        # Do not rearrange the mean's conversion: a kernel given its mean distance must build the
        # positional form's lookup to the bit, or the canonical references move.
        scale = isnothing(sd) ?
                (2 * _kilometres(meandist, :meandist)) / sqrt(pi) :
                sqrt(2) * _kilometres(sd, :sd)
        return new(scale, Float64(thresh))
    end
end

function Base.show(io::IO, k::GaussianKernel)
    return print(io, "GaussianKernel(meandist = ", _kernelphrase(_meandist(k)),
                 ", thresh = ", k.thresh, ")")
end

"""
    LongTailKernel(; shape, scale, meandist, thresh = 1.0e-9)

A dispersal kernel following the two-dimensional Student's t distribution, the "2Dt" of Clark et
al. (1999, Ecology 80: 1475-1494), whose tail falls off as a power of distance rather than
exponentially, so that rare long-distance dispersal is not effectively impossible. The density of
landing at distance `r` is `(shape - 1) / (pi * scale^2) * (1 + (r / scale)^2)^-shape` per unit
area.

Give `shape` and exactly one of `scale` and `meandist`, which are not equivalent: the mean distance
travelled is `scale * sqrt(pi) / 2 * gamma(shape - 3/2) / gamma(shape - 1)`, which is finite only for
`shape > 3/2`, and its variance is finite only for `shape > 2`. The same `scale` with a smaller
`shape` reaches further.

# Keywords

  - `shape`: how heavy the tail is, greater than 1, without which the density does not integrate.
  - `scale`: the distance scale in the density above.
  - `meandist`: the mean distance travelled, refused for `shape <= 3/2`, where it is infinite.
  - `thresh`: the probability below which a destination is dropped, which is what bounds the lookup
    table to a finite neighbourhood.

# Fields

  - `scale`: as above, in kilometres.
  - `shape`, `thresh`: as above.
"""
struct LongTailKernel <: AbstractKernel
    scale::typeof(1.0km)
    shape::Float64
    thresh::Float64

    function LongTailKernel(; shape::Real,
                            scale::Union{Nothing, Unitful.Length} = nothing,
                            meandist::Union{Nothing, Unitful.Length} = nothing,
                            thresh::Real = 1.0e-9)
        shape > 1 ||
            throw(ArgumentError("a `LongTailKernel`'s `shape` must exceed 1, or its density " *
                                "does not integrate; got $shape."))
        _onedistance(LongTailKernel, (scale = scale, meandist = meandist))
        if isnothing(meandist)
            a = _kilometres(scale, :scale)
        else
            shape > 3 / 2 ||
                throw(ArgumentError("a `LongTailKernel` of `shape` $shape has no finite mean " *
                                    "distance, which needs `shape` above 3/2; give its " *
                                    "`scale` instead."))
            a = 2 * _kilometres(meandist, :meandist) /
                (sqrt(pi) * _gammaratio(shape))
        end
        return new(a, Float64(shape), Float64(thresh))
    end
end

function Base.show(io::IO, k::LongTailKernel)
    return print(io, "LongTailKernel(shape = ", k.shape, ", scale = ",
                 _kernelphrase(k.scale), ", thresh = ", k.thresh, ")")
end

"""
    AbstractMovement

Which individuals disperse at all: [`BirthOnlyMovement`](@ref) (only the newly born, plant-like),
[`AlwaysMovement`](@ref) (any individual, animal-like) or [`NoMovement`](@ref) (none).

A movement carries the per-species part of dispersal - one [`AbstractKernel`](@ref) per species, and
one `disperse_safely` flag per species. The grid's own edges are not part of it: they belong to the
[`GridHabitat`](@ref), because two species on one grid cannot be on different topologies.

`disperse_safely` says what becomes of an individual dispersing into a **dead cell** - one off the
grid, or an inactive one. `true` redistributes it among the reachable destinations as though it had
never aimed there; `false` loses it, so the share of the kernel pointing at dead cells becomes
mortality.

**It is per species because it is a fact about the disperser**, not about the map. A wind-dispersed
seed blown out to sea is gone; an animal-dispersed one is carried to somewhere the animal can
actually reach, so two species on one grid may legitimately differ. It cannot empty a cell outright:
staying put is one of the kernel's own destinations, and a cell only disperses while it is active, so
some share always survives.
"""
abstract type AbstractMovement end

"""
    BirthOnlyMovement{K <: AbstractKernel} <: AbstractMovement

Only individuals that have just been born disperse - plant-like.

# Fields

  - `kernels`: one [`AbstractKernel`](@ref) per species.
  - `disperse_safely`: one flag per species, as described on [`AbstractMovement`](@ref). Defaults to
    `true` for every species.
"""
struct BirthOnlyMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
    disperse_safely::Vector{Bool}

    function BirthOnlyMovement(kernels::Vector{K},
                               disperse_safely::AbstractVector{Bool} = fill(true,
                                                                            length(kernels))) where {K <:
                                                                                                     AbstractKernel}
        _checkdispersesafely(kernels, disperse_safely)
        return new{K}(kernels, collect(disperse_safely))
    end
end

"""
    AlwaysMovement{K <: AbstractKernel} <: AbstractMovement

Any individual may disperse, not only the newly born - animal-like.

# Fields

  - `kernels`: one [`AbstractKernel`](@ref) per species.
  - `disperse_safely`: one flag per species, as described on [`AbstractMovement`](@ref). Defaults to
    `true` for every species.
"""
struct AlwaysMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
    disperse_safely::Vector{Bool}

    function AlwaysMovement(kernels::Vector{K},
                            disperse_safely::AbstractVector{Bool} = fill(true,
                                                                         length(kernels))) where {K <:
                                                                                                  AbstractKernel}
        _checkdispersesafely(kernels, disperse_safely)
        return new{K}(kernels, collect(disperse_safely))
    end
end

"""
    NoMovement{K <: AbstractKernel} <: AbstractMovement

Nothing disperses: every individual stays in the cell it was born in.

# Fields

  - `kernels`: one [`AbstractKernel`](@ref) per species, carried but never drawn from, so that a
    model can be switched between movements without rebuilding its species.
"""
struct NoMovement{K <: AbstractKernel} <: AbstractMovement
    kernels::Vector{K}
end

# Accepts (and ignores) a `disperse_safely` vector, so `NoMovement` can be constructed the same way
# as `BirthOnlyMovement`/`AlwaysMovement` when the movement type is chosen as a value and called
# uniformly (`movement(kernels, disperse_safely)`). Nothing disperses, so nothing can be lost.
function NoMovement(kernels::Vector{K},
                    ::AbstractVector{Bool}) where {K <: AbstractKernel}
    return NoMovement(kernels)
end
@doc (@doc NoMovement) NoMovement(::Vector{K} where {K <: AbstractKernel},
                                  ::AbstractVector{Bool})

"""
    Lookup

One species' dispersal neighbourhood, precomputed once and read by the hot loop: every destination
its kernel can reach, as an offset from the source cell, with the probability of landing there.

The last two fields are scratch, rewritten on every cell of every timestep. Holding them here rather
than allocating them per call is what keeps the movement step allocation-free.

# Fields

  - `y`, `x`: the destination offsets, relative to the cell dispersing. In that order, as everywhere
    else in the package.
  - `p`: the kernel's probability for each, fixed for the run.
  - `pnew`: those probabilities renormalised over the destinations actually available from this cell.
  - `moves`: how many individuals go to each destination this step.
"""
struct Lookup
    y::Vector{Int64}
    x::Vector{Int64}
    p::Vector{Float64}
    pnew::Vector{Float64}
    moves::Vector{Int64}
end

# The only place `Lookup`'s fields are filled positionally, so the only place their `(y, x)` order
# has to be got right by hand.
function Lookup(df::DataFrame)
    return Lookup(df[!, :Y],
                  df[!, :X],
                  df[!, :Prob],
                  zeros(Float64, nrow(df)),
                  zeros(Int64, nrow(df)))
end

# One line, because the default prints all five vectors - measured at 230 320 characters for a
# single species' neighbourhood. What identifies a lookup is how far it reaches, not the numbers in
# it: the last three fields are scratch rewritten on every cell of every timestep.
function Base.show(io::IO, l::Lookup)
    return print(io, "Lookup($(length(l.y)) destinations)")
end

# One line, because the default prints every kernel, tens of thousands of characters at a thousand
# species. What identifies a movement is its kind, how many species it covers and how far they reach.
# Written once on the supertype, since `getkernels` is the whole interface.
function Base.show(io::IO, m::AbstractMovement)
    kernels = getkernels(m)
    return print(io, nameof(typeof(m)), "(", length(kernels), " species, ",
                 isempty(kernels) ? "no kernels" :
                 "$(nameof(eltype(kernels))) $(_dispersalphrase(kernels))",
                 _unsafephrase(m), ")")
end

function Base.show(io::IO, ::MIME"text/plain", m::AbstractMovement)
    kernels = getkernels(m)
    println(io, nameof(typeof(m)))
    println(io, "  species          ", length(kernels))
    println(io, "  kernel           ",
            isempty(kernels) ? "none" : nameof(eltype(kernels)))
    println(io, "  dispersal        ", _dispersalphrase(kernels))
    print(io, "  disperse_safely  ", _safephrase(m))
    return nothing
end

# == Functions ==================================================================================

"""
    getkernels(m::AbstractMovement)

Return the vector of dispersal kernels a movement holds, one per species.
"""
getkernels(m::BirthOnlyMovement) = m.kernels
getkernels(m::AlwaysMovement) = m.kernels
getkernels(m::NoMovement) = m.kernels

"""
    dispersesafely(m::AbstractMovement, sp::Integer)

Return whether species `sp` is redistributed (`true`) or lost (`false`) when its dispersal is aimed
at a dead cell. `NoMovement` answers `true`, since nothing disperses and so nothing can be lost.
"""
dispersesafely(m::BirthOnlyMovement, sp::Integer) = m.disperse_safely[sp]
dispersesafely(m::AlwaysMovement, sp::Integer) = m.disperse_safely[sp]
dispersesafely(::NoMovement, ::Integer) = true

# A vector's values as one phrase: the single value where all agree, else `lo to hi`. Quantities
# and floats are rounded to three digits, so a 2.4 km kernel reads as `2.4 km`.
function _rangephrase(v)
    isempty(v) && return "none"
    lo, hi = extrema(v)
    lo == hi && return _roundphrase(lo)
    return _roundphrase(lo) * " to " * _roundphrase(hi)
end
function _roundphrase(x::Unitful.Quantity)
    return string(round(typeof(1.0 * unit(x)), x, digits = 3))
end
_roundphrase(x::Real) = string(round(x, digits = 3))

# How far a movement's kernels reach, for its display: a Gaussian kernel's mean distance, and a
# long-tailed kernel's scale, since its mean may be infinite.
function _dispersalphrase(kernels::AbstractVector)
    return isempty(kernels) ? "none" : _reachphrase(kernels)
end
function _reachphrase(kernels::AbstractVector{GaussianKernel})
    return "mean " * _rangephrase(map(_meandist, kernels))
end
function _reachphrase(kernels::AbstractVector)
    return "scale " * _rangephrase(map(k -> k.scale, kernels))
end

# A kernel's distance in its own display, to six significant figures, since the mean is recovered
# from the stored scale and may differ from the distance given in its last bit.
function _kernelphrase(d::Unitful.Length)
    return string(round(typeof(1.0km), d, sigdigits = 6))
end

# A Gaussian kernel's mean distance, from its stored scale.
_meandist(k::GaussianKernel) = k.scale * sqrt(pi) / 2

# The ratio `gamma(shape - 3/2) / gamma(shape - 1)` in a long-tailed kernel's mean distance, through
# the log-gamma function so that a large shape cannot overflow.
_gammaratio(shape::Real) = exp(loggamma(shape - 3 / 2) - loggamma(shape - 1))

# Refuse a kernel given other than exactly one of its alternative distances, naming them.
function _onedistance(::Type{K}, given::NamedTuple) where {K <: AbstractKernel}
    count(!isnothing, values(given)) == 1 && return nothing
    names = join(("`$k`" for k in keys(given)), " and ")
    got = all(isnothing, values(given)) ? "neither" : "both"
    return throw(ArgumentError("a `$(nameof(K))` takes exactly one of $names; got $got."))
end

# A kernel's distance keyword in kilometres, refused unless positive.
function _kilometres(d::Unitful.Length, name::Symbol)
    d > zero(d) ||
        throw(ArgumentError("a kernel's `$name` must be positive; got $d."))
    return uconvert(km, d)
end

# How many species disperse unsafely, for the compact form: nothing when none do. `NoMovement`
# carries no flags, so it has nothing to say.
_unsafephrase(::NoMovement) = ""
function _unsafephrase(m::AbstractMovement)
    n = count(!, m.disperse_safely)
    return n == 0 ? "" : ", $n not disperse_safely"
end

# The `disperse_safely` line of the display form.
_safephrase(::NoMovement) = "not applicable, nothing disperses"
function _safephrase(m::AbstractMovement)
    n = count(m.disperse_safely)
    return n == length(m.disperse_safely) ? "all" :
           "$n of $(length(m.disperse_safely)) species"
end

# One flag per kernel, or the two vectors silently describe different species - `zip` would truncate
# and the last species would take a neighbour's setting.
function _checkdispersesafely(kernels, disperse_safely)
    length(disperse_safely) == length(kernels) ||
        error("`disperse_safely` has $(length(disperse_safely)) entries but there are " *
              "$(length(kernels)) dispersal kernels: it is a per-species setting, so it needs one " *
              "entry per species (or a single value applied to all).")
    return nothing
end
