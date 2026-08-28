# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Everything the simulation knows about the species themselves: what they tolerate, what they
# demand, how they move, their demographic rates, and the phylogeny relating them.

using Diversity
using Missings
using Distributions
using DataFrames

"""
    SpeciesList{TL <: AbstractTolerance, DM <: AbstractDemand,
                MO <: AbstractMovement, T <: AbstractTypes,
                P <: AbstractParams} <: AbstractTypes

Everything the simulation knows about the species themselves — the species-side counterpart of a
[`GridHabitat`](@ref), and one of the two halves [`build_ecosystem`](@ref) checks against each other
name for name.

It is itself a `Diversity.AbstractTypes`, but it **holds** the real types in `types` rather than
being them, so every `AbstractTypes` question is forwarded there; see `src/DiversityInterface.jl`.

# Fields

  - `names`: one name per species.
  - `tolerance`: the niche of each species, an [`AbstractTolerance`](@ref) matching the habitat's
    regime axis for axis.
  - `demand`: what each species takes from the shared pool, an [`AbstractDemand`](@ref) matching the
    habitat's supply the same way.
  - `abun`: the abundances a run starts from, before `populate!` scatters them.
  - `types`: the similarity structure — a `UniqueTypes` where species are maximally distinct, or a
    `PhyloBranches` derived from a phylogeny.
  - `movement`: which individuals disperse and how far, an [`AbstractMovement`](@ref).
  - `params`: the demographic rates, an [`AbstractParams`](@ref).
  - `native`: whether each species is native, one flag each.
  - `susceptible`: optional per-species disease susceptibility, `missing` where unset.

# Type parameters

  - `TL`, `DM`, `MO`, `T`, `P`: the types of `tolerance`, `demand`, `movement`, `types` and `params`.
    `T` is what a method keys on to require a phylogeny, since `Diversity`'s own `Sim` slot cannot
    carry that constraint.
"""
mutable struct SpeciesList{TL <: AbstractTolerance,
                           DM <: AbstractDemand,
                           MO <: AbstractMovement,
                           T <: AbstractTypes,
                           P <: AbstractParams} <: AbstractSpecies
    names::Vector{String}
    tolerance::TL
    abun::Vector{Int64}
    demand::DM
    types::T
    movement::MO
    params::P
    native::Vector{Bool}
    susceptible::Vector{Union{Missing, Float64}}

    function SpeciesList{TL, DM, MO, T, P}(names::Vector{String},
                                           tolerance::TL,
                                           abun::Vector{Int64},
                                           demand::DM,
                                           types::T,
                                           movement::MO,
                                           params::P,
                                           native::Vector{Bool}) where {TL <:
                                                                        AbstractTolerance,
                                                                        DM <:
                                                                        AbstractDemand,
                                                                        MO <:
                                                                        AbstractMovement,
                                                                        T <:
                                                                        AbstractTypes,
                                                                        P <:
                                                                        AbstractParams}
        # Check dimensions
        equal_param = equalpop(params, length(names))
        sus = Vector{Union{Missing, Float64}}(undef, length(names))
        return new{TL, DM, MO, T, typeof(equal_param)}(names,
                                                       tolerance,
                                                       abun,
                                                       demand,
                                                       types,
                                                       movement,
                                                       equal_param,
                                                       native,
                                                       sus)
    end
end
# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# The species-side counterpart of `GridHabitat`'s display, and the asymmetry of having one without
# the other is the argument for it: an ecosystem is species in an environment, and both halves
# should answer "what have I got here?" the same way.
#
# The similarity structure is named because it is what decides whether a diversity measure means
# anything — `UniqueTypes` and a phylogeny give different answers to the same question.
function Base.show(io::IO, sl::SpeciesList)
    return print(io,
                 "SpeciesList($(length(sl.names)) species, ",
                 sprint(show, sl.tolerance), " | ", sprint(show, sl.demand),
                 ")")
end

function Base.show(io::IO, ::MIME"text/plain", sl::SpeciesList)
    println(io, "SpeciesList")
    println(io, "  species   ", length(sl.names))
    println(io, "  tolerance ", sprint(show, sl.tolerance))
    println(io, "  demand    ", sprint(show, sl.demand))
    println(io, "  movement  ", nameof(typeof(sl.movement)))
    println(io, "  params    ", nameof(typeof(sl.params)))
    print(io, "  types     ", nameof(typeof(sl.types)))
    return nothing
end

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

# Outer rather than inner, and that is the point: a constructor belongs inside a `struct` only to
# reach `new`, and this one derives the default names and delegates. Leaving it inside would imply a
# second route to `new` when there is only one.
function SpeciesList{TL, DM, MO, T, P}(tolerance::TL,
                                       abun::Vector{Int64},
                                       demand::DM,
                                       types::T,
                                       movement::MO,
                                       params::P,
                                       native::Vector{Bool}) where {TL <:
                                                                    AbstractTolerance,
                                                                    DM <:
                                                                    AbstractDemand,
                                                                    MO <:
                                                                    AbstractMovement,
                                                                    T <:
                                                                    AbstractTypes,
                                                                    P <:
                                                                    AbstractParams}
    return SpeciesList{TL, DM, MO, T, P}(map(x -> "$x", eachindex(abun)),
                                         tolerance,
                                         abun,
                                         demand,
                                         types,
                                         movement,
                                         params,
                                         native)
end

# The two tree-generating constructors live in `EcoSISTEMPhyloExt`. They do not merely accept a
# phylogeny — they build one and derive a `PhyloBranches` similarity from it — so they cannot work
# without `Phylo`. Their docstrings stay here, immediately above the stubs, because `api.md`'s
# `@autodocs` cannot see into an extension.
#
# Every other `SpeciesList` constructor is unaffected: they take the tolerances, and the similarity,
# explicitly, which is the path `build_species` uses.
"""
    SpeciesList(numspecies::Int64, numtraits::Int64, abun::Vector{Int64},
      demand::DM, movement::MO, params::P, native::Vector{Bool})
    SpeciesList(numspecies::Int64, numtraits::Int64, abun::Vector{Int64},
      demand::DM, movement::MO, params::P, native::Vector{Bool},
      switch::Vector{Float64})
    SpeciesList(numspecies::Int64, numtraits::Int64, abun::Vector{Int64},
      demand::DM, movement::MO, phy::T, params::P, native::Vector{Bool})

Create a `SpeciesList` whose categorical niche tolerances are **evolved along a random ultrametric
phylogeny**, with a `PhyloBranches` similarity structure computed from that tree.

The third form takes a similarity structure explicitly instead of computing one from the tree; the
tolerances are still evolved along a random ultrametric tree.

**All three forms need `Phylo`** — each builds the tree, or forwards to one that does — so
`using Phylo` is required and without it the call is a `MethodError`. The implementation is in
`EcoSISTEMPhyloExt`.

# Arguments

  - `numspecies`: how many species.
  - `numtraits`: how many categorical niche classes the tolerances are drawn from.
  - `abun`: the starting abundances.
  - `demand`: what each species takes from the shared pool.
  - `movement`: which individuals disperse and how far.
  - `params`: the demographic rates.
  - `native`: whether each species is native.
  - `switch`: the rate of trait change along branches. Defaults to `[0.5]`.
  - `phy`: a similarity structure to use in place of the one the tree would give.
"""
function SpeciesList(numspecies::Int64,
                     numtraits::Int64,
                     abun::Vector{Int64},
                     demand::DM,
                     movement::MO,
                     params::P,
                     native::Vector{Bool}) where {DM <: AbstractDemand,
                                                  MO <: AbstractMovement,
                                                  P <: AbstractParams}
    return SpeciesList(numspecies, numtraits, abun, demand, movement, params,
                       native, [0.5])
end
# ---------------------------------------------------------------------------
# The size-based (`SizeDemand`) constructor is parked — still commented out — in
# `ext/EcoSISTEMPhyloExt/`. It evolves body mass along a phylogeny, so were it ever uncommented it
# would need `Phylo` like the two constructors above, and keeping it beside them keeps the open
# decision about `SizeDemand` in one place rather than two. See `src/Demand.jl`.

"""
    SpeciesList(numspecies::Int64, tolerance::TL, abun::Vector{Int64}, demand::DM,
      movement::MO, params::P, native::Vector{Bool})

Create a `SpeciesList` from tolerances supplied directly, rather than evolved along a phylogeny.
Species are treated as maximally distinct, with a `UniqueTypes` similarity structure.

# Arguments

  - `numspecies`: how many species.
  - `tolerance`: their niches, as an [`AbstractTolerance`](@ref).
  - `abun`: the starting abundances. Padded with zeros if shorter than `numspecies`.
  - `demand`: what each species takes from the shared pool.
  - `movement`: which individuals disperse and how far.
  - `params`: the demographic rates.
  - `native`: whether each species is native.
"""
function SpeciesList(numspecies::Int64,
                     tolerance::TL,
                     abun::Vector{Int64},
                     demand::DM,
                     movement::MO,
                     params::P,
                     native::Vector{Bool}) where {TL <: AbstractTolerance,
                                                  DM <: AbstractDemand,
                                                  MO <: AbstractMovement,
                                                  P <: AbstractParams}
    names = map(x -> "$x", 1:numspecies)
    # Create similarity matrix (for now identity)
    ty = UniqueTypes(numspecies)
    # Draw random set of abundances from distribution
    if length(abun) < numspecies
        abun = vcat(abun, fill(0, numspecies - length(abun)))
    end
    # error out when abun dist and NumberSpecies are not the same (same for resource dist)
    length(abun) == numspecies || throw(DimensionMismatch("Abundance vector
                                            doesn't match number species"))
    counttypes(demand) == numspecies || throw(DimensionMismatch("Demand vector
                                            doesn't match number species"))
    return SpeciesList{typeof(tolerance),
                       typeof(demand),
                       typeof(movement),
                       typeof(ty),
                       typeof(params)}(names,
                                       tolerance,
                                       abun,
                                       demand,
                                       ty,
                                       movement,
                                       params,
                                       native)
end

# The species-side half of the four-way accessor family; the docstring covering both methods lives
# with the family in `Ecosystem.jl`. This returns the **demand** itself. The total resource usage is
# `_getdemand(abun, demand)` in `Demand.jl`, which is what the hot loop reads.
getdemand(sppl::SpeciesList) = sppl.demand

# The only route from this package to a `UniqueTypes`, so that a species' name cannot be lost on the
# way. Diversity offers `UniqueTypes(count)` as well as `UniqueTypes(names)`, and the count form
# generates `"1"`, `"2"`, ... — right for a genuinely anonymous community, wrong wherever names are
# held, as they always are here. Leaving nothing that can reach the lossy form is what keeps that from
# happening twice.
#
# The renaming is invisible on the synthetic route, where `build_species` already yields numeric
# names. It bites a phylogeny-built `SpeciesList`, whose types carry tip labels, so a test of this has
# to use real names to prove anything.
function _uniquetypes(names::AbstractVector{<:AbstractString})
    return UniqueTypes(collect(String, names))
end

# --- Per-species fields from what a caller wrote ------------------------------

# Turn a per-species argument into a length-`n` vector: a scalar is filled to
# every species, a vector is passed through after checking it has exactly `n`
# entries (erroring clearly, and naming the argument, if it does not).
function _tofield(x::AbstractVector, n::Integer, name::AbstractString)
    length(x) == n ||
        throw(DimensionMismatch("`$name` has $(length(x)) entries but there are $n species"))
    return collect(x)
end

_tofield(x, n::Integer, name::AbstractString) = fill(x, n)

# A per-species abundance vector: an integer total is split at random across
# species (seedable); an explicit vector is validated and used as-is.
function _abundances(abundance::AbstractVector, n::Integer, seed)
    length(abundance) == n ||
        throw(DimensionMismatch("`abundance` has $(length(abundance)) entries but there are $n species"))
    return Vector{Int64}(abundance)
end

function _abundances(abundance, n::Integer, seed)
    isnothing(seed) || Random.seed!(seed)
    return rand(Multinomial(Int(abundance), n))
end
