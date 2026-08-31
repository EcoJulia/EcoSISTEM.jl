# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a species takes from the shared pool - the `Resource` mirror of a tolerance.

using RecipesBase
using Unitful
using Unitful.DefaultSymbols
using DimensionalData
using DimensionalData.Lookups: NoLookup

"""
    Demand{A <: NicheAxis, V, X} <: AbstractDemand{A, V}

What each species takes from the shared pool on one niche axis - the `Resource` mirror of a
tolerance, and one of the two halves that must line up with the environment.

Name the axis to build one, `Demand{SolarRadiation}(values)`. The values are converted to that axis's
own canonical resource unit - `kJ/day` for [`SolarRadiation`](@ref), `L/day` for
[`Precipitation`](@ref), `g/day` for [`CarbonFlux`](@ref) - so any scale of the right dimension is
accepted, while a wrong dimension, or an axis that declares no resource at all, is refused at the
call rather than deep in the simulation loop.

# Fields

  - `resource`: one demand per species, in the axis's canonical unit.
  - `exchange_rate`: what a demand is measured against. Defaults to `1/mean(resource)`, so a species'
    demand is read relative to the assemblage it is part of rather than as an absolute. It scales
    birth and death equally, so it sets the tempo of turnover without moving any equilibrium.

# Type parameters

  - `A`: the niche axis these demands are stated on.
  - `V`: the value type, which follows from `canonicalunit(Resource, A)`.
  - `X`: the exchange rate's type, the inverse of `V`.
"""
struct Demand{A <: NicheAxis, V, X} <: AbstractDemand{A, V}
    resource::Vector{V}
    exchange_rate::X

    # The axis is the only declaration of the unit, exactly as on the supply side: the value type
    # follows from `canonicalunit(Resource, A)` rather than being pinned in the type, so the axis and
    # the unit cannot be two statements that disagree.
    #
    # The cost is that dispatch can no longer refuse a bad unit on its own, so the check is made here:
    # `_resourceunit` refuses a wrong dimension **and** an axis that is not a resource, naming which
    # of the two mistakes it is.
    function Demand{A}(resource::AbstractVector{<:Unitful.Quantity{Float64}},
                       exchange_rate::Unitful.Quantity{Float64} = 1.0 /
                                                                  mean(resource)) where {A <:
                                                                                         NicheAxis}
        u = _resourceunit(eltype(resource), A, "demand")
        r = collect(uconvert.(u, resource))
        x = uconvert(inv(u), exchange_rate)
        return new{A, eltype(r), typeof(x)}(r, x)
    end
end
# The species-side mirror of `Supply{A}`'s display in `Layer.jl`, with species where a supply has
# cells. `exchange_rate` is deliberately absent: it is community-relative machinery rather than
# something a reader identifies a demand by.
function Base.show(io::IO, d::Demand)
    return print(io,
                 "Demand{$(nameof(axisof(d)))}($(length(d.resource)) species, $(unit(eltype(d))))")
end

# == Functions ==================================================================================

# A bare-number vector is refused here rather than by the inner constructor's signature. Left to
# dispatch it gives a `MethodError` naming `Vector{Float64}` - true, but saying nothing about what is
# wrong. Routing it through the same check a united vector takes gets it the same sentence, naming
# the unit a demand on this axis is measured in.
#
# The demand-side mirror of `_tocanonu(x::Real, u)` on the regime side: a bare number carries no
# statement of what it measures, so it is refused wherever it appears rather than being stamped with
# whatever unit happened to be expected.
function Demand{A}(resource::AbstractVector{<:Real}) where {A <:
                                                            NicheAxis}
    _resourceunit(eltype(resource), A, "demand")
    return error("a demand on axis $A cannot be built from bare numbers.")   # unreachable
end

# How many resources a demand asks for - the width of the `totaldemand` cache, one column per resource.
# A single demand asks for one; a collection asks for one per member (a `SpeciesRequirementCollection`
# below), read straight off the backing's type so the cache is sized without an instance.
numdemands(::Type{<:Demand}) = 1

# A cell's total demand for one resource: each species' per-individual need times how many are
# there. This is the `E` the hot loop divides supply by, and the reason `getdemand` returns the
# demand itself rather than a number - the total is a different question, asked here.
_getdemand(abun::Vector{Int64}, demand::Demand) = sum(abun .* demand.resource)

# The niche axis a demand is stated on - the counterpart of `axisof` for a supply, and what
# `_checkaligned` would use to compare the two sides by axis rather than by stored unit.
axisof(::Type{<:Demand{A}}) where {A <: NicheAxis} = A
axisof(demand::Demand) = axisof(typeof(demand))

# ---------------------------------------------------------------------------
# `SizeDemand` - commented out, deliberately, pending a decision
# ---------------------------------------------------------------------------
# This is parked code, not dead code, and it must not be quietly deleted or quietly revived.
#
# `SizeDemand` was a demand whose per-species values are body size evolved along the phylogeny,
# used directly as a rate in the free/dimensionless currency. When that currency's *supply* was
# removed it became unusable: a species' demands must line up with the environment's supplies
# (`_checkaligned`), and no `𝐓^-1` supply exists any more - measured, it constructed happily and
# every `build_ecosystem` then refused it.
#
# Why it is kept at all: its two extra fields, `pop_mass_rel` and `area`, belong to a **metabolic**
# size-scaling reading - Kleiber-style - which is a real modelling idea nobody has implemented.
# Deleting the code would throw that away along with the allometry already worked out in the
# `SpeciesList` constructor commented out beside it (`Species.jl`).
#
# The open decision is to implement the metabolic reading, which needs an axis of its own so that it
# can be a `Demand{...}` like the rest, or to delete this outright. Until then it is neither exported,
# nor deprecated, nor reachable.
#
# The `SpeciesList(numspecies, numtraits, pop_mass, mean, var, area, ...)` constructor is commented
# out with it, being its only caller and unable to work without it.
#
# """
#     SizeDemand <: AbstractDemand{typeof(1.0/day)}
#
# A resource demand scaled by species size - a per-species free-resource rate for each
# species, `/day` (`pop_mass_rel`/`area` are used only once, at
# initial-abundance construction time in `SpeciesList`, not by the ongoing demand:supply
# balance).
# """
# struct SizeDemand <: AbstractDemand{_SimpleRate}
#     resource::Vector{_SimpleRate}
#     pop_mass_rel::Float64
#     area::Unitful.Area
#     exchange_rate::typeof(1.0 * day)
#
#     function SizeDemand(resource::Vector{<:Real},
#                         pop_mass_rel::Float64,
#                         area::Unitful.Area,
#                         exchange_rate::Unitful.Quantity{Float64} = 1.0 * day)
#         return new(resource .* unit(_SimpleRate), pop_mass_rel, area,
#                    uconvert(day, exchange_rate))
#     end
# end
#
# Base.length(demand::SizeDemand) = length(demand.resource)
# function _getdemand(abun::Vector{Int64}, demand::SizeDemand)
#     return sum(abun .* demand.resource)
# end

# Worth knowing when demanding `CarbonFlux`: NPP is *potential* productivity estimated from climate,
# by the Miami model, so carbon, solar radiation and water are not independent - an NPP supply
# already embeds the light and water a plant could use, and demanding all three counts one limitation
# twice. Nothing forbids it, since supplies combine by Liebig's minimum and a generous resource never
# binds, but it should be a deliberate choice.

# Several demands are a `SpeciesRequirementCollection` whose members happen to be demands, exactly as
# several supplies are a `LayerCollection{Resource}`; the role is read off the members. A collection
# implements the standard container interface (`src/collections.jl`), so write `values(x)`, `keys(x)`
# and `NamedTuple(x)`, which work identically on a leaf as a one-member container. There is
# deliberately no `Base.length` on the collection: on a demand `length` already means the species
# count, not the arity.

# The demand as a side of a pairing check - see `_checkaligned` (`collections.jl`). No `kinds`, for
# the same reason as the supply it is matched against.
_demandside(d::AbstractDemand) = _side(d, "species demand", false)

function numdemands(::Type{SpeciesRequirementCollection{Resource, A, C}}) where {A,
                                                                                 C}
    return fieldcount(C)
end

function _getdemand(abun::Vector{Int64},
                    demand::SpeciesRequirementCollection{Resource})
    return [_getdemand(abun, d) for d in values(demand)]
end

# A plot title for a supply, keyed on its **axis** - the same rule as everything else on this path,
# and no longer on the value's dimension, which cannot tell water from wind. The three shipped
# resources get the wording a reader expects; the fallback names any other axis from its own
# declaration, so a new resource plots sensibly with no entry added here.
_resourcetitle(::Type{<:SolarRadiation}) = "Solar Radiation (kJ/day)"
_resourcetitle(::Type{<:Precipitation}) = "Available water (L/day)"
_resourcetitle(::Type{<:CarbonFlux}) = "Net primary productivity (g/day)"
function _resourcetitle(::Type{A}) where {A <: NicheAxis}
    return "$(nameof(A)) ($(canonicalunit(Resource, A)))"
end

# --- Building a demand from what a caller wrote -------------------------------

# Ask the axis for its demand type, refusing an axis that is not a resource in the one wording
# `_notaresource` owns. Without this `demandtype` returns `nothing` and the caller invokes it, giving
# `MethodError: objects of type Nothing are not callable` - a sentinel every call site would have to
# remember to test.
function _demandtype(axis)
    T = demandtype(axis)
    isnothing(T) && _notaresource(axis, "demand")
    return T
end

# Build a demand on the **named** axis. The axis is passed in rather than inferred from the values'
# unit: `demandtype` is asked of the axis, which is the direction that cannot be wrong, since two
# axes may legitimately share a unit.
function _demand(axis, resource::AbstractVector)
    return _demandtype(axis)(float.(resource))
end
