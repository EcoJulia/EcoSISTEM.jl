# SPDX-License-Identifier: LGPL-3.0-or-later

import Base.eltype

"""
    AbstractTolerance{T}

Abstract supertype for all tolerance types, parameterised by their response type of any type `T`.
"""
abstract type AbstractTolerance{T} end

eltype(::AbstractTolerance{D}) where {D} = D

"""
    BasicTrait{T} <: AbstractTolerance{T}

Basic trait type that holds information on a single trait for each species, of
any type `T`.
"""
struct DiscreteTolerance{D} <: AbstractTolerance{D}
    val::Vector{D}
end

iscontinuous(trait::DiscreteTolerance) = false

struct LandCoverTolerance{D <: Number} <: AbstractTolerance{D}
    vals::Vector{Vector{D}}
end

iscontinuous(trait::LandCoverTolerance) = false

function LandCoverTolerance(vals::Vector{Vector{<:AbstractFloat}})
    return LandCoverTolerance{typeof(1.0)}(vals)
end

"""
    DiscreteEvolve(numTraits::Int64, tree::BinaryTree)

Evolve a discrete switching trait along a BinaryTree, `tree`. Takes in a number
of traits, `numTraits` to be switched between and rate to switch between traits,
`switch_rate` with default value of 0.5.
"""
function DiscreteEvolve(numTraits::Int64, tree::BinaryTree, switch_rate = 0.5)
    # Create traits and assign to tips
    traits = DataFrame(trait1 = 1:numTraits)
    assign_traits!(tree, 0.5, traits)
    # Get traits from tree
    return DiscreteTolerance(Array(get_traits(tree, true)[:, 1]))
end

"""
    ContinuousEvolve(val::Union{Float64, Unitful.Quantity{Float64}}, var::Union{Float64, Unitful.Quantity{Float64}}, tree::BinaryTree)

Evolve a continuous trait along a BinaryTree, `tree` via Brownian motion. Takes
in a starting value, `val` and a variance, `var`.
"""
function ContinuousEvolve(val::Union{Float64, Unitful.Quantity{Float64}},
                          var::Union{Float64, Unitful.Quantity{Float64}},
                          tree::BinaryTree)
    # Create traits and assign to tips
    numspecies = length(getleafnames(tree))
    traits = DataFrame(start = ustrip(val), σ² = ustrip(var))
    assign_traits!(tree, traits)
    # Get traits from tree
    newtraits = get_traits(tree, true)
    newtraits[!, :start] = newtraits[!, :start] .* unit(val)
    return NicheTolerance(Unclassified, Normal, newtraits[!, :start],
                          newtraits[!, :σ²])
end

"""
    ContinuousTolerance{C <: Number} <: AbstractTolerance{T}

Abstract trait type that holds information on a single continuous trait for each
species, of any Number type `C`.
"""
abstract type ContinuousTolerance{C <: Number} <: AbstractTolerance{C} end

"""
    NicheTolerance{A <: NicheAxis, C, D} <: ContinuousTolerance{C}

Trait type holding a **binned** regime preference for each species on niche axis `A`, as one **built**
continuous response distribution per species — `dists::Vector{D}`, where `D` is a concrete
`Distributions.ContinuousUnivariateDistribution` (e.g. `Trapezoid{Float64}` or `Uniform{Float64}`). The
distributions are built once in the **support frame** `C` — the unit their support/domain is measured in,
which is also the trait's `eltype` and the unit the matching regime must be in. In the hot loop the
[`NicheSuitability`](@ref) nichefit fetches `dists[sp]` — see [`getdist`](@ref) — and evaluates the density at
`ustrip(regime)`, already a bare number in frame `C`, so there is no per-call conversion or allocation.
"""
struct NicheTolerance{A <: NicheAxis, C, D} <: ContinuousTolerance{C}
    dists::Vector{D}
end

# Build a `NicheTolerance` whose distributions live in the `frame` unit (= its `eltype` / the required regime unit),
# reading each row's bare parameters as being in `input_unit` and converting to `frame` *per role* (a
# location properly/affine, a scale as an interval, a shape dimensionless — see `read_distribution`; a
# shape-only family is placed via `offset`/`scale`). The single place the storage frame is fixed.
function _buildniche(::Type{A}, ::Type{D}, mat::AbstractMatrix, input_unit,
                     frame;
                     offset, scale, probes) where {A <: NicheAxis, D}
    roles = param_roles_resolved(D; probes = probes)
    dists = [read_distribution(D, input_unit, mat[sp, :]; canonical = frame,
                               offset = offset, scale = scale, roles = roles)
             for sp in Base.axes(mat, 1)]
    return NicheTolerance{A, typeof(1.0 * frame), eltype(dists)}(dists)
end

"""
    NicheTolerance(::Type{A}, ::Type{D}, dist::Matrix; support = canonicalunit(A()),
        offset = nothing, scale = nothing, probes = …)

Build a [`NicheTolerance`](@ref) on axis `A` with response distribution `D` from a **bare** parameter matrix (one
species per row). `support` is the **frame** the distribution is built in — the unit its support/domain is
measured in, the trait's `eltype`, and the unit the matching regime must be in (defaults to the axis's
canonical unit). The bare numbers are taken to be in that frame. `offset`/`scale` (references in `support`)
place a *shape-only* distribution (`Beta`, `LogNormal`, …) on the dimensioned axis via a `LocationScale`.
"""
function NicheTolerance(::Type{A}, ::Type{D}, dist::Matrix;
                        support = canonicalunit(A()), offset = nothing,
                        scale = nothing,
                        probes = _default_probes(D)) where {A <: NicheAxis, D}
    dimension(support) == dimension(canonicalunit(A())) ||
        error("NicheTolerance support unit $support and axis $A's canonical unit " *
              "$(canonicalunit(A())) have different dimensions.")
    return _buildniche(A, D, dist, support, support; offset = offset,
                       scale = scale, probes = probes)
end

# Impute the *input* unit of a set of parameter vectors — the unit their (bare) magnitudes are read in:
# their single shared *dimensioned* unit (dimensionless/shape vectors ignored), or the `support` frame when
# all bare (so bare vectors are read in the frame, like the matrix constructor). Mixed units error.
function _impute_input_unit(params, support)
    us = unique(unit(eltype(v)) for v in params if unit(eltype(v)) !== NoUnits)
    isempty(us) && return support
    length(us) == 1 ||
        error("NicheTolerance parameter vectors carry differing units $(us); pass values in one unit.")
    return only(us)
end

"""
    NicheTolerance(::Type{A}, ::Type{D}, params::AbstractVector...; support = canonicalunit(A()),
        offset = nothing, scale = nothing, probes = …)

Build a [`NicheTolerance`](@ref) on axis `A` with response distribution `D` from **one (unitful) vector per
parameter** of `D` — the ergonomic, programmatic counterpart of the `dist::Matrix` constructor (and of the
old `GaussTrait(A, mean, sd)`). Each `params` vector holds one distribution parameter across species, e.g.
`NicheTolerance(MeanTemperature, Normal, opts, vars)` (a `μ` and a `σ` vector) or
`NicheTolerance(Precipitation, Gamma, shape, scale_vec)`. The vectors' own units are respected (read per role); the
distribution is built in the `support` **frame** (the unit its support/domain — and the matching regime —
is measured in, default canonical), converting the inputs to it. So `support = K` and `support = u"°C"`
build the *same* preference in the K and °C frames respectively (`pdf(getdist(bin, sp), ustrip(support, x))`
is correct in each); the frame is the trait's `eltype`, and a NicheTolerance only matches a regime in that same unit.
All parameter vectors must share a unit (a mixed set errors) and have one entry per species.
"""
function NicheTolerance(::Type{A}, ::Type{D}, params::AbstractVector...;
                        support = canonicalunit(A()), offset = nothing,
                        scale = nothing,
                        probes = _default_probes(D)) where {A <: NicheAxis, D}
    isempty(params) &&
        error("NicheTolerance needs at least one parameter vector for $D.")
    n = length(first(params))
    all(v -> length(v) == n, params) ||
        error("NicheTolerance parameter vectors must all have the same length (one entry per species).")
    dimension(support) == dimension(canonicalunit(A())) ||
        error("NicheTolerance support unit $support and axis $A's canonical unit " *
              "$(canonicalunit(A())) have different dimensions.")
    input_unit = _impute_input_unit(params, support)
    # A unitful vector → its magnitude in `input_unit` (its role, applied in `_buildniche`, decides position
    # vs interval); a bare vector passes straight through.
    cols = map(v -> unit(eltype(v)) === NoUnits ? float.(v) :
                    ustrip.(input_unit, v),
               params)
    return _buildniche(A, D, reduce(hcat, cols), input_unit, support;
                       offset = offset, scale = scale, probes = probes)
end

iscontinuous(::NicheTolerance) = true

"""    TempTolerance — a binned temperature preference (a [`NicheTolerance`](@ref) on `MeanTemperature`, K frame, with a `Trapezoid` response). """
const TempTolerance = NicheTolerance{MeanTemperature,
                                     typeof(1.0 *
                                            canonicalunit(MeanTemperature())),
                                     Trapezoid{Float64}}
"""    RainTolerance — a binned rainfall preference (a [`NicheTolerance`](@ref) on `Precipitation`, mm frame, with a `Uniform` response). """
const RainTolerance = NicheTolerance{Precipitation,
                                     typeof(1.0 *
                                            canonicalunit(Precipitation())),
                                     Uniform{Float64}}

# Back-compat constructors preserving the old `TempTolerance(matrix)` / `RainTolerance(matrix)` call form.
TempTolerance(dist::Matrix) = NicheTolerance(MeanTemperature, Trapezoid, dist)
RainTolerance(dist::Matrix) = NicheTolerance(Precipitation, Uniform, dist)

# Convenience: the matching [`NicheSuitability`](@ref) nichefit for a `NicheTolerance`, taking its `TR` (unit) from the
# trait's axis — so callers building an ecosystem by hand need not re-type the unit.
NicheSuitability(t::NicheTolerance) = NicheSuitability{eltype(t)}()

"""
    ToleranceCollection2{T1, T2} <: AbstractTolerance{Tuple{T1, T2}}

Trait collection that holds two trait types, `TR1` and `TR2`.
"""
struct ToleranceCollection2{T1, T2} <: AbstractTolerance{Tuple{T1, T2}}
    one::T1
    two::T2
end

function iscontinuous(trait::ToleranceCollection2)
    return [iscontinuous(trait.one), iscontinuous(trait.two)]
end

eltype(trait::ToleranceCollection2) = [eltype(trait.one), eltype(trait.two)]

"""
    ToleranceCollection3{T1, T2, T3} <: AbstractTolerance{Tuple{T1, T2, T3}}

Trait collection that holds three trait types, `TR1`, `TR2` and `TR3`.
"""
struct ToleranceCollection3{T1, T2, T3} <: AbstractTolerance{Tuple{T1, T2, T3}}
    one::T1
    two::T2
    three::T3
end

function iscontinuous(trait::ToleranceCollection3)
    return [
        iscontinuous(trait.one),
        iscontinuous(trait.two),
        iscontinuous(trait.three)
    ]
end

function Base.eltype(trait::ToleranceCollection3)
    return [eltype(trait.one), eltype(trait.two), eltype(trait.three)]
end
