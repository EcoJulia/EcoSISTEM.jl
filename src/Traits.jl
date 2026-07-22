# SPDX-License-Identifier: LGPL-3.0-or-later

import Base.eltype

"""
    AbstractTraits{T}

Abstract supertype for all trait types, parameterised by traits of any type `T`.
"""
abstract type AbstractTraits{T} end

eltype(::AbstractTraits{D}) where {D} = D

"""
    BasicTrait{T} <: AbstractTraits{T}

Basic trait type that holds information on a single trait for each species, of
any type `T`.
"""
struct DiscreteTrait{D} <: AbstractTraits{D}
    val::Vector{D}
end

iscontinuous(trait::DiscreteTrait) = false

struct LCtrait{D <: Number} <: AbstractTraits{D}
    vals::Vector{Vector{D}}
end

iscontinuous(trait::LCtrait) = false

function LCtrait(vals::Vector{Vector{<:AbstractFloat}})
    return LCtrait{typeof(1.0)}(vals)
end

"""
    DiscreteEvolve(numTraits::Int64, tree::BinaryTree)

Evolve a discrete switching trait along a BinaryTree, `tree`. Takes in a number
of traits, `numTraits` to be switched between and rate to switch between traits,
`switch_rate` with default value of 0.5.
"""
function DiscreteEvolve(numTraits::Int64, tree::BinaryTree, switch_rate = 0.5)
    # Create traits and assign to tips
    trts = DataFrame(trait1 = 1:numTraits)
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    return DiscreteTrait(Array(get_traits(tree, true)[:, 1]))
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
    trts = DataFrame(start = ustrip(val), σ² = ustrip(var))
    assign_traits!(tree, trts)
    # Get traits from tree
    newtrts = get_traits(tree, true)
    newtrts[!, :start] = newtrts[!, :start] .* unit(val)
    return Bin(Unclassified, Normal, newtrts[!, :start], newtrts[!, :σ²])
end

"""
    ContinuousTrait{C <: Number} <: AbstractTraits{T}

Abstract trait type that holds information on a single continuous trait for each
species, of any Number type `C`.
"""
abstract type ContinuousTrait{C <: Number} <: AbstractTraits{C} end


"""
    Bin{A <: NicheAxis, C, D} <: ContinuousTrait{C}

Trait type holding a **binned** habitat preference for each species on niche axis `A`, as one **built**
continuous response distribution per species — `dists::Vector{D}`, where `D` is a concrete
`Distributions.ContinuousUnivariateDistribution` (e.g. `Trapezoid{Float64}` or `Uniform{Float64}`). The
distributions are built once in the **support frame** `C` — the unit their support/domain is measured in,
which is also the trait's [`eltype`](@ref) and the unit the matching habitat must be in. In the hot loop the
[`DistRel`](@ref) relationship fetches `dists[sp]` — see [`getdist`](@ref) — and evaluates the density at
`ustrip(habitat)`, already a bare number in frame `C`, so there is no per-call conversion or allocation.
"""
struct Bin{A <: NicheAxis, C, D} <: ContinuousTrait{C}
    dists::Vector{D}
end

# Build a `Bin` whose distributions live in the `frame` unit (= its `eltype` / the required habitat unit),
# reading each row's bare parameters as being in `input_unit` and converting to `frame` *per role* (a
# location properly/affine, a scale as an interval, a shape dimensionless — see `read_distribution`; a
# shape-only family is placed via `offset`/`scale`). The single place the storage frame is fixed.
function _buildbin(::Type{A}, ::Type{D}, mat::AbstractMatrix, input_unit, frame;
                   offset, scale, probes) where {A <: NicheAxis, D}
    roles = param_roles_resolved(D; probes = probes)
    dists = [read_distribution(D, input_unit, mat[sp, :]; canonical = frame,
                               offset = offset, scale = scale, roles = roles)
             for sp in Base.axes(mat, 1)]
    return Bin{A, typeof(1.0 * frame), eltype(dists)}(dists)
end

"""
    Bin(::Type{A}, ::Type{D}, dist::Matrix; support = canonicalunit(A()),
        offset = nothing, scale = nothing, probes = …)

Build a [`Bin`](@ref) on axis `A` with response distribution `D` from a **bare** parameter matrix (one
species per row). `support` is the **frame** the distribution is built in — the unit its support/domain is
measured in, the trait's `eltype`, and the unit the matching habitat must be in (defaults to the axis's
canonical unit). The bare numbers are taken to be in that frame. `offset`/`scale` (references in `support`)
place a *shape-only* distribution (`Beta`, `LogNormal`, …) on the dimensioned axis via a `LocationScale`.
"""
function Bin(::Type{A}, ::Type{D}, dist::Matrix;
             support = canonicalunit(A()), offset = nothing, scale = nothing,
             probes = _default_probes(D)) where {A <: NicheAxis, D}
    dimension(support) == dimension(canonicalunit(A())) ||
        error("Bin support unit $support and axis $A's canonical unit " *
              "$(canonicalunit(A())) have different dimensions.")
    return _buildbin(A, D, dist, support, support; offset = offset,
                     scale = scale, probes = probes)
end

# Impute the *input* unit of a set of parameter vectors — the unit their (bare) magnitudes are read in:
# their single shared *dimensioned* unit (dimensionless/shape vectors ignored), or the `support` frame when
# all bare (so bare vectors are read in the frame, like the matrix constructor). Mixed units error.
function _impute_input_unit(params, support)
    us = unique(unit(eltype(v)) for v in params if unit(eltype(v)) !== NoUnits)
    isempty(us) && return support
    length(us) == 1 ||
        error("Bin parameter vectors carry differing units $(us); pass values in one unit.")
    return only(us)
end

"""
    Bin(::Type{A}, ::Type{D}, params::AbstractVector...; support = canonicalunit(A()),
        offset = nothing, scale = nothing, probes = …)

Build a [`Bin`](@ref) on axis `A` with response distribution `D` from **one (unitful) vector per
parameter** of `D` — the ergonomic, programmatic counterpart of the `dist::Matrix` constructor (and of the
old `GaussTrait(A, mean, sd)`). Each `params` vector holds one distribution parameter across species, e.g.
`Bin(MeanTemperature, Normal, opts, vars)` (a `μ` and a `σ` vector) or
`Bin(Precipitation, Gamma, shape, scale_vec)`. The vectors' own units are respected (read per role); the
distribution is built in the `support` **frame** (the unit its support/domain — and the matching habitat —
is measured in, default canonical), converting the inputs to it. So `support = K` and `support = u"°C"`
build the *same* preference in the K and °C frames respectively (`pdf(getdist(bin, sp), ustrip(support, x))`
is correct in each); the frame is the trait's `eltype`, and a Bin only matches a habitat in that same unit.
All parameter vectors must share a unit (a mixed set errors) and have one entry per species.
"""
function Bin(::Type{A}, ::Type{D}, params::AbstractVector...;
             support = canonicalunit(A()), offset = nothing, scale = nothing,
             probes = _default_probes(D)) where {A <: NicheAxis, D}
    isempty(params) && error("Bin needs at least one parameter vector for $D.")
    n = length(first(params))
    all(v -> length(v) == n, params) ||
        error("Bin parameter vectors must all have the same length (one entry per species).")
    dimension(support) == dimension(canonicalunit(A())) ||
        error("Bin support unit $support and axis $A's canonical unit " *
              "$(canonicalunit(A())) have different dimensions.")
    input_unit = _impute_input_unit(params, support)
    # A unitful vector → its magnitude in `input_unit` (its role, applied in `_buildbin`, decides position
    # vs interval); a bare vector passes straight through.
    cols = map(v -> unit(eltype(v)) === NoUnits ? float.(v) :
                    ustrip.(input_unit, v),
               params)
    return _buildbin(A, D, reduce(hcat, cols), input_unit, support;
                     offset = offset, scale = scale, probes = probes)
end

iscontinuous(::Bin) = true

"""    TempBin — a binned temperature preference (a [`Bin`](@ref) on `MeanTemperature`, K frame, with a `Trapezoid` response). """
const TempBin = Bin{MeanTemperature,
                    typeof(1.0 * canonicalunit(MeanTemperature())),
                    Trapezoid{Float64}}
"""    RainBin — a binned rainfall preference (a [`Bin`](@ref) on `Precipitation`, mm frame, with a `Uniform` response). """
const RainBin = Bin{Precipitation, typeof(1.0 * canonicalunit(Precipitation())),
                    Uniform{Float64}}

# Back-compat constructors preserving the old `TempBin(matrix)` / `RainBin(matrix)` call form.
TempBin(dist::Matrix) = Bin(MeanTemperature, Trapezoid, dist)
RainBin(dist::Matrix) = Bin(Precipitation, Uniform, dist)

# Convenience: the matching [`DistRel`](@ref) relationship for a `Bin`, taking its `TR` (unit) from the
# trait's axis — so callers building an ecosystem by hand need not re-type the unit.
DistRel(t::Bin) = DistRel{eltype(t)}()

"""
    TraitCollection2{T1, T2} <: AbstractTraits{Tuple{T1, T2}}

Trait collection that holds two trait types, `TR1` and `TR2`.
"""
struct TraitCollection2{T1, T2} <: AbstractTraits{Tuple{T1, T2}}
    one::T1
    two::T2
end

function iscontinuous(trait::TraitCollection2)
    return [iscontinuous(trait.one), iscontinuous(trait.two)]
end

eltype(trait::TraitCollection2) = [eltype(trait.one), eltype(trait.two)]

"""
    TraitCollection3{T1, T2, T3} <: AbstractTraits{Tuple{T1, T2, T3}}

Trait collection that holds three trait types, `TR1`, `TR2` and `TR3`.
"""
struct TraitCollection3{T1, T2, T3} <: AbstractTraits{Tuple{T1, T2, T3}}
    one::T1
    two::T2
    three::T3
end

function iscontinuous(trait::TraitCollection3)
    return [
        iscontinuous(trait.one),
        iscontinuous(trait.two),
        iscontinuous(trait.three)
    ]
end

function Base.eltype(trait::TraitCollection3)
    return [eltype(trait.one), eltype(trait.two), eltype(trait.three)]
end
