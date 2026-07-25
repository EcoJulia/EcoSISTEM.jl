# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoBase
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

# ---------------------------------------------------------------------------
# Layer dynamics
# ---------------------------------------------------------------------------
# The per-layer update rule. Named `LayerUpdate` (was `HabitatUpdate`) — it drives any layer, regime or
# supply, not just a habitat.

"""
    LayerUpdate{F <: Function, DT}

Stores the update rule for a layer. `changefun` is applied to the layer matrix at each
timestep, and `rate` is the rate of change with appropriate units.
"""
struct LayerUpdate{F <: Function, DT}
    changefun::F
    rate::DT
end

"""
    LayerUpdate(changefun::F, rate::DT, ::Type{D})

Construct a `LayerUpdate` with a type-checked rate. Errors if `rate * 1month` does not
have dimensions `D`.
"""
function LayerUpdate(changefun::F, rate::DT,
                     ::Type{D}) where {F <: Function, DT,
                                       D <: Unitful.Dimensions}
    dimension(rate * 1month) isa D ||
        error("Failed to match types $(rate * 1month) vs $D")
    return LayerUpdate(changefun, rate)
end

# ---------------------------------------------------------------------------
# AbstractLayer — the materialised, hot-loop grid-layer family
# ---------------------------------------------------------------------------
# A materialised layer is `matrix (+ time) + size + dynamics`, tagged with a `Role`
# (`Condition` condition vs `Resource` resource) and a `NicheAxis` (what it measures). It
# unifies today's `ContinuousRegime`/`ContinuousTimeRegime`/`DiscreteRegime` (and, from sub-step 3,
# the supply structs). During the fold the old `*Hab` names below are **aliases** over these
# types, so existing methods/constructors keep working (defaulting to the `Unclassified`
# axis) until each path is threaded with its real axis.

"""
    AbstractLayer{R <: Role}

Abstract supertype of the materialised, hot-loop grid layers, tagged with their
[`Role`](@ref) (`Condition` condition or `Resource` resource). Concrete kinds are
[`ContinuousLayer`](@ref), [`DiscreteLayer`](@ref) and the collections
[`LayerCollection2`](@ref)/`LayerCollection3`.
"""
abstract type AbstractLayer{R <: Role} <: EcoBase.AbstractGrid end

"""
    Unclassified <: NicheAxis

The default niche axis for a layer built without a declared one (an unclassified layer —
it still has data and a unit, but is not modelled on a named axis).
"""
struct Unclassified <: NicheAxis end

# An unclassified axis carries no modelled unit — treat it as dimensionless so axis-less continuous
# tolerances (e.g. a `NicheTolerance` built from bare data) resolve `canonicalunit`/`eltype` without erroring.
canonicalunit(::Unclassified) = NoUnits

# The third Resource-role canonical unit (Precipitation/SolarRadiation are in NicheInfo.jl): the
# free/simple per-individual demand rate. `day^-1`, not `1/day` — the latter is a Quantity, not a
# Units, breaking every other canonicalunit's uconvert-target contract; matches the literal
# already used in AbioticEnv.jl's `cancel` third branch.
canonicalunit(::Type{Resource}, ::Unclassified) = day^-1

"""
    ContinuousLayer{R <: Role, A <: NicheAxis, V <: Number, Arr <: AbstractArray{V}}

A continuous (numeric) grid layer of role `R` on niche axis `A`, holding a value `matrix`
(eltype `V`). The array type `Arr` carries the **rank**: a `Matrix{V}` is a static layer
(`time` unused, `= 1`), an `AbstractArray{V, 3}` a monthly time series indexed by `time`.
"""
mutable struct ContinuousLayer{R <: Role, A <: NicheAxis, V <: Number,
                               Arr <: AbstractArray{V}} <: AbstractLayer{R}
    matrix::Arr
    time::Int64
    size::Unitful.Length
    dynamics::LayerUpdate
end

"""
    DiscreteLayer{A <: NicheAxis, V, Arr <: AbstractArray{V}}

A discrete (categorical) grid layer on niche axis `A` — always a `Condition` (there is no
discrete supply).
"""
mutable struct DiscreteLayer{A <: NicheAxis, V, Arr <: AbstractArray{V}} <:
               AbstractLayer{Condition}
    matrix::Arr
    size::Unitful.Length
    dynamics::LayerUpdate
end

"""
    LayerCollection2{R <: Role, L1, L2}
    LayerCollection3{R <: Role, L1, L2, L3}

A collection of 2 or 3 layers of the same role `R` over one grid (e.g. temperature +
rainfall). Sub-layer types are kept concrete in named fields for hot-loop type stability.
"""
struct LayerCollection2{R <: Role, L1 <: AbstractLayer{R},
                        L2 <: AbstractLayer{R}} <: AbstractLayer{R}
    one::L1
    two::L2
end
struct LayerCollection3{R <: Role, L1 <: AbstractLayer{R},
                        L2 <: AbstractLayer{R},
                        L3 <: AbstractLayer{R}} <: AbstractLayer{R}
    one::L1
    two::L2
    three::L3
end

# ---------------------------------------------------------------------------
# Back-compat aliases + constructors (regime role)
# ---------------------------------------------------------------------------

"""
    AbstractRegime

An environmental condition matched to species tolerances — any `Condition`-role
[`AbstractLayer`](@ref).
"""
const AbstractRegime = AbstractLayer{Condition}

"""    ContinuousRegime{C} — a static continuous regime (a `Condition`-role [`ContinuousLayer`](@ref) over a `Matrix{C}`). """
const ContinuousRegime{C} = ContinuousLayer{Condition, A, C,
                                            Matrix{C}} where {A}
"""    ContinuousTimeRegime{C, M} — a monthly time-varying continuous regime (a `Condition`-role [`ContinuousLayer`](@ref) over a 3-D array). """
const ContinuousTimeRegime{C,
                           M <: AbstractArray{C, 3}} = ContinuousLayer{Condition,
                                                                       A, C,
                                                                       M} where {A}
"""    DiscreteRegime{D} — a categorical regime (a [`DiscreteLayer`](@ref), e.g. land cover). """
const DiscreteRegime{D} = DiscreteLayer{A, D, Matrix{D}} where {A}
"""    RegimeCollection2{H1, H2} — two regimes over one grid (e.g. temperature + rainfall). """
const RegimeCollection2{H1, H2} = LayerCollection2{Condition, H1, H2}
"""    RegimeCollection3{H1, H2, H3} — three regimes over one grid. """
const RegimeCollection3{H1, H2, H3} = LayerCollection3{Condition, H1, H2, H3}

# Old positional constructors → new layer types, defaulting to the `Unclassified` axis and
# time = 1 (the axis-aware `materialise` path constructs `ContinuousLayer{Condition, A}` with
# the real axis directly).
function ContinuousRegime(matrix::Matrix{C}, size::Unitful.Length,
                          dynamics::LayerUpdate) where {C}
    return ContinuousLayer{Condition, Unclassified, C, Matrix{C}}(matrix, 1,
                                                                  size,
                                                                  dynamics)
end
function ContinuousTimeRegime(matrix::AbstractArray{C, 3}, time::Integer,
                              size::Unitful.Length,
                              dynamics::LayerUpdate) where {C}
    return ContinuousLayer{Condition, Unclassified, C, typeof(matrix)}(matrix,
                                                                       time,
                                                                       size,
                                                                       dynamics)
end
function DiscreteRegime(matrix::Matrix{D}, size::Unitful.Length,
                        dynamics::LayerUpdate) where {D}
    return DiscreteLayer{Unclassified, D, Matrix{D}}(matrix, size, dynamics)
end

# The collection aliases fix `R = Condition` but leave the sub-layer params free, so the
# auto-generated constructor doesn't cover `RegimeCollection2(l1, l2)`; forward to the bare
# `LayerCollection` constructor (which infers the role from the sub-layers).
function RegimeCollection2(l1::AbstractLayer{Condition},
                           l2::AbstractLayer{Condition})
    return LayerCollection2(l1, l2)
end
function RegimeCollection3(l1::AbstractLayer{Condition},
                           l2::AbstractLayer{Condition},
                           l3::AbstractLayer{Condition})
    return LayerCollection3(l1, l2, l3)
end

# ---------------------------------------------------------------------------
# Back-compat aliases (supply role)
# ---------------------------------------------------------------------------
# A supply is a `Resource`-role layer: the old supply structs are aliases over
# `ContinuousLayer{Resource, axis, V, Arr}` — static supplies are 2-D (`Matrix`), the
# time-varying `*TimeSupply`s 3-D (`Array`, indexed by `time`). The axis records the
# resource measured; the (unused) `size` and the `dynamics` rule are filled by the
# constructors in Energy.jl. `SimpleSupply` (free energy) has no resource axis.

"""
    AbstractSupply

A resource layer consumed by species — any `Resource`-role [`AbstractLayer`](@ref).
"""
const AbstractSupply = AbstractLayer{Resource}

# A time-varying supply: a `Resource`-role continuous layer whose array is 3-D (monthly).
const AbstractTimeSupply = ContinuousLayer{Resource, A, V,
                                           Arr} where {A, V,
                                                       Arr <:
                                                       AbstractArray{V, 3}}

# Single source of truth for each Resource rate's element type, derived from `canonicalunit`
# (NicheInfo.jl / above) rather than re-spelling the literal unit at every alias/constructor site.
const _SimpleRate = typeof(1.0 * canonicalunit(Resource, Unclassified()))
const _SolarRate = typeof(1.0 * canonicalunit(Resource, SolarRadiation()))
const _WaterRate = typeof(1.0 * canonicalunit(Resource, Precipitation()))

"""    SimpleSupply — a free-resource supply rate (`/day`), one float per cell. """
const SimpleSupply = ContinuousLayer{Resource, Unclassified, _SimpleRate,
                                     Matrix{_SimpleRate}}
"""    SolarSupply — a static solar-energy supply rate (`kJ/day`). """
const SolarSupply = ContinuousLayer{Resource, SolarRadiation, _SolarRate,
                                    Matrix{_SolarRate}}
"""    WaterSupply — a static available-water supply rate (`L/day`). """
const WaterSupply = ContinuousLayer{Resource, Precipitation, _WaterRate,
                                    Matrix{_WaterRate}}
"""    SolarTimeSupply — a monthly time-varying solar-energy supply rate (`kJ/day`). """
const SolarTimeSupply = ContinuousLayer{Resource, SolarRadiation, _SolarRate,
                                        Array{_SolarRate, 3}}
"""    WaterTimeSupply — a monthly time-varying available-water supply rate (`L/day`). """
const WaterTimeSupply = ContinuousLayer{Resource, Precipitation, _WaterRate,
                                        Array{_WaterRate, 3}}
"""    SupplyCollection2{B1, B2} — two supplies over one grid (e.g. solar + water). """
const SupplyCollection2{B1, B2} = LayerCollection2{Resource, B1, B2}

# As for `RegimeCollection2`, the partially-applied alias has no auto-generated
# constructor; forward to the bare `LayerCollection2` (which infers the role).
function SupplyCollection2(b1::AbstractLayer{Resource},
                           b2::AbstractLayer{Resource})
    return LayerCollection2(b1, b2)
end

# ---------------------------------------------------------------------------
# Axis-driven canonicalisation + axis re-tagging (shared by the hand `*AE`
# constructors in `AbioticEnv.jl` and the `build_*`/`materialise` path in `Simplify.jl`)
# ---------------------------------------------------------------------------

# Canonicalise a regime *value* (a position) to its layer axis's unit, `canonicalunit(A())`: a unitful
# value converts (proper affine — canonical units are absolute, so no interval subtlety); a bare value
# attaches the unit. An axis with no canonical unit (`Unclassified`/dimensionless, `NoUnits`/`nothing`)
# keeps the value's magnitude but still **absolutises** an affine unit (°C→K) — regimes are always in an
# absolute unit, which the downstream dynamics/rate machinery assumes. The single, axis-driven replacement
# for the old dimension-sniffing conversions.
function _canonical(x, ::Type{A}) where {A <: NicheAxis}
    return _tocanon(x, canonicalunit(A()))
end
function _tocanon(x, u)
    return (u === NoUnits || isnothing(u)) ? _absolutise(float(x)) :
           _tocanon_u(float(x), u)
end
_tocanon_u(x::Unitful.AbstractQuantity, u) = uconvert(u, x)
_tocanon_u(x::Real, u) = x * u
function _absolutise(x::Unitful.AbstractQuantity)
    return uconvert(Unitful.absoluteunit(unit(x)), x)
end
_absolutise(x::Real) = x

# Re-tag a materialised layer with its niche axis `A` — a phantom type parameter, so this shares the
# arrays. The low-level constructors build an `Unclassified` layer that this narrows to the real axis.
function _reaxis(l::ContinuousLayer{Condition, A0, V, Arr},
                 ::Type{A}) where {A0, V, Arr, A <: NicheAxis}
    return ContinuousLayer{Condition, A, V, Arr}(l.matrix, l.time, l.size,
                                                 l.dynamics)
end
function _reaxis(l::DiscreteLayer{A0, V, Arr},
                 ::Type{A}) where {A0, V, Arr, A <: NicheAxis}
    return DiscreteLayer{A, V, Arr}(l.matrix, l.size, l.dynamics)
end

# Symmetric value accessors for the environment's two layers: a `Regime` supplies a `condition` grid, a
# `Supply` supplies a `resource` grid (mirroring the species side's `.tolerance` / `.demand`).
condition(regime::AbstractRegime) = regime.matrix
resource(supply::AbstractSupply) = supply.matrix
