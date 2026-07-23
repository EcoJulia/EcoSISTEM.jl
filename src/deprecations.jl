# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations — main `EcoSISTEM` module
#
# Every deprecated public trait/nichefit API is collected here (sorted into
# sections by context) and included late in `EcoSISTEM.jl`, after all the types
# it shims. Each shim warns — so downstream code gets a migration message rather
# than a silent `MethodError` — and forwards to the current API. Mirrored by
# `test/test_deprecations.jl`. The `ClimatePref` submodule keeps its own
# deprecations in `src/ClimatePref/deprecations.jl` (a single file cannot span
# two modules).
# ===========================================================================

# ---------------------------------------------------------------------------
# Trait line: `GaussTrait` → `NicheTolerance(A, Normal, …)`
#
# A Gaussian preference is just the `Normal` case of a `NicheTolerance`. The redesign takes
# a trait's unit from its niche axis (not its data), so a *unitful* preference
# must name an axis; only dimensionless (bare) data may stay axis-less.
# ---------------------------------------------------------------------------
"""
    GaussTrait(mean, sd)
    GaussTrait(::Type{A}, mean, sd)

!!! warning "Deprecated"
    `GaussTrait` is deprecated and will be removed. A Gaussian preference is just the `Normal` case of a
    [`NicheTolerance`](@ref), so use one directly — `NicheTolerance(A, Normal, mean, sd)` (one unitful vector per parameter,
    `support` imputed). This shim forwards to that. `GaussTrait(::Type{A}, mean, sd)` names the niche axis
    `A` explicitly; the axis-less `GaussTrait(mean, sd)` form takes dimensionless (bare) data as
    `Unclassified`, and — for back-compatibility only — *infers* the axis from a unitful vector's unit
    (temperature → `MeanTemperature`, length → `Precipitation`) with an extra warning, because a
    [`NicheTolerance`](@ref)'s unit must come from its axis, not its data.
"""
function GaussTrait(::Type{A}, mean::AbstractVector,
                    sd::AbstractVector) where {A <: NicheAxis}
    Base.depwarn("`GaussTrait` is deprecated; use `NicheTolerance($A, Normal, mean, sd)` instead.",
                 :GaussTrait)
    return NicheTolerance(A, Normal, mean, sd)
end

# Axis-less bare (dimensionless) form: an `Unclassified` NicheTolerance, as before.
function GaussTrait(mean::AbstractVector{<:Real}, sd::AbstractVector)
    return GaussTrait(Unclassified, mean, sd)
end

# Map a bare unit back to the niche axis it canonically measures. This is a *deprecated* fallback for the
# axis-less unitful `GaussTrait(mean, sd)` form only: inferring a layer's meaning from its unit is exactly
# the anti-pattern the axis redesign removes, so it is intentionally minimal (temperature/length) and errors
# for anything else, telling the caller to name the axis.
function _axis_from_unit(u)
    dimension(u) == dimension(u"K") && return MeanTemperature
    dimension(u) == dimension(u"mm") && return Precipitation
    return error("Cannot infer a niche axis from unit `$u`; name the axis explicitly, e.g. " *
                 "`GaussTrait(MeanTemperature, mean, sd)` or `NicheTolerance(axis, Normal, mean, sd)`.")
end

# Axis-less *unitful* form (doubly deprecated): kept working via unit→axis inference so downstream code gets
# a clear migration message instead of a `MethodError`, but the meaning-from-unit inference will be dropped.
function GaussTrait(mean::AbstractVector{<:Unitful.AbstractQuantity},
                    sd::AbstractVector)
    u = unit(eltype(mean))
    axis = _axis_from_unit(u)
    Base.depwarn("`GaussTrait(mean, sd)` with unitful data infers the niche axis from its unit " *
                 "($u → $axis); this meaning-from-unit fallback is deprecated. Name the axis: " *
                 "`GaussTrait($axis, mean, sd)` or `NicheTolerance($axis, Normal, mean, sd)`.",
                 :GaussTrait)
    return NicheTolerance(axis, Normal, mean, sd)
end

# ---------------------------------------------------------------------------
# Trait line: `Gauss` / `Trapeze` / `Unif` relationships → `NicheSuitability`
#
# A Gaussian/trapezoidal/uniform preference is just the `Normal`/`Trapezoid`/`Uniform` case of a `NicheTolerance`, so
# all three relationships are `NicheSuitability` now. Each is kept as a distinct type that warns on construction but
# shares `NicheSuitability`'s 2-arg density functor, so old hand-built ecosystems still build and evaluate.
# ---------------------------------------------------------------------------
"""
    Gauss{TR} <: AbstractNicheFit{TR}

!!! warning "Deprecated"
    `Gauss` is deprecated and will be removed; use [`NicheSuitability`](@ref) instead (a Gaussian preference is the
    `Normal` case of a [`NicheTolerance`](@ref)). This shim shares `NicheSuitability`'s 2-argument density functor and also
    retains the legacy 3-argument `(current, opt, sd)` Gaussian call for back-compatibility.
"""
mutable struct Gauss{TR} <: AbstractNicheFit{TR}
    function Gauss{TR}() where {TR}
        Base.depwarn("`Gauss` is deprecated; use `NicheSuitability{TR}()` instead.",
                     :Gauss)
        return new{TR}()
    end
end

"""
    Trapeze{TR} <: AbstractNicheFit{TR}

!!! warning "Deprecated"
    `Trapeze` is deprecated and will be removed; use [`NicheSuitability`](@ref) instead (it pairs with a `Trapezoid`
    [`NicheTolerance`](@ref)).
"""
mutable struct Trapeze{TR} <: AbstractNicheFit{TR}
    function Trapeze{TR}() where {TR}
        Base.depwarn("`Trapeze` is deprecated; use `NicheSuitability{TR}()` instead.",
                     :Trapeze)
        return new{TR}()
    end
end

"""
    Unif{TR} <: AbstractNicheFit{TR}

!!! warning "Deprecated"
    `Unif` is deprecated and will be removed; use [`NicheSuitability`](@ref) instead (it pairs with a `Uniform`
    [`NicheTolerance`](@ref)).
"""
mutable struct Unif{TR} <: AbstractNicheFit{TR}
    function Unif{TR}() where {TR}
        Base.depwarn("`Unif` is deprecated; use `NicheSuitability{TR}()` instead.",
                     :Unif)
        return new{TR}()
    end
end

# The deprecated shims share `NicheSuitability`'s density functor + trait interface.
for Old in (:Gauss, :Trapeze, :Unif)
    @eval begin
        function (::$Old{TR})(dist::ContinuousUnivariateDistribution,
                              current::TR) where {TR}
            return pdf(dist, ustrip(current))
        end
        iscontinuous(::$Old) = true
        Base.eltype(::$Old{TR}) where {TR} = TR
    end
end

# Legacy 3-argument Gaussian functor `Gauss{TR}()(current, opt, sd)`, restored for back-compatibility: the
# hand-written Gaussian density (a dimensionless `Quantity`) that `Gauss` evaluated before it became a
# `NicheSuitability` shim. New code should build a `Normal` `NicheTolerance` and use `NicheSuitability`'s 2-argument functor instead.
function (::Gauss{TR})(current::TR, opt::TR, var::TR) where {TR}
    pref = 1.0 / sqrt(2 * π * var^2) *
           exp(-abs(current - opt)^2 / (2 * var^2))
    return pref * unit(current)
end

# ---------------------------------------------------------------------------
# Resource line: `Resource` → `Supply` (v0.4.0 rename; the environment's resource layer)
# ---------------------------------------------------------------------------
# `AbstractSupply`/`AbstractTimeSupply` are unexported, so their shims don't export either (the `false`).
Base.@deprecate_binding AbstractBudget AbstractSupply false
Base.@deprecate_binding AbstractTimeBudget AbstractTimeSupply false
Base.@deprecate_binding SimpleBudget SimpleSupply
Base.@deprecate_binding SolarBudget SolarSupply
Base.@deprecate_binding SolarTimeBudget SolarTimeSupply
Base.@deprecate_binding WaterBudget WaterSupply
Base.@deprecate_binding WaterTimeBudget WaterTimeSupply
Base.@deprecate_binding VolWaterBudget VolWaterSupply
Base.@deprecate_binding VolWaterTimeBudget VolWaterTimeSupply
Base.@deprecate_binding BudgetCollection2 SupplyCollection2
@deprecate getbudget(eco) getsupply(eco)

# ---------------------------------------------------------------------------
# Resource line: `Requirement` → `Demand` (v0.4.0 rename; the species' resource need)
# ---------------------------------------------------------------------------
Base.@deprecate_binding AbstractRequirement AbstractDemand false
Base.@deprecate_binding Abstract1Requirement Abstract1Demand false
Base.@deprecate_binding Abstract2Requirements Abstract2Demands false
Base.@deprecate_binding SimpleRequirement SimpleDemand
Base.@deprecate_binding SizeRequirement SizeDemand
Base.@deprecate_binding SolarRequirement SolarDemand
Base.@deprecate_binding WaterRequirement WaterDemand
Base.@deprecate_binding VolWaterRequirement VolWaterDemand
Base.@deprecate_binding ReqCollection2 DemandCollection2

# ---------------------------------------------------------------------------
# Condition line: `Condition`(-role layer) → `Regime` (v0.4.0 rename; the environment's condition layer)
# ---------------------------------------------------------------------------
Base.@deprecate_binding ContinuousHab ContinuousRegime
Base.@deprecate_binding ContinuousTimeHab ContinuousTimeRegime
Base.@deprecate_binding DiscreteHab DiscreteRegime
Base.@deprecate_binding HabitatCollection2 RegimeCollection2
Base.@deprecate_binding HabitatCollection3 RegimeCollection3
@deprecate gethabitat getregime

# ---------------------------------------------------------------------------
# Environment container: `AbioticEnv`/`GridAbioticEnv` → `Condition` (v0.4.0 rename). NB the *condition
# layer* `AbstractHabitat` was renamed to `AbstractRegime`, freeing `AbstractHabitat` for the environment;
# v0.4.0's (unexported) `AbstractHabitat` therefore changes meaning — a NEWS breaking note, not a shim.
# ---------------------------------------------------------------------------
Base.@deprecate_binding AbstractAbiotic AbstractHabitat false
Base.@deprecate_binding GridAbioticEnv GridHabitat

# `reenergise!` → `resupply!` (v0.4.0 rename; "energy" is a misnomer — the resource isn't always energy)
@deprecate reenergise! resupply!

# ---------------------------------------------------------------------------
# Condition line: `Trait`/`Bin` → `Tolerance`/`NicheTolerance` (v0.4.0 rename; the species' condition response)
# ---------------------------------------------------------------------------
Base.@deprecate_binding AbstractTraits AbstractTolerance false
Base.@deprecate_binding ContinuousTrait ContinuousTolerance false
Base.@deprecate_binding DiscreteTrait DiscreteTolerance
Base.@deprecate_binding LCtrait LCtolerance
Base.@deprecate_binding TraitCollection2 ToleranceCollection2
Base.@deprecate_binding TraitCollection3 ToleranceCollection3
Base.@deprecate_binding TempBin TempTolerance
Base.@deprecate_binding RainBin RainTolerance
@deprecate traitpopulate!(args...) tolerancepopulate!(args...)
@deprecate traitrepopulate!(args...) tolerancerepopulate!(args...)

# ---------------------------------------------------------------------------
# Condition line: the matcher `TraitRelationship`/`Match`/… → `NicheFit`/`Suitability` (v0.4.0 rename)
# ---------------------------------------------------------------------------
Base.@deprecate_binding AbstractTraitRelationship AbstractNicheFit false
Base.@deprecate_binding Match MatchSuitability
Base.@deprecate_binding LCmatch LCsuitability
Base.@deprecate_binding NoRelContinuous NoFitContinuous
Base.@deprecate_binding NoRelDiscrete NoFitDiscrete
Base.@deprecate_binding multiplicativeTR2 multiplicativeFit2
Base.@deprecate_binding multiplicativeTR3 multiplicativeFit3
Base.@deprecate_binding additiveTR2 additiveFit2
Base.@deprecate_binding additiveTR3 additiveFit3
@deprecate gettraitrel getnichefit
