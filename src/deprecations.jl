# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations — main `EcoSISTEM` module
#
# Every deprecated public trait/relationship API is collected here (sorted into
# sections by context) and included late in `EcoSISTEM.jl`, after all the types
# it shims. Each shim warns — so downstream code gets a migration message rather
# than a silent `MethodError` — and forwards to the current API. Mirrored by
# `test/test_deprecations.jl`. The `ClimatePref` submodule keeps its own
# deprecations in `src/ClimatePref/deprecations.jl` (a single file cannot span
# two modules).
# ===========================================================================

# ---------------------------------------------------------------------------
# Trait line: `GaussTrait` → `Bin(A, Normal, …)`
#
# A Gaussian preference is just the `Normal` case of a `Bin`. The redesign takes
# a trait's unit from its niche axis (not its data), so a *unitful* preference
# must name an axis; only dimensionless (bare) data may stay axis-less.
# ---------------------------------------------------------------------------
"""
    GaussTrait(mean, sd)
    GaussTrait(::Type{A}, mean, sd)

!!! warning "Deprecated"
    `GaussTrait` is deprecated and will be removed. A Gaussian preference is just the `Normal` case of a
    [`Bin`](@ref), so use one directly — `Bin(A, Normal, mean, sd)` (one unitful vector per parameter,
    `support` imputed). This shim forwards to that. `GaussTrait(::Type{A}, mean, sd)` names the niche axis
    `A` explicitly; the axis-less `GaussTrait(mean, sd)` form takes dimensionless (bare) data as
    `Unclassified`, and — for back-compatibility only — *infers* the axis from a unitful vector's unit
    (temperature → `MeanTemperature`, length → `Precipitation`) with an extra warning, because a
    [`Bin`](@ref)'s unit must come from its axis, not its data.
"""
function GaussTrait(::Type{A}, mean::AbstractVector,
                    sd::AbstractVector) where {A <: NicheAxis}
    Base.depwarn("`GaussTrait` is deprecated; use `Bin($A, Normal, mean, sd)` instead.",
                 :GaussTrait)
    return Bin(A, Normal, mean, sd)
end

# Axis-less bare (dimensionless) form: an `Unclassified` Bin, as before.
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
                 "`GaussTrait(MeanTemperature, mean, sd)` or `Bin(axis, Normal, mean, sd)`.")
end

# Axis-less *unitful* form (doubly deprecated): kept working via unit→axis inference so downstream code gets
# a clear migration message instead of a `MethodError`, but the meaning-from-unit inference will be dropped.
function GaussTrait(mean::AbstractVector{<:Unitful.AbstractQuantity},
                    sd::AbstractVector)
    u = unit(eltype(mean))
    axis = _axis_from_unit(u)
    Base.depwarn("`GaussTrait(mean, sd)` with unitful data infers the niche axis from its unit " *
                 "($u → $axis); this meaning-from-unit fallback is deprecated. Name the axis: " *
                 "`GaussTrait($axis, mean, sd)` or `Bin($axis, Normal, mean, sd)`.",
                 :GaussTrait)
    return Bin(axis, Normal, mean, sd)
end

# ---------------------------------------------------------------------------
# Trait line: `Gauss` / `Trapeze` / `Unif` relationships → `DistRel`
#
# A Gaussian/trapezoidal/uniform preference is just the `Normal`/`Trapezoid`/`Uniform` case of a `Bin`, so
# all three relationships are `DistRel` now. Each is kept as a distinct type that warns on construction but
# shares `DistRel`'s 2-arg density functor, so old hand-built ecosystems still build and evaluate.
# ---------------------------------------------------------------------------
"""
    Gauss{TR} <: AbstractTraitRelationship{TR}

!!! warning "Deprecated"
    `Gauss` is deprecated and will be removed; use [`DistRel`](@ref) instead (a Gaussian preference is the
    `Normal` case of a [`Bin`](@ref)). This shim shares `DistRel`'s 2-argument density functor and also
    retains the legacy 3-argument `(current, opt, sd)` Gaussian call for back-compatibility.
"""
mutable struct Gauss{TR} <: AbstractTraitRelationship{TR}
    function Gauss{TR}() where {TR}
        Base.depwarn("`Gauss` is deprecated; use `DistRel{TR}()` instead.",
                     :Gauss)
        return new{TR}()
    end
end

"""
    Trapeze{TR} <: AbstractTraitRelationship{TR}

!!! warning "Deprecated"
    `Trapeze` is deprecated and will be removed; use [`DistRel`](@ref) instead (it pairs with a `Trapezoid`
    [`Bin`](@ref)).
"""
mutable struct Trapeze{TR} <: AbstractTraitRelationship{TR}
    function Trapeze{TR}() where {TR}
        Base.depwarn("`Trapeze` is deprecated; use `DistRel{TR}()` instead.",
                     :Trapeze)
        return new{TR}()
    end
end

"""
    Unif{TR} <: AbstractTraitRelationship{TR}

!!! warning "Deprecated"
    `Unif` is deprecated and will be removed; use [`DistRel`](@ref) instead (it pairs with a `Uniform`
    [`Bin`](@ref)).
"""
mutable struct Unif{TR} <: AbstractTraitRelationship{TR}
    function Unif{TR}() where {TR}
        Base.depwarn("`Unif` is deprecated; use `DistRel{TR}()` instead.",
                     :Unif)
        return new{TR}()
    end
end

# The deprecated shims share `DistRel`'s density functor + trait interface.
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
# `DistRel` shim. New code should build a `Normal` `Bin` and use `DistRel`'s 2-argument functor instead.
function (::Gauss{TR})(current::TR, opt::TR, var::TR) where {TR}
    pref = 1.0 / sqrt(2 * π * var^2) *
           exp(-abs(current - opt)^2 / (2 * var^2))
    return pref * unit(current)
end
