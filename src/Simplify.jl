# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Distributions
using Random
using AxisArrays
using RasterDataSources

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# The niche axis of a threaded trait / layer (`Unclassified` when none was declared).
axisof(::NicheTolerance{A}) where {A} = A
axisof(::ContinuousLayer{R, A}) where {R, A} = A
axisof(::DiscreteLayer{A}) where {A} = A
axisof(_) = Unclassified

# A trait may only be matched to a layer on the same niche axis (`Unclassified` matches
# anything, for the un-threaded/back-compat paths).
function _checkaxis(t, h)
    at, ah = axisof(t), axisof(h)
    (at === Unclassified || ah === Unclassified || at === ah) ||
        error("A trait on the `$at` axis cannot be matched to a layer on the `$ah` axis — check the order (or axes) of the species niches against the environment layers.")
    return nothing
end

# Infer the trait–environment nichefit from the trait type, checking the niche axes.
# A `NicheTolerance`'s response distribution is built per-species in `_suitability`, so the nichefit
# (`NicheSuitability`) is distribution-agnostic — it only needs the regime's unit for its `TR`.
function _default_suitability(t::NicheTolerance, regime)
    _checkaxis(t, regime)
    return NicheSuitability{eltype(regime)}()
end
function _default_suitability(::DiscreteTolerance, regime)
    return MatchSuitability{eltype(regime)}()
end
function _default_suitability(::LandCoverTolerance, regime)
    return LandCoverSuitability{eltype(regime)}()
end
# A collection of tolerances over a matching collection of regimes infers each sub-tolerance's
# nichefit and combines them multiplicatively.
function _default_suitability(t::ToleranceCollection2, h::RegimeCollection2)
    return multiplicativeFit2(_default_suitability(t.one, h.one),
                              _default_suitability(t.two, h.two))
end
function _default_suitability(t::ToleranceCollection3, h::RegimeCollection3)
    return multiplicativeFit3(_default_suitability(t.one, h.one),
                              _default_suitability(t.two, h.two),
                              _default_suitability(t.three, h.three))
end
function _default_suitability(traits::AbstractTolerance, regime)
    return error("Cannot infer a trait nichefit for traits of type $(typeof(traits)); pass one explicitly with `nichefit = …`.")
end
