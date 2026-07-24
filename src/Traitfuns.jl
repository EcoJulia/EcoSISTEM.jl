# SPDX-License-Identifier: LGPL-3.0-or-later

"""
    suitability(eco::AbstractEcosystem, pos::Int64, sp::Int64)

Calculate nichefit between the current environment and a species' particular
trait.

"""
function suitability(eco::AbstractEcosystem, pos::Int64, sp::Int64)
    regime = eco.habitat.regime
    tolerance = eco.spplist.tolerance
    nichefit = eco.nichefit
    return _suitability(regime, tolerance, nichefit, pos, sp)
end

function _suitability(regime::RegimeCollection2,
                      tolerance::ToleranceCollection2,
                      nichefit::R,
                      pos::Int64,
                      sp::Int64) where {R <: AbstractNicheFit}
    res1 = _suitability(regime.one, tolerance.one, nichefit.one, pos, sp)
    res2 = _suitability(regime.two, tolerance.two, nichefit.two, pos, sp)
    return combinefit(nichefit)(res1, res2)
end
function _suitability(regime::ContinuousRegime,
                      tolerance::NicheTolerance,
                      nichefit::R,
                      pos::Int64,
                      sp::Int64) where {R <: AbstractNicheFit}
    h = getregime(regime, pos)
    # Static-regime counterpart of the `ContinuousTimeRegime` method below: fetch the species'
    # pre-built response distribution and evaluate it via `nichefit` (no per-call construction/allocation).
    return nichefit(getdist(tolerance, sp), h)
end
function _suitability(regime::ContinuousTimeRegime,
                      tolerance::NicheTolerance,
                      nichefit::R,
                      pos::Int64,
                      sp::Int64) where {R <: AbstractNicheFit}
    h = getregime(regime, pos)
    # Fetch the species' pre-built response distribution and evaluate it via `nichefit` — no per-call
    # construction or allocation (the distributions were built once when the `NicheTolerance` was made).
    return nichefit(getdist(tolerance, sp), h)
end
function _suitability(regime::DiscreteRegime,
                      tolerance::DiscreteTolerance,
                      nichefit::R,
                      pos::Int64,
                      sp::Int64) where {R <: AbstractNicheFit}
    currentniche = getregime(regime, pos)
    preference = getpref(tolerance, sp)
    return nichefit(currentniche, preference)
end
function _suitability(regime::RegimeCollection3,
                      tolerance::ToleranceCollection3,
                      nichefit::R,
                      pos::Int64,
                      sp::Int64) where {R <: AbstractNicheFit}
    res1 = _suitability(regime.one, tolerance.one, nichefit.one, pos, sp)
    res2 = _suitability(regime.two, tolerance.two, nichefit.two, pos, sp)
    res3 = _suitability(regime.three, tolerance.three, nichefit.three, pos, sp)
    return combinefit(nichefit)(res1, res2, res3)
end

function _suitability(regime::DiscreteRegime,
                      tolerance::LandCoverTolerance,
                      nichefit::R,
                      pos::Int64,
                      sp::Int64) where {R <: AbstractNicheFit}
    h = getregime(regime, pos)
    vals = getpref(tolerance, sp)
    return nichefit(h, vals)
end

"""
    getpref(tolerance::LandCoverTolerance, sp::Int64)

Extract the land cover preference values for species `sp` from an `LandCoverTolerance`.
"""
function getpref(tolerance::LandCoverTolerance, sp::Int64)
    return tolerance.vals[sp]
end

"""
    getdist(tolerance::NicheTolerance, sp::Int64)

Return the pre-built response distribution for species `sp` from a [`NicheTolerance`](@ref) trait. This is what
`_suitability` uses in the hot loop — a plain vector fetch, no per-call construction or allocation.
"""
function getdist(tolerance::NicheTolerance, sp::Int64)
    return tolerance.dists[sp]
end

# Per-species niche optima (distribution means) of a `NicheTolerance`, as a bare vector in the axis's canonical
# frame — the accessor the old `GaussTrait.mean` field readers migrate to (e.g. size tolerance evolved by
# `ContinuousEvolve`).
_nichemeans(tolerance::NicheTolerance) = Distributions.mean.(tolerance.dists)

"""
    getpref(tolerance::NicheTolerance, sp::Int64)

Return the parameters of species `sp`'s response distribution from a [`NicheTolerance`](@ref) trait (a
back-compat shim over [`getdist`](@ref); the parameters are no longer stored separately).
"""
function getpref(tolerance::NicheTolerance, sp::Int64)
    return Distributions.params(getdist(tolerance, sp))
end

"""
    getpref(tolerance::DiscreteTolerance, sp::Int64)

Extract the discrete niche preference value for species `sp` from a
[`DiscreteTolerance`](@ref).
"""
function getpref(tolerance::DiscreteTolerance, sp::Int64)
    return tolerance.val[sp]
end

"""
    getpref(tolerance::T, field::Symbol) where T <: AbstractTolerance

Extract trait preferences for all species in the ecosystem.

"""
function getpref(tolerance::T, field::Symbol) where {T <: AbstractTolerance}
    return getfield(tolerance, field)
end

"""
    getrelationship(nichefit::R, field::Symbol) where R <: AbstractNicheFit

Extract the trait nichefit of all species in the ecosystem.

"""
function getrelationship(nichefit::R,
                         field::Symbol) where {R <: AbstractNicheFit}
    return getfield(nichefit, field)
end
