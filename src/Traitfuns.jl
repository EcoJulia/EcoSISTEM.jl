# SPDX-License-Identifier: LGPL-3.0-or-later

"""
    traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64)

Calculate relationship between the current environment and a species' particular
trait.

"""
function traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64)
    regime = eco.habitat.regime
    tolerance = eco.spplist.tolerance
    rel = eco.relationship
    return _traitfun(regime, tolerance, rel, pos, sp)
end

function _traitfun(regime::RegimeCollection2,
                   tolerance::ToleranceCollection2,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    res1 = _traitfun(regime.one, tolerance.one, rel.one, pos, sp)
    res2 = _traitfun(regime.two, tolerance.two, rel.two, pos, sp)
    return combineTR(rel)(res1, res2)
end
function _traitfun(regime::ContinuousRegime,
                   tolerance::NicheTolerance,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = getregime(regime, pos)
    # Static-regime counterpart of the `ContinuousTimeRegime` method below: fetch the species'
    # pre-built response distribution and evaluate it via `rel` (no per-call construction/allocation).
    return rel(getdist(tolerance, sp), h)
end
function _traitfun(regime::ContinuousTimeRegime,
                   tolerance::NicheTolerance,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = getregime(regime, pos)
    # Fetch the species' pre-built response distribution and evaluate it via `rel` — no per-call
    # construction or allocation (the distributions were built once when the `NicheTolerance` was made).
    return rel(getdist(tolerance, sp), h)
end
function _traitfun(regime::DiscreteRegime,
                   tolerance::DiscreteTolerance,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    currentniche = getregime(regime, pos)
    preference = getpref(tolerance, sp)
    return rel(currentniche, preference)
end
function _traitfun(regime::RegimeCollection3,
                   tolerance::ToleranceCollection3,
                   rel::R,
                   pos::Int64,
                   spp::Int64) where {R <: AbstractTraitRelationship}
    res1 = _traitfun(regime.one, tolerance.one, rel.one, pos, spp)
    res2 = _traitfun(regime.two, tolerance.two, rel.two, pos, spp)
    res3 = _traitfun(regime.three, tolerance.three, rel.three, pos, spp)
    return combineTR(rel)(res1, res2, res3)
end

function _traitfun(regime::DiscreteRegime,
                   tolerance::LCtolerance,
                   rel::R,
                   pos::Int64,
                   spp::Int64) where {R <: AbstractTraitRelationship}
    h = getregime(regime, pos)
    vals = getpref(tolerance, spp)
    return rel(h, vals)
end

"""
    getpref(tolerance::LCtolerance, spp::Int64)

Extract the land cover preference values for species `spp` from an `LCtolerance`.
"""
function getpref(tolerance::LCtolerance, spp::Int64)
    return tolerance.vals[spp]
end

"""
    getdist(tolerance::NicheTolerance, sp::Int64)

Return the pre-built response distribution for species `sp` from a [`NicheTolerance`](@ref) trait. This is what
`_traitfun` uses in the hot loop — a plain vector fetch, no per-call construction or allocation.
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
    getrelationship(rel::R, field::Symbol) where R <: AbstractTraitRelationship

Extract the trait relationship of all species in the ecosystem.

"""
function getrelationship(rel::R,
                         field::Symbol) where {R <: AbstractTraitRelationship}
    return getfield(rel, field)
end
