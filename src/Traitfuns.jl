# SPDX-License-Identifier: LGPL-3.0-or-later

"""
    traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64)

Calculate relationship between the current environment and a species' particular
trait.

"""
function traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64)
    hab = eco.abenv.habitat
    trts = eco.spplist.traits
    rel = eco.relationship
    return _traitfun(hab, trts, rel, pos, sp)
end

function _traitfun(hab::HabitatCollection2,
                   trts::TraitCollection2,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    res1 = _traitfun(hab.one, trts.one, rel.one, pos, sp)
    res2 = _traitfun(hab.two, trts.two, rel.two, pos, sp)
    return combineTR(rel)(res1, res2)
end
function _traitfun(hab::ContinuousHab,
                   trts::GaussTrait,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = gethabitat(hab, pos)
    mean, var = getpref(trts, sp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab,
                   trts::GaussTrait,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = gethabitat(hab, pos)
    mean, var = getpref(trts, sp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab,
                   trts::TempBin,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = gethabitat(hab, pos)
    (a, b, c, d) = getpref(trts, sp)
    return rel(Trapezoid(a, b, c, d), h)
end
function _traitfun(hab::ContinuousTimeHab,
                   trts::RainBin,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = gethabitat(hab, pos)
    (a, b) = getpref(trts, sp)
    return rel(Uniform(a, b), h)
end
function _traitfun(hab::DiscreteHab,
                   trts::DiscreteTrait,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    currentniche = gethabitat(hab, pos)
    preference = getpref(trts, sp)
    return rel(currentniche, preference)
end
function _traitfun(hab::HabitatCollection3,
                   trts::TraitCollection3,
                   rel::R,
                   pos::Int64,
                   spp::Int64) where {R <: AbstractTraitRelationship}
    res1 = _traitfun(hab.one, trts.one, rel.one, pos, spp)
    res2 = _traitfun(hab.two, trts.two, rel.two, pos, spp)
    res3 = _traitfun(hab.three, trts.three, rel.three, pos, spp)
    return combineTR(rel)(res1, res2, res3)
end

function _traitfun(hab::DiscreteHab,
                   trts::LCtrait,
                   rel::R,
                   pos::Int64,
                   spp::Int64) where {R <: AbstractTraitRelationship}
    h = gethabitat(hab, pos)
    vals = getpref(trts, spp)
    return rel(h, vals)
end

"""
    getpref(traits::LCtrait, spp::Int64)

Extract the land cover preference values for species `spp` from an `LCtrait`.
"""
function getpref(traits::LCtrait, spp::Int64)
    return traits.vals[spp]
end

"""
    getpref(traits::GaussTrait, sp::Int64)

Extract the Gaussian habitat preference optimum (`mean`) and standard deviation
(`sd`) for species `sp` from a [`GaussTrait`](@ref). Returns a tuple `(mean, sd)`.
"""
function getpref(traits::GaussTrait, sp::Int64)
    return traits.mean[sp], traits.sd[sp]
end

"""
    getpref(traits::TempBin, sp::Int64)

Extract the trapezoid distribution parameters `(a, b, c, d)` for species `sp`
from a [`TempBin`](@ref) trait.
"""
function getpref(traits::TempBin, sp::Int64)
    return traits.dist[sp, :]
end

"""
    getpref(traits::RainBin, sp::Int64)

Extract the uniform distribution parameters `(a, b)` for species `sp` from a
[`RainBin`](@ref) trait.
"""
function getpref(traits::RainBin, sp::Int64)
    return traits.dist[sp, :]
end

"""
    getpref(traits::DiscreteTrait, sp::Int64)

Extract the discrete niche preference value for species `sp` from a
[`DiscreteTrait`](@ref).
"""
function getpref(traits::DiscreteTrait, sp::Int64)
    return traits.val[sp]
end

"""
    getpref(traits::T, field::Symbol) where T <: AbstractTraits

Extract trait preferences for all species in the ecosystem.

"""
function getpref(traits::T, field::Symbol) where {T <: AbstractTraits}
    return getfield(traits, field)
end

"""
    getrelationship(rel::R, field::Symbol) where R <: AbstractTraitRelationship

Extract the trait relationship of all species in the ecosystem.

"""
function getrelationship(rel::R,
                         field::Symbol) where {R <: AbstractTraitRelationship}
    return getfield(rel, field)
end
