"""
    traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64, ::S) where S <: AbstractSpeciesTypes

Function to calculate relationship between the current environment and a species' particular trait.

"""
function traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64, ::S) where S <: AbstractSpeciesTypes
    hab = eco.abenv.habitat
    trts = getspeciestraits(eco)
    rel = eco.relationship
    _traitfun(hab, trts, rel, pos, sp)
end

"""
    traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64, ::P) where P <: AbstractPathogenTypes

Function to calculate relationship between the current environment and a pathogen's particular trait.

"""
function traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64, ::P) where P <: AbstractPathogenTypes
    hab = eco.abenv.habitat
    trts = getpathogentraits(eco)
    rel = eco.relationship
    _traitfun(hab, trts, rel, pos, sp)
end

function _traitfun(hab::HabitatCollection2, trts::TraitCollection2, rel::R, pos::Int64, sp::Int64) where R <: AbstractTraitRelationship
    res1 = _traitfun(hab.h1, trts.t1, rel.tr1, pos, sp)
    res2 = _traitfun(hab.h2, trts.t2, rel.tr2, pos, sp)
    return combineTR(rel)(res1, res2)
end

function _traitfun(hab::HabitatCollection3, trts::TraitCollection3, rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
    res1 = _traitfun(hab.h1, trts.t1, rel.tr1, pos, spp)
    res2 = _traitfun(hab.h2, trts.t2, rel.tr2, pos, spp)
    res3 = _traitfun(hab.h3, trts.t3, rel.tr3, pos, spp)
    return combineTR(rel)(res1, res2, res3)
end

function _traitfun(hab::ContinuousHab, trts::GaussTrait,
    rel::R, pos::Int64, sp::Int64) where R <: AbstractTraitRelationship
    h = hab.matrix[pos]
    mean, var = getpref(trts, sp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab, trts::GaussTrait,
    rel::R, pos::Int64, sp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        mean, var = getpref(trts, sp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab, trts::TempBin,
    rel::R, pos::Int64, sp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        (a, b, c, d) = getpref(trts, sp)
    return rel(Trapezoid(a, b, c, d), h)
end
function _traitfun(hab::ContinuousTimeHab, trts::RainBin,
    rel::R, pos::Int64, sp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        (a, b) = getpref(trts, sp)
    return rel(Uniform(a, b), h)
end
function _traitfun(hab::DiscreteHab, trts::DiscreteTrait,
    rel::R, pos::Int64, sp::Int64) where R <: AbstractTraitRelationship
        currentniche = gethabitat(hab, pos)
        preference = getpref(trts, sp)
    return rel(currentniche, preference)
end

function _traitfun(hab::DiscreteHab, trts::LCtrait,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        vals = getpref(trts, spp)
    return rel(h, vals)
end

function getpref(traits::LCtrait, spp::Int64)
  return traits.vals[spp]
end

function getpref(traits::GaussTrait, sp::Int64)
  return traits.mean[sp], traits.var[sp]
end
function getpref(traits::TempBin, sp::Int64)
  return traits.dist[sp, :]
end
function getpref(traits::RainBin, sp::Int64)
  return traits.dist[sp, :]
end
function getpref(traits::DiscreteTrait, sp::Int64)
  return traits.val[sp]
end

"""
    getpref(traits::T, field::Symbol) where T <: AbstractTraits

Function to extract trait preferences for all species in the ecosystem.

"""
function getpref(traits::T, field::Symbol) where T <: AbstractTraits
  return getfield(traits, field)
end

"""
    getpref(traits::T, field::Symbol) where T <: AbstractTraits

Function to extract the trait relationship of all species in the ecosystem.

"""
function getrelationship(rel::R, field::Symbol) where R <: AbstractTraitRelationship
  return getfield(rel, field)
end
