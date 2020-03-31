"""
    traitfun(eco::Ecosystem, pos::Int64, sp::Int64)

Function to calculate relationship between the current environment and a species' particular trait.

"""
function traitfun(eco::Ecosystem, pos::Int64, sp::Int64)
    hab = eco.abenv.habitat
    trts = eco.spplist.traits
    rel = eco.relationship
  _traitfun(hab, trts, rel, pos, spp)
end

function _traitfun(hab::HabitatCollection2, trts::TraitCollection2, rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
    res1 = _traitfun(hab.h1, trts.t1, rel.tr1, pos, spp)
    res2 = _traitfun(hab.h2, trts.t2, rel.tr2, pos, spp)
    return combineTR(rel)(res1, res2)
end
function _traitfun(hab::ContinuousHab, trts::GaussTrait,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        mean, var = getpref(trts, spp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab, trts::GaussTrait,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        mean, var = getpref(trts, spp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab, trts::TempBin,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        (a, b, c, d) = getpref(trts, spp)
    return rel(Trapezoid(a, b, c, d), h)
end
function _traitfun(hab::ContinuousTimeHab, trts::RainBin,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        (a, b) = getpref(trts, spp)
    return rel(Uniform(a, b), h)
end
function _traitfun(hab::DiscreteHab, trts::DiscreteTrait,
    rel::R, pos::Int64, spp::Int64) where R <: AbstractTraitRelationship
        currentniche = gethabitat(hab, pos)
        preference = getpref(trts, spp)
    return rel(currentniche, preference)
end

function getpref(traits::GaussTrait, spp::Int64)
  return traits.mean[spp], traits.var[spp]
end
function getpref(traits::TempBin, spp::Int64)
  return traits.dist[spp, :]
end
function getpref(traits::RainBin, spp::Int64)
  return traits.dist[spp, :]
end
function getpref(traits::DiscreteTrait, spp::Int64)
  return traits.val[spp]
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
