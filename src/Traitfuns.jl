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
                   trts::Bin,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = gethabitat(hab, pos)
    # Static-habitat counterpart of the `ContinuousTimeHab` method below: fetch the species'
    # pre-built response distribution and evaluate it via `rel` (no per-call construction/allocation).
    return rel(getdist(trts, sp), h)
end
function _traitfun(hab::ContinuousTimeHab,
                   trts::Bin,
                   rel::R,
                   pos::Int64,
                   sp::Int64) where {R <: AbstractTraitRelationship}
    h = gethabitat(hab, pos)
    # Fetch the species' pre-built response distribution and evaluate it via `rel` — no per-call
    # construction or allocation (the distributions were built once when the `Bin` was made).
    return rel(getdist(trts, sp), h)
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
    getdist(traits::Bin, sp::Int64)

Return the pre-built response distribution for species `sp` from a [`Bin`](@ref) trait. This is what
`_traitfun` uses in the hot loop — a plain vector fetch, no per-call construction or allocation.
"""
function getdist(traits::Bin, sp::Int64)
    return traits.dists[sp]
end

# Per-species niche optima (distribution means) of a `Bin`, as a bare vector in the axis's canonical
# frame — the accessor the old `GaussTrait.mean` field readers migrate to (e.g. size traits evolved by
# `ContinuousEvolve`).
_nichemeans(traits::Bin) = Distributions.mean.(traits.dists)

"""
    getpref(traits::Bin, sp::Int64)

Return the parameters of species `sp`'s response distribution from a [`Bin`](@ref) trait (a
back-compat shim over [`getdist`](@ref); the parameters are no longer stored separately).
"""
function getpref(traits::Bin, sp::Int64)
    return Distributions.params(getdist(traits, sp))
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
