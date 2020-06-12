"""
    traitfun(epi::AbstractEpiSystem, pos::Int64, sp::Int64)

Function to calculate relationship between the current environment and a particular trait of a disease class.

"""
function traitfun(epi::AbstractEpiSystem, pos::Int64, sp::Int64)
    hab = epi.epienv.habitat
    trts = gettraits(epi)
    rel = epi.relationship
    return _traitfun(hab, trts, rel, pos, sp)
end

gettraits(epi::EpiSystem) = epi.epilist.virus.traits
