"""
    traitfun(eco::AbstractEcosystem, pos::Int64, sp::Int64)

Function to calculate relationship between the current environment and a species' particular trait.

"""
function traitfun(epi::AbstractEpiSystem, pos::Int64, sp::Int64)
    hab = epi.epienv.habitat
    trts = epi.epilist.traits
    rel = epi.relationship
  _traitfun(hab, trts, rel, pos, sp)
end
