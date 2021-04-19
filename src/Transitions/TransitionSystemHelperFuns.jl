"""
    gettraitrel(eco::Ecosystem)

Function to extract trait relationships.
"""
function gettraitrel(eco::AbstractEcosystem)
  return eco.relationship
end

"""
    gethabitat(eco::Ecosystem)

Function to extract habitat from Ecosystem object.
"""
function gethabitat(eco::AbstractEcosystem)
  return eco.abenv.habitat
end
"""
    getbudget(eco::Ecosystem)

Function to extract budget from Ecosystem object.
"""
function getbudget(eco::AbstractEcosystem)
    return _getbudget(eco.abenv.budget)
end

"""
    getsize(eco::Ecosystem)

Function to extract size of habitat from Ecosystem object.
"""
function getsize(eco::AbstractEcosystem)
  return _getsize(eco.abenv.habitat)
end

"""
    getgridsize(eco::Ecosystem)

Function to extract grid cell size of habitat from Ecosystem object.
"""
function getgridsize(eco::AbstractEcosystem)
  return _getgridsize(eco.abenv.habitat)
end

"""
    getdimension(eco::Ecosystem)

Function to extract dimension of habitat from Ecosystem object.
"""
function getdimension(eco::AbstractEcosystem)
    return _getdimension(eco.abenv.habitat)
end

"""
    getdispersaldist(eco::Ecosystem)

Function to extract average dispersal distance of species from Ecosystem object.
Returns a vector of distances, unless a specific species is provided as a String
or Integer.
"""
function getdispersaldist(eco::A, sp::Int64) where A <: AbstractEcosystem
    return getdispersaldist(eco.spplist.species.movement, sp)
end
function getdispersaldist(eco::A, sp::String) where A <: AbstractEcosystem
    sp = findfirst(eco.spplist.species.names .== sp)
    return getdispersaldist(eco.spplist.species.movement, sp)
end

"""
    getdispersalvar(eco::Ecosystem)

Function to extract dispersal varaince of species from Ecosystem object.
Returns a vector of distances, unless a specific species is provided as a String
or Integer.
"""
function getdispersalvar(eco::A, sp::Int64) where A <: AbstractEcosystem
    return getdispersalvar(eco.spplist.species.movement, sp)
end
function getdispersalvar(eco::A, sp::String) where A <: AbstractEcosystem
    sp = findfirst(eco.spplist.species.names .== sp)
    return getdispersalvar(eco.spplist.species.movement, sp)
end

"""
    getlookup(eco::Ecosystem)

Function to extract movement lookup table of species from Ecosystem object.
"""
function getlookup(eco::A, sp::Int64) where A <: AbstractEcosystem
    return _getlookup(eco.lookup, sp)
end

function getlookup(eco::Ecosystem, id::Int64, movetype::MovementType)
    if movetype == homeMovement
        return eco.lookup.homelookup[id, :]
    elseif movetype == workMovement
        return eco.lookup.worklookup[id, :]
    else
        return error("No other movement types currently implemented")
    end
end

"""
    resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝐓^-1)})

Function to reset the rate of habitat change for a species.
"""
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, typeof(𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, Unitful.Dimensions{()})
end
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, typeof(dimension(1K)))
end
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, 𝐓^-1})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, Unitful.Dimensions{()})
end
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, 𝚯*𝐓^-1})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, typeof(dimension(1K)))
end

"""
    resettime!(eco::AbstractEcosystem)

Function to reset the time of habitat change for a species.
"""
function resettime!(eco::AbstractEcosystem)
    _resettime!(eco.abenv.habitat)
end

"""
    addspecies!(eco::Ecosystem, abun::Int64)

Function to add a species to the Ecosystem.
"""
function addspecies!(eco::Ecosystem, abun::Int64)
    eco.abundances.matrix = vcat(eco.abundances.matrix, zeros(1, size(eco.abundances.matrix, 2)))
    eco.abundances.grid = reshape(eco.abundances.matrix, (counttypes(eco.spplist, true)+1, _getdimension(eco.abenv.habitat)...))
    repopulate!(eco, abun)
    push!(eco.spplist.species.names, string.(counttypes(eco.spplist, true)+1))
    append!(eco.spplist.species.abun, abun)
    append!(eco.spplist.species.native, true)
    addtraits!(eco.spplist.species.traits)
    addmovement!(eco.spplist.species.movement)
    addparams!(eco.spplist.params)
    addrequirement!(eco.spplist.species.requirement)
    addtypes!(eco.spplist.species.types)
end
function addtraits!(tr::GaussTrait)
    append!(tr.mean, tr.mean[end])
    append!(tr.var, tr.var[end])
end

function addtraits!(tr::DiscreteTrait)
    append!(tr.val, rand(tr.val))
end

addmovement!(mv::AbstractMovement) = push!(mv.kernels, mv.kernels[end])

function addparams!(pr::AbstractParams)
    append!(pr.birth, pr.birth[end])
    append!(pr.death, pr.death[end])
end

addrequirement!(rq::AbstractRequirement) = append!(rq.energy, rq.energy[end])

function addtypes!(ut::UniqueTypes)
    ut = UniqueTypes(ut.num+1)
end

"""
    invalidatecaches!(eco::AbstractEcosystem)

Function to invalidate Ecosystem caches at the end of a timestep.
"""
function invalidatecaches!(eco::AbstractEcosystem)
    _invalidatecaches!(eco, eco.cache)
end
function _invalidatecaches!(eco::A, cache::Cache) where A <: AbstractEcosystem
    eco.ordinariness = missing
    eco.cache.netmigration .= 0
    eco.cache.valid = false
end
function _invalidatecaches!(eco::A, cache::PlantCache) where A <: AbstractEcosystem
    eco.ordinariness = missing
    eco.cache.netmigration .= 0
    eco.cache.seedbank .= 0
    eco.cache.valid = false
end
function _invalidatecaches!(eco::Ecosystem, cache::EpiCache)
    eco.ordinariness = missing
    eco.cache.virusmigration .= 0
    eco.cache.valid = false
end
