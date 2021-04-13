const TimeType = typeof(1.0/year)
mutable struct BirthProcess <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeType
end

mutable struct DeathProcess <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeType
end

mutable struct AllDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

mutable struct UpdateEnergy <: AbstractSetUp
    update_fun::Function
end

mutable struct UpdateEnvironment <: AbstractWindDown
    update_fun::Function
end
function update_environment!(eco::Ecosystem, timestep::Unitful.Time)
    eco.abundances.matrix .+= eco.cache.netmigration

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
end

function create_transitions(spplist::SpeciesList, abenv::GridAbioticEnv)

    births = [BirthProcess(spp, loc, spplist.params.birth[spp])  for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]
    deaths = [DeathProcess(spp, loc, spplist.params.death[spp]) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]
    state_list = [births; deaths]
    place_list = [AllDisperse(spp, loc) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]
    before = UpdateEnergy(update_energy_usage!)
    after = UpdateEnvironment(update_environment!)
    return TransitionList(before, state_list, place_list, after)
end
