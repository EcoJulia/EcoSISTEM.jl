mutable struct BirthProcess{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct DeathProcess{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct BirthDeathProcess{U <: Unitful.Units} <: AbstractStateTransition
    birth::BirthProcess{U}
    death::DeathProcess{U}
end

mutable struct AllDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end


function create_transitions(spplist::SpeciesList, abenv::GridAbioticEnv)

    state_list = [BirthDeathProcess(BirthProcess(spp, loc, spplist.params.birth[spp]), DeathProcess(spp, loc, spplist.params.death[spp])) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]

    place_list = [AllDisperse(spp, loc) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]

    return TransitionList(state_list, place_list)
end
