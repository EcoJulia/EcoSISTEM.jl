mutable struct ViralLoad{U <: Unitful.Units} <: AbstractStateTransition
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct Exposure{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    force_prob::TimeUnitType{U}
    virus_prob::TimeUnitType{U}
end

function getprob(rule::Exposure)
    return rule.force_prob, rule.virus_prob
end

mutable struct Infection{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::TimeUnitType{U}
end

mutable struct Recovery{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::TimeUnitType{U}
end

mutable struct SEIR{U <: Unitful.Units} <: AbstractStateTransition
    exposure::Exposure{U}
    infection::Infection{U}
    recovery::Recovery{U}
end

mutable struct ForceProduce{U <: Unitful.Units} <: AbstractPlaceTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

function getprob(rule::ForceProduce)
    return rule.prob
end

mutable struct ForceDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

mutable struct Force{U <: Unitful.Units}
    forceprod::ForceProduce{U}
    forcedisp::ForceDisperse
end

function create_transition_list(epilist::EpiList, epienv::GridEpiEnv)
    params = epilist.params

    state_list = [SEIR(Exposure(1, loc, 2, params.transition_force[2, 1], params.transition_virus[2, 1]), Infection(2, loc, 3, params.transition[3, 2]), Recovery(3, loc, 4, params.transition[4, 3])) for loc in eachindex(epienv.habitat.matrix)]

    virus_list = [ViralLoad(loc, params.virus_decay) for loc in eachindex(epienv.habitat.matrix)]

    state_list = [virus_list; state_list]

    place_list = [ForceDisperse(spp, loc) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]
    force_list = [ForceProduce(spp, loc, params.virus_growth[spp]) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]
    place_list = [force_list; place_list]

    return TransitionList(state_list, place_list)
end
