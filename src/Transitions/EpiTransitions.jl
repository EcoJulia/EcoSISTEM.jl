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

mutable struct EnvExposure{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    force_prob::TimeUnitType{U}
    virus_prob::TimeUnitType{U}
    env_field::Symbol
    env_param::Unitful.Quantity{Float64}
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
    paramDat = params.transition_dat
    num_cats = length(paramDat[1, :from_id])

    state_list = [SEIR(Exposure(paramDat[1, :from_id][n], loc, paramDat[1, :to_id][n], paramDat[1, :prob].force[n], paramDat[1, :prob].env[n]), Infection(paramDat[2, :from_id][n], loc, paramDat[2, :to_id][n], paramDat[2, :prob][n]), Recovery(paramDat[3, :from_id][n], loc, paramDat[3, :to_id][n], paramDat[3, :prob][n])) for n in 1:num_cats for loc in eachindex(epienv.habitat.matrix)]

    virus_list = [ViralLoad(loc, params.virus_decay) for loc in eachindex(epienv.habitat.matrix)]

    state_list = [virus_list; state_list]

    place_list = [ForceDisperse(spp, loc) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]
    force_list = [ForceProduce(spp, loc, params.virus_growth[spp]) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]
    place_list = [force_list; place_list]

    return TransitionList(state_list, place_list)
end
