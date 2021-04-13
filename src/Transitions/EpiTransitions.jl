const DayType = typeof(1.0/day)
mutable struct ViralLoad <: AbstractStateTransition
    location::Int64
    prob::DayType
end

mutable struct Exposure <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    force_prob::DayType
    virus_prob::DayType
end

mutable struct EnvExposure <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    force_prob::DayType
    virus_prob::DayType
    env_param::Float64
    get_env::Function
end

function getprob(rule::Exposure)
    return rule.force_prob, rule.virus_prob
end

function getprob(rule::EnvExposure)
    return rule.force_prob, rule.virus_prob
end

mutable struct Infection <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

mutable struct Recovery <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

mutable struct ForceProduce <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::DayType
end

function getprob(rule::ForceProduce)
    return rule.prob
end

mutable struct ForceDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

function get_env(habitat::A) where A <: AbstractHabitat
    return habitat.matrix
end


mutable struct UpdateEpiEnvironment <: AbstractWindDown
    update_fun::Function
end
function update_virus_cache!(epi::Ecosystem)
    force_cats = epi.spplist.virus.force_cats
    human_to_force = epi.spplist.human.human_to_force
    locs = size(virus(epi.abundances), 2)
    vm = zeros(eltype(epi.cache.virusmigration), length(force_cats), locs)
    classes = length(epi.spplist.human.names)
    Threads.@threads for i in 1:classes
        for j in 1:locs
            vm[human_to_force[i], j] += epi.cache.virusmigration[i, j]
        end
    end
    virus(epi.abundances)[force_cats, :] .= vm
end

function update_epi_environment!(epi::Ecosystem, timestep::Unitful.Time)
    update_virus_cache!(epi)
    # Invalidate all caches for next update
    invalidatecaches!(epi)
end

function create_transition_list_SEIR(epilist::EpiList, epienv::GridEpiEnv)
    params = epilist.params
    paramDat = params.transition_dat
    num_cats = length(paramDat[1, :from_id])

    exposure = [Exposure(paramDat[1, :from_id][n], loc, paramDat[1, :to_id][n], paramDat[1, :prob].force[n], paramDat[1, :prob].env[n]) for n in 1:num_cats for loc in eachindex(epienv.habitat.matrix)]
    infection = [Infection(paramDat[2, :from_id][n], loc, paramDat[2, :to_id][n], paramDat[2, :prob][n]) for n in 1:num_cats for loc in eachindex(epienv.habitat.matrix)]
    recovery = [Recovery(paramDat[3, :from_id][n], loc, paramDat[3, :to_id][n], paramDat[3, :prob][n]) for n in 1:num_cats for loc in eachindex(epienv.habitat.matrix)]

    virus_list = [ViralLoad(loc, params.virus_decay) for loc in eachindex(epienv.habitat.matrix)]
    force_list = [ForceProduce(spp, loc, params.virus_growth[spp]) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]

    state_list = [force_list; virus_list; exposure; infection; recovery]

    place_list = [ForceDisperse(spp, loc) for spp in eachindex(epilist.human.names) for loc in eachindex(epienv.habitat.matrix)]

    before = missing
    after = UpdateEpiEnvironment(update_epi_environment!)
    return TransitionList(before, state_list, place_list, after)
end

function create_transition_list_env_SEIR(epilist::EpiList, epienv::GridEpiEnv, env_params::NamedTuple)
    params = epilist.params
    paramDat = params.transition_dat
    num_cats = length(paramDat[1, :from_id])

    exposure = [EnvExposure(paramDat[1, :from_id][n], loc,
                paramDat[1, :to_id][n], paramDat[1, :prob].force[n],
                paramDat[1, :prob].env[n], env_params.env_exposure, get_env)
                for n in 1:num_cats for loc in eachindex(epienv.habitat.matrix)]
    infection = [Infection(paramDat[2, :from_id][n], loc,
                paramDat[2, :to_id][n],
                paramDat[2, :prob][n])
                for n in 1:num_cats for loc in eachindex(epienv.habitat.matrix)]
    recovery = [Recovery(paramDat[3, :from_id][n], loc,
                paramDat[3, :to_id][n],
                paramDat[3, :prob][n])
                for n in 1:num_cats for loc in eachindex(epienv.habitat.matrix)]

    virus_list = [ViralLoad(loc, params.virus_decay)
                    for loc in eachindex(epienv.habitat.matrix)]
    force_list = [ForceProduce(spp, loc, params.virus_growth[spp])
                    for spp in eachindex(epilist.human.names)
                    for loc in eachindex(epienv.habitat.matrix)]

    state_list = [force_list; virus_list; exposure; infection; recovery]

    place_list = [ForceDisperse(spp, loc)
                    for spp in eachindex(epilist.human.names)
                    for loc in eachindex(epienv.habitat.matrix)]

    before = missing
    after = UpdateEpiEnvironment(update_epi_environment!)
    return TransitionList(before, state_list, place_list, after)
end
