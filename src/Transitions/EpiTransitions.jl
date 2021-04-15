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

mutable struct DevelopSymptoms <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

mutable struct Hospitalise <: AbstractStateTransition
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

mutable struct DeathFromInfection <: AbstractStateTransition
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
    force_cats = epi.spplist.pathogens.force_cats
    human_to_force = epi.spplist.species.human_to_force
    locs = size(virus(epi.abundances), 2)
    vm = zeros(eltype(epi.cache.virusmigration), length(force_cats), locs)
    classes = length(epi.spplist.species.names)
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
