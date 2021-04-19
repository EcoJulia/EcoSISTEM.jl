const DayType = typeof(1.0/day)
"""
    ViralLoad <: AbstractStateTransition

Transition for force of infection to settle into
    environmental reservoir.
"""
mutable struct ViralLoad <: AbstractStateTransition
    location::Int64
    prob::DayType
end

"""
    Exposure <: AbstractStateTransition

Transition for exposure of susceptibles to force of infection
    and environmental reservoir at set probabilities,
    `force_prob` and `virus_prob` respectively.
"""
mutable struct Exposure <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    force_prob::DayType
    virus_prob::DayType
end

"""
    EnvExposure <: AbstractStateTransition

Transition for exposure of susceptibles to force of infection
    and environmental reservoir at set probabilities,
    `force_prob` and `virus_prob` respectively. This exposure
    mechanism is influenced by the local environment (fetched
    through the `get_env` function), controlled by the strength
    of `env_param`.
"""
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

"""
    Infection <: AbstractStateTransition

Transition from exposed to infectious categories at a
set probability, `prob`. Any infectious category can
be designated through `destination`.
"""
mutable struct Infection <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

"""
    DevelopSymptoms <: AbstractStateTransition

Transition from pre-infectious to infectious categories at a
set probability, `prob`. Any infectious category can
be designated through `destination`.
"""
mutable struct DevelopSymptoms <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

"""
    Hospitalise <: AbstractStateTransition

Transition from infectious to hospitalised categories at a
set probability, `prob`.
"""
mutable struct Hospitalise <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

"""
    Recovery <: AbstractStateTransition

Transition from infected to recovered category at a
set probability, `prob`.
"""
mutable struct Recovery <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

"""
    DeathFromInfection <: AbstractStateTransition

Transition from infected to dead category at a
set probability, `prob`.
"""
mutable struct DeathFromInfection <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::DayType
end

"""
    ForceProduce <: AbstractStateTransition

Transition to produce force of infection from an infectious category at a
set probability, `prob`.
"""
mutable struct ForceProduce <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::DayType
end

function getprob(rule::ForceProduce)
    return rule.prob
end

"""
    ForceDisperse <: AbstractStateTransition

Transition to disperse force of infection from a `location`.
"""
mutable struct ForceDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

"""
    get_env(habitat::A) where A <: AbstractHabitat

Generic function to extract a habitat matrix from any
habitat type.
"""
function get_env(habitat::A) where A <: AbstractHabitat
    return habitat.matrix
end

"""
    UpdateEpiEnvironment <: AbstractWindDown

Transition to update caches at end of update timestep,
described in `update_fun`.
"""
mutable struct UpdateEpiEnvironment <: AbstractWindDown
    update_fun::Function
end

"""
    update_virus_cache!(epi::Ecosystem)

Function to update the virus abundances with cached
migration moves.
"""
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

"""
    update_epi_environment!(epi::Ecosystem, timestep::Unitful.Time)

Function to update the virus abundances with cached
migration moves and invalidate all caches.
"""
function update_epi_environment!(epi::Ecosystem, timestep::Unitful.Time)
    update_virus_cache!(epi)
    # Invalidate all caches for next update
    invalidatecaches!(epi)
end
