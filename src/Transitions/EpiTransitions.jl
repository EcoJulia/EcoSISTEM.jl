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
    EnvViralLoad <: AbstractStateTransition

Transition for force of infection to settle into environmental reservoir,
influenced by the local environment, accessed through `get_env`.

"""
mutable struct EnvViralLoad <: AbstractStateTransition
    location::Int64
    prob::DayType
    get_env::Function
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

function getprob(rule::Exposure)
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
function get_env(habitat::A, loc::Int64) where A <: AbstractHabitat
    return habitat.matrix[loc]
end

"""
    SeedInfection <: AbstractWindDown

Transition to update caches at end of update timestep,
described in `update_fun`.
"""
mutable struct SeedInfection <: AbstractSetUp
    update_fun::Function
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
    classes = length(epi.spplist.species.names)
    Threads.@threads for i in 1:classes
        for j in 1:locs
            epi.cache.forcemigration[human_to_force[i], j] += epi.cache.virusmigration[i, j]
        end
    end
    virus(epi.abundances)[force_cats, :] .= epi.cache.forcemigration
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
    # Update environment - habitat and energy budgets
    habitatupdate!(epi, timestep)
    applycontrols!(epi, timestep)
end


function deterministic_seed!(epi::Ecosystem, controls::Lockdown, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    if (epi.cache.initial_infected > 0) && (controls.current_date < controls.lockdown_date)
        inf = rand(rng, Poisson(epi.cache.initial_infected * timestep /controls.lockdown_date))
        sus_id = epi.spplist.species.susceptible
        exp_id = sus_id .+ length(epi.spplist.species.susceptible)
        pos = epi.cache.ordered_active[1:inf]
        for i in sus_id 
            for j in eachindex(pos)
                if (human(epi.abundances)[sus_id[i], pos[j]] >= 1)
                    human(epi.abundances)[sus_id[i], pos[j]] -= 1
                    human(epi.abundances)[exp_id[i], pos[j]] += 1
                end
            end
        end
    elseif controls.current_date == controls.lockdown_date
        @info "Lockdown initiated - $(sum(human(epi.abundances)[epi.spplist.species.susceptible .+ length(epi.spplist.species.susceptible), :])) individuals infected"
    end
    return controls
end