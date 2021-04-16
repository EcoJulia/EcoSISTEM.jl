const TimeType = typeof(1.0/year)

"""
    _run_rule!(eco::Ecosystem, rule::BirthProcess, timestep::Unitful.Time)

Stochastic birth process for a location and species, with an associated probability.
"""
mutable struct BirthProcess <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeType
end

"""
    _run_rule!(eco::Ecosystem, rule::GenerateSeed, timestep::Unitful.Time)

Stochastic seed generation process for a location and species, with an associated probability.
"""
mutable struct GenerateSeed <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeType
end

"""
    _run_rule!(eco::Ecosystem, rule::DeathProcess, timestep::Unitful.Time)

Stochastic death process for a location and species, with an associated probability.
"""
mutable struct DeathProcess <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeType
end

"""
    _run_rule!(eco::Ecosystem, rule::AllDisperse, timestep::Unitful.Time)

Stochastic dispersal for a species at a location.
"""
mutable struct AllDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

"""
    _run_rule!(eco::Ecosystem, rule::SeedDisperse, timestep::Unitful.Time)

Stochastic seed dispersal for a species at a location.
"""
mutable struct SeedDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

"""
    _run_rule!(eco::Ecosystem, rule::UpdateEnergy, timestep::Unitful.Time)

Houses a function to update energy usage across the Ecosystem.
"""
mutable struct UpdateEnergy <: AbstractSetUp
    update_fun::Function
end

"""
    _run_rule!(eco::Ecosystem, rule::UpdateEnvironment, timestep::Unitful.Time)

Houses a function to update the habitat and resources across the Ecosystem.
"""
mutable struct UpdateEnvironment <: AbstractWindDown
    update_fun::Function
end

"""
    update_environment!(eco::Ecosystem, timestep::Unitful.Time)

Function to update the habitat and resources across the Ecosystem.
"""
function update_environment!(eco::Ecosystem, timestep::Unitful.Time)
    eco.abundances.matrix .+= eco.cache.netmigration

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
end
