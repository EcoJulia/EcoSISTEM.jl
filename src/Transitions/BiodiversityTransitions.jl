const TimeType = typeof(1.0/year)

"""
    BirthProcess <: AbstractStateTransition

Stochastic birth process for a location and species, with an associated probability. An optional destination 
    argument can be added if the birth should be added to a different category (i.e. age structured).
"""
mutable struct BirthProcess <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::TimeType
    function BirthProcess(species::Int64, location::Int64, prob::T) where T
        prob = uconvert(unit(TimeType), prob)
        new(species, location, species, prob)
    end
    function BirthProcess(species::Int64, location::Int64, destination::Int64, prob::T) where T
        prob = uconvert(unit(TimeType), prob)
        new(species, location, destination, prob)
    end
end

"""
    GenerateSeed <: AbstractStateTransition

Stochastic seed generation process for a location and species, with an associated probability. An optional destination 
    argument can be added if the seed produced should be added to a different category (i.e. age structured).
"""
mutable struct GenerateSeed <: AbstractStateTransition
    species::Int64
    location::Int64
    destination::Int64
    prob::TimeType
    function GenerateSeed(species::Int64, location::Int64, prob::T) where T
        prob = uconvert(unit(TimeType), prob)
        new(species, location, species, prob)
    end
    function GenerateSeed(species::Int64, location::Int64, destination::Int64, prob::T) where T
        prob = uconvert(unit(TimeType), prob)
        new(species, location, destination, prob)
    end
end

"""
    DeathProcess <: AbstractStateTransition

Stochastic death process for a location and species, with an associated probability.
"""
mutable struct DeathProcess <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeType
    function DeathProcess(species::Int64, location::Int64, prob::T) where T
        prob = uconvert(unit(TimeType), prob)
        new(species, location, prob)
    end
end

"""
    AllDisperse <: AbstractPlaceTransition

Stochastic dispersal for a species at a location.
"""
mutable struct AllDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

"""
    SeedDisperse <: AbstractPlaceTransition

Stochastic seed dispersal for a species at a location.
"""
mutable struct SeedDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

"""
    UpdateEnergy <: AbstractSetUp

Houses a function to update energy usage across the Ecosystem.
"""
mutable struct UpdateEnergy <: AbstractSetUp
    update_fun::Function
end

"""
    UpdateEnvironment <: AbstractWindDown

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
