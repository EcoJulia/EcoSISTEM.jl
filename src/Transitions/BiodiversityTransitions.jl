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
