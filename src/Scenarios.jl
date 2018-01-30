
using Unitful
using Unitful.DefaultSymbols
using myunitful

"""
    AbstractScenario

Abstract supertype for all whole ecosystem change scenarios
"""
abstract type AbstractScenario end

"""
    SimpleScenario <: AbstractScenario

This scenario type holds a function that acts to change the entire ecosystem.
"""
mutable struct SimpleScenario <: AbstractScenario
    fun::Function
    rate::Quantity{Float64, typeof(𝐓^-1)}
end
"""
    RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that randomly removes a portion of habitat at a certain rate.
"""
function RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, typeof(𝐓^-1)})
    v = uconvert(unit(timestep)^-1, rate)
    pos = find(eco.abenv.budget.matrix .> 0.0)
    howmany = jbinom(1, length(pos), ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] = 0.0
end
"""
    ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that removes a portion of habitat at a certain rate, which spreads through
joined pieces of habitat.
"""
function ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, typeof(𝐓^-1)})
    v = uconvert(unit(timestep)^-1, rate)
    if all(eco.abenv.budget.matrix .> 0.0)
        pos = find(eco.abenv.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.abenv.budget.matrix[smp] = 0.0
    else
        pos = find(eco.abenv.budget.matrix .== 0.0)
        howmany = jbinom(1, length(find(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.matrix,1)
        x, y = convert_coords(pos, width)
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.matrix[i]=0.0
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] = 0.0
end
function HabitatReplacement(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, typeof(𝐓^-1)})
    v = uconvert(unit(timestep)^-1, rate)
    pos = length(eco.abenv.budget.matrix)
    howmany = jbinom(1, pos, ustrip(v))[1]
    smp = sample(1:pos, howmany)
    eco.abenv.habitat.matrix[smp] = maximum(eco.spplist.traits.val) + 1
end

function UniformDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     spp = 1:size(eco.abundances.matrix, 1)
     meanabun = mean(eco.abundances.matrix)
     avlost = rate * timestep * meanabun
     for i in spp
        if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos, Int(round(avlost)),
          replace=true)
         eco.abundances.matrix[i, smp] .-= 1
        end
     end
end

function ProportionalDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     spp = 1:size(eco.abundances.matrix, 1)
     currentabun = mapslices(sum, eco.abundances.matrix, 2)
     avlost = rate * timestep .* currentabun
     for i in spp
         if any(eco.abundances.matrix[i, :] .> 0)
             pos = find(eco.abundances.matrix[i, :] .> 0)
             smp = sample(pos, Int(round(avlost[i])),
              replace=true)
             eco.abundances.matrix[i, smp] .-= 1
         end
     end
end

function LargeDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     !all(eco.spplist.requirement.energy .== maximum(eco.spplist.requirement.energy)) ||
     error("All species have the same requirement")
     largest = find(eco.spplist.requirement.energy .>
     (0.5 * maximum(eco.spplist.requirement.energy)))
     currentabun = mapslices(sum, eco.abundances.matrix, 2)
     avlost = rate * timestep .* currentabun
     for i in largest
         if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos, Int(round(avlost[i])),
          replace=true)
         eco.abundances.matrix[i, smp] .-= 1
     end
     end
end

function RareDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     currentabun = mapslices(sum, eco.abundances.matrix, 2)
     rarest = find(currentabun .< (0.5 * maximum(currentabun)))
     avlost = rate * timestep .* currentabun
     for i in rarest
         if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos, Int(round(avlost[i])),
          replace=true)
         eco.abundances.matrix[i, smp] .-= 1
     end
     end
end
function CommonDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     currentabun = mapslices(sum, eco.abundances.matrix, 2)
     common = find(currentabun .> (0.5 * maximum(currentabun)))
     avlost = rate * timestep .* currentabun
     for i in common
         if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos, Int(round(avlost[i])),
          replace=true)
         eco.abundances.matrix[i, smp] .-= 1
     end
     end
end

function Invasive(eco::Ecosystem, timestep::Unitful.Time,
    rate::Quantity{Float64, typeof(𝐓^-1)})
    invasive = find(eco.spplist.native .== false)
    numspecies = length(eco.spplist.names)
    currentabun = mapslices(sum, eco.abundances.matrix, 2)
    if all(eco.spplist.names .!= "sensitive")
        sensitive = sample(1:numspecies, 5)
        eco.spplist.names[sensitive] = "sensitive"
    end
    avgain = rate * timestep .* currentabun
    for i in invasive
        pos = find(eco.abundances.matrix[i, :])
        smp = sample(pos, Int(round(avgain[i])),
         replace=true)
        eco.abundances.matrix[i, smp] .+= 1
        end
    sensitive = find(eco.spplist.names .== "sensitive")
    for j in sensitive
        pos = find(eco.abundances.matrix[j, :])
        smp = sample(pos, Int(round(avgain[j])),
         replace=true)
        eco.abundances.matrix[j, smp] .-= 1
    end
end



mutable struct DisturbanceScenario <: AbstractScenario
    fun::Function
    loss::Quantity{Float64, typeof(𝐓^-1)}
    level::Float64
    recovery::Quantity{Float64, typeof(𝐓^-1)}
    lag::Unitful.Time
end
"""

A function to create varying degrees of habitat disturbance
"""
function HabitatDisturbance!(eco::Ecosystem, timestep::Unitful.Time, loss::Quantity{Float64, typeof(𝐓^-1)},
    level::Float64, recovery::Quantity{Float64, typeof(𝐓^-1)}, delay::Unitful.Time,
    currentstep::Unitful.Time)
    level < maximum(eco.abenv.budget.matrix) || error("Energy level must be lower than current")
    if ustrip(recovery) == 0.0
        HabitatDisturbance!(eco, timestep, loss, level)
    elseif (ustrip(recovery) > 0.0) & (level == 0.0)
        HabitatDisturbance!(eco, timestep, loss, recovery)
    elseif (ustrip(recovery)> 0.0) & (level > 0.0) & lag == 0.0s
        HabitatDisturbance!(eco, timestep, loss, level, recovery)
    else
        HabitatDisturbance!(eco, timestep, loss, delay, recovery, currentstep)
    end
end

function HabitatDisturbance!(eco::Ecosystem, timestep::Unitful.Time, loss::Quantity{Float64, typeof(𝐓^-1)},
    level::Float64)
    v = uconvert(unit(timestep)^-1, loss)
    if all(eco.abenv.budget.matrix .> level)
        pos = find(eco.abenv.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.abenv.budget.matrix[smp] = level
    else
        pos = find(eco.abenv.budget.matrix .== level)
        howmany = jbinom(1, length(find(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.matrix,1)
        x, y = convert_coords(pos, width)
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.matrix[i] = level
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] = level
end

function HabitatDisturbance!(eco::Ecosystem, timestep::Unitful.Time, loss::Quantity{Float64, typeof(𝐓^-1)},
    recovery::Quantity{Float64, typeof(𝐓^-1)})
    # Recovery step
    # Add recovered energy per timestep
    eco.abenv.budget.matrix[eco.abenv.budget.matrix .<
    maximum(eco.abenv.budget.matrix)] =
    eco.abenv.budget.matrix[eco.abenv.budget.matrix .<
      maximum(eco.abenv.budget.matrix)] .+ (recovery * timestep)
    # Make sure to cap any squares that are now above original starting energy
    eco.abenv.budget.matrix[eco.abenv.budget.matrix .>
      maximum(eco.abenv.budget.matrix)] = maximum(eco.abenv.budget.matrix)

    # Disturbance step
    v = uconvert(unit(timestep)^-1, loss)
    # When starting, pick some grid squares at random
    if all(eco.abenv.budget.matrix .> 0)
        pos = find(eco.abenv.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.abenv.budget.matrix[smp] = 0.0
    # For subsequent steps, only choose from neighbours
    else
        pos = find(eco.abenv.budget.matrix .== 0.0)
        howmany = jbinom(1, length(find(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.matrix,1)
        x, y = convert_coords(pos, width)
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.matrix[i]=0.0
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] = 0.0
end

function HabitatDisturbance!(eco::Ecosystem, timestep::Unitful.Time, loss::Quantity{Float64, typeof(𝐓^-1)},
    level::Float64, recovery::Quantity{Float64, typeof(𝐓^-1)})
    # Recovery step
    # Add recovered energy per timestep
    eco.abenv.budget.matrix[eco.abenv.budget.matrix .<
    maximum(eco.abenv.budget.matrix)] =
    eco.abenv.budget.matrix[eco.abenv.budget.matrix .<
      maximum(eco.abenv.budget.matrix)] .+ (recovery * timestep)
    # Make sure to cap any squares that are now above original starting energy
    eco.abenv.budget.matrix[eco.abenv.budget.matrix .>
      maximum(eco.abenv.budget.matrix)] = maximum(eco.abenv.budget.matrix)

    # Disturbance step
    v = uconvert(unit(timestep)^-1, loss)
    # When starting, pick some grid squares at random
    if all(eco.abenv.budget.matrix .> level)
        pos = find(eco.abenv.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.abenv.budget.matrix[smp] = level
    # For subsequent steps, only choose from neighbours
    else
        pos = find(eco.abenv.budget.matrix .== level)
        howmany = jbinom(1, length(find(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.matrix,1)
        x, y = convert_coords(pos, width)
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.matrix[i] = level
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] = level
end

function HabitatDisturbance!(eco::Ecosystem, timestep::Unitful.Time, loss::Quantity{Float64, typeof(𝐓^-1)},
    delay::Unitful.Time, recovery::Quantity{Float64, typeof(𝐓^-1)},
    currentstep::Unitful.Time)
    # Delayed recovery step
    # Add recovered energy per timestep
    # Only after delay period
    if currentstep > delay
        eco.abenv.budget.matrix[eco.abenv.budget.matrix .<
        maximum(eco.abenv.budget.matrix)] =
        eco.abenv.budget.matrix[eco.abenv.budget.matrix .<
          maximum(eco.abenv.budget.matrix)] .+ (recovery * timestep)
        # Make sure to cap any squares that are now above original starting energy
        eco.abenv.budget.matrix[eco.abenv.budget.matrix .>
          maximum(eco.abenv.budget.matrix)] = maximum(eco.abenv.budget.matrix)
     end

    # Disturbance step
    v = uconvert(unit(timestep)^-1, loss)
    # When starting, pick some grid squares at random
    if all(eco.abenv.budget.matrix .> level)
        pos = find(eco.abenv.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.abenv.budget.matrix[smp] = level
    # For subsequent steps, only choose from neighbours
    else
        pos = find(eco.abenv.budget.matrix .== level)
        howmany = jbinom(1, length(find(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.matrix,1)
        x, y = convert_coords(pos, width)
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.matrix[i] = level
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] = level
end


function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::SimpleScenario,
    currentstep::Unitful.Time)
    scenario.fun(eco, timestep, scenario.rate)
end
function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::DisturbanceScenario,
    currentstep::Unitful.Time)
    scenario.fun(eco, timestep, scenario.loss, scenario.level, scenario.recovery,
    scenario.lag, currentstep)
end
