
using Unitful
using Unitful.DefaultSymbols
using MyUnitful

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
    eco.abundances.matrix[:, smp] .= 0.0
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
        eco.abundances.matrix[:, smp] .= 0.0
    else
        pos = find(eco.abenv.budget.matrix .== 0.0)
        howmany = jbinom(1, length(find(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.matrix,1)
        x, y = convert_coords(pos, width)
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.matrix[i]=0.0
        eco.abundances.matrix[:, i] .= 0.0
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] = 0.0
    eco.abundances.matrix[:, smp] .= 0.0
end

"""
     HabitatReplacement(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that randomly replaces a portion of habitat with another.
"""
function HabitatReplacement(eco::Ecosystem, timestep::Unitful.Time, rate::Quantity{Float64, typeof(𝐓^-1)})
    v = uconvert(unit(timestep)^-1, rate)
    pos = length(eco.abenv.budget.matrix)
    howmany = jbinom(1, pos, ustrip(v))[1]
    smp = sample(1:pos, howmany)
    eco.abenv.habitat.matrix[smp] = maximum(eco.spplist.traits.val) + 1
end

function SusceptibleDecline(eco::Ecosystem, timestep::Unitful.Time,
            rate::Quantity{Float64, typeof(𝐓^-1)})
     eco.spplist.susceptible *= (1 + rate * timestep)
     currentabun = mapslices(sum, eco.abundances.matrix, 2)
     spp = 1:size(eco.abundances.matrix, 1)
     prob  = 1 ./ (1 .+ exp.(-eco.spplist.susceptible))
     howmany =
         map(currentabun, prob) do x, y
             jbinom(1,x, y)[1]
         end
    for i in spp
        for j in 1:howmany[i]
            if any(eco.abundances.matrix[i, :] .> 0)
            pos = find(eco.abundances.matrix[i, :] .> 0)
            smp = sample(pos)
            eco.abundances.matrix[i, smp] .-= 1
        end
    end
    end
end

"""
    UniformDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces each species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`.
"""
function UniformDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     spp = 1:size(eco.abundances.matrix, 1)
     meanabun = mean(eco.spplist.abun)
     avlost = rate * timestep * meanabun
     for i in spp
         for j in 1:round(Int64, avlost)
             if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos)
         eco.abundances.matrix[i, smp] .-= 1
     end
        end
     end
end

"""
    ProportionalDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces each species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size.
"""
function ProportionalDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     spp = 1:size(eco.abundances.matrix, 1)
     avlost = rate * timestep .* eco.spplist.abun
     for i in spp
         if any(eco.abundances.matrix[i, :] .> 0)
             pos = find(eco.abundances.matrix[i, :] .> 0)
             smp = sample(pos, round(Int64, avlost[i]),
              replace = true)
             eco.abundances.matrix[i, smp] .-= 1
         end
     end
end

"""
    LargeDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces the largest species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size. The
    largest species are those that have greater than half of the maximum energy
    in the system.
"""
function LargeDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     !all(eco.spplist.requirement.energy .==
        maximum(eco.spplist.requirement.energy)) ||
        error("All species have the same requirement")
     largest = find(eco.spplist.requirement.energy .>
        quantile(eco.spplist.requirement.energy)[3])
     avlost = rate * timestep .* eco.spplist.abun
     for i in largest
         if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos, round(Int64, avlost[i]),
          replace = true)
         eco.abundances.matrix[i, smp] .-= 1
     end
     end
end

"""
    RareDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces the rarest species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size. The
    rarest species are those that have abundances lower than the 25th percentile.
"""
function RareDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     originalabun = eco.spplist.abun
     rarest = find(originalabun .< quantile(originalabun)[2])
     avlost = rate * timestep .* eco.spplist.abun
     for i in rarest
         if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos, round(Int64, avlost[i]),
          replace = true)
         eco.abundances.matrix[i, smp] .-= 1
     end
     end
end

"""
    CommonDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces the most common species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size. The
    rarest species are those that have abundances greater than the 75th percentile.
"""
function CommonDecline(eco::Ecosystem, timestep::Unitful.Time,
     rate::Quantity{Float64, typeof(𝐓^-1)})
     originalabun = eco.spplist.abun
     common = find(originalabun .> quantile(originalabun)[4])
     avlost = rate * timestep .* eco.spplist.abun
     for i in common
         if any(eco.abundances.matrix[i, :] .> 0)
         pos = find(eco.abundances.matrix[i, :] .> 0)
         smp = sample(pos, round(Int64, avlost[i]),
          replace = true)
         eco.abundances.matrix[i, smp] .-= 1
     end
     end
end

"""
    Invasive(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that introduces an invasive species into the ecosystem, `eco`,
    that gains abundance at each `timestep` at a particular rate, `rate` and
    reduces 5 designated 'sensitive' species by an equivalent amount.
"""

function Invasive(eco::Ecosystem, timestep::Unitful.Time,
    rate::Quantity{Float64, typeof(𝐓^-1)})
    qual = eco.spplist.native .== false
    invasive = find(qual)
    natives = find(!qual)
    invasive_abun = eco.spplist.abun[invasive]
    avgain = rate * timestep .* invasive_abun
    for i in eachindex(invasive)
        for j in 1:round(Int64, avgain[i])
            pos = find(mapslices(sum, eco.abundances.matrix[natives, :], 1) .> 0)
            smp = sample(pos)
            eco.abundances.matrix[invasive[i], smp] .+= 1
            nonzero = find(eco.abundances.matrix[natives, smp] .> 0)
            eco.abundances.matrix[sample(nonzero), smp] .-= 1
        end
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
