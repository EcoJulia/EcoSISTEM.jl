
using Unitful
using Unitful.DefaultSymbols
using MyUnitful

"""
    AbstractScenario

Abstract supertype for all whole ecosystem change scenarios
"""
abstract type AbstractScenario end

RateType = typeof(1.0/year)

if VERSION >= v"0.7"
    """
        SimpleScenario <: AbstractScenario

    This scenario type holds a function that acts to change the entire ecosystem.
    """
    mutable struct SimpleScenario <: AbstractScenario
        fun::Function
        rate::Union{Quantity{Float64, 𝐓^-1}, Quantity{Float64, 𝚯*𝐓^-1}, Quantity{Float64, 𝐋*𝐓^-1}}
    end
else
    """
        SimpleScenario <: AbstractScenario

    This scenario type holds a function that acts to change the entire ecosystem.
    """
    mutable struct SimpleScenario <: AbstractScenario
        fun::Function
        rate::Union{Quantity{Float64, typeof(𝐓^-1)}, Quantity{Float64, typeof(𝚯*𝐓^-1)}, Quantity{Float64, typeof(𝐋*𝐓^-1)}}
    end
end
GLOBAL_typedict["SimpleScenario"] = SimpleScenario


function TempIncrease(eco::Ecosystem, timestep::Unitful.Time, rate::typeof(1.0K/year))
    resetrate!(eco, rate)
end
"""
    RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that randomly removes a portion of habitat at a certain rate.
"""
function RandHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    v = uconvert(unit(timestep)^-1, rate)
    pos = findall(eco.abenv.budget.matrix .> 0.0)
    howmany = jbinom(1, length(pos), ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] .= 0.0
    eco.abundances.grid[:, smp] .= 0.0
end
GLOBAL_funcdict["RandHabitatLoss!"] = RandHabitatLoss!
"""
    ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that removes a portion of habitat at a certain rate, which spreads through
joined pieces of habitat.
"""
function ClustHabitatLoss!(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    v = uconvert(unit(timestep)^-1, rate)
    if all(eco.abenv.budget.matrix .> 0.0)
        pos = findall(eco.abenv.active)
        howmany = jbinom(1, length(pos), ustrip(v))[1]
        smp = sample(pos, howmany)
        eco.abenv.budget.matrix[smp] .= 0.0
        eco.abundances.grid[:, smp] .= 0.0
    else
        pos = findall(eco.abenv.budget.matrix .== 0.0)
        howmany = jbinom(1, length(findall(eco.abenv.active)), ustrip(v))[1]
        width = size(eco.abenv.budget.matrix,1)
        x, y = [ x[1] for x in pos ], [ x[2] for x in pos ]
        neighbours = get_neighbours(eco.abenv.budget.matrix, x, y, 8)
        smp = sample(1:size(neighbours,1), howmany)
        i = convert_coords(neighbours[smp, 1], neighbours[smp, 2], width)
        eco.abenv.budget.matrix[i] .= 0.0
        eco.abundances.matrix[:, i] .= 0.0
    end
    # Add in additional start points
    howmany = jbinom(1, 1, ustrip(v))[1]
    smp = sample(pos, howmany)
    eco.abenv.budget.matrix[smp] .= 0.0
    eco.abundances.grid[:, smp] .= 0.0
end
GLOBAL_funcdict["ClustHabitatLoss!"] = ClustHabitatLoss!

"""
     HabitatReplacement(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that randomly replaces a portion of habitat with another.
"""
function HabitatReplacement(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    v = uconvert(unit(timestep)^-1, rate)
    pos = length(eco.abenv.budget.matrix)
    howmany = jbinom(1, pos, ustrip(v))[1]
    smp = sample(1:pos, howmany)
    eco.abenv.habitat.matrix[smp] = maximum(eco.spplist.traits.val) + 1
end

GLOBAL_funcdict["HabitatReplacement"] = HabitatReplacement

function SusceptibleDecline(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
     eco.spplist.susceptible *= (1 + rate * timestep)
     currentabun = mapslices(sum, eco.abundances.matrix, dims = 2)
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
GLOBAL_funcdict["SusceptibleDecline"] = SusceptibleDecline
"""
    UniformDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces each species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`.
"""
function UniformDecline(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
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
GLOBAL_funcdict["UniformDecline"] = UniformDecline
"""
    ProportionalDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces each species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size.
"""
function ProportionalDecline(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
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
GLOBAL_funcdict["ProportionalDecline"] = ProportionalDecline
"""
    LargeDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces the largest species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size. The
    largest species are those that have greater than half of the maximum energy
    in the system.
"""
function LargeDecline(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
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
GLOBAL_funcdict["LargeDecline"] = LargeDecline
"""
    RareDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces the rarest species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size. The
    rarest species are those that have abundances lower than the 25th percentile.
"""
function RareDecline(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
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
GLOBAL_funcdict["RareDecline"] = RareDecline
"""
    CommonDecline(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that reduces the most common species in the Ecosystem, `eco`, per `timestep`,
    at a particular `rate`, proportional to their starting population size. The
    rarest species are those that have abundances greater than the 75th percentile.
"""
function CommonDecline(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
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
GLOBAL_funcdict["CommonDecline"] = CommonDecline
"""
    Invasive(eco::Ecosystem, timestep::Unitful.Time,
        rate::Quantity{Float64, typeof(𝐓^-1)})

A function that introduces an invasive species into the ecosystem, `eco`,
    that gains abundance at each `timestep` at a particular rate, `rate` and
    reduces 5 designated 'sensitive' species by an equivalent amount.
"""
function Invasive(eco::Ecosystem, timestep::Unitful.Time, rate::RateType)
    qual = eco.spplist.native .== false
    invasive = findall(qual)
    natives = findall(.!qual)
    invasive_abun = eco.spplist.abun[invasive]
    avgain = rate * timestep .* invasive_abun
    for i in eachindex(invasive)
        gains = round(Int64, avgain[i]) * eco.spplist.requirement.energy[invasive[i]]
        avail_space = mapslices(sum, eco.abundances.matrix .* eco.spplist.requirement.energy, dims = 1)
        pos = find(avail_space .> gains)
        smp = sample(pos)
        eco.abundances.matrix[invasive[i], smp] .+= 1
        while gains .> 0
            spp = sample(find(eco.abundances.matrix[:, smp].>0))
            eco.abundances.matrix[spp, smp] .-= 1
            gains -= eco.spplist.requirement.energy[spp]
        end
    end
end
GLOBAL_funcdict["Invasive"] = Invasive

function runscenario!(eco::Ecosystem, timestep::Unitful.Time, scenario::SimpleScenario, currentstep::Unitful.Time)
    scenario.fun(eco, timestep, scenario.rate)
end
