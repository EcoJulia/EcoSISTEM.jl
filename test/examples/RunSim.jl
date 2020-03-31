using Simulation
using Unitful
using Simulation.Units
using JLD
using Printf

import Simulation.runscenario!
function simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  scenario::SC, divfuns::Array{Function, 1}, q::Float64, cacheInterval::Unitful.Time, cacheFolder::String, scenario_name::String, rep::Int64) where SC <: Simulation.AbstractScenario
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
      update!(eco, timestep);
      runscenario!(eco, timestep, scenario, time_seq[i]);
      # Record diversity profiles for each measure
      if time_seq[i] in record_seq
          counting = counting + 1
          for j in eachindex(divfuns)
              storage[j, counting] = divfuns[j](eco, q)[!, :diversity][1]
          end
      end
      # Save cache of abundances
      if mod(time_seq[i], cacheInterval) == 0year
          JLD.save(joinpath(cacheFolder, scenario_name * (@sprintf "%02d" uconvert(NoUnits,time_seq[i]/cacheInterval)) * (@sprintf "%03d.jld" rep)), "abun", eco.abundances.matrix)
      end
  end
  storage
end

function runsim!(div::Array{Float64, 3}, abenv::AB, paramDict::Dict, simDict::Dict, rep::Int64, folder::String = "./", recreate = true) where AB <: Simulation.AbstractAbiotic
    cacheFolder = joinpath(folder, "cache")
    isdir(cacheFolder) || mkdir(cacheFolder)
    eco = create_eco(paramDict, abenv, bound = paramDict["bound"])
    totalK = (sum(eco.abenv.budget.b1.matrix), sum(eco.abenv.budget.b2.matrix))
    scenario = simDict["scenarios"]
    scenario_names = simDict["scenario_names"]
    divfuns = simDict["divfuns"]
    q = simDict["q"]
    for i in 1:length(scenario)
        if recreate
            recreate_eco!(eco, totalK, sum(eco.spplist.abun))
        end
        thisstore = view(div, :, :, i)
        simulate!(eco, simDict["burnin"], simDict["timestep"])
        simulate_record_diversity!(thisstore, eco, simDict["times"], simDict["interval"], simDict["timestep"],
        scenario[i], divfuns, q, simDict["cacheInterval"], cacheFolder, scenario_names[i], rep)
        adjusttemp!(eco, scenario[i], simDict["times"])
    end
    JLD.save((folder * @sprintf "%03d.jld" rep), "div", div)
    print(rep, "\n")
end

function adjusttemp!(eco::Ecosystem, scenario::Union{SimpleScenario, FluctScenario}, times::Unitful.Time)
    if scenario.fun == TempIncrease! && scenario.rate > 0.0K/year
        resetrate!(eco, 0.0K/year)
        eco.abenv.habitat.matrix .-= (scenario.rate * times)
    end
end

function adjusttemp!(eco::Ecosystem, scenario::MultiScenario, times::Unitful.Time)
    if scenario.sc1.fun == TempIncrease! && scenario.sc1.rate > 0.0K/year
        resetrate!(eco, 0.0K/year)
        eco.abenv.habitat.matrix .-= (scenario.sc1.rate * times)
    end
end

function dispersalrun!(div::Array{Float64, 3}, abenv::AB, paramDict::Dict, simDict::Dict, rep::Int64, folder::String = "./") where AB <: Simulation.AbstractAbiotic
    cacheFolder = joinpath(folder, "cache")
    isdir(cacheFolder) || mkdir(cacheFolder)
    eco = create_eco(paramDict, abenv, bound = paramDict["bound"])
    totalK = (sum(eco.abenv.budget.b1.matrix), sum(eco.abenv.budget.b2.matrix))
    scenario = simDict["scenarios"]
    scenario_names = simDict["scenario_names"]
    divfuns = simDict["divfuns"]
    q = simDict["q"]
    eco.abundances.matrix .= 0
    eco.abundances.grid[1, :, 1] .= rand(Multinomial(eco.spplist.abun[1], 10))
    eco.abundances.grid[2, :, end] .= rand(Multinomial(eco.spplist.abun[2], 10))
    for i in 1:length(scenario)
        thisstore = view(div, :, :, i)
        simulate!(eco, simDict["burnin"], simDict["timestep"])
        simulate_record_diversity!(thisstore, eco, simDict["times"], simDict["interval"], simDict["timestep"],
        scenario[i], divfuns, q, simDict["cacheInterval"], cacheFolder, scenario_names[i], rep)
        adjusttemp!(eco, scenario[i], simDict["times"])
    end
    JLD.save((folder * @sprintf "%03d.jld" rep), "div", div)
    print(rep, "\n")
end
