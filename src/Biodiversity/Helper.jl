using Unitful
using Unitful.DefaultSymbols
using Diversity
using EcoSISTEM.Units
import Diversity.Gamma
using Compat

"""
    simulate!(eco::Ecosystem, duration::Unitful.Time, interval::Unitful.Time,
         timestep::Unitful.Time)

Function to run an ecosystem, `eco` for specified length of times, `duration`,
for a particular timestep, 'timestep'.
"""
function biodiversity_simulate!(eco::AbstractEcosystem, duration::Unitful.Time, timestep::Unitful.Time)
  times = length(0s:timestep:duration)
  for i in 1:times
    biodiversity_update!(eco, timestep)
  end
end
function biodiversity_simulate!(cache::CachedEcosystem, srt::Unitful.Time, timestep::Unitful.Time)
  eco = Ecosystem{typeof(cache.abenv), typeof(cache.spplist),
  typeof(cache.relationship)}(copy(cache.abundances.matrix[srt]),
   cache.spplist, cache.abenv,
  cache.ordinariness, cache.relationship, cache.lookup, cache.cache)
  biodiversity_update!(eco, timestep)
  cache.abundances.matrix[srt + timestep] = eco.abundances
end

function generate_storage(eco::Ecosystem, times::Int64, reps::Int64)
  numSpecies = length(eco.spplist.species.abun)
  gridSize = _countsubcommunities(eco.abenv.habitat)
  abun = Array{Int64, 4}(Compat.undef, numSpecies, gridSize, times, reps)
end
function generate_storage(eco::Ecosystem, qs::Int64, times::Int64, reps::Int64)
  gridSize = _countsubcommunities(eco.abenv.habitat)
  abun = Array{Float64, 4}(Compat.undef, gridSize, qs, times, reps)
end

"""
    simulate_record!(eco::Ecosystem, duration::Unitful.Time, interval::Unitful.Time,
         timestep::Unitful.Time)

Function to run an ecosystem, `eco` for specified length of times, `duration`,
for a particular timestep, 'timestep', and time interval for abundances to be
recorded, `interval`. Optionally, there may also be a scenario by which the
whole ecosystem is updated, such as removal of habitat patches.
"""
function biodiversity_simulate_record!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = eco.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    biodiversity_update!(eco, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = eco.abundances.matrix
    end
  end
  storage
end

function biodiversity_simulate_record!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time,
  scenario::AbstractScenario)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = eco.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    biodiversity_update!(eco, timestep);
    runscenario!(eco, timestep, scenario, time_seq[i]);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = eco.abundances.matrix
    end
  end
  storage
end

"""
    simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
      times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
      scenario::SimpleScenario, divfun::Function, qs::Float64)

Function to run an ecosystem, `eco` for specified length of times, `duration`,
for a particular timestep, 'timestep', and time interval for a diversity to be
calculated and recorded, `interval`. Optionally, there may also be a scenario by which the
whole ecosystem is updated, such as removal of habitat patches.
"""
function biodiversity_simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  scenario::SimpleScenario, divfun::F, qs::Vector{Float64}) where {F<:Function}
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
    biodiversity_update!(eco, timestep);
    runscenario!(eco, timestep, scenario, time_seq[i]);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = reshape(divfun(eco, qs)[!, :diversity],
      countsubcommunities(eco), length(qs))
    end
  end
  storage
end

function biodiversity_simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  divfun::F, qs::Vector{Float64}) where {F<:Function}
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = ifelse(iseven(size(storage, 3)),timestep:interval:times, 0s:interval:times)
  time_seq = ifelse(iseven(size(storage, 3)), timestep:timestep:times, 0s:timestep:times)
  counting = 0
  for i in 1:length(time_seq)
    biodiversity_update!(eco, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      diversity = divfun(eco, qs)[!, :diversity]
      storage[:, :, counting] = reshape(diversity,
      Int(length(diversity)/ length(qs)), length(qs))
    end
  end
  storage
end
function biodiversity_simulate_record_diversity!(storage::AbstractArray,
    storage2::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  qs::Vector{Float64})
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = ifelse(iseven(size(storage, 3)),timestep:interval:times, 0s:interval:times)
  time_seq = ifelse(iseven(size(storage, 3)), timestep:timestep:times, 0s:timestep:times)
  counting = 0
  for i in 1:length(time_seq)
    biodiversity_update!(eco, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      measures = [NormalisedAlpha, NormalisedBeta, Gamma]
      for k in 1:3
      dm = measures[k](eco)
      diversity = subdiv(dm, qs)[!, :diversity]
      diversity2 = metadiv(dm, qs)[!, :diversity]
      storage[:, :, k, counting] = reshape(diversity,
      Int(length(diversity)/ length(qs)), length(qs))
      storage2[:, k, counting] = diversity2
    end
    end
  end
  return storage, storage2
end
function biodiversity_simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  divfuns::Array{Function, 1}, q::Float64)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
    biodiversity_update!(eco, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      for j in eachindex(divfuns)
          storage[:, j, counting] .= divfuns[j](eco, q)[!, :diversity][1]
      end
    end
  end
  storage
end
function biodiversity_simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  scenario::SimpleScenario, divfuns::Array{Function, 1}, q::Float64)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
      biodiversity_update!(eco, timestep);
      runscenario!(eco, timestep, scenario, time_seq[i]);
      if time_seq[i] in record_seq
          counting = counting + 1
          for j in eachindex(divfuns)
              storage[:, j, counting] .= divfuns[j](eco, q)[!, :diversity][1]
          end
      end
  end
  storage
end
