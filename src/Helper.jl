using Unitful
using Unitful.DefaultSymbols
using Diversity
using MyUnitful
"""
    simulate!(eco::Ecosystem, duration::Unitful.Time, interval::Unitful.Time,
         timestep::Unitful.Time)

Function to run an ecosystem, `eco` for specified length of times, `duration`,
for a particular timestep, 'timestep'.
"""
function simulate!(eco::Ecosystem, duration::Unitful.Time, timestep::Unitful.Time)
  times = length(0s:timestep:duration)
  for i in 1:times
    update!(eco, timestep)
  end
end

function generate_storage(eco::Ecosystem, times::Int64, reps::Int64)
  numSpecies = length(eco.spplist.abun)
  gridSize = _countsubcommunities(eco.abenv.habitat)
  abun = Array{Int64, 4}(numSpecies, gridSize, times, reps)
end
function generate_storage(eco::Ecosystem, qs::Int64, times::Int64, reps::Int64)
  gridSize = _countsubcommunities(eco.abenv.habitat)
  abun = Array{Float64, 4}(gridSize, qs, times, reps)
end
"""
    simulate!(eco::Ecosystem, duration::Unitful.Time, interval::Unitful.Time,
         timestep::Unitful.Time)

Function to run an ecosystem, `eco` for specified length of times, `duration`,
for a particular timestep, 'timestep', and time interval for abundances to be
recorded, `interval`. Optionally, there may also be a scenario by which the
whole ecosystem is updated, such as removal of habitat patches.
"""
function simulate_record!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = eco.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    update!(eco, timestep);
    if any(time_seq[i].==record_seq)
      counting = counting + 1
      storage[:, :, counting] = eco.abundances.matrix
    end
  end
  storage
end

function simulate_record!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time, timestep::Unitful.Time,
  scenario::AbstractScenario)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = eco.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    update!(eco, timestep);
    runscenario!(eco, timestep, scenario, time_seq[i]);
    if any(time_seq[i].==record_seq)
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
function simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  scenario::SimpleScenario, divfun::Function, qs::Vector{Float64})
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
    update!(eco, timestep);
    runscenario!(eco, timestep, scenario, time_seq[i]);
    if any(time_seq[i].==record_seq)
      counting = counting + 1
      storage[:, :, counting] = reshape(divfun(eco, qs)[:diversity],
      countsubcommunities(eco), length(qs))
    end
  end
  storage
end

function simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  divfun::Function, qs::Vector{Float64})
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
    update!(eco, timestep);
    if any(time_seq[i].==record_seq)
      counting = counting + 1
      diversity = divfun(eco, qs)[:diversity]
      storage[:, :, counting] = reshape(diversity,
      Int(length(diversity)/ length(qs)), length(qs))
    end
  end
  storage
end
function simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  divfuns::Array{Function, 1}, q::Float64)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
    update!(eco, timestep);
    if any(time_seq[i].==record_seq)
      counting = counting + 1
      for j in eachindex(divfuns)
          storage[:, j, counting] = divfuns[j](eco, q)[:diversity][1]
      end
    end
  end
  storage
end
function simulate_record_diversity!(storage::AbstractArray, eco::Ecosystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time,
  scenario::SimpleScenario, divfuns::Array{Function, 1}, q::Float64)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
      update!(eco, timestep);
      runscenario!(eco, timestep, scenario, time_seq[i]);
      if any(time_seq[i].==record_seq)
          counting = counting + 1
          for j in eachindex(divfuns)
              storage[:, j, counting] = divfuns[j](eco, q)[:diversity][1]
          end
      end
  end
  storage
end

function cleanup!(abun::Array{Int64, 4})
    zeroabun = mapslices(x -> all(x.!=0), abun, [1, 2, 4])[1, 1, :, 1]
    abun = abun[:, :, zeroabun, :]
end

function cleanup!(div::Array{Float64, 4})
    nadiv = mapslices(x -> all(!isnan.(x)), div, [1, 2, 4])[1, 1, :, 1]
    div = div[:, :, nadiv, :]
end


function expected_counts(grd::Array{Float64, 3}, sq::Int64)
  grd = convert(Array{Int64}, grd)
  total = mapslices(sum, grd , length(size(grd)))[:, :,  1]
  grd = grd[:, :, sq]
  _expected_counts(total, grd, sq)
end


function expected_counts(grd::Array{Float64, 4}, sq::Int64)
  grd = convert(Array{Int64}, grd)
  total = mapslices(sum, grd , length(size(grd)))[:, :, :,  1]
  grd = grd[:, :, :, sq]
  _expected_counts(total, grd, sq)
end

function _expected_counts(total::Array{Int64}, grd::Array{Int64}, sq::Int64)
  grd = grd[reshape(total, size(grd)).>0]
  total = total[total.>0]

  actual = counts(grd+1, maximum(grd+1))
  actual = convert(Array{Float64,1}, actual)

  expected_dist = zeros(Float64, (length(total), maximum(total)+1))
  for i in 1:length(total)
    expected_dist[i, 1:(total[i]+1)] = repmat([1/(total[i]+1)], total[i]+1)
  end
  expected = mapslices(sum, expected_dist, 1)

  # Cut expected values to length of actual
  expected = expected[1:length(actual)]

  return [expected, actual]
end


function expected_counts(grd::Array{Float64}, sq::Int64, spp::Int64)
  spp_grd = grd[:, spp, :, :]
  expected_counts(spp_grd, sq)
end
