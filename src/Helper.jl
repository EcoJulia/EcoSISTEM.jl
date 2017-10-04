using Unitful
using Unitful.DefaultSymbols
"""
run_sim_spatial(eco::Ecosystem, param::AbstractVector,
   times::Int64, burnin::Int64, interval::Int64, reps::Int64, birth_move::Bool)

Function to run an ecosystem, `eco`, through a simulation for a set of parameters,
`param`, specified number of times, `times` and certain number of repetitions,
`reps`, with a burnin period, `burnin`, time interval for abundances to be recorded,
`interval`. There is also the option to run migration over all abundances or only
those in the birth pulse, `birth_move`.
"""
function simulate!(eco::Ecosystem, duration::Unitful.Time,
  interval::Unitful.Time, timestep::Unitful.Time)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
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
