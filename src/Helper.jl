using ProgressMeter

"""
run_sim_spatial(eco::Ecosystem, param::AbstractVector,
   times::Int64, burnin::Int64, interval::Int64, reps::Int64, birth_move::Bool)

Function to run an ecosystem, `eco`, through a simulation for a set of parameters,
`param`, specified number of times, `times` and certain number of repetitions,
`reps`, with a burnin period, `burnin`, time interval for abundances to be recorded,
`interval`. There is also the option to run migration over all abundances or only
those in the birth pulse, `birth_move`.
"""
function run_sim_spatial(eco::Ecosystem, param::AbstractVector,
   times::Int64, burnin::Int64, interval::Int64, reps::Int64, birth_move::Bool)

  birth = param[1]
  death = param[2]
  timestep = param[3]
  l = param[4]
  s = param[5]
  numSpecies = length(eco.spplist.abun)
  time_seq = collect(burnin:interval:times)
  gridSize = length(eco.abenv.habitat.matrix)
  abun = zeros(length(time_seq)+1, numSpecies, reps, gridSize);
  #ener = zeros(length(time_seq)+1, reps)

  if birth_move
    update_fun=update_birth_move!
  else
    update_fun=update!
  end

  @showprogress 1 "Computing..." for j in 1:reps
    repopulate!(eco, false)

    abun[1, :, j, :] = eco.abundances.matrix
    counting = 1
    for i in 1:times
        update_fun(eco, birth, death, l, s, timestep);
        if any(i.==time_seq)
            counting = counting+1
            abun[counting, :, j, :] = eco.abundances.matrix
      end
    end
  end
  abun
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
