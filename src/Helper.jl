using ProgressMeter

function run_sim(eco, params::AbstractVector, times::Int64, reps::Int64)

  birth = param[1]
  death = param[2]
  timestep = param[3]
  l = param[4]
  s = param[5]

  gridSize = grid[1] *  grid[2]
  abun = zeros(times+1, numSpecies, gridSize, reps); ener = zeros(times+1, gridSize, reps)

  for j in 1:reps

    repopulate!(eco, false)

    for k in 1:gridSize
      abun[1,:, k,  j] = eco.abundances.matrix[:, 1]
      ener[1, k,  j] = sum(eco.spplist.abun .* eco.spplist.energy.energy)
    end

    for i in 1:times
        update!(eco, birth, death, move, l, s, timestep)
        for g in 1:gridSize
          abun[i+1, :, g, j] = eco.abundances.matrix[: , g]
          ener[i+1, g, j] = sum(eco.abundances.matrix[: , g] .* eco.spplist.energy.energy)
        end
    end
    map
  end
  mean_abun = mapslices(mean, abun, 4)
  sd_abun = mapslices(std, abun, 4)
  mean_ener = mapslices(mean, ener, 3)
  sd_ener = mapslices(std, ener, 3)
  [abun, mean_abun, sd_abun, mean_ener, sd_ener]
end

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
  ener = zeros(length(time_seq)+1, reps)

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
        update_fun(eco, birth, death, l, s, timestep); #print(eco.abundances)
        if any(i.==time_seq)
            counting = counting+1
            abun[counting, :, j, :] = eco.abundances.matrix
      end
    end
  end
  abun
end

function expected_counts(grd::Array{Int64, 3}, sq::Int64)

  total = mapslices(sum, grd , length(size(grd)))[:, :,  1]
  grd = grd[:, :, sq]

  grd = grd[total.>0]
  total = total[total.>0]

  actual = counts(grd+1, maximum(grd+1))
  actual = convert(Array{Float64,1}, actual)

  # Calculate
  expected_dist = zeros(Float64, (length(total), maximum(total)+1))
  for i in 1:length(total)
    expected_dist[i, 1:(total[i]+1)] = repmat([1/(total[i]+1)], total[i]+1)
  end
  expected = mapslices(sum, expected_dist, 1)

  # Cut expected values to length of actual
  expected = expected[1:length(actual)]

  return [expected, actual]
end

function expected_counts(grd::Array{Int64, 4}, sq::Int64)

  total = mapslices(sum, grd , length(size(grd)))[:, :, :,  1]
  grd = grd[:, :, :, sq]

  grd = grd[total.>0]
  total = total[total.>0]

  actual = counts(grd+1, maximum(grd+1))
  actual = convert(Array{Float64,1}, actual)

  # Calculate
  expected_dist = zeros(Float64, (length(total), maximum(total)+1))
  for i in 1:length(total)
    expected_dist[i, 1:(total[i]+1)] = repmat([1/(total[i]+1)], total[i]+1)
  end
  expected = mapslices(sum, expected_dist, 1)

  # Cut expected values to length of actual
  expected = expected[1:length(actual)]

  return [expected, actual]
end
