using ProgressMeter

function run_sim(eco, params::AbstractVector, times::Int64, reps::Int64)

  birth = param[1]
  death = param[2]
  move = param[3]
  timestep = param[4]
  l = param[5]
  s = param[6]

  gridSize = grid[1] *  grid[2]
  abun = zeros(times+1, numSpecies, gridSize, reps); ener = zeros(times+1, gridSize, reps)

  for j in 1:reps

    populate!(eco, false)

    for k in 1:gridSize
      abun[1,:, k,  j] = eco.partition.abundances[:, 1]
      ener[1, k,  j] = sum(eco.spplist.abun .* eco.spplist.energy.energy)
    end

    for i in 1:times
        update!(eco, birth, death, move, l, s, timestep)
        for g in 1:gridSize
          abun[i+1, :,g, j] = eco.partition.abundances[: , g]
          ener[i+1, g, j] = sum(eco.partition.abundances[: , g] .* eco.spplist.energy.energy)
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

function run_sim_spatial(eco::Ecosystem, params::AbstractVector,
   times::Int64, burnin::Int64, interval::Int64, reps::Int64, birth_move::Bool)

  birth = param[1]
  death = param[2]
  move = param[3]
  timestep = param[4]
  l = param[5]
  s = param[6]
  time_seq = collect(burnin:interval:times)
  gridSize = grid[1] *  grid[2]
  abun = zeros(length(time_seq), numSpecies, reps, gridSize); ener = zeros(length(time_seq)+1, reps)

  if birth_move
    update_fun=update_birth_move!
  else
    update_fun=update!
  end

  @showprogress 1 "Computing..." for j in 1:reps
    populate!(eco, false)

    abun[1, :, j, :] = eco.partition.abundances[:, :]
    counting = 1
    for i in 1:times
        update_fun(eco, birth, death, move, l, s, timestep); #print(eco.partition.abundances)
        if any(i.==time_seq)
            counting = counting+1
            abun[counting, :, j, :] = eco.partition.abundances[: , :]
      end
    end
  end
  abun
end
