using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

function simulate!(epi::AbstractEpiSystem, duration::Unitful.Time, timestep::Unitful.Time)
  times = length(0s:timestep:duration)
  for i in 1:times
    update!(epi, timestep)
  end
end

function simulate_record!(storage::AbstractArray, epi::EpiSystem,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = epi.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    update!(epi, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = epi.abundances.matrix
    end
    print(".")
  end
  storage
end
