using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

"""
    simulate!(
        epi::AbstractEpiSystem,
        duration::Unitful.Time,
        interval::Unitful.Time,
        timestep::Unitful.Time;
        save=false,
        save_path=pwd(),
    )

Run an epidemiological system, `epi`, for specified length of times, `duration`, for a
particular timestep, `timestep`. If `save=true`, inputs and outputs are saved as JLSO files
at `save_path`.
"""
function simulate!(
    epi::AbstractEpiSystem,
    duration::Unitful.Time,
    timestep::Unitful.Time;
    save=false,
    save_path=pwd(),
)
  times = length(0s:timestep:duration)
  if save && !isdir(save_path) # Create the directory if it doesn't already exist.
      mkdir(save_path)
  end
  save && Simulation.save(joinpath(save_path, "initial_system.jlso"), epi)
  save && JLSO.save(
    joinpath(save_path, "configuration.jlso"),
    :duration => duration,
    :timestep => timestep,
  )
  for i in 1:times
    update!(epi, timestep)
  end
  save && Simulation.save(joinpath(save_path, "final_system.jlso"), epi)
end

"""
    simulate_record!(
        storage::AbstractArray,
        epi::EpiSystem,
        times::Unitful.Time,
        interval::Unitful.Time,
        timestep::Unitful.Time;
        save=false,
        save_path=pwd(),
    )

Run an epidemiological system, `epi`, for specified length of times, `times`,
for a particular timestep, `timestep`, and time interval for abundances to be
recorded, `interval`. Optionally, there may also be a scenario by which the
whole ecosystem is updated, such as removal of habitat patches.
"""
function simulate_record!(
    storage::AbstractArray,
    epi::EpiSystem,
    times::Unitful.Time,
    interval::Unitful.Time,
    timestep::Unitful.Time;
    save=false,
    save_path=pwd(),
)
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = epi.abundances.matrix
  counting = 1
  if save && !isdir(save_path) # Create the directory if it doesn't already exist.
      mkdir(save_path)
  end
  save && Simulation.save(joinpath(save_path, "initial_system.jlso"), epi)
  save && JLSO.save(
    joinpath(save_path, "configuration.jlso"),
    :times => times,
    :timestep => timestep,
    :interval => interval,
  )
  for i in 2:length(time_seq)
    update!(epi, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = epi.abundances.matrix
    end
    print(".")
  end
  save && Simulation.save(joinpath(save_path, "final_system.jlso"), epi)
  save && JLSO.save(joinpath(save_path, "abundances.jlso"), :abundances => storage)
  return storage
end
