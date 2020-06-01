using AxisArrays
using HDF5
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units

"""
    simulate!(
        epi::AbstractEpiSystem,
        duration::Unitful.Time,
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
  # save pre-simulation inputs
  if save && !isdir(save_path) # Create the directory if it doesn't already exist.
      mkpath(save_path)
  end
  save && Simulation.save(joinpath(save_path, "initial_system.jlso"), epi)
  save && JLSO.save(
    joinpath(save_path, "configuration.jlso"),
    :duration => duration,
    :timestep => timestep,
  )

  times = length(0s:timestep:duration)

  for i in 1:times
    update!(epi, timestep)
  end

  # save simulation results
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

  # save pre-simulation inputs
  if save && !isdir(save_path) # Create the directory if it doesn't already exist.
      mkpath(save_path)
  end
  save && Simulation.save(joinpath(save_path, "initial_system.jlso"), epi)
  save && JLSO.save(
    joinpath(save_path, "configuration.jlso"),
    :storage => storage,
    :times => times,
    :timestep => timestep,
    :interval => interval,
  )

  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = epi.abundances.matrix
  counting = 1

  # - initialise and save the first timestep abuns/storage to HDF5
  # construct axes for abuns/storage matrix
  grid_id = map(Iterators.product(axisvalues(epi.epienv.initial_population)...)) do (x,y)
      return string.(Int(x), "-", Int(y))
  end
  # TODO: confirm converting `grid_id` from matrix to vector in the way below gives the
  # correct order assumed in the model
  grid_id = vec(grid_id)
  axes = (;
      compartment = epi.epilist.human.names,
      grid_id = grid_id,
      times = string.(record_seq)
  )
  if save
      initialise_output_abuns(
          storage,
          axes;
          h5fn=joinpath(save_path, "abundances.h5")
      )
      update_output_abuns(
          epi.abundances.matrix,
          counting;
          h5fn=joinpath(save_path, "abundances.h5")
      )
  end

  # - simulate each timestep
  for i in 2:length(time_seq)
    update!(epi, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = epi.abundances.matrix
      if save
          update_output_abuns(
              epi.abundances.matrix,
              counting;
              h5fn=joinpath(save_path, "abundances.h5")
          )
      end
    end
  end

  # save simulation results
  save && Simulation.save(joinpath(save_path, "final_system.jlso"), epi)
  return storage
end

"""
    initialise_output_abuns(
        abuns::Array,
        axes::NamedTuple;
        h5fn=joinpath(pwd(),"abundances.h5")
    )

Create an HDF5 file `h5fn` to store abundance. Preallocate a fix-sized array to store
abundance for each compartment, grid location and timestep. Fill in values for compartments,
grid locations and timesteps.
"""
function initialise_output_abuns(
    abuns::Array,
    axes::NamedTuple;
    h5fn=joinpath(pwd(),"abundances.h5")
)
    # - sanity check
    if size(abuns) != length.(values(axes))
        throw(DimensionMismatch(
            "abundances array size $(size(abuns)) does not match axes size " *
            "axes: $(keys(axes)). vector length: $(length.(values(axes)))"
        ))
    end

    # - initialise HDF5 file
    h5open(h5fn, "w") do fid
        # - create group
        group = g_create(fid, "abundances")
        attrs(group)["Description"] = string(
            "Contains the abundances for each compartment and geographic location ",
            "through the simulation duration."
        )
        # - add data to the group
        # initialise abuns 3D array
        dset_abuns = d_create(
            group,
            "abuns",
            datatype(Int),
            dataspace(size(abuns)),
            "chunk",
            (size.(Ref(abuns), [1,2])...,1)
        )
        # fill in axes information
        for k in keys(axes)
            group[string(k)] = axes[k]
        end
    end
end

"""
    update_output_abuns(
        abuns_t::Matrix,
        timestep::Int;
        h5fn=joinpath(pwd(),"abundances.h5")
    )

Update the existing HDF5 file `h5fn` with the abundance matrix at a certain timestep.
"""
function update_output_abuns(
    abuns_t::Matrix,
    timestep::Int;
    h5fn=joinpath(pwd(),"abundances.h5")
)
    if !(isfile(h5fn) && ishdf5(h5fn))
        throw(ArgumentError(
            "$(h5fn) does not exist or is not a valid HDF5 file"
        ))
    end
    h5open(h5fn, "cw") do fid
        fid["abundances"]["abuns"][:,:,timestep] = abuns_t
    end
end
