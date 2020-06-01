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

  for i in 2:length(time_seq)
    update!(epi, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = epi.abundances.matrix
    end
  end

  # save simulation results
  save && Simulation.save(joinpath(save_path, "final_system.jlso"), epi)
  save && JLSO.save(joinpath(save_path, "abundances.jlso"), :abundances => storage)
  return storage
end

function initialise_output_abuns(
    abuns::Array,
    axes::NamedTuple;
    h5fn=joinpath(pwd(),"abundances.h5")
)
    # - sanity check
    if size(abuns) != length.(values(axes))
        throw(DimensionMismatch(
            "abundances array size $(size(abuns)) does not match axes size ", *
            "axes: $(keys(axes)). vector length: $(length.(values(axes)))"
        ))
    end

    # - initialise HDF5 file
    h5open(h5fn, "w") do fid
        # - create group
        group = g_create(fid, "abundances")
        attrs(group)["Description"] = String(
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

function update_output_abuns(
    abuns::Array,
    axes::NamedTuple,
    timestep::Int;
    h5fn=joinpath(pwd(),"abundances.h5")
)
    h5open(h5fn, "cw") do fid
        fid["abundances"]["abuns"][:,:,timestep] = epi.abundances.matrix
    end
end



# TEMP
using AxisArrays
using HDF5


grid_id = map(Iterators.product(axisvalues(epi.epienv.initial_population)...)) do (x,y)
    return string.(Int(x), "-", Int(y))
end
grid_id = vec(grid_id) # TODO: construct `grid-id` correctly


axes = (;
    compartment = epi.epilist.human.names,
    grid_id = grid_id,
    times = string.(0s:interval:times) # TODO, use `record_seq`, better to convert `s` to `day`
)

# before iteration


# for each iteration
fid = h5open("test.h5", "cw")
fid["abundances"]["abuns"][:,:,1] = epi.abundances.matrix # TODO 1 to `counting`
close(fid)

# - to check what's saved
f = h5open("test.h5", "r")
abuns_data = read(f["abundances"]["abuns"])
close(f)
