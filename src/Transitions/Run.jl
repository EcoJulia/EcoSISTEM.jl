using AxisArrays
using HDF5
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

import HDF5: ishdf5

function run_rule!(eco::Ecosystem, rule::R, timestep::Unitful.Time) where R <: AbstractStateTransition
    if typeof(rule) == BirthProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == DeathProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Exposure
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == EnvExposure
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Infection
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Recovery
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == ViralLoad
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == ForceProduce
        _run_rule!(eco, rule, timestep)
    end
end

function run_rule!(eco::Ecosystem, rule::R, timestep::Unitful.Time) where R <: AbstractPlaceTransition
    if typeof(rule) == AllDisperse
        _run_rule!(eco, rule)
    elseif typeof(rule) == ForceDisperse
        _run_rule!(eco, rule)
    end
end

function run_rule!(eco::Ecosystem, rule::R, timestep::Unitful.Time) where R <: AbstractSetUp
    if typeof(rule) == UpdateEnergy
        _run_rule!(eco, rule, timestep)
    end
end

function run_rule!(eco::Ecosystem, rule::Missing, timestep::Unitful.Time)
    return @warn "No setup"
end

function run_rule!(eco::Ecosystem, rule::R, timestep::Unitful.Time) where R <: AbstractWindDown
    if typeof(rule) == UpdateEnvironment
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == UpdateEpiEnvironment
        _run_rule!(eco, rule, timestep)
    end
end

function update!(eco::Ecosystem, timestep::Unitful.Time)

    run_rule!(eco, eco.transitions.setup, timestep)

    Threads.@threads for st in eco.transitions.state
        run_rule!(eco, st, timestep)
    end

    Threads.@threads for pl in eco.transitions.place
        run_rule!(eco, pl, timestep)
    end

    run_rule!(eco, eco.transitions.winddown, timestep)

end

"""
    function new_simulate!(
        eco::AbstractEcosystem,
        duration::Unitful.Time,
        timestep::Unitful.Time;
        save=false,
        save_path=pwd(),
    )

Run an ecosystem, `eco`, for specified length of times, `duration`, for a
particular timestep, `timestep`. If `save=true`, inputs and outputs are saved as JLSO files
at `save_path`.
"""
function new_simulate!(
    eco::AbstractEcosystem,
    duration::Unitful.Time,
    timestep::Unitful.Time;
    save=false,
    save_path=pwd(),
)
  # save pre-simulation inputs
  if save && !isdir(save_path) # Create the directory if it doesn't already exist.
      mkpath(save_path)
  end
  save && EcoSISTEM.save(joinpath(save_path, "initial_system.jlso"), eco)
  save && JLSO.save(
    joinpath(save_path, "configuration.jlso"),
    :duration => duration,
    :timestep => timestep,
  )

  times = length(0s:timestep:duration)

  for i in 1:times
    update!(eco, timestep)
  end

  # save simulation results
  save && EcoSISTEM.save(joinpath(save_path, "final_system.jlso"), eco)
end

"""
    new_simulate_record!(
        storage::AbstractArray,
        eco::Ecosystem,
        times::Unitful.Time,
        interval::Unitful.Time,
        timestep::Unitful.Time;
        save=false,
        save_path=pwd(),
    )

Run an ecosystem, `ecoepi`, for specified length of times, `times`,
for a particular timestep, `timestep`, and time interval for abundances to be
recorded, `interval`. Optionally, there may also be a scenario by which the
whole ecosystem is updated, such as removal of habitat patches.
"""
function new_simulate_record!(
    storage::AbstractArray,
    eco::Ecosystem,
    times::Unitful.Time,
    interval::Unitful.Time,
    timestep::Unitful.Time;
    save=false,
    save_path=pwd(),
)
  mod(interval,timestep) == zero(mod(interval,timestep)) ||
    error("Interval must be a multiple of timestep")

  # save pre-simulation inputs
  if save && !isdir(save_path) # Create the directory if it doesn't already exist.
      mkpath(save_path)
  end
  save && EcoSISTEM.save(joinpath(save_path, "initial_system.jlso"), eco)
  save && JLSO.save(
    joinpath(save_path, "configuration.jlso"),
    :storage => storage,
    :times => times,
    :timestep => timestep,
    :interval => interval,
  )

  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = eco.abundances.matrix
  counting = 1

  # - initialise and save the first timestep abuns/storage to HDF5
  # construct axes for abuns/storage matrix
  ax = AxisArrays.axes(eco.abundances.grid)[end-1:end]
  grid_id = map(Iterators.product(ax...)) do (x,y)
      return string.(x, "-", y)
  end
  # TODO: confirm converting `grid_id` from matrix to vector in the way below gives the
  # correct order assumed in the model
  grid_id = vec(grid_id)
  axes = (;
      compartment = getnames(eco.spplist),
      grid_id = grid_id,
      times = string.(uconvert.(day, 1.0 .* record_seq))
  )
  if save
      initialise_output_abuns(
          storage,
          axes;
          h5fn=joinpath(save_path, "abundances.h5")
      )
      update_output_abuns(
          eco.abundances.matrix,
          counting;
          h5fn=joinpath(save_path, "abundances.h5")
      )
  end

  # - simulate each timestep
  for i in 2:length(time_seq)
    update!(eco, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = eco.abundances.matrix
      if save
          update_output_abuns(
              eco.abundances.matrix,
              counting;
              h5fn=joinpath(save_path, "abundances.h5")
          )
      end
    end
  end

  # save simulation results
  save && EcoSISTEM.save(joinpath(save_path, "final_system.jlso"), eco)
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
        group = create_group(fid, "abundances")
        attributes(group)["Description"] = string(
            "Contains the abundances for each compartment and geographic location ",
            "through the simulation duration."
        )
        # - add data to the group
        # initialise abuns 3D array
        dset_abuns = create_dataset(
            group,
            "abuns",
            datatype(eltype(abuns)),
            dataspace(size(abuns)),
            chunk=(size.(Ref(abuns), [1,2])...,1)
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
    abuns_t::Matrix{U},
    timestep::Int;
    h5fn=joinpath(pwd(),"abundances.h5")
) where U <: Integer
    if !(isfile(h5fn) && ishdf5(h5fn))
        throw(ArgumentError(
            "$(h5fn) does not exist or is not a valid HDF5 file"
        ))
    end
    h5open(h5fn, "cw") do fid
        fid["abundances"]["abuns"][:,:,timestep] = abuns_t
    end
end
