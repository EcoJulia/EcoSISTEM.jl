using AxisArrays
using HDF5
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Printf

import HDF5: ishdf5

"""
    run_rule!(eco::Ecosystem, rule::AbstractTransition, timestep::Unitful.Time)

Implement `_run_rule!` function for a particular rule type, `R`, for one timestep.
"""
function run_rule!(eco::Ecosystem, rule::AbstractStateTransition, timestep::Unitful.Time)
    if typeof(rule) == BirthProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == GenerateSeed
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == DeathProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Exposure
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == EnvExposure
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Infection
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == DevelopSymptoms
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Hospitalise
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == DeathFromInfection
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == Recovery
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == ViralLoad
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == ForceProduce
        _run_rule!(eco, rule, timestep)
    else
        _run_rule!(eco, rule, timestep)
    end
end

function run_rule!(eco::Ecosystem, rule::AbstractPlaceTransition, timestep::Unitful.Time)
    if typeof(rule) == AllDisperse
        _run_rule!(eco, rule)
    elseif typeof(rule) == SeedDisperse
        _run_rule!(eco, rule)
    elseif typeof(rule) == ForceDisperse
        _run_rule!(eco, rule)
    else
        _run_rule!(eco, rule)
    end
end

function run_rule!(eco::Ecosystem, rule::AbstractSetUp, timestep::Unitful.Time)
    if typeof(rule) == UpdateEnergy
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == SeedInfection
        _run_rule!(eco, rule, timestep)
    else
        _run_rule!(eco, rule, timestep)
    end
end

function run_rule!(eco::Ecosystem, rule::Missing, timestep::Unitful.Time)
    return @warn "No setup"
end

function run_rule!(eco::Ecosystem, rule::AbstractWindDown, timestep::Unitful.Time)
    if typeof(rule) == UpdateEnvironment
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == UpdateEpiEnvironment
        _run_rule!(eco, rule, timestep)
    else
        _run_rule!(eco, rule, timestep)
    end
end


function getsetup(rule::R) where R <: AbstractSetUp
    return rule
end
function getstate(rule::R) where R <: AbstractStateTransition
    return rule
end
function getplace(rule::R) where R <: AbstractPlaceTransition
    return rule
end
function getwinddown(rule::R) where R <: AbstractWindDown
    return rule
end

"""
    update!(eco::Ecosystem, timestep::Unitful.Time)

Update an Ecosystem by one timestep, running through different
transitions, including set up, state and place transitions, and
winddown.
"""
function update!(eco::Ecosystem, timestep::Unitful.Time, ::TransitionList)

    Threads.@threads for su in eachindex(eco.transitions.setup)
        su = getsetup(getindex(eco.transitions.setup, su))
        run_rule!(eco, su, timestep)
    end

    Threads.@threads for st in eachindex(eco.transitions.state)
        st = getstate(getindex(eco.transitions.state, st))
        run_rule!(eco, st, timestep)
    end

    Threads.@threads for pl in eachindex(eco.transitions.place)
        pl = getplace(getindex(eco.transitions.place, pl))
        run_rule!(eco, pl, timestep)
    end

    Threads.@threads for wd in eachindex(eco.transitions.winddown)
        wd = getwinddown(getindex(eco.transitions.winddown, wd))
        run_rule!(eco, wd, timestep)
    end

end

"""
    update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing) where L <: EpiLandscape

Update an Ecosystem by one timestep, running through the older version
of the epidemiology code.
"""
function update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing) where L <: EpiLandscape
    epi_update!(eco, timestep)
end

"""
    update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing) where L <: GridLandscape

Update an Ecosystem by one timestep, running through the older version
of the biodiversity code.
"""
function update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing) where L <: GridLandscape
    biodiversity_update!(eco, timestep)
end

function generate_storage(eco::Ecosystem, times::Int64, reps::Int64)
  numSpecies = counttypes(eco.spplist)
  gridSize = _countsubcommunities(eco.abenv.habitat)
  abun = Array{Int64, 4}(undef, numSpecies, gridSize, times, reps)
end

"""
    function simulate!(
        eco::Ecosystem,
        duration::Unitful.Time,
        timestep::Unitful.Time;
        save=false,
        save_path=pwd(),
    )

Run an ecosystem, `eco`, for specified length of times, `duration`, for a
particular timestep, `timestep`. If `save=true`, inputs and outputs are saved as JLSO files
at `save_path`.
"""
function simulate!(
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
    update!(eco, timestep, gettransitions(eco))
  end

  # save simulation results
  save && EcoSISTEM.save(joinpath(save_path, "final_system.jlso"), eco)
end

"""
    simulate_record!(
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
function simulate_record!(
    storage::AbstractArray,
    eco::AbstractEcosystem,
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
      initialise_output(
          storage,
          axes;
          h5fn=joinpath(save_path, "abundances.h5")
      )
      update_output(
          eco.abundances.matrix,
          counting;
          h5fn=joinpath(save_path, "abundances.h5")
      )
  end

  # - simulate each timestep
  for i in 2:length(time_seq)
    update!(eco, timestep, gettransitions(eco));
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = eco.abundances.matrix
      if save
          update_output(
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
    initialise_output(
        abuns::Array,
        axes::NamedTuple;
        h5fn=joinpath(pwd(),"abundances.h5")
    )

Create an HDF5 file `h5fn` to store abundance. Preallocate a fix-sized array to store
abundance for each compartment, grid location and timestep. Fill in values for compartments,
grid locations and timesteps.
"""
function initialise_output(
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
    update_output(
        abuns_t::Matrix,
        timestep::Int;
        h5fn=joinpath(pwd(),"abundances.h5")
    )

Update the existing HDF5 file `h5fn` with the abundance matrix at a certain timestep.
"""
function update_output(
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

function simulate!(eco::AbstractEcosystem, times::Unitful.Time, timestep::Unitful.Time,
  cacheInterval::Unitful.Time, cacheFolder::String, scenario_name::String)
  time_seq = 0s:timestep:times
  for i in 1:length(time_seq)
      update!(eco, timestep, gettransitions(eco));
      # Save cache of abundances
      if mod(time_seq[i], cacheInterval) == 0year
          @save joinpath(cacheFolder, scenario_name * (@sprintf "%02d.jld" uconvert(NoUnits,time_seq[i]/cacheInterval)))  abun=eco.abundances.matrix
      end
  end
end
