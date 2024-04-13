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
function run_rule!(eco::Ecosystem, rule::S, timestep::Unitful.Time, specialise = false) where S <: AbstractStateTransition
    if specialise
        _run_rule!(eco, rule, timestep)
    else
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
end

function run_rule!(eco::Ecosystem, rule::P, timestep::Unitful.Time, specialise = false) where P <: AbstractPlaceTransition
    if specialise

        _run_rule!(eco, rule)
    else
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
end

function run_rule!(eco::Ecosystem, rule::SU, timestep::Unitful.Time, specialise = false) where SU <: AbstractSetUp
    if specialise
        _run_rule!(eco, rule, timestep)
    else
        if typeof(rule) == UpdateEnergy
            _run_rule!(eco, rule, timestep)
        elseif typeof(rule) == SeedInfection
            _run_rule!(eco, rule, timestep)
        else
            _run_rule!(eco, rule, timestep)
        end
    end
end

function run_rule!(eco::Ecosystem, rule::Missing, timestep::Unitful.Time, specialise = false)
    return @warn "No setup"
end

function run_rule!(eco::Ecosystem, rule::WD, timestep::Unitful.Time, specialise = false) where WD <: AbstractWindDown
    if specialise
        _run_rule!(eco, rule, timestep)
    else
        if typeof(rule) == UpdateEnvironment
            _run_rule!(eco, rule, timestep)
        elseif typeof(rule) == UpdateEpiEnvironment
            _run_rule!(eco, rule, timestep)
        else
            _run_rule!(eco, rule, timestep)
        end
    end
end

@generated function run_generated!(eco::E, rule::AbstractStateTransition, timestep::Unitful.Time) where {StateTypes <: AbstractStateTransition, L, Part, SL,
    TR, LU, C, SU, PL, WD, TL <: TransitionList{SU, StateTypes, PL, WD}, E <: Ecosystem{L, Part, SL, TR, LU, C, TL}}
    generated = quote end
    if isabstracttype(StateTypes)
        types = rsubtypes(StateTypes)
    else
        types = Base.uniontypes(StateTypes)
    end
    for type in types
        push!(generated.args, :(if typeof(rule) == $type
                                    return _run_rule!(eco, rule, timestep)
                                end)
             )
    end
    push!(generated.args, :(@error "Reached an unmatched rule type ($(typeof(rule)))"))
    return generated
end

@generated function run_generated!(eco::E, rule::AbstractSetUp, timestep::Unitful.Time) where {SetUpTypes <: AbstractSetUp, L, Part, SL,
    TR, LU, C, ST, PL, WD, TL <: TransitionList{SetUpTypes, ST, PL, WD}, E <: Ecosystem{L, Part, SL, TR, LU, C, TL}}
    generated = quote end
    if isabstracttype(SetUpTypes)
        types = rsubtypes(SetUpTypes)
    else
        types = Base.uniontypes(SetUpTypes)
    end
    for type in types
        push!(generated.args, :(if typeof(rule) == $type
                                    return _run_rule!(eco, rule, timestep)
                                end)
             )
    end
    push!(generated.args, :(@error "Reached an unmatched rule type ($(typeof(rule)))"))
    return generated
end

@generated function run_generated!(eco::E, rule::AbstractPlaceTransition) where {PlaceTypes <: AbstractPlaceTransition, L, Part, SL,
    TR, LU, C, SU, ST, WD, TL <: TransitionList{SU, ST, PlaceTypes, WD}, E <: Ecosystem{L, Part, SL, TR, LU, C, TL}}
    generated = quote end
    if isabstracttype(PlaceTypes)
        types = rsubtypes(PlaceTypes)
    else
        types = Base.uniontypes(PlaceTypes)
    end
    for type in types
        push!(generated.args, :(if typeof(rule) == $type
                                    return _run_rule!(eco, rule)
                                end)
             )
    end
    push!(generated.args, :(@error "Reached an unmatched rule type ($(typeof(rule)))"))
    return generated
end

@generated function run_generated!(eco::E, rule::AbstractWindDown, timestep::Unitful.Time) where {WindDownTypes <: AbstractWindDown, L, Part, SL,
    TR, LU, C, SU, ST, PL, TL <: TransitionList{SU, ST, PL, WindDownTypes}, E <: Ecosystem{L, Part, SL, TR, LU, C, TL}}
    generated = quote end
    if isabstracttype(WindDownTypes)
        types = rsubtypes(WindDownTypes)
    else
        types = Base.uniontypes(WindDownTypes)
    end
    for type in types
        push!(generated.args, :(if typeof(rule) == $type
                                    return _run_rule!(eco, rule, timestep)
                                end)
             )
    end
    push!(generated.args, :(@error "Reached an unmatched rule type ($(typeof(rule)))"))
    return generated
end

"""
    update!(eco::Ecosystem, timestep::Unitful.Time)

Update an Ecosystem by one timestep, running through different
transitions, including set up, state and place transitions, and
winddown.
"""
function update!(eco::Ecosystem, timestep::Unitful.Time, ::TransitionList, specialise = false)
    if specialise
        run! = _run_rule!
    else
        run! = run_generated!
    end

    Threads.@threads for su in eco.transitions.setup
        run!(eco, su, timestep)
    end

    Threads.@threads for st in eco.transitions.state
        run!(eco, st, timestep)
    end

    Threads.@threads for pl in eco.transitions.placelist
        for p in pl
            run!(eco, eco.transitions.place[p])
        end
    end

    Threads.@threads for wd in eco.transitions.winddown
        run!(eco, wd, timestep)
    end

end

"""
    update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing) where L <: EpiLandscape

Update an Ecosystem by one timestep, running through the older version
of the epidemiology code.
"""
function update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing, specialise = false) where L <: EpiLandscape
    epi_update!(eco, timestep)
end

"""
    update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing) where L <: GridLandscape

Update an Ecosystem by one timestep, running through the older version
of the biodiversity code.
"""
function update!(eco::AbstractEcosystem{L}, timestep::Unitful.Time, ::Nothing, specialise = false) where L <: GridLandscape
    biodiversity_update!(eco, timestep)
end

function generate_storage(eco::Ecosystem, times::Int64, reps::Int64)
  numSpecies = counttypes(eco.spplist)
  gridSize = _countsubcommunities(eco.abenv.habitat)
  return Array{Int64, 4}(undef, numSpecies, gridSize, times, reps)
end

function generate_storage(eco::Ecosystem, qs::Int64, times::Int64, reps::Int64)
    gridSize = _countsubcommunities(eco.abenv.habitat)
    return Array{Float64, 4}(undef, gridSize, qs, times, reps)
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
    specialise = false,
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
    update!(eco, timestep, gettransitions(eco), specialise)
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
    specialise = false,
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
    update!(eco, timestep, gettransitions(eco), specialise);
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
