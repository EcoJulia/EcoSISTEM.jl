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
    elseif typeof(rule) == Missing
        _run_rule!(eco, rule, timestep)
    end
end

function run_rule!(eco::Ecosystem, rule::R, timestep::Unitful.Time) where R <: AbstractWindDown
    if typeof(rule) == UpdateEnvironment
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) == UpdateEpiEnvironment
        _run_rule!(eco, rule, timestep)
    end
end

function new_update!(eco::Ecosystem, timestep::Unitful.Time)

    run_rule!(eco, eco.transitions.setup, timestep)

    Threads.@threads for st in eco.transitions.state
        run_rule!(eco, st, timestep)
    end

    Threads.@threads for pl in eco.transitions.place
        run_rule!(eco, pl, timestep)
    end

    run_rule!(eco, eco.transitions.winddown, timestep)

end

function new_simulate!(eco::E, duration::Unitful.Time, timestep::Unitful.Time) where E <: AbstractEcosystem
  times = length(0s:timestep:duration)
  for i in 1:times
    new_update!(eco, timestep)
  end
end

function new_simulate_record!(storage::AbstractArray, eco::E,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time) where E <: AbstractEcosystem
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = eco.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    new_update!(eco, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = eco.abundances.matrix
    end
  end
  storage
end
