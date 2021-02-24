function _run_rule!(eco::Ecosystem, rule::BirthProcess{U}, timestep) where U <: Unitful.Units
    rng = eco.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, loc, spp)
        birthprob = getprob(rule) * timestep * adjusted_birth
        newbirthprob = 1.0 - exp(-birthprob)
        births = rand(rng, Poisson(eco.abundances.matrix[spp, loc] * newbirthprob))
        eco.abundances.matrix[spp, loc] += births
    end
end

function _run_rule!(eco::Ecosystem, rule::DeathProcess{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = eco.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, loc, spp)
        deathprob = getprob(rule) * timestep * adjusted_death
        newdeathprob = 1.0 - exp(-deathprob)
        deaths = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newdeathprob))
        eco.abundances.matrix[spp, loc] -= deaths
    end
end

function _run_rule!(eco::Ecosystem, rule::BirthDeathProcess{U}, timestep::Unitful.Time) where U <: Unitful.Units
    _run_rule!(eco, rule.birth, timestep)
    _run_rule!(eco, rule.death, timestep)
end

function _run_rule!(eco::Ecosystem, rule::AllDisperse)
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        move!(eco, eco.spplist.movement, loc, spp, eco.cache.netmigration, eco.abundances.matrix[spp, loc])
    end
end

function run_rule!(eco::Ecosystem, rule::R, timestep::Unitful.Time) where R <: AbstractTransition
    if typeof(rule) <: BirthProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) <: DeathProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) <: BirthDeathProcess
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) <: AllDisperse
        _run_rule!(eco, rule)
    end
end

function new_update!(eco::Ecosystem, timestep::Unitful.Time)

    # Set the overall energy budget of that square
    update_energy_usage!(eco)

    Threads.@threads for st in eco.transitions.state
        run_rule!(eco, st, timestep)
    end
    Threads.@threads for pl in eco.transitions.place
        run_rule!(eco, pl, timestep)
    end

    eco.abundances.matrix .+= eco.cache.netmigration

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
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
