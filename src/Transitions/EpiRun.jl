
function _run_rule!(epi::EpiSystem, rule::Exposure, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if epi.epienv.active[loc]
        params = epi.epilist.params
        force_cats = epi.epilist.virus.force_cats
        age_cat = epi.epilist.human.human_to_force
        N = sum_pop(epi.abundances.matrix, loc)
        env_inf = virus(epi.abundances)[1, loc] /
            (N^params.freq_vs_density_env)
        force_inf = (params.age_mixing[age_cat[spp], :] ⋅ virus(epi.abundances)[force_cats, loc]) /
            (N^params.freq_vs_density_force)
        expprob = (getprob(rule)[1] * force_inf + getprob(rule)[2] * env_inf) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exposures = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newexpprob))
        human(epi.abundances)[spp, loc] -= exposures
        human(epi.abundances)[dest, loc] += exposures
    end
end

function _run_rule!(epi::EpiSystem, rule::EnvExposure, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    get_env = rule.get_env
    env_exposure = rule.env_param
    if epi.epienv.active[loc]
        params = epi.epilist.params
        force_cats = epi.epilist.virus.force_cats
        age_cat = epi.epilist.human.human_to_force
        N = sum_pop(epi.abundances.matrix, loc)
        env = @view get_env(epi.epienv.habitat)[loc]
        env_transition = env[1] * env_exposure/mean(epi.epienv.habitat.matrix)
        env_inf = virus(epi.abundances)[1, loc] /
            (N^params.freq_vs_density_env)
        force_inf = env_transition * (params.age_mixing[age_cat[spp], :] ⋅ virus(epi.abundances)[force_cats, loc]) /
            (N^params.freq_vs_density_force)
        expprob = (getprob(rule)[1] * force_inf + getprob(rule)[2] * env_inf) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exposures = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newexpprob))
        human(epi.abundances)[spp, loc] -= exposures
        human(epi.abundances)[dest, loc] += exposures
    end
end

function _run_rule!(epi::EpiSystem, rule::Infection, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if epi.epienv.active[loc]
        infprob = getprob(rule) * timestep
        newinfprob = 1.0 - exp(-infprob)
        infs = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newinfprob))
        human(epi.abundances)[spp, loc] -= infs
        human(epi.abundances)[dest, loc] += infs
    end
end

function _run_rule!(epi::EpiSystem, rule::Recovery, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if epi.epienv.active[loc]
        recprob = getprob(rule) * timestep
        newrecprob = 1.0 - exp(-recprob)
        recs = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newrecprob))
        human(epi.abundances)[spp, loc] -= recs
        human(epi.abundances)[dest, loc] += recs
    end
end

function _run_rule!(epi::EpiSystem, rule::ForceProduce, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    # Calculate effective rates
    birthrate = getprob(rule) * timestep * human(epi.abundances)[spp, loc]
    births = rand(rng, Poisson(birthrate))
    # Spread force of infection over space
    if !iszero(births)
        virusmove!(epi, spp, loc, epi.cache.virusmigration, births)
    end
end

function _run_rule!(epi::EpiSystem, rule::ViralLoad, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    params = epi.epilist.params
    loc = getlocation(rule)
    traitmatch = traitfun(epi, loc, 1)
    deathrate = getprob(rule) * timestep * traitmatch^-1
    # Convert death rate into 0 - 1 probability
    deathprob = 1.0 - exp(-deathrate)

    # Calculate how much virus degrades in the environment
    deaths = rand(rng, Binomial(virus(epi.abundances)[1, loc], deathprob))
    # Force of infection on average around half of timestep in environment
    survivalprob = exp(-deathrate/2.0)

    # So this much force of infection survives in the environment
    force_cats = epi.epilist.virus.force_cats
    env_virus = rand(rng, Binomial(sum(virus(epi.abundances)[force_cats, loc], dims=1)[1], survivalprob * params.env_virus_scale))

    # Now update virus in environment and force of infection
    virus(epi.abundances)[1, loc] += env_virus - deaths
end

function _run_rule!(epi::EpiSystem, rule::ForceDisperse)
    rng = epi.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    dist = Poisson(epi.cache.virusmigration[spp, loc])
    epi.cache.virusmigration[spp, loc] = rand(rng, dist)
end

function run_rule!(epi::EpiSystem, rule::R, timestep::Unitful.Time) where R <: AbstractStateTransition
    if typeof(rule) == Exposure
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) == EnvExposure
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) == Infection
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) == Recovery
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) == ViralLoad
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) == ForceProduce
        _run_rule!(epi, rule, timestep)
    end
end

function run_rule!(epi::EpiSystem, rule::R, timestep::Unitful.Time) where R <: AbstractPlaceTransition
    if typeof(rule) == ForceDisperse
        _run_rule!(epi, rule)
    end
end

function update_virus_cache!(epi::EpiSystem)
    force_cats = epi.epilist.virus.force_cats
    human_to_force = epi.epilist.human.human_to_force
    locs = size(virus(epi.abundances), 2)
    vm = zeros(eltype(epi.cache.virusmigration), length(force_cats), locs)
    classes = length(epi.epilist.human.names)
    Threads.@threads for i in 1:classes
        for j in 1:locs
            vm[human_to_force[i], j] += epi.cache.virusmigration[i, j]
        end
    end
    virus(epi.abundances)[force_cats, :] .= vm
end

function new_update!(epi::EpiSystem, timestep::Unitful.Time)


    Threads.@threads for st in epi.transitions.state
        run_rule!(epi, st, timestep)
    end

    Threads.@threads for pl in epi.transitions.place
        run_rule!(epi, pl, timestep)
    end

    update_virus_cache!(epi)

    # Invalidate all caches for next update
    invalidatecaches!(epi)

end


function new_simulate!(epi::E, duration::Unitful.Time, timestep::Unitful.Time) where E <: AbstractEpiSystem
  times = length(0s:timestep:duration)
  for i in 1:times
    new_update!(epi, timestep)
  end
end

function new_simulate_record!(storage::AbstractArray, epi::E,
  times::Unitful.Time, interval::Unitful.Time,timestep::Unitful.Time) where E <: AbstractEpiSystem
  ustrip(mod(interval,timestep)) == 0.0 || error("Interval must be a multiple of timestep")
  record_seq = 0s:interval:times
  time_seq = 0s:timestep:times
  storage[:, :, 1] = epi.abundances.matrix
  counting = 1
  for i in 2:length(time_seq)
    new_update!(epi, timestep);
    if time_seq[i] in record_seq
      counting = counting + 1
      storage[:, :, counting] = epi.abundances.matrix
    end
  end
  storage
end
