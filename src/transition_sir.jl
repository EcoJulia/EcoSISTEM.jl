mutable struct ViralLoad{U <: Unitful.Units} <: AbstractStateTransition
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct Exposure{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    force_prob::TimeUnitType{U}
    virus_prob::TimeUnitType{U}
end

function get_prob(rule::Exposure)
    return rule.force_prob, rule.virus_prob
end

mutable struct Infection{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct Recovery{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

mutable struct SEIR{U <: Unitful.Units} <: AbstractStateTransition
    exposure::Exposure{U}
    infection::Infection{U}
    recovery::Recovery{U}
end

mutable struct Force{U <: Unitful.Units} <: AbstractPlaceTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
end

function get_prob(rule::Force)
    return rule.prob
end

mutable struct ForceDisperse <: AbstractPlaceTransition
    species::Int64
    location::Int64
end

function run_rule!(epi::EpiSystem, rule::Exposure{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    params = epi.epilist.params
    if epi.abenv.active[loc]
        N = sum_pop(epi.abundances.matrix, loc)
        env_inf = virus(epi.abundances)[1, loc] /
            (N^params.freq_vs_density_env)
        force_inf = virus(epi.abundances)[spp, loc] /
            (N^params.freq_vs_density_force)
        expprob = (getprob(rule)[1] * force_inf + getprob(rule)[2] * env_inf) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exposures = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newexpprob))
        human(epi.abundances)[spp, loc] -= exposures
        human(epi.abundances)[spp + 1, loc] += exposures
    end
end

function run_rule!(epi::EpiSystem, rule::Infection{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if epi.abenv.active[loc]
        infprob = getprob(rule) * timestep
        newinfprob = 1.0 - exp(-infprob)
        infs = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newinfprob))
        human(epi.abundances)[spp, loc] -= infs
        human(epi.abundances)[spp + 1, loc] += infs
    end
end

function run_rule!(epi::EpiSystem, rule::Recovery{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = epi.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if epi.abenv.active[loc]
        recprob = getprob(rule) * timestep
        newrecprob = 1.0 - exp(-recprob)
        recs = rand(rng, Binomial(epi.abundances.matrix[spp, loc], newrecprob))
        human(epi.abundances)[spp, loc] -= recs
        human(epi.abundances)[spp + 1, loc] += recs
    end
end


function run_rule!(epi::EpiSystem, rule::SEIR{U}, timestep::Unitful.Time) where U <: Unitful.Units
    run_rule!(epi, rule.exposure, timestep)
    run_rule!(epi, rule.infection, timestep)
    run_rule!(epi, rule.recovery, timestep)
end

function _run_rule!(epi::EpiSystem, rule::Force{U}, timestep::Unitful.Time) where U <: Unitful.Units
    spp = getspecies(rule)
    loc = getlocation(rule)
    # Calculate effective rates
    birthrate = get_prob(rule) * timestep * human(epi.abundances)[spp, loc]
    births = rand(rng, Poisson(birthrate))
    # Spread force of infection over space
    if !iszero(births)
        virusmove!(epi, spp, loc, epi.cache.virusmigration, births)
    end
end

function _run_rule!(epi::EpiSystem, rule::ViralLoad{U}, timestep::Unitful.Time) where U <: Unitful.Units
    params = epi.epilist.params
    loc = getlocation(rule)
    traitmatch = traitfun(epi, loc, 1)
    deathrate = get_prob(rule) * timestep * traitmatch^-1
    # Convert death rate into 0 - 1 probability
    deathprob = 1.0 - exp(-deathrate)

    # Calculate how much virus degrades in the environment
    deaths = rand(rng, Binomial(virus(epi.abundances)[1, loc], deathprob))

    # Force of infection on average around half of timestep in environment
    survivalprob = exp(-deathrate/2.0)

    # So this much force of infection survives in the environment
    env_virus = rand(rng, Binomial(virus(epi.abundances)[1, loc], survivalprob * params.env_virus_scale))

    # Now update virus in environment and force of infection
    virus(epi.abundances)[1, loc] += env_virus - deaths
end

function _run_rule!(epi::EpiSystem, rule::ForceDisperse)
    spp = getspecies(rule)
    loc = getlocation(rule)
    dist = Poisson(epi.cache.virusmigration[spp, loc])
    epi.cache.virusmigration[spp, loc] = rand(rng, dist)
end

function run_rule!(epi::EpiSystem, rule::R, timestep::Unitful.Time) where R <: AbstractTransition
    if typeof(rule) <: Exposure
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: Infection
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: Recovery
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: ForceDisperse
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: ViralLoad
        _run_rule!(epi, rule, timestep)
    elseif typeof(rule) <: Force
        _run_rule!(epi, rule)
    end
end

function create_transitions(epilist::EpiList, epienv::GridEpiEnv)
    params = epilist.params

    state_list = [SEIR(Exposure(1, loc, params.transition_force[1, 2], params.transition_virus[1, 2]), Infection(2, loc, params.transition[2, 3]), Recovery(3, loc, params.transition[3, 4])) for loc in eachindex(epienv.habitat.matrix)]

    virus_list = [ViralLoad(loc, params.virus_decay) for loc in eachindex(epienv.habitat.matrix)]

    state_list = [virus_list; state_list]

    place_list = [ForceDisperse(spp, loc) for spp in eachindex(epilist.virus.names) for loc in eachindex(epienv.habitat.matrix)]

    return TransitionList(state_list, place_list)
end

function new_update!(epi::EpiSystem, timestep::Unitful.Time)

    Threads.@threads for st in epi.transitions.state
        run_rule!(epi, st, timestep)
    end
    Threads.@threads for pl in epi.transitions.place
        run_rule!(epi, pl, timestep)
    end

    epi.abundances.matrix .+= epi.cache.netmigration

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
