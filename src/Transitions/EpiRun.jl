
"""
    _run_rule!(eco::Ecosystem, rule::Exposure, timestep::Unitful.Time)

Stochastic exposure process for a location, housed inside `rule`,
 for one timestep. Moves from susceptible category, `species` to
exposed category, `destination`.
"""
function _run_rule!(eco::Ecosystem, rule::Exposure, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        params = eco.spplist.params
        force_cats = eco.spplist.pathogens.force_cats
        age_cat = eco.spplist.species.human_to_force
        N = sum_pop(eco.abundances.matrix, loc)
        env_inf = virus(eco.abundances)[1, loc] /
            (N^params.freq_vs_density_env)
        force_inf = (params.age_mixing[age_cat[spp], :] ⋅ virus(eco.abundances)[force_cats, loc]) /
            (N^params.freq_vs_density_force)
        expprob = (getprob(rule)[1] * force_inf + getprob(rule)[2] * env_inf) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exposures = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newexpprob))
        human(eco.abundances)[spp, loc] -= exposures
        human(eco.abundances)[dest, loc] += exposures
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::EnvExposure, timestep::Unitful.Time)

Stochastic environmentally-driven exposure process for a location, housed inside `rule`,
 for one timestep. Moves from susceptible category, `species` to
exposed category, `destination`.
"""
function _run_rule!(eco::Ecosystem, rule::EnvExposure, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    get_env = rule.get_env
    env_exposure = rule.env_param
    if eco.abenv.active[loc]
        params = eco.spplist.params
        force_cats = eco.spplist.pathogens.force_cats
        age_cat = eco.spplist.species.human_to_force
        N = sum_pop(eco.abundances.matrix, loc)
        env = @view get_env(eco.abenv.habitat)[loc]
        env_transition = env[1] * env_exposure/mean(eco.abenv.habitat.matrix)
        env_inf = virus(eco.abundances)[1, loc] /
            (N^params.freq_vs_density_env)
        force_inf = env_transition * (params.age_mixing[age_cat[spp], :] ⋅ virus(eco.abundances)[force_cats, loc]) /
            (N^params.freq_vs_density_force)
        expprob = (getprob(rule)[1] * force_inf + getprob(rule)[2] * env_inf) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exposures = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newexpprob))
        human(eco.abundances)[spp, loc] -= exposures
        human(eco.abundances)[dest, loc] += exposures
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::Union{Infection, DevelopSymptoms, Hospitalise,
        DeathFromInfection, Recovery}, timestep::Unitful.Time)

Stochastic epi transition process for a location, housed inside `rule`,
 for one timestep. Moves from source category, `species` to
destination category, `destination`.
"""
function _run_rule!(eco::Ecosystem, rule::Union{Infection, DevelopSymptoms, Hospitalise,
    DeathFromInfection, Recovery}, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        prob = getprob(rule) * timestep
        newprob = 1.0 - exp(-prob)
        changes = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newprob))
        human(eco.abundances)[spp, loc] -= changes
        human(eco.abundances)[dest, loc] += changes
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::ForceProduce, timestep::Unitful.Time)

Generation of force of infection for a location.
"""
function _run_rule!(eco::Ecosystem, rule::ForceProduce, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    # Calculate effective rates
    birthrate = getprob(rule) * timestep * human(eco.abundances)[spp, loc]
    births = rand(rng, Poisson(birthrate))
    # Spread force of infection over space
    if !iszero(births)
        virusmove!(eco, spp, loc, eco.cache.virusmigration, births)
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::ViralLoad, timestep::Unitful.Time)

Settling of force of infection into environmental reservoir for a location.
"""
function _run_rule!(eco::Ecosystem, rule::ViralLoad, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    params = eco.spplist.params
    loc = getlocation(rule)
    traitmatch = traitfun(eco, loc, 1)
    deathrate = getprob(rule) * timestep * traitmatch^-1
    # Convert death rate into 0 - 1 probability
    deathprob = 1.0 - exp(-deathrate)

    # Calculate how much virus degrades in the environment
    deaths = rand(rng, Binomial(virus(eco.abundances)[1, loc], deathprob))
    # Force of infection on average around half of timestep in environment
    survivalprob = exp(-deathrate/2.0)

    # So this much force of infection survives in the environment
    force_cats = eco.spplist.pathogens.force_cats
    env_virus = rand(rng, Binomial(sum(virus(eco.abundances)[force_cats, loc], dims=1)[1], survivalprob * params.env_virus_scale))

    # Now update virus in environment and force of infection
    virus(eco.abundances)[1, loc] += env_virus - deaths
end

"""
    _run_rule!(eco::Ecosystem, rule::ForceDisperse)

Random dispersal of force of infection from a location.
"""
function _run_rule!(eco::Ecosystem, rule::ForceDisperse)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    dist = Poisson(eco.cache.virusmigration[spp, loc])
    eco.cache.virusmigration[spp, loc] = rand(rng, dist)
end

"""
    _run_rule!(eco::Ecosystem, rule::SeedInfection)

Seed initial infected individuals.
"""
function _run_rule!(eco::Ecosystem, rule::SeedInfection, timestep::Unitful.Time)
    rule.update_fun(eco, eco.abenv.control, timestep)
end

"""
    _run_rule!(eco::Ecosystem, rule::UpdateEpiEnvironment)

Update the ecosystem caches.
"""
function _run_rule!(eco::Ecosystem, rule::UpdateEpiEnvironment, timestep::Unitful.Time)
    rule.update_fun(eco, timestep)
end
