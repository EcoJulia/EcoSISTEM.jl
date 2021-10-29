"""
    _run_rule!(eco::Ecosystem, rule::BirthProcess, timestep::Unitful.Time)

Stochastic birth process for a location and species, house inside `rule`,
 for one timestep.
"""
function _run_rule!(eco::Ecosystem, rule::BirthProcess, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if (eco.abenv.active[loc]) & (eco.cache.totalE[loc, 1] > 0)
        adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, loc, spp)
        birthrate = getprob(rule) * timestep * adjusted_birth |> NoUnits
        births = rand(rng, Poisson(eco.abundances.matrix[spp, loc] * birthrate))
        eco.abundances.matrix[dest, loc] += births
    end
end
"""
    _run_rule!(eco::Ecosystem, rule::GenerateSeed, timestep::Unitful.Time)

Stochastic seeding process for a location and species, house inside `rule`,
 for one timestep.
"""
function _run_rule!(eco::Ecosystem, rule::GenerateSeed, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    dest = getdestination(rule)
    loc = getlocation(rule)
    if (eco.abenv.active[loc]) && (eco.cache.totalE[loc, 1] > 0)
        adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, loc, spp)
        birthrate = getprob(rule) * timestep * adjusted_birth |> NoUnits
        births = rand(rng, Poisson(eco.abundances.matrix[spp, loc] * birthrate))
        eco.cache.seedbank[dest, loc] = births
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::DeathProcess, timestep::Unitful.Time)

Stochastic death process for a location and species, house inside `rule`,
 for one timestep.
"""
function _run_rule!(eco::Ecosystem, rule::DeathProcess, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if (eco.abenv.active[loc]) && (eco.cache.totalE[loc, 1] > 0)
        adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, loc, spp)
        deathprob = getprob(rule) * timestep * adjusted_death
        newdeathprob = 1.0 - exp(-deathprob)
        deaths = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newdeathprob))
        eco.abundances.matrix[spp, loc] -= deaths
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::AllDisperse, timestep::Unitful.Time)

Stochastic dispersal process across the ecosystem for a species from a location, house inside `rule`,
 for one timestep.
"""
function _run_rule!(eco::Ecosystem, rule::AllDisperse)
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        move!(eco, eco.spplist.species.movement, loc, spp, eco.cache.netmigration, eco.abundances.matrix[spp, loc])
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::SeedDisperse, timestep::Unitful.Time)

Stochastic seed dispersal process across the ecosystem for a species from a location, house inside `rule`,
 for one timestep.
"""
function _run_rule!(eco::Ecosystem, rule::SeedDisperse)
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        move!(eco, eco.spplist.species.movement, loc, spp, eco.cache.netmigration, eco.cache.seedbank[spp, loc])
    end
end

"""
    _run_rule!(eco::Ecosystem, rule::UpdateEnergy, timestep::Unitful.Time)

Calculate energy usage across the Ecosystem for one timestep.
"""
function _run_rule!(eco::Ecosystem, rule::UpdateEnergy, timestep::Unitful.Time)
    rule.update_fun(eco)
end

"""
    _run_rule!(eco::Ecosystem, rule::UpdateEnvironment, timestep::Unitful.Time)

Update habitat and resource budget across the Ecosystem for one timestep.
"""
function _run_rule!(eco::Ecosystem, rule::UpdateEnvironment, timestep::Unitful.Time)
    rule.update_fun(eco, timestep)
end
