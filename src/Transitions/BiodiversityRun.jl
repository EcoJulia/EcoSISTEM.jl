function _run_rule!(eco::Ecosystem, rule::BirthProcess, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
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

function _run_rule!(eco::Ecosystem, rule::DeathProcess, timestep::Unitful.Time)
    rng = eco.abundances.rngs[Threads.threadid()]
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

function _run_rule!(eco::Ecosystem, rule::AllDisperse)
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        move!(eco, eco.spplist.species.movement, loc, spp, eco.cache.netmigration, eco.abundances.matrix[spp, loc])
    end
end

function _run_rule!(eco::Ecosystem, rule::UpdateEnergy, timestep::Unitful.Time)
    rule.update_fun(eco)
end

function _run_rule!(eco::Ecosystem, rule::UpdateEnvironment, timestep::Unitful.Time)
    rule.update_fun(eco, timestep)
end
