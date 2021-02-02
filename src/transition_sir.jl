mutable struct Exposure{U <: Unitful.Units} <: AbstractStateTransition
    species::Int64
    location::Int64
    prob::TimeUnitType{U}
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

function run_rule!(eco::Episystem, rule::Exposure{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = eco.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        expprob = getprob(rule) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exposures = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newexpprob))
        eco.abundances.matrix[spp, loc] -= exposures
        eco.abundances.matrix[spp + 1, loc] += exposures
    end
end

function run_rule!(eco::Episystem, rule::Infection{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = eco.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        infprob = getprob(rule) * timestep
        newinfprob = 1.0 - exp(-infprob)
        infs = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newinfprob))
        eco.abundances.matrix[spp, loc] -= infs
        eco.abundances.matrix[spp + 1, loc] += infs
    end
end

function run_rule!(eco::Episystem, rule::Recovery{U}, timestep::Unitful.Time) where U <: Unitful.Units
    rng = eco.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        recprob = getprob(rule) * timestep
        newrecprob = 1.0 - exp(-recprob)
        recs = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newrecprob))
        eco.abundances.matrix[spp, loc] -= recs
        eco.abundances.matrix[spp + 1, loc] += recs
    end
end


function run_rule!(eco::Episystem, rule::SEIR{U}, timestep::Unitful.Time) where U <: Unitful.Units
    run_rule!(eco, rule.exposure, timestep)
    run_rule!(eco, rule.infection, timestep)
    run_rule!(eco, rule.recovery, timestep)
end

function _run_rule!(eco::Episystem, rule::AllDisperse)
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        move!(eco, eco.spplist.movement, loc, spp, eco.cache.netmigration, eco.abundances.matrix[spp, loc])
    end
end

function run_rule!(eco::Episystem, rule::R, timestep::Unitful.Time) where R <: AbstractTransition
    if typeof(rule) <: Exposure
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) <: Infection
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) <: Recovery
        _run_rule!(eco, rule, timestep)
    elseif typeof(rule) <: AllDisperse
        _run_rule!(eco, rule)
    end
end

function create_epi_transitions(spplist::SpeciesList, abenv::GridAbioticEnv)

    state_list = [SEIR(Exposure(1, loc, spplist.params.birth[1]), Infection(2, loc, spplist.params.death[2]), Recovery(3, loc, spplist.params.death[3])) for loc in eachindex(abenv.habitat.matrix)]

    place_list = [AllDisperse(spp, loc) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]

    return TransitionList(state_list, place_list)
end

function new_update!(epi::Episystem, timestep::Unitful.Time)

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
