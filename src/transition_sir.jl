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

function run_rule!(eco::Ecosystem, rule::Exposure, timestep::TimeUnitType{U}) where U <: Unitful.Units
    rng = eco.abundances.seed[Threads.threadid()]
    spp = getspecies(rule)
    loc = getlocation(rule)
    if eco.abenv.active[loc]
        expprob = getprob(rule) * timestep
        newexpprob = 1.0 - exp(-expprob)
        exp = rand(rng, Binomial(eco.abundances.matrix[spp, loc], newexpprob))
        eco.abundances.matrix[spp, loc] -= exp
        eco.abundances.matrix[spp + 1, loc] += exp
    end
end

function run_rule!(eco::Ecosystem, rule::Infection, timestep::TimeUnitType{U}) where U <: Unitful.Units
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

function run_rule!(eco::Ecosystem, rule::Recovery, timestep::TimeUnitType{U}) where U <: Unitful.Units
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



function create_epi_transitions(spplist::SpeciesList, abenv::GridAbioticEnv)

    state_list = [SEIR(Exposure(spp, loc, spplist.params.birth[spp]), Infection(spp, loc, spplist.params.death[spp]), Recovery(spp, loc, spplist.params.death[spp])) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]

    place_list = [AllDisperse(spp, loc) for spp in eachindex(spplist.names) for loc in eachindex(abenv.habitat.matrix)]

    return TransitionList(state_list, place_list)
end
