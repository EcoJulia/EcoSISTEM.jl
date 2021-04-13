using StatsBase
using Compat
using LinearAlgebra

"""
    update!(epi::Ecosystem, time::Unitful.Time)
Function to update disease and virus class abundances and environment for one timestep.
"""
function epi_update!(epi::Ecosystem, timestep::Unitful.Time)

    # Seed initial infecteds
    seedinfected!(epi, epi.abenv.control, timestep)

    # Virus movement loop
    virusupdate!(epi, timestep)

    # Birth/death/infection/recovery loop of each class
    classupdate!(epi, timestep)

    # Invalidate all caches for next update
    invalidatecaches!(epi)

    # Update environment - habitat and energy budgets
    habitatupdate!(epi, timestep)
    applycontrols!(epi, timestep)
end

function seedinfected!(epi::Ecosystem, controls::NoControl, timestep::Unitful.Time)
    return controls
end

function seedinfected!(epi::Ecosystem, controls::Lockdown, timestep::Unitful.Time)
    rng = epi.abundances.rngs[Threads.threadid()]
    if (epi.cache.initial_infected > 0) && (controls.current_date < controls.lockdown_date)
        inf = rand(rng, Poisson(epi.cache.initial_infected * timestep /controls.lockdown_date))
        sus_id = sample(epi.spplist.human.susceptible, inf)
        exp_id = sus_id .+ length(epi.spplist.human.susceptible)
        summed_exp = sum(human(epi.abundances)[exp_id, :], dims = 1)[1, :]
        pos = sample(rng, epi.cache.ordered_active, weights(summed_exp[epi.cache.ordered_active]), inf)
        for i in 1:inf
            if (human(epi.abundances)[sus_id[i], pos[i]] > 0)
                human(epi.abundances)[sus_id[i], pos[i]] -= 1
                human(epi.abundances)[exp_id[i], pos[i]] += 1
            end
        end
    elseif controls.current_date == controls.lockdown_date
        @info "Lockdown initiated - $(sum(human(epi.abundances)[epi.spplist.human.susceptible .+ length(epi.spplist.human.susceptible), :])) individuals infected"
    end
    return controls
end

"""
    virusupdate!(epi::Ecosystem, time::Unitful.Time)
Function to update virus abundances and disperse for one timestep.
"""
function virusupdate!(epi::Ecosystem, timestep::Unitful.Time)
    dims = _countsubcommunities(epi.abenv.habitat)
    width = getdimension(epi)[1]
    params = epi.spplist.params
    rng = epi.abundances.rngs[Threads.threadid()]
    classes = findall((params.virus_growth .* timestep) .> 0)

    # Convert 1D dimension to 2D coordinates with convert_coords(epi, j, width)
    # Check which grid cells are active, only iterate along those
    activejindices = findall(j -> epi.abenv.active[convert_coords(epi, j, width)...], 1:dims)
    # Loop through grid squares
    function firstloop(i)
        for j in activejindices
            # Calculate how much birth and death should be adjusted

            # Calculate effective rates
            birthrate = params.virus_growth[i] * timestep * human(epi.abundances)[i, j]
            births = rand(rng, Poisson(birthrate))

            # Spread force of infection over space
            if !iszero(births)
                virusmove!(epi, i, j, epi.cache.virusmigration, births)
            end
        end
    end

    firstlooptasks = Dict{Int, Task}()
    # spawn tasks on threads through the first loop
    for i in classes
        firstlooptasks[i] = Threads.@spawn firstloop(i)
    end

    force_cats = epi.spplist.virus.force_cats
    human_to_force = epi.spplist.human.human_to_force
    ages = length(unique(human_to_force))

    function secondloop(j)
        vm = zeros(eltype(epi.cache.virusmigration), ages)
        for i in classes
            # after wait virusmigration[i, j] will be up to date
            haskey(firstlooptasks, i) && Threads.wait(firstlooptasks[i])
            iszero(epi.cache.virusmigration[i, j]) && continue
            dist = Poisson(epi.cache.virusmigration[i, j])
            epi.cache.virusmigration[i, j] = rand(rng, dist)
            # Put force of infection from this class into right group
            vm[human_to_force[i]] += epi.cache.virusmigration[i, j]
        end

        # viral species 1 is assumed to be the environmental virus
        traitmatch = traitfun(epi, j, 1)
        deathrate = params.virus_decay * timestep * traitmatch^-1
        # Convert death rate into 0 - 1 probability
        deathprob = 1.0 - exp(-deathrate)

        # Calculate how much virus degrades in the environment
        deaths = rand(rng, Binomial(virus(epi.abundances)[1, j], deathprob))

        # Force of infection on average around half of timestep in environment
        survivalprob = exp(-deathrate/2.0)

        # So this much force of infection survives in the environment
        env_virus = rand(rng, Binomial(Int(sum(vm)), survivalprob * params.env_virus_scale))

        # Now update virus in environment and force of infection
        virus(epi.abundances)[1, j] += env_virus - deaths
        virus(epi.abundances)[force_cats, j] .= vm
    end


    jindices = 1:size(epi.cache.virusmigration, 2)
    threadedjindices = [jindices[j:Threads.nthreads():end]
                        for j in 1:Threads.nthreads()]
    @assert all(sort(unique(vcat(threadedjindices...))) .== jindices)
    Threads.@threads for jrange in threadedjindices
       for j in jrange
           secondloop(j)
       end
    end
end

"""
    sum_pop(m::Matrix{Int64}, i::Int64)
Function to sum a population matrix, `m`, without memory allocation, at a grid location `i`.
"""
function sum_pop(m::Matrix{R}, i::Int64) where R <: Real
    n = zero(R)
    @inbounds for j in 1:size(m, 1)
        n += m[j, i]
    end
    return n
end

function sum_pop(v::Vector{R}) where R <: Real
    n = zero(R)
    @inbounds for i in 1:length(v)
        n += v[i]
    end
    return n
end


"""
    classupdate!(epi::Ecosystem, timestep::Unitful.Time)
Function to update disease class abundances for one timestep.
"""
function classupdate!(epi::Ecosystem, timestep::Unitful.Time)
    # Calculate dimenions of habitat and number of classes
    dims = _countsubcommunities(epi.abenv.habitat)
    params = epi.spplist.params
    width = getdimension(epi)[1]
    classes = 1:size(human(epi.abundances), 1)

    human_to_force = epi.spplist.human.human_to_force
    force_cats = epi.spplist.virus.force_cats

    # Convert 1D dimension to 2D coordinates with convert_coords(epi, j, width)
    # Check which grid cells are active, only iterate along those
    activejindices = findall(j->epi.abenv.active[convert_coords(epi, j, width)...], 1:dims)
    threadedjindices = [activejindices[j:Threads.nthreads():end] for j in 1:Threads.nthreads()]
    @assert all(sort(unique(vcat(threadedjindices...))) .== activejindices)
    Threads.@threads for jrange in threadedjindices; for j in jrange
        # will result +=/-= 0 at end of inner loop, so safe to skip
        iszero(sum_pop(human(epi.abundances), j)) && continue

        rng = epi.abundances.rngs[Threads.threadid()]
        N = sum_pop(epi.abundances.matrix, j)
        # Loop through classes in chosen square
        for i in classes
            # Births
            births = rand(rng, Binomial(human(epi.abundances)[i, j],
                                        params.births[i] * timestep))
            human(epi.abundances)[1, j] += births

            # Note - transpose transition matrices makes iteration over k faster
            # Calculate force of inf and env inf
            for k in 1:size(params.transition_virus, 1)
                # Skip if there are no people in k at location j
                iszero(human(epi.abundances)[k, j]) && continue
                # Skip if there are no transitions from k to i
                params.transition[i, k] +
                    params.transition_virus[i, k] +
                    params.transition_force[i, k] > zero(inv(timestep)) ||
                    continue

                k_age_cat = human_to_force[k]

                # Environmental infection rate from k to i
                env_inf = virus(epi.abundances)[1, j] /
                    (N^params.freq_vs_density_env)

                # Direct transmission infection rate from k to i
                force_inf = (params.age_mixing[k_age_cat, :] ⋅
                             virus(epi.abundances)[force_cats, j]) /
                    (N^params.freq_vs_density_force)

                # Add to baseline transitional probabilities from k to i
                trans_val = params.transition[i, k] +
                    params.transition_virus[i, k] * env_inf +
                    params.transition_force[i, k] * force_inf
                trans_prob = 1.0 - exp(-trans_val * timestep)

                # Skip if probability is zero
                iszero(trans_prob) && continue

                # Make transitions
                trans = rand(rng,
                             Binomial(human(epi.abundances)[k, j], trans_prob))
                human(epi.abundances)[i, j] += trans
                human(epi.abundances)[k, j] -= trans
            end
        end
    end; end

end

"""
    populate!(ml::EpiLandscape, spplist::EpiList, abenv::EE, rel::R)
Function to populate an EpiLandscape with information on each disease class in the EpiList.
"""
function populate!(ml::EpiLandscape, spplist::EpiList, abenv::EE, rel::R) where {EE <: AbstractEpiEnv, R <: AbstractTraitRelationship}
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]

    rng = ml.rngs[Threads.threadid()]
    # Loop through classes
    for i in eachindex(spplist.human.abun)
        rand!(rng, Multinomial(spplist.human.abun[i], len), (@view ml.matrix[i, :]))
    end
    for i in eachindex(spplist.virus.abun)
        rand!(rng, Multinomial(spplist.virus.abun[i], len), (@view ml.matrix_v[i, :]))
    end
end

"""
    applycontrols!(epi::Ecosystem, timestep::Unitful.Time)
Function to apply control strategies to an Ecosystem for one timestep.
"""
function applycontrols!(epi::AbstractEcosystem, timestep::Unitful.Time)
    _applycontrols!(epi, epi.abenv.control, timestep)
end

function _applycontrols!(epi::AbstractEcosystem, controls::NoControl, timestep::Unitful.Time)
    return controls
end

function _applycontrols!(epi::AbstractEcosystem, controls::Lockdown, timestep::Unitful.Time)
    if controls.current_date >= controls.lockdown_date
        epi.spplist.human.work_balance .= 0.0
        controls.current_date += timestep
    else
        controls.current_date += timestep
    end
    return controls
end

"""
    virusmove!(epi::AbstractEcosystem, id::Int64, pos::Int64, grd::Array{Int64, 2}, newvirus::Int64)

Function to calculate the movement of force of infection `id` from a given position in the landscape `pos`, using the lookup table found in the Ecosystem and updating the movement patterns on a cached grid, `grd`. The number of new virus is provided, so that movement only takes place as part of the generation process.
"""
function virusmove!(epi::AbstractEcosystem, id::Int64, pos::Int64, grd::Array{Float64, 2}, newvirus::Int64)
    # Add in home movements
    home = epi.lookup.homelookup
    home_scale = newvirus * epi.spplist.human.home_balance[id]
    if home_scale > zero(home_scale)
        for nzi in home.colptr[pos]:(home.colptr[pos+1]-1)
            grd[id, home.rowval[nzi]] += home_scale * home.nzval[nzi]
        end
    end

    # Add in work movements
    work = epi.lookup.worklookup
    work_scale = newvirus * epi.spplist.human.work_balance[id]
    if work_scale > zero(work_scale)
        for nzi in work.colptr[pos]:(work.colptr[pos+1]-1)
            grd[id, work.rowval[nzi]] += work_scale * work.nzval[nzi]
        end
    end

    return epi
end


"""
    ukChange(epi::Ecosystem, hab::ContinuousHab, timestep::Unitful.Time)

Function to step the uk climate forward by one timestep. Will repeat if time counter becomes greater than the number of dimensions in the habitat.
"""
function ukChange(epi::Ecosystem, hab::ContinuousTimeHab, timestep::Unitful.Time)
    daystep = uconvert(day, timestep)
    hab.time += round(Int64, daystep/day)
    if hab.time > size(hab.matrix, 3)
        hab.time = 1
        Compat.@warn "More timesteps than available, have repeated"
    end
end

function quickmod(x::T, h::T) where {T<:Integer}
    while !(-1 < x < h)
        if x <= -1
            x += h
        elseif x >= h
            x -= h
        end
    end
    return x
end
