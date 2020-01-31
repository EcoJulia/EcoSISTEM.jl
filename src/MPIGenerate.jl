
"""
    update!(eco::Ecosystem, time::Unitful.Time)
Function to update a ecosystem abundances and environment for one timestep.
"""
function update!(eco::MPIEcosystem, timestep::Unitful.Time, ::Val{N}) where N

    # Calculate dimenions of habitat and number of species
    dims = _countsubcommunities(eco.abenv.habitat)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall energy budget of that square
    if eco.rank == 0
        update_energy_usage!(eco, Val{N}())
    end
    MPI.Bcast!(eco.cache.totalE, 0, eco.comm)
    eco.cache.valid = true

    # Loop through species in chosen square
    Threads.@threads for j in eco.firstspecies:(eco.firstspecies + eco.counts[eco.rank + 1] - 1)
        rng = eco.abundances.seed[Threads.threadid()]
        # Loop through grid squares
        for i in 1:dims
            # Calculate how much birth and death should be adjusted
            adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, i, j)

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(i, width)
            # Check if grid cell currently active
            if eco.abenv.active[x, y] && (eco.cache.totalE[i, 1] > 0)

                currentabun = @view eco.abundances.matrix[:, i]

                # Calculate effective rates
                birthprob = params.birth[j] * timestep * adjusted_birth
                deathprob = params.death[j] * timestep * adjusted_death

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                (newbirthprob >= 0) & (newdeathprob >= 0) || error("Birth: $newbirthprob \n Death: $newdeathprob \n \n i: $i \n j: $j")
                # Calculate how many births and deaths
                births = rand(rng, Poisson(currentabun[j] * newbirthprob))
                deaths = rand(rng, Binomial(currentabun[j], newdeathprob))

                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Calculate moves and write to cache
                move!(eco, eco.spplist.movement, i, j, eco.cache.netmigration, births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.matrix .+= eco.cache.netmigration
    gatherings = copy(eco.counts)
    gatherings .*= dims
    m = collect(eco.abundances.matrix')
    MPI.Allgatherv!(MPI.IN_PLACE, m, gatherings, eco.comm)
    eco.abundances.matrix .= m'

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
end

function getlookup(eco::MPIEcosystem, spp::Int64)
    return eco.lookup[spp - eco.firstspecies + 1]
end
