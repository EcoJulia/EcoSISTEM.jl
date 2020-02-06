"""
    update!(eco::MPIEcosystem, timestep::Unitful.Time) where N
Function to update an MPIEcosystem abundances and environment for one timestep.
"""
function update!(eco::MPIEcosystem, timestep::Unitful.Time)

    # Calculate dimenions of habitat and number of species
    nsc = countsubcommunities(eco)
    nsp = counttypes(eco)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall energy budget of that square
    update_energy_usage!(eco)
    MPI.Allgatherv!(MPI.IN_PLACE, eco.cache.totalE,
                    eco.sccounts, eco.comm)
    eco.cache.valid = true

    # Loop through species in chosen square
    Threads.@threads for j in eco.firstspecies:(eco.firstspecies + eco.speciescounts[eco.rank + 1] - 1)
        rng = eco.abundances.seed[Threads.threadid()]
        # Loop through grid squares
        for i in 1:nsc
            # Calculate how much birth and death should be adjusted
            adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, i, j)

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(i, width)
            # Check if grid cell currently active
            if eco.abenv.active[x, y] && (eco.cache.totalE[i, 1] > 0)

                # Calculate effective rates
                birthprob = params.birth[j] * timestep * adjusted_birth
                deathprob = params.death[j] * timestep * adjusted_death

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                (newbirthprob >= 0) & (newdeathprob >= 0) || error("Birth: $newbirthprob \n Death: $newdeathprob \n \n i: $i \n j: $j")
                # Calculate how many births and deaths
                births = rand(rng, Poisson(eco.abundances.matrix[j, i] * newbirthprob))
                deaths = rand(rng, Binomial(eco.abundances.matrix[j, i], newdeathprob))

                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Calculate moves and write to cache
                move!(eco, eco.spplist.movement, i, j, eco.cache.netmigration, births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.matrix .+= eco.cache.netmigration
    gatherings = copy(eco.speciescounts)
    gatherings .*= nsc
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

function update_energy_usage!(eco::MPIEcosystem{A, SpeciesList{Tr,  Req, B, C, D}, E}) where {A, B, C, D, E, Tr, Req <: Abstract1Requirement}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.requirement.energy

    # Loop through grid squares
    Threads.@threads for i in eco.firstsc:(eco.firstsc + eco.sccounts[eco.rank + 1] - 1)
        eco.cache.totalE[i, 1] = ((@view eco.abundances.matrix[:, i]) ⋅ ϵ̄) * eco.spplist.requirement.exchange_rate
    end
    eco.cache.valid = true
end

function update_energy_usage!(eco::MPIEcosystem{A, SpeciesList{Tr,  Req, B, C, D}, E}) where {A, B, C, D, E, Tr, Req <: Abstract2Requirements}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.r1.energy
    ϵ̄2 = eco.spplist.requirement.r2.energy

    # Loop through grid squares
    Threads.@threads for i in eco.firstsc:(eco.firstsc + eco.sccounts[eco.rank + 1] - 1)
        currentabun = @view eco.abundances.matrix[:, i]
        eco.cache.totalE[i, 1] = (currentabun ⋅ ϵ̄1) * eco.spplist.requirement.r1.exchange_rate
        eco.cache.totalE[i, 2] = (currentabun ⋅ ϵ̄2) * eco.spplist.requirement.r2.exchange_rate
    end
    eco.cache.valid = true
end
