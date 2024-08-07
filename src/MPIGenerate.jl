# SPDX-License-Identifier: LGPL-3.0-or-later

import EcoSISTEM
using MPI
using LinearAlgebra
using Distributions

using EcoSISTEM: AbstractAbiotic, Abstract1Requirement, Abstract2Requirements
using EcoSISTEM: AbstractHabitat, AbstractBudget, AbstractTraitRelationship
using EcoSISTEM: energy_adjustment, invalidatecaches!, habitatupdate!,
                 budgetupdate!, BirthOnlyMovement, BudgetCollection2

"""
    update!(eco::MPIEcosystem, timestep::Unitful.Time) where N
Function to update an MPIEcosystem abundances and environment for one timestep.
"""
function EcoSISTEM.update!(eco::MPIEcosystem, timestep::Unitful.Time)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    # Calculate dimenions of habitat and number of species
    numsc = countsubcommunities(eco)
    params = eco.spplist.params
    # Set the overall energy budget of that square
    EcoSISTEM.update_energy_usage!(eco)
    MPI.Allgatherv!(MPI.VBuffer(eco.cache.totalE, eco.sccounts), comm)
    eco.cache.valid = true

    # Loop through species in chosen square
    Threads.@threads for mpisp in 1:eco.sppcounts[rank + 1]
        truesp = eco.firstsp + mpisp - 1
        rng = eco.abundances.rngs[Threads.threadid()]
        # Loop through grid squares
        for sc in 1:numsc
            # Calculate how much birth and death should be adjusted
            adjusted_birth, adjusted_death = energy_adjustment(eco,
                                                               eco.abenv.budget,
                                                               sc, truesp)

            # Convert 1D dimension to 2D coordinates
            (x, y) = EcoSISTEM.convert_coords(eco, sc)
            # Check if grid cell currently active
            if eco.abenv.active[x, y] && (eco.cache.totalE[sc, 1] > 0)

                # Calculate effective rates
                birthprob = params.birth[truesp] * timestep * adjusted_birth
                deathprob = params.death[truesp] * timestep * adjusted_death

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                (newbirthprob >= 0) & (newdeathprob >= 0) ||
                    error("Birth: $newbirthprob \n Death: $newdeathprob \n \n sc: $sc \n sp: $truesp")
                # Calculate how many births and deaths
                births = rand(rng,
                              Poisson(eco.abundances.rows_matrix[mpisp, sc] *
                                      newbirthprob))
                deaths = rand(rng,
                              Binomial(eco.abundances.rows_matrix[mpisp, sc],
                                       newdeathprob))

                # Update population
                eco.abundances.rows_matrix[mpisp, sc] += (births - deaths)

                # Calculate moves and write to cache
                EcoSISTEM.move!(eco, eco.spplist.movement, sc, truesp,
                                eco.cache.netmigration, births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.rows_matrix .+= eco.cache.netmigration
    EcoSISTEM.synchronise_from_rows!(eco.abundances)

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    return budgetupdate!(eco, timestep)
end

function EcoSISTEM.getlookup(eco::MPIEcosystem, sp::Int64)
    return eco.lookup[sp - eco.firstsp + 1]
end

function EcoSISTEM.update_energy_usage!(eco::MPIEcosystem{MPIGL, A,
                                                          EcoSISTEM.SpeciesList{Tr,
                                                                                Req,
                                                                                B,
                                                                                C,
                                                                                D},
                                                          E}) where
         {MPIGL <: MPIGridLandscape, A, B, C, D, E, Tr,
          Req <: Abstract1Requirement}
    !eco.cache.valid || return true

    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.requirement.energy
    mats = eco.abundances.reshaped_cols

    # Loop through grid squares
    Threads.@threads for sc in 1:eco.sccounts[rank + 1]
        truesc = eco.firstsc + sc - 1
        eco.cache.totalE[truesc, 1] = 0.0
        spindex = 1
        for block in eachindex(mats)
            nextsp = spindex + eco.sppcounts[block] - 1
            currentabun = @view mats[block][:, sc]
            e1 = @view ϵ̄[spindex:nextsp]
            eco.cache.totalE[truesc, 1] += (currentabun ⋅ e1) *
                                           eco.spplist.requirement.exchange_rate
            spindex = nextsp + 1
        end
    end
    return eco.cache.valid = true
end

function EcoSISTEM.update_energy_usage!(eco::MPIEcosystem{MPIGL, A,
                                                          EcoSISTEM.SpeciesList{Tr,
                                                                                Req,
                                                                                B,
                                                                                C,
                                                                                D},
                                                          E}) where {
                                                                     MPIGL <:
                                                                     MPIGridLandscape,
                                                                     A,
                                                                     B, C, D, E,
                                                                     Tr,
                                                                     Req <:
                                                                     Abstract2Requirements
                                                                     }
    !eco.cache.valid || return true

    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.r1.energy
    ϵ̄2 = eco.spplist.requirement.r2.energy
    mats = eco.abundances.reshaped_cols

    # Loop through grid squares
    Threads.@threads for sc in 1:eco.sccounts[rank + 1]
        truesc = eco.firstsc + sc - 1
        eco.cache.totalE[truesc, 1] = 0.0
        eco.cache.totalE[truesc, 2] = 0.0
        spindex = 1
        for block in eachindex(mats)
            nextsp = spindex + eco.sppcounts[block] - 1
            currentabun = @view mats[block][:, sc]
            e1 = @view ϵ̄1[spindex:nextsp]
            eco.cache.totalE[truesc, 1] += (currentabun ⋅ e1) *
                                           eco.spplist.requirement.r1.exchange_rate
            e2 = @view ϵ̄2[spindex:nextsp]
            eco.cache.totalE[truesc, 2] += (currentabun ⋅ e2) *
                                           eco.spplist.requirement.r2.exchange_rate
            spindex = nextsp + 1
        end
    end
    return eco.cache.valid = true
end

using EcoSISTEM: getdimension, getboundary, calc_lookup_moves!
function EcoSISTEM.move!(eco::MPIEcosystem, ::BirthOnlyMovement, sc::Int64,
                         truesp::Int64,
                         grd::Matrix{Int64}, births::Int64)
    width, height = getdimension(eco)
    (x, y) = EcoSISTEM.convert_coords(eco, sc, width)
    lookup = EcoSISTEM.getlookup(eco, truesp)
    calc_lookup_moves!(getboundary(eco.spplist.movement), x, y, truesp, eco,
                       births)
    # Lose moves from current grid square
    mpisp = truesp - eco.firstsp + 1
    grd[mpisp, sc] -= births
    # Map moves to location in grid
    mov = lookup.moves
    for i in eachindex(lookup.x)
        newx = mod(lookup.x[i] + x - 1, width) + 1
        newy = mod(lookup.y[i] + y - 1, height) + 1
        loc = EcoSISTEM.convert_coords(eco, (newx, newy), width)
        grd[mpisp, loc] += mov[i]
    end
    return eco
end

using EcoSISTEM: _getdimension, _getbudget
function EcoSISTEM.populate!(ml::MPIGridLandscape,
                             spplist::EcoSISTEM.SpeciesList, abenv::AB,
                             rel::R) where {AB <: AbstractAbiotic,
                                            R <: AbstractTraitRelationship}
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b = reshape(ustrip.(_getbudget(abenv.budget)), size(grid))
    activity = reshape(copy(abenv.active), size(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    B = b ./ sum(b)
    # Loop through species
    abundances = @view spplist.abun[(ml.rows_tuple.first):(ml.rows_tuple.last)]
    for mpisp in eachindex(abundances)
        rand!(Multinomial(abundances[mpisp], B),
              (@view ml.rows_matrix[mpisp, :]))
    end
    return EcoSISTEM.synchronise_from_rows!(ml)
end

function EcoSISTEM.populate!(ml::MPIGridLandscape,
                             spplist::EcoSISTEM.SpeciesList,
                             abenv::EcoSISTEM.GridAbioticEnv{H,
                                                             BudgetCollection2{B1,
                                                                               B2}},
                             rel::R) where {H <: AbstractHabitat,
                                            B1 <: AbstractBudget,
                                            B2 <: AbstractBudget,
                                            R <: AbstractTraitRelationship}
    # Calculate size of habitat
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b1 = reshape(copy(_getbudget(abenv.budget, :b1)), size(grid))
    b2 = reshape(copy(_getbudget(abenv.budget, :b2)), size(grid))
    units1 = unit(b1[1])
    units2 = unit(b2[1])
    activity = reshape(copy(abenv.active), size(grid))
    b1[.!activity] .= 0.0 * units1
    b2[.!activity] .= 0.0 * units2
    B = (b1 ./ sum(b1)) .* (b2 ./ sum(b2))
    # Loop through species
    abundances = @view spplist.abun[(ml.rows_tuple.first):(ml.rows_tuple.last)]
    for mpisp in eachindex(abundances)
        rand!(Multinomial(abundances[mpisp], B ./ sum(B)),
              (@view ml.rows_matrix[mpisp, :]))
    end
    return EcoSISTEM.synchronise_from_rows!(ml)
end
