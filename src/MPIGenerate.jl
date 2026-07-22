# SPDX-License-Identifier: LGPL-3.0-or-later

import EcoSISTEM
using MPI
using LinearAlgebra
using Distributions
using Random

using EcoSISTEM: AbstractAbiotic, Abstract1Demand, Abstract2Demands
using EcoSISTEM: AbstractHabitat, AbstractSupply, AbstractTraitRelationship
using EcoSISTEM:
                 resource_adjustment,
                 invalidatecaches!,
                 habitatupdate!,
                 supplyupdate!,
                 BirthOnlyMovement,
                 SupplyCollection2

"""
    update!(eco::MPIEcosystem, timestep::Unitful.Time)

Update an `MPIEcosystem`'s abundances and environment for one timestep,
computing births, deaths, and dispersal in parallel across threads and MPI
nodes.
"""
function EcoSISTEM.update!(eco::MPIEcosystem, timestep::Unitful.Time)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    # Calculate dimenions of habitat and number of species
    numsc = countsubcommunities(eco)
    params = eco.spplist.params
    # Set the overall resource supply of that square
    EcoSISTEM.update_resource_usage!(eco)
    # Share per-cell resource usage across ranks. `totalE` is (numsc, numdemands)
    # and each rank owns a contiguous block of cells (rows); gather one demand
    # column at a time so that multi-demand environments (where the columns
    # are not contiguous in the flat buffer) are combined correctly.
    for r in axes(eco.cache.totalE, 2)
        MPI.Allgatherv!(MPI.VBuffer(view(eco.cache.totalE, :, r), eco.sccounts),
                        comm)
    end
    eco.cache.valid = true

    # Loop through this rank's species in cache-line-sized contiguous blocks
    # (see `EcoSISTEM.species_blocksize`), cells outside the inner species loop,
    # so a block's species (adjacent rows of the column-major rows_matrix) are
    # touched as one cache line. The active/resource gate is per-cell. Each species
    # is still drawn only by its owning thread, in ascending-cell order, so
    # per-species RNG streams stay race-free and reproducible.
    nlocal = eco.sppcounts[rank + 1]
    block = EcoSISTEM.species_blocksize()
    nblocks = cld(nlocal, block)
    # :greedy hands the cache-line-sized species blocks to cores as they free up
    # (dynamic load balancing); blocks are independent so results are unchanged.
    Threads.@threads :greedy for b in 1:nblocks
        mpistart = (b - 1) * block + 1
        mpiend = min(b * block, nlocal)
        # Loop through grid squares
        for sc in 1:numsc
            # Convert 1D dimension to 2D coordinates
            (x, y) = EcoSISTEM.convert_coords(eco, sc)
            # Check if grid cell currently active
            (eco.abenv.active[x, y] && (eco.cache.totalE[sc, 1] > 0)) ||
                continue
            for mpisp in mpistart:mpiend
                truesp = eco.firstsp + mpisp - 1
                rng = EcoSISTEM.getrng(eco, truesp)
                # Calculate how much birth and death should be adjusted
                adjusted_birth, adjusted_death = resource_adjustment(eco,
                                                                     eco.abenv.supply,
                                                                     sc, truesp)

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
                EcoSISTEM.move!(eco,
                                eco.spplist.movement,
                                sc,
                                truesp,
                                eco.cache.netmigration,
                                births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.rows_matrix .+= eco.cache.netmigration
    EcoSISTEM.synchronise_from_rows!(eco.abundances)

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and resource supplies
    habitatupdate!(eco, timestep)
    return supplyupdate!(eco, timestep)
end

"""
    getlookup(eco::MPIEcosystem, sp::Int64)

Return the movement lookup table for species `sp` from an `MPIEcosystem`,
adjusting the species index by the node's `firstsp` offset.
"""
function EcoSISTEM.getlookup(eco::MPIEcosystem, sp::Int64)
    return eco.lookup[sp - eco.firstsp + 1]
end

"""
    update_resource_usage!(eco::MPIEcosystem)

Update the total resource usage cache for a single-resource `MPIEcosystem`,
summing each species' abundance × resource demand across all MPI blocks and
writing results into `eco.cache.totalE`.
"""
function EcoSISTEM.update_resource_usage!(eco::MPIEcosystem{MPIGL, A,
                                                            EcoSISTEM.SpeciesList{Tr,
                                                                                  Req,
                                                                                  B,
                                                                                  C,
                                                                                  D},
                                                            E}) where {MPIGL <:
                                                                       MPIGridLandscape,
                                                                       A, B, C,
                                                                       D,
                                                                       E, Tr,
                                                                       Req <:
                                                                       Abstract1Demand}
    !eco.cache.valid || return true

    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    # Get resource supplies of species in square
    ϵ̄ = eco.spplist.demand.resource
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
                                           eco.spplist.demand.exchange_rate
            spindex = nextsp + 1
        end
    end
    return eco.cache.valid = true
end

"""
    update_resource_usage!(eco::MPIEcosystem)

Two-resource variant of `update_resource_usage!`; updates both columns of
`eco.cache.totalE` for environments with `Abstract2Demands`.
"""
function EcoSISTEM.update_resource_usage!(eco::MPIEcosystem{MPIGL, A,
                                                            EcoSISTEM.SpeciesList{Tr,
                                                                                  Req,
                                                                                  B,
                                                                                  C,
                                                                                  D},
                                                            E}) where {MPIGL <:
                                                                       MPIGridLandscape,
                                                                       A, B, C,
                                                                       D,
                                                                       E, Tr,
                                                                       Req <:
                                                                       Abstract2Demands}
    !eco.cache.valid || return true

    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    # Get resource supplies of species in square
    ϵ̄1 = eco.spplist.demand.one.resource
    ϵ̄2 = eco.spplist.demand.two.resource
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
                                           eco.spplist.demand.one.exchange_rate
            e2 = @view ϵ̄2[spindex:nextsp]
            eco.cache.totalE[truesc, 2] += (currentabun ⋅ e2) *
                                           eco.spplist.demand.two.exchange_rate
            spindex = nextsp + 1
        end
    end
    return eco.cache.valid = true
end

using EcoSISTEM: getdimension, getboundary, calc_lookup_moves!
"""
    move!(eco::MPIEcosystem, ::BirthOnlyMovement, sc::Int64, truesp::Int64,
        grd::Matrix{Int64}, births::Int64)

Apply dispersal for `births` new individuals of species `truesp` from grid cell
`sc` using the [`BirthOnlyMovement`](@ref) kernel, writing net moves into `grd`.
"""
function EcoSISTEM.move!(eco::MPIEcosystem,
                         ::BirthOnlyMovement,
                         sc::Int64,
                         truesp::Int64,
                         grd::Matrix{Int64},
                         births::Int64)
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

using EcoSISTEM: _getdimension, _getsupply
"""
    populate!(ml::MPIGridLandscape, spplist::SpeciesList, abenv::AB, rel::R)

Populate an `MPIGridLandscape` by distributing each species' abundance across
active grid cells proportionally to the available supply, then synchronising
from rows to columns across all MPI nodes.
"""
function EcoSISTEM.populate!(ml::MPIGridLandscape,
                             spplist::EcoSISTEM.SpeciesList,
                             abenv::AB,
                             rel::R,
                             rngs::Vector{Random.Xoshiro}) where {AB <:
                                                                  AbstractAbiotic,
                                                                  R <:
                                                                  AbstractTraitRelationship}
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of supply
    b = reshape(ustrip.(_getsupply(abenv.supply)), size(grid))
    activity = reshape(copy(abenv.active), size(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    B = b ./ sum(b)
    # Loop through owned species, drawing from each species' global RNG stream
    abundances = @view spplist.abun[(ml.rows_tuple.first):(ml.rows_tuple.last)]
    for mpisp in eachindex(abundances)
        truesp = ml.rows_tuple.first + mpisp - 1
        rand!(rngs[truesp], Multinomial(abundances[mpisp], B),
              (@view ml.rows_matrix[mpisp, :]))
    end
    return EcoSISTEM.synchronise_from_rows!(ml)
end

"""
    populate!(ml::MPIGridLandscape, spplist::SpeciesList,
        abenv::GridAbioticEnv{H, SupplyCollection2{B1, B2}}, rel::R)

Two-supply variant of `populate!`; distributes abundances proportionally to the
product of the two normalised supply fractions.
"""
function EcoSISTEM.populate!(ml::MPIGridLandscape,
                             spplist::EcoSISTEM.SpeciesList,
                             abenv::EcoSISTEM.GridAbioticEnv{H,
                                                             SupplyCollection2{B1,
                                                                               B2}},
                             rel::R,
                             rngs::Vector{Random.Xoshiro}) where {H <:
                                                                  AbstractHabitat,
                                                                  B1 <:
                                                                  AbstractSupply,
                                                                  B2 <:
                                                                  AbstractSupply,
                                                                  R <:
                                                                  AbstractTraitRelationship}
    # Calculate size of habitat
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of supply
    b1 = reshape(copy(_getsupply(abenv.supply, :one)), size(grid))
    b2 = reshape(copy(_getsupply(abenv.supply, :two)), size(grid))
    units1 = unit(b1[1])
    units2 = unit(b2[1])
    activity = reshape(copy(abenv.active), size(grid))
    b1[.!activity] .= 0.0 * units1
    b2[.!activity] .= 0.0 * units2
    B = (b1 ./ sum(b1)) .* (b2 ./ sum(b2))
    # Loop through owned species, drawing from each species' global RNG stream
    abundances = @view spplist.abun[(ml.rows_tuple.first):(ml.rows_tuple.last)]
    for mpisp in eachindex(abundances)
        truesp = ml.rows_tuple.first + mpisp - 1
        rand!(rngs[truesp], Multinomial(abundances[mpisp], B ./ sum(B)),
              (@view ml.rows_matrix[mpisp, :]))
    end
    return EcoSISTEM.synchronise_from_rows!(ml)
end
