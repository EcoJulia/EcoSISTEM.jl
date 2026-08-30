# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The distributed hot loop: `update!`, `move!`, `update_resource_usage!` and `populate!`.
#
# **These DUPLICATE the serial implementations in `src/dynamics.jl` rather than calling them**, and
# that is the file's defining hazard: a rule phrased as "the hot loop" has two places to apply, and a
# `grep` over `src/` finds only one of them. A change to the demographics, the dispersal step or the
# resource accounting must be made here as well, or the two diverge silently — a distributed run
# would then disagree with a serial one on the same seed, which is the property the whole design
# exists to guarantee.
#
# Everything here is written to keep a result independent of how the work was divided: the per-species
# RNG streams are addressed by global species index, never by rank-local position.

import EcoSISTEM
using MPI
using LinearAlgebra
using Distributions
using Random

using EcoSISTEM: AbstractHabitat, Demand, SpeciesRequirementCollection
using EcoSISTEM: AbstractRegime, AbstractSupply, AbstractNicheFit
using EcoSISTEM: Resource, LayerCollection
using EcoSISTEM:
                 resource_adjustment,
                 invalidatecaches!,
                 regimeupdate!,
                 supplyupdate!,
                 AlwaysMovement,
                 BirthOnlyMovement

"""
    update!(eco::MPIEcosystem, timestep::Unitful.Time)

Update an `MPIEcosystem`'s abundances and environment for one timestep,
computing births, deaths, and dispersal in parallel across threads and MPI
nodes.
"""
function EcoSISTEM.update!(eco::MPIEcosystem, timestep::Unitful.Time)
    return EcoSISTEM.update!(eco, timestep, nothing)
end

"""
    update!(eco::MPIEcosystem, timestep::Unitful.Time, intervention)

Update a distributed ecosystem for one timestep, applying any scheduled [`Intervention`](@ref).

 **The schedule/region machinery is already rank-safe**: selections come from the counter-based
`hash((seed, :intervention, k, step))` stream and the `active` mask and layers are replicated on
every rank, so every rank computes the same cells and makes the same edit without communicating.

 **All six operations work**, including the abundance ones: by the time interventions run the
landscape is row-partitioned, so a rank holds its own species across *all* cells and every abundance
write is rank-local.
"""
function EcoSISTEM.update!(eco::MPIEcosystem, timestep::Unitful.Time,
                           intervention)
    _checkmpiinterventions(intervention)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    # Calculate dimenions of regime and number of species
    numsc = countsubcommunities(eco)
    params = eco.spplist.params
    # Set the overall resource supply of that square
    EcoSISTEM.update_resource_usage!(eco)
    # Share per-cell resource usage across ranks. `totaldemand` is (numsc, numdemands)
    # and each rank owns a contiguous block of cells (rows); gather one demand
    # column at a time so that multi-demand environments (where the columns
    # are not contiguous in the flat buffer) are combined correctly.
    for r in axes(eco.cache.totaldemand, 2)
        MPI.Allgatherv!(MPI.VBuffer(view(eco.cache.totaldemand, :, r),
                                    eco.sccounts),
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
        spstart = (b - 1) * block + 1
        spend = min(b * block, nlocal)
        # Loop through grid squares
        for sc in 1:numsc
            # Convert 1D dimension to 2D coordinates
            (y, x) = EcoSISTEM.convert_coords(eco, sc)
            # Check if grid cell currently active
            (eco.habitat.active[y, x] && (eco.cache.totaldemand[sc, 1] > 0)) ||
                continue
            for mpisp in spstart:spend
                truesp = eco.firstsp + mpisp - 1
                rng = EcoSISTEM.getrng(eco, truesp)
                # Calculate how much birth and death should be adjusted
                adjusted_birth, adjusted_death = resource_adjustment(eco,
                                                                     eco.habitat.supply,
                                                                     sc, truesp)

                # Both are per-individual rates over the timestep. Only the death one becomes
                # a probability: deaths are drawn per individual (Binomial), while births are a
                # count (Poisson) whose mean is the rate itself. `NoUnits` is needed on the birth
                # rate because it reaches `Poisson` as a bare number; the death rate is made
                # dimensionless by `exp`.
                birthrate = params.birth[truesp] * timestep * adjusted_birth |>
                            NoUnits
                deathrate = params.death[truesp] * timestep * adjusted_death

                deathprob = 1.0 - exp(-deathrate)

                (birthrate >= 0) & (deathprob >= 0) ||
                    error("Birth: $birthrate \n Death: $deathprob \n \n sc: $sc \n sp: $truesp")
                # Calculate how many births and deaths
                births = rand(rng,
                              Poisson(eco.abundances.rows_matrix[mpisp, sc] *
                                      birthrate))
                deaths = rand(rng,
                              Binomial(eco.abundances.rows_matrix[mpisp, sc],
                                       deathprob))

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

    # Advance the simulation clock before the layers change, so a layer sees the time it is
    # changing *to*. Every rank advances identically.
    EcoSISTEM._advanceclock!(eco, timestep)

    # Same position as the serial loop: after the dynamics, after the clock, before the layers —
    # so a `SetChange` bites this step. Every rank runs this identically and independently.
    EcoSISTEM.applyinterventions!(eco, intervention,
                                  EcoSISTEM.simulationtime(eco), timestep,
                                  EcoSISTEM._stepnumber(eco, timestep))

    # Update environment - regime and resource supplies
    regimeupdate!(eco, timestep)
    return supplyupdate!(eco, timestep)
end

# `AlwaysMovement` is not supported under MPI, and this method exists to say so.
#
# The serial `move!(::AbstractEcosystem, ::AlwaysMovement, ...)` in `src/dynamics.jl` is typed on
# `AbstractEcosystem`, so it catches an
# `MPIEcosystem` too, and disperses the *whole current population* by reading
# `eco.abundances.matrix[sp, i]` — a field an `MPIGridLandscape` has not got, since it holds
# `rows_matrix`/`cols_vector` instead. Without this method that surfaces as a bare `FieldError`
# from inside a `Threads.@threads` task on the very first timestep, arriving as a
# `TaskFailedException` with no hint that the movement type is what chose it.
#
# **`BirthOnlyMovement` is unaffected** and is why nothing caught this: it disperses only the
# newborns, which the MPI loop passes in explicitly as `births`, so it never reads the landscape —
# and `test/SmallMPItest.jl:63` uses it.
#
# A real implementation needs more than the field rename: the count would come from
# `rows_matrix[mpisp, sc]`, indexed by this rank's **local** species row rather than the global
# `truesp` the loop passes on.
function EcoSISTEM.move!(::MPIEcosystem, ::AlwaysMovement, ::Int64, ::Int64,
                         ::Matrix{Int64}, ::Int64)
    return error("`AlwaysMovement` is not supported in a distributed (MPI) run: it disperses " *
                 "each species' whole population, which is not available rank-locally. Use " *
                 "`movement = BirthOnlyMovement` (dispersal at birth only), or run serially with " *
                 "`build_ecosystem(...; distributed = false)`.")
end

# **Every operation is rank-safe, including the abundance ones**, because of *when* interventions
# run: `synchronise_from_rows!` has already put the landscape in its row-partitioned phase, so this
# rank's `rows_matrix` holds **its own species across all cells**. An abundance write is therefore
# entirely local — the rank owning that species applies the whole edit, the others do nothing — and
# `Deactivate`'s clearance is local too, since each rank clears its own species at the same cells.
# Simpler than the plan's "draw globally, then apply only the entries this rank owns", which was
# written for a harder problem than the phase ordering actually presents.
function EcoSISTEM._ownedabundances(eco::MPIEcosystem)
    return (rows = eco.abundances.rows_matrix, firstspecies = eco.firstsp)
end

# **`AddSpecies` is the one operation that is genuinely not rank-local**, and the contrast with the
# abundance operations is the point: those need no partition work at all, because interventions run
# after `synchronise_from_rows!` when a rank owns whole species across every cell. Adding a species
# changes *which species each rank owns* — the partition itself — so it needs the landscape
# reallocated and redistributed on every rank in step. Refused for now rather than half-done.
function _checkmpiinterventions(intervention)
    for iv in EcoSISTEM._interventions(intervention)
        any(op -> op isa EcoSISTEM.AddSpecies, iv.operations) || continue
        error("`AddSpecies` is not yet supported under MPI. Species are partitioned across ranks, " *
              "so adding one changes the partition itself — unlike every other operation, which is " *
              "rank-local. Add the species before the run, or run this scenario serially.")
    end
    return nothing
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
writing results into `eco.cache.totaldemand`.
"""
function EcoSISTEM.update_resource_usage!(eco::MPIEcosystem{MPIGL, A,
                                                            EcoSISTEM.SpeciesList{TL,
                                                                                  DM,
                                                                                  B,
                                                                                  C,
                                                                                  D},
                                                            E}) where {MPIGL <:
                                                                       MPIGridLandscape,
                                                                       A, B, C,
                                                                       D,
                                                                       E, TL,
                                                                       DM <:
                                                                       Demand}
    !eco.cache.valid || return true

    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    # Get resource supplies of species in square
    ϵ̄ = eco.spplist.demand.resource
    mats = eco.abundances.reshaped_cols

    # Loop through grid squares
    Threads.@threads for sc in 1:eco.sccounts[rank + 1]
        truesc = eco.firstsc + sc - 1
        eco.cache.totaldemand[truesc, 1] = 0.0
        spindex = 1
        for block in eachindex(mats)
            nextsp = spindex + eco.sppcounts[block] - 1
            currentabun = @view mats[block][:, sc]
            e1 = @view ϵ̄[spindex:nextsp]
            eco.cache.totaldemand[truesc, 1] += (currentabun ⋅ e1) *
                                                eco.spplist.demand.exchange_rate
            spindex = nextsp + 1
        end
    end
    return eco.cache.valid = true
end

"""
    update_resource_usage!(eco::MPIEcosystem)

Multi-resource variant of `update_resource_usage!`; updates one column of
`eco.cache.totaldemand` per member of a [`SpeciesRequirementCollection`](@ref).
"""
function EcoSISTEM.update_resource_usage!(eco::MPIEcosystem{MPIGL, A,
                                                            EcoSISTEM.SpeciesList{TL,
                                                                                  DM,
                                                                                  B,
                                                                                  C,
                                                                                  D},
                                                            E}) where {MPIGL <:
                                                                       MPIGridLandscape,
                                                                       A, B, C,
                                                                       D,
                                                                       E, TL,
                                                                       DM <:
                                                                       EcoSISTEM.SpeciesRequirementCollection}
    !eco.cache.valid || return true

    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    # Get resource supplies of species in square
    ds = values(eco.spplist.demand)
    mats = eco.abundances.reshaped_cols

    # Loop through grid squares
    Threads.@threads for sc in 1:eco.sccounts[rank + 1]
        truesc = eco.firstsc + sc - 1
        EcoSISTEM._zerototaldemand!(eco.cache.totaldemand, truesc, ds, 1)
        spindex = 1
        for block in eachindex(mats)
            nextsp = spindex + eco.sppcounts[block] - 1
            currentabun = @view mats[block][:, sc]
            EcoSISTEM._addtotaldemand!(eco.cache.totaldemand, truesc,
                                       currentabun,
                                       spindex, nextsp, ds, 1)
            spindex = nextsp + 1
        end
    end
    return eco.cache.valid = true
end

using EcoSISTEM: getgridshape, calc_lookup_moves!
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
    height, width = getgridshape(eco)
    (y, x) = EcoSISTEM.convert_coords(eco, sc, height)
    lookup = EcoSISTEM.getlookup(eco, truesp)
    calc_lookup_moves!(eco.habitat.topology, y, x, truesp, eco,
                       births)
    # Lose moves from current grid square
    mpisp = truesp - eco.firstsp + 1
    grd[mpisp, sc] -= births
    # Map moves to location in grid
    moves = lookup.moves
    for i in eachindex(lookup.y)
        newy = mod(lookup.y[i] + y - 1, height) + 1
        newx = mod(lookup.x[i] + x - 1, width) + 1
        loc = EcoSISTEM.convert_coords(eco, (newy, newx), height)
        grd[mpisp, loc] += moves[i]
    end
    return eco
end

using EcoSISTEM: getgridshape, _getsupply
"""
    populate!(ml::MPIGridLandscape, spplist::SpeciesList, habitat::AB, nichefit::NF)

Populate an `MPIGridLandscape` by distributing each species' abundance across
active grid cells proportionally to the available supply, then synchronising
from rows to columns across all MPI nodes.
"""
function EcoSISTEM.populate!(ml::MPIGridLandscape,
                             spplist::EcoSISTEM.SpeciesList,
                             habitat::AB,
                             nichefit::NF,
                             rngs::Vector{Random.Xoshiro}) where {AB <:
                                                                  AbstractHabitat,
                                                                  NF <:
                                                                  AbstractNicheFit}
    dim = getgridshape(habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of supply
    b = reshape(parent(ustrip.(_getsupply(habitat.supply))), size(grid))
    activity = reshape(copy(parent(habitat.active)), size(grid))
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
        habitat::GridHabitat{H, <:LayerCollection{Resource}}, nichefit::NF)

Multi-supply variant of `populate!`; distributes abundances proportionally to the
product of the normalised supply fractions.
"""
function EcoSISTEM.populate!(ml::MPIGridLandscape,
                             spplist::EcoSISTEM.SpeciesList,
                             habitat::EcoSISTEM.GridHabitat{H,
                                                            <:LayerCollection{Resource}},
                             nichefit::NF,
                             rngs::Vector{Random.Xoshiro}) where {H <:
                                                                  AbstractRegime,
                                                                  NF <:
                                                                  AbstractNicheFit}
    # Calculate size of regime
    dim = getgridshape(habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    activity = reshape(copy(parent(habitat.active)), size(grid))
    # Set up a copy of each supply, zeroed outside the active cells and normalised to its own
    # total, then multiply them together so a cell needs every resource to be worth populating.
    fractions = EcoSISTEM._zipmap(values(habitat.supply)) do supply
        b = reshape(parent(copy(_getsupply(supply))), size(grid))
        b[.!activity] .= zero(eltype(b))
        return b ./ sum(b)
    end
    B = EcoSISTEM._fold(fractions) do f1, f2
        return f1 .* f2
    end
    # Loop through owned species, drawing from each species' global RNG stream
    abundances = @view spplist.abun[(ml.rows_tuple.first):(ml.rows_tuple.last)]
    for mpisp in eachindex(abundances)
        truesp = ml.rows_tuple.first + mpisp - 1
        rand!(rngs[truesp], Multinomial(abundances[mpisp], B ./ sum(B)),
              (@view ml.rows_matrix[mpisp, :]))
    end
    return EcoSISTEM.synchronise_from_rows!(ml)
end
