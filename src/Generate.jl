# SPDX-License-Identifier: LGPL-3.0-or-later

using StatsBase
using LinearAlgebra
using Random

"""
    get_neighbours(mat::Matrix, x_coord::Int64, y_coord::Int64, chess::Int64=4)

Get the neighbours of a grid square in a matrix in 4 or 8 directions
"""
function get_neighbours(mat::Matrix, x_coord::Int64, y_coord::Int64,
                        chess::Int64 = 4)
    # Calculate dimensions
    dims = size(mat)
    x_coord <= dims[1] && y_coord <= dims[2] ||
        error("Coordinates outside grid")
    # Include 4 directions
    if chess == 4
        neighbour_vec = [x_coord y_coord-1
                         x_coord y_coord+1
                         x_coord-1 y_coord
                         x_coord+1 y_coord]
        # Include 8 directions
    elseif chess == 8
        neighbour_vec = [x_coord y_coord-1
                         x_coord y_coord+1
                         x_coord-1 y_coord
                         x_coord+1 y_coord
                         x_coord-1 y_coord-1
                         x_coord-1 y_coord+1
                         x_coord+1 y_coord-1
                         x_coord+1 y_coord+1]
    else
        # Give error if other number chosen than 4 or 8
        error("Can only calculate neighbours in 4 or 8 directions")
    end
    # Remove answers outside of the dimensions of the matrix
    remove = vcat(mapslices(all,
                            [
                             neighbour_vec .>= 1 neighbour_vec[:, 1] .<=
                                                 dims[1] neighbour_vec[:, 2] .<=
                                                         dims[2]],
                            dims = 2)...)
    neighbour_vec = neighbour_vec[remove, :]
    return neighbour_vec
end
"""
    get_neighbours(mat::Matrix, x_coord::Vector{Int64}, y_coord::Vector{Int64},
        chess::Int64=4)

As [`get_neighbours`](@ref) but accepts vectors of coordinates and returns the
combined neighbours for all positions.
"""
function get_neighbours(mat::Matrix,
                        x_coord::Vector{Int64},
                        y_coord::Vector{Int64},
                        chess::Int64 = 4)
    neighbours = map(n -> get_neighbours(mat, x_coord[n], y_coord[n], chess),
                     eachindex(x_coord))
    return vcat(neighbours...)
end

"""
    update!(eco::AbstractEcosystem, timestep::Unitful.Time)

Dispatch function that selects the appropriate threaded implementation of
[`update!`](@ref) based on the number of available threads.
"""
function update!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return update!(eco, timestep, Val{Threads.nthreads()}())
end

"""
    update!(eco::Ecosystem, timestep::Unitful.Time)

Update an ecosystem's abundances and environment for one timestep.
"""
function update!(eco::Ecosystem, timestep::Unitful.Time)

    # Calculate dimenions of habitat and number of species
    dims = _countsubcommunities(eco.abenv.habitat)
    spp = size(eco.abundances.grid, 1)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall energy budget of that square
    update_energy_usage!(eco)

    # Loop through species in cache-line-sized contiguous blocks (see
    # `species_blocksize`): each thread owns whole blocks, and the cell loop sits
    # outside the inner species loop so a block's species — adjacent rows of the
    # column-major (species, cells) matrix — are touched as one cache line. The
    # active/energy gate is per-cell, so it lifts outside the species loop. Each
    # species is still drawn only by its owning thread, in ascending-cell order,
    # so per-species RNG streams stay race-free and reproducible.
    block = species_blocksize()
    nblocks = cld(spp, block)
    # :greedy hands the cache-line-sized species blocks to cores as they free up
    # (dynamic load balancing); blocks are independent so results are unchanged.
    Threads.@threads :greedy for b in 1:nblocks
        jstart = (b - 1) * block + 1
        jend = min(b * block, spp)
        # Loop through grid squares
        for i in 1:dims
            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(eco, i, width)
            # Check if grid cell currently active
            (eco.abenv.active[x, y] && (eco.cache.totalE[i, 1] > 0)) || continue
            for j in jstart:jend
                rng = getrng(eco, j)
                # Calculate how much birth and death should be adjusted
                adjusted_birth, adjusted_death = energy_adjustment(eco,
                                                                   eco.abenv.budget,
                                                                   i, j)

                # Calculate effective rates
                birthrate = params.birth[j] * timestep * adjusted_birth |>
                            NoUnits
                deathparam = params.death[j] * timestep * adjusted_death

                # Turn deathparam into probability and cancel units of birthrate
                deathprob = 1.0 - exp(-deathparam)

                (birthrate >= 0) & (deathprob >= 0) ||
                    error("Birth: $birthrate \n Death: $deathprob \n \n i: $i \n j: $j")
                # Calculate how many births and deaths
                births = rand(rng,
                              Poisson(eco.abundances.matrix[j, i] * birthrate))
                deaths = rand(rng,
                              Binomial(eco.abundances.matrix[j, i], deathprob))

                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Calculate moves and write to cache
                move!(eco, eco.spplist.movement, i, j, eco.cache.netmigration,
                      births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.matrix .+= eco.cache.netmigration

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    return budgetupdate!(eco, timestep)
end

"""
    update_energy_usage!(eco::Ecosystem)
Calculate how much energy has been used up by the current species in each grid
square in the ecosystem, `eco`. This function is parameterised on whether the
species have one type of energy requirement or two.
"""
function update_energy_usage!(eco::AbstractEcosystem{A,
                                                     SpeciesList{Tr, Req, B, C,
                                                                 D}, E}) where {A,
                                                                                B,
                                                                                C,
                                                                                D,
                                                                                E,
                                                                                Tr,
                                                                                Req <:
                                                                                Abstract1Requirement}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.requirement.energy

    # Loop through grid squares
    Threads.@threads for i in Base.axes(eco.abundances.matrix, 2)
        eco.cache.totalE[i, 1] = ((@view eco.abundances.matrix[:, i]) ⋅ ϵ̄) *
                                 eco.spplist.requirement.exchange_rate
    end
    return eco.cache.valid = true
end

function update_energy_usage!(eco::AbstractEcosystem{A,
                                                     SpeciesList{Tr, Req, B, C,
                                                                 D}, E}) where {A,
                                                                                B,
                                                                                C,
                                                                                D,
                                                                                E,
                                                                                Tr,
                                                                                Req <:
                                                                                Abstract2Requirements}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.one.energy
    ϵ̄2 = eco.spplist.requirement.two.energy

    # Loop through grid squares
    Threads.@threads for i in Base.axes(eco.abundances.matrix, 2)
        currentabun = @view eco.abundances.matrix[:, i]
        eco.cache.totalE[i, 1] = (currentabun ⋅ ϵ̄1) *
                                 eco.spplist.requirement.one.exchange_rate
        eco.cache.totalE[i, 2] = (currentabun ⋅ ϵ̄2) *
                                 eco.spplist.requirement.two.exchange_rate
    end
    return eco.cache.valid = true
end

"""
    energy_adjustment(eco::Ecosystem, bud::AbstractBudget, i::Int64, sp::Int64)

Calculate how much birth and death rates should be adjusted by, according to how
much energy is available, `bud`, in the grid square, `i`, and how much energy
the species, `sp`, requires.
"""
function energy_adjustment(eco::AbstractEcosystem, bud::AbstractBudget,
                           i::Int64, sp::Int64)
    # NoGrowth freezes the population
    eco.spplist.params isa NoGrowth && return (0.0, 0.0)

    # Otherwise adjust birth/death rates by the available energy.
    return _energy_adjustment(eco, bud, i, sp)
end

# Birth and death rate multipliers for a single-requirement environment. Weighs
# the species' own energy requirement (`ϵ̄`) and how well its traits match the cell
# (`ϵ̄real`) against the energy available in the cell (`K`) relative to the total
# demand there (`E`): births are boosted when energy is plentiful (`K/E`, capped at
# `params.boost`) and deaths rise as demand approaches the budget (`E/K`). Called
# only for growing populations — [`energy_adjustment`](@ref) short-circuits NoGrowth.
function _energy_adjustment(eco::AbstractEcosystem, bud::AbstractBudget,
                            i::Int64, sp::Int64)
    params = eco.spplist.params
    width = getdimension(eco)[1]
    (x, y) = convert_coords(eco, i, width)
    K = getbudget(eco)[x, y] * eco.spplist.requirement.exchange_rate
    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.requirement.energy[sp] *
        eco.spplist.requirement.exchange_rate
    E = eco.cache.totalE[i, 1]
    # Traits
    ϵ̄real = 1 / traitfun(eco, i, sp)
    # Alter rates by energy available in current pop & own requirements
    birth_energy = ϵ̄^-params.longevity * ϵ̄real^-params.survival *
                   min(K / E, params.boost)
    death_energy = ϵ̄^-params.longevity * ϵ̄real^params.survival * (E / K)
    return birth_energy, death_energy
end

# As above but for a two-requirement environment (e.g. solar energy and water),
# combining the two budgets. The species is limited by whichever resource is
# scarcest: births use the `min` of the two availability ratios (`K1/E1`, `K2/E2`,
# still capped at `params.boost`) and deaths the `max` of the two demand ratios
# (`E1/K1`, `E2/K2`), so both requirements must be met for the population to grow.
function _energy_adjustment(eco::AbstractEcosystem,
                            bud::BudgetCollection2,
                            i::Int64,
                            sp::Int64)
    width = getdimension(eco)[1]
    (x, y) = convert_coords(eco, i, width)
    params = eco.spplist.params
    K1 = _getbudget(eco.abenv.budget.one)[x, y] *
         eco.spplist.requirement.one.exchange_rate
    K2 = _getbudget(eco.abenv.budget.two)[x, y] *
         eco.spplist.requirement.two.exchange_rate
    # Get abundances of square we are interested in
    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.one.energy[sp] *
         eco.spplist.requirement.one.exchange_rate
    ϵ̄2 = eco.spplist.requirement.two.energy[sp] *
         eco.spplist.requirement.two.exchange_rate
    E1 = eco.cache.totalE[i, 1]
    E2 = eco.cache.totalE[i, 2]
    ϵ̄real1 = 1 / traitfun(eco, i, sp)
    ϵ̄real2 = 1 / traitfun(eco, i, sp)
    # Alter rates by energy available in current pop & own requirements
    birth_energy = (ϵ̄1 * ϵ̄2)^-params.longevity *
                   (ϵ̄real1 * ϵ̄real2)^-params.survival *
                   min(K1 / E1, K2 / E2, params.boost)
    death_energy = (ϵ̄1 * ϵ̄2)^-params.longevity *
                   (ϵ̄real1 * ϵ̄real2)^params.survival *
                   max(E1 / K1, E2 / K2)
    return birth_energy, death_energy
end

"""
    convert_coords(eco, i::Int64, width::Int64)
    convert_coords(eco, x::Int64, y::Int64, width::Int64)
Convert coordinates from two-dimensional (`x`,`y`) format to one dimension
(`i`), or vice versa, using the `width` of the grid. This function can also be
applied to arrays of coordinates.
"""
function convert_coords(eco::AbstractEcosystem,
                        i::Int64,
                        width::Int64 = getdimension(eco)[1])
    x = ((i - 1) % width) + 1
    y = div((i - 1), width) + 1
    return (x, y)
end
function convert_coords(eco::AbstractEcosystem,
                        pos::Tuple{Int64, Int64},
                        width::Int64 = getdimension(eco)[1])
    i = pos[1] + width * (pos[2] - 1)
    return i
end

function convert_coords(i::Int64, width::Int64)
    x = ((i - 1) % width) + 1
    y = div((i - 1), width) + 1
    return (x, y)
end
function convert_coords(x::Int64, y::Int64, width::Int64)
    i = x + width * (y - 1)
    return i
end
"""
    calc_lookup_moves!(bound, x::Int64, y::Int64, sp::Int64, eco::Ecosystem, abun::Int64)

Calculate the number of moves taken by a species, `sp`, from a specific grid
square location (`x`, `y`). There is a boundary condition, `bound`, which
determines how the species can move across space (see AbstractBoundary). The
total abundance of individuals is given in `abun`, which may be the number of
births in the timestep, or total individuals.
"""
function calc_lookup_moves!(bound::NoBoundary,
                            x::Int64,
                            y::Int64,
                            sp::Int64,
                            eco::AbstractEcosystem,
                            abun::Int64)
    lookup = getlookup(eco, sp)
    maxX = getdimension(eco)[1] - x
    maxY = getdimension(eco)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        valid = (-x < lookup.x[i] <= maxX) &&
                (-y < lookup.y[i] <= maxY) &&
                (eco.abenv.active[lookup.x[i] + x, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    return rand!(getrng(eco, sp), dist, lookup.moves)
end

function calc_lookup_moves!(bound::Cylinder,
                            x::Int64,
                            y::Int64,
                            sp::Int64,
                            eco::AbstractEcosystem,
                            abun::Int64)
    lookup = getlookup(eco, sp)
    maxX = getdimension(eco)[1] - x
    maxY = getdimension(eco)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x :
               mod(lookup.x[i] + x - 1, getdimension(eco)[1]) + 1

        valid = (-y < lookup.y[i] <= maxY) &&
                (eco.abenv.active[newx, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    return rand!(getrng(eco, sp), dist, lookup.moves)
end

function calc_lookup_moves!(bound::Torus,
                            x::Int64,
                            y::Int64,
                            sp::Int64,
                            eco::AbstractEcosystem,
                            abun::Int64)
    lookup = getlookup(eco, sp)
    maxX = getdimension(eco)[1] - x
    maxY = getdimension(eco)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x :
               mod(lookup.x[i] + x - 1, getdimension(eco)[1]) + 1
        newy = -y < lookup.y[i] <= maxY ? lookup.y[i] + y :
               mod(lookup.y[i] + y - 1, getdimension(eco)[2]) + 1
        valid = eco.abenv.active[newx, newy]

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    return rand!(getrng(eco, sp), dist, lookup.moves)
end

"""
    move!(eco::Ecosystem, ::AbstractMovement, i::Int64, sp::Int64, grd::Matrix{Int64}, abun::Int64)

Calculate the movement of species `sp` from a given position in the landscape
`i`, using the lookup table found in the [`Ecosystem`](@ref) and updating the
movement patterns on a cached grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process,
instead of the entire population
"""
function move!(eco::AbstractEcosystem,
               ::AlwaysMovement,
               i::Int64,
               sp::Int64,
               grd::Matrix{Int64},
               ::Int64)
    # "Animal-like": the whole current population disperses
    return _move!(eco, i, sp, grd, eco.abundances.matrix[sp, i])
end

function move!(eco::AbstractEcosystem,
               ::NoMovement,
               ::Int64,
               ::Int64,
               ::Matrix{Int64},
               ::Int64)
    return eco
end

function move!(eco::AbstractEcosystem,
               ::BirthOnlyMovement,
               i::Int64,
               sp::Int64,
               grd::Matrix{Int64},
               births::Int64)
    # "Plant-like": only the newly born individuals disperse
    return _move!(eco, i, sp, grd, births)
end

# Common dispersal code: move `amount` individuals of species `sp` from
# position `i`, scattering them across the landscape according to the
# species' lookup table.
function _move!(eco::AbstractEcosystem,
                i::Int64,
                sp::Int64,
                grd::Matrix{Int64},
                amount::Int64)
    width, height = getdimension(eco)
    (x, y) = convert_coords(eco, i, width)
    lookup = getlookup(eco, sp)
    calc_lookup_moves!(getboundary(eco.spplist.movement), x, y, sp, eco, amount)
    # Lose moves from current grid square
    grd[sp, i] -= amount
    # Map moves to location in grid
    mov = lookup.moves
    for j in eachindex(lookup.x)
        newx = mod(lookup.x[j] + x - 1, width) + 1
        newy = mod(lookup.y[j] + y - 1, height) + 1
        loc = convert_coords(eco, (newx, newy), width)
        grd[sp, loc] += mov[j]
    end
    return eco
end

# Return the two ingredients the population routines share when spreading
# individuals across the grid:
#   `grid`     - a vector of the linear indices `1:ncells` of every grid cell,
#                used both to size/flatten the budget and to sample cell locations;
#   `activity` - a flattened copy of `abenv.active`, the boolean mask of which cells
#                are habitable. Callers zero the budget of inactive (`false`) cells so
#                that no individuals are ever placed outside the active region.
function _gridactivity(abenv::AbstractAbiotic)
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    activity = reshape(copy(abenv.active), len)
    return grid, activity
end

"""
    populate!(ml::GridLandscape, spplist::SpeciesList, abenv::AbstractAbiotic,
              rel::AbstractTraitRelationship, rngs::Vector{Random.Xoshiro})
    populate!(ml::GridLandscape, spplist::SpeciesList,
              abenv::GridAbioticEnv{H, BudgetCollection2{B1, B2}}, rel, rngs)

Populate the grid landscape `ml` by randomly scattering each species' total
abundance (taken from `spplist.abun`) across the grid cells, choosing each cell
with probability proportional to its available energy budget. Inactive cells are
given zero probability, so no individuals are placed outside the habitable
region. Each species is drawn from its own generator in `rngs`, so the result is
reproducible and independent of the number of threads or MPI processes.

`rel` is unused by these resource-based methods; it is accepted only so that they
share a signature with [`traitpopulate!`](@ref) and can be passed
interchangeably as the population function when constructing an
[`Ecosystem`](@ref). For a two-budget environment (`BudgetCollection2`) the
sampling weight of a cell is the product of its two separately normalised
budgets.
"""
function populate!(ml::GridLandscape,
                   spplist::SpeciesList,
                   abenv::AB,
                   rel::R,
                   rngs::Vector{Random.Xoshiro}) where {AB <: AbstractAbiotic,
                                                        R <:
                                                        AbstractTraitRelationship}
    grid, activity = _gridactivity(abenv)
    # Set up copy of budget
    b = reshape(ustrip.(_getbudget(abenv.budget)), length(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    B = b ./ sum(b)
    # Loop through species, drawing from each species' own RNG stream
    for i in eachindex(spplist.abun)
        rand!(rngs[i], Multinomial(spplist.abun[i], B), (@view ml.matrix[i, :]))
    end
end

function populate!(ml::GridLandscape,
                   spplist::SpeciesList,
                   abenv::GridAbioticEnv{H, BudgetCollection2{B1, B2}},
                   rel::R,
                   rngs::Vector{Random.Xoshiro}) where {H <: AbstractHabitat,
                                                        B1 <: AbstractBudget,
                                                        B2 <: AbstractBudget,
                                                        R <:
                                                        AbstractTraitRelationship}
    # Calculate size of habitat
    grid, activity = _gridactivity(abenv)
    # Set up copy of budget
    b1 = reshape(copy(_getbudget(abenv.budget, :one)), length(grid))
    b2 = reshape(copy(_getbudget(abenv.budget, :two)), length(grid))
    units1 = unit(b1[1])
    units2 = unit(b2[1])
    b1[.!activity] .= 0.0 * units1
    b2[.!activity] .= 0.0 * units2
    B = (b1 ./ sum(b1)) .* (b2 ./ sum(b2))
    # Loop through species, drawing from each species' own RNG stream
    for i in eachindex(spplist.abun)
        rand!(rngs[i], Multinomial(spplist.abun[i], B ./ sum(B)),
              (@view ml.matrix[i, :]))
    end
end

"""
    repopulate!(eco::Ecosystem)
    repopulate!(eco::Ecosystem, abun::Int64)

Repopulate an ecosystem `eco` by redistributing abundances according to resource
availability. If an `abun` parameter is given, that number of individuals of the
final species is added at randomly sampled locations instead.
"""
function repopulate!(eco::Ecosystem)
    eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
    eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun),
                                        length(eco.spplist.abun)))
    return populate!(eco.abundances, eco.spplist, eco.abenv, eco.relationship,
                     eco.rngs)
end

function repopulate!(eco::Ecosystem, abun::Int64)
    grid, activity = _gridactivity(eco.abenv)
    # Set up copy of budget
    b = reshape(copy(_getbudget(eco.abenv.budget)), length(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    # Draw locations from the last species' own RNG stream
    pos = sample(getrng(eco, lastindex(eco.rngs)), grid[b .> (0 * units)], abun)
    # Add individual to this location
    map(pos) do p
        return eco.abundances.matrix[end, p] += 1
    end
end

"""
    traitpopulate!(ml::GridLandscape, spplist::SpeciesList, abenv::AbstractAbiotic,
                   rel::AbstractTraitRelationship, rngs::Vector{Random.Xoshiro})

Populate the grid landscape `ml` by scattering each species' total abundance
(taken from `spplist.abun`) across the grid cells with probability proportional
to how well the species' traits match each cell's environment, as scored by the
trait relationship `rel` applied to `spplist.traits` and `abenv.habitat`. Where a
species matches no cell the distribution falls back to uniform. Only native
species (those flagged in `spplist.native`) are placed; non-native species are
left empty.

This is the trait-based counterpart of [`populate!`](@ref), which instead weights
cells by their available energy budget.
"""
function traitpopulate!(ml::GridLandscape,
                        spplist::SpeciesList,
                        abenv::AB,
                        rel::R,
                        rngs::Vector{Random.Xoshiro}) where {AB <:
                                                             AbstractAbiotic,
                                                             R <:
                                                             AbstractTraitRelationship}
    # Calculate size of habitat
    dim = _getdimension(abenv.habitat)
    numsquares = dim[1] * dim[2]
    numspp = length(spplist.names)
    maxrng = spplist.traits.mean .+ spplist.traits.sd
    minrng = spplist.traits.mean .- spplist.traits.sd
    hab = reshape(abenv.habitat.matrix, numsquares)
    probabilities = [_traitfun(abenv.habitat, spplist.traits, rel, i, sp)
                     for i in 1:numsquares,
                         sp in 1:numspp]
    # Loop through species, drawing from each species' own RNG stream
    for i in eachindex(spplist.abun)
        if spplist.native[i]
            # Get abundance of species
            probs = probabilities[:, i] ./ sum(probabilities[:, i])
            probs[isnan.(probs)] .= 1 / numsquares
            abun = rand(rngs[i], Multinomial(spplist.abun[i], probs))
            # Add individual to this location
            ml.matrix[i, :] .+= abun
        end
    end
end

"""
    traitrepopulate!(eco::Ecosystem)

Repopulate an ecosystem `eco` according to how well species traits match their
environment, redistributing the total abundance across species at random.
"""
function traitrepopulate!(eco::Ecosystem)
    eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
    eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun),
                                        length(eco.spplist.abun)))
    return traitpopulate!(eco.abundances, eco.spplist, eco.abenv,
                          eco.relationship, eco.rngs)
end

"""
    emptypopulate!(ml::GridLandscape, spplist::SpeciesList, abenv::AB, rel::R,
                   rngs::Vector{Random.Xoshiro}) where {AB <: EcoSISTEM.AbstractAbiotic, R <: EcoSISTEM.AbstractTraitRelationship}

Placeholder population function that leaves the landscape empty and warns.
"""
function emptypopulate!(ml::GridLandscape,
                        spplist::SpeciesList,
                        abenv::AB,
                        rel::R,
                        rngs::Vector{Random.Xoshiro}) where {AB <:
                                                             EcoSISTEM.AbstractAbiotic,
                                                             R <:
                                                             EcoSISTEM.AbstractTraitRelationship}
    @warn "Ecosystem not populated!"
end
"""
    reenergise!(eco::Ecosystem, budget::Union{Float64, Unitful.Quantity{Float64}}, grid::Tuple{Int64, Int64})
Refill an ecosystem `eco`, with energy from a budget value, `budget` and a grid
size.
"""
function reenergise!(eco::Ecosystem,
                     budget::Union{Float64, Unitful.Quantity{Float64}},
                     grid::Tuple{Int64, Int64})
    return fill!(eco.abenv.budget.matrix, budget / (grid[1] * grid[2]))
end
function reenergise!(eco::Ecosystem,
                     budget::Tuple{Unitful.Quantity{Float64},
                                   Unitful.Quantity{Float64}},
                     grid::Tuple{Int64, Int64})
    fill!(eco.abenv.budget.one.matrix, budget[1] / (grid[1] * grid[2]))
    return fill!(eco.abenv.budget.two.matrix, budget[2] / (grid[1] * grid[2]))
end
