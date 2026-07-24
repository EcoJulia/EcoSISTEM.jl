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

    # Calculate dimenions of regime and number of species
    dims = _countsubcommunities(eco.habitat.regime)
    nspp = size(eco.abundances.grid, 1)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall resource supply of that square
    update_resource_usage!(eco)

    # Loop through species in cache-line-sized contiguous blocks (see
    # `species_blocksize`): each thread owns whole blocks, and the cell loop sits
    # outside the inner species loop so a block's species — adjacent rows of the
    # column-major (species, cells) matrix — are touched as one cache line. The
    # active/resource gate is per-cell, so it lifts outside the species loop. Each
    # species is still drawn only by its owning thread, in ascending-cell order,
    # so per-species RNG streams stay race-free and reproducible.
    block = species_blocksize()
    nblocks = cld(nspp, block)
    # :greedy hands the cache-line-sized species blocks to cores as they free up
    # (dynamic load balancing); blocks are independent so results are unchanged.
    Threads.@threads :greedy for b in 1:nblocks
        jstart = (b - 1) * block + 1
        jend = min(b * block, nspp)
        # Loop through grid squares
        for i in 1:dims
            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(eco, i, width)
            # Check if grid cell currently active
            (eco.habitat.active[x, y] && (eco.cache.totalE[i, 1] > 0)) ||
                continue
            for j in jstart:jend
                rng = getrng(eco, j)
                # Calculate how much birth and death should be adjusted
                adjusted_birth, adjusted_death = resource_adjustment(eco,
                                                                     eco.habitat.supply,
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

    # Update environment - regime and resource supplies
    regimeupdate!(eco, timestep)
    return supplyupdate!(eco, timestep)
end

"""
    update_resource_usage!(eco::Ecosystem)
Calculate how much resource has been used up by the current species in each grid
square in the ecosystem, `eco`. This function is parameterised on whether the
species have one type of resource demand or two.
"""
function update_resource_usage!(eco::AbstractEcosystem{A,
                                                       SpeciesList{TL, DM, B,
                                                                   C,
                                                                   D}, E}) where {A,
                                                                                  B,
                                                                                  C,
                                                                                  D,
                                                                                  E,
                                                                                  TL,
                                                                                  DM <:
                                                                                  Abstract1Demand}
    !eco.cache.valid || return true

    # Get resource supplies of species in square
    ϵ̄ = eco.spplist.demand.resource

    # Loop through grid squares
    Threads.@threads for i in Base.axes(eco.abundances.matrix, 2)
        eco.cache.totalE[i, 1] = ((@view eco.abundances.matrix[:, i]) ⋅ ϵ̄) *
                                 eco.spplist.demand.exchange_rate
    end
    return eco.cache.valid = true
end

function update_resource_usage!(eco::AbstractEcosystem{A,
                                                       SpeciesList{TL, DM, B,
                                                                   C,
                                                                   D}, E}) where {A,
                                                                                  B,
                                                                                  C,
                                                                                  D,
                                                                                  E,
                                                                                  TL,
                                                                                  DM <:
                                                                                  Abstract2Demands}
    !eco.cache.valid || return true

    # Get resource supplies of species in square
    ϵ̄1 = eco.spplist.demand.one.resource
    ϵ̄2 = eco.spplist.demand.two.resource

    # Loop through grid squares
    Threads.@threads for i in Base.axes(eco.abundances.matrix, 2)
        currentabun = @view eco.abundances.matrix[:, i]
        eco.cache.totalE[i, 1] = (currentabun ⋅ ϵ̄1) *
                                 eco.spplist.demand.one.exchange_rate
        eco.cache.totalE[i, 2] = (currentabun ⋅ ϵ̄2) *
                                 eco.spplist.demand.two.exchange_rate
    end
    return eco.cache.valid = true
end

"""
    resource_adjustment(eco::Ecosystem, supply::AbstractSupply, i::Int64, sp::Int64)

Calculate how much birth and death rates should be adjusted by, according to how
much resource is available, `supply`, in the grid square, `i`, and how much resource
the species, `sp`, requires.
"""
function resource_adjustment(eco::AbstractEcosystem, supply::AbstractSupply,
                             i::Int64, sp::Int64)
    # NoGrowth freezes the population
    eco.spplist.params isa NoGrowth && return (0.0, 0.0)

    # Otherwise adjust birth/death rates by the available resource.
    return _resource_adjustment(eco, supply, i, sp)
end

# Birth and death rate multipliers for a single-demand environment. Weighs
# the species' own resource demand (`ϵ̄`) and how well its tolerances match the cell
# (`ϵ̄real`) against the resource available in the cell (`K`) relative to the total
# demand there (`E`): births are boosted when resource is plentiful (`K/E`, capped at
# `params.boost`) and deaths rise as demand approaches the supply (`E/K`). Called
# only for growing populations — [`resource_adjustment`](@ref) short-circuits NoGrowth.
function _resource_adjustment(eco::AbstractEcosystem, supply::AbstractSupply,
                              i::Int64, sp::Int64)
    params = eco.spplist.params
    width = getdimension(eco)[1]
    (x, y) = convert_coords(eco, i, width)
    K = getsupply(eco)[x, y] * eco.spplist.demand.exchange_rate
    # Get resource supplies of species in square
    ϵ̄ = eco.spplist.demand.resource[sp] *
        eco.spplist.demand.exchange_rate
    E = eco.cache.totalE[i, 1]
    # Traits
    ϵ̄real = 1 / suitability(eco, i, sp)
    # Alter rates by resource available in current pop & own demands
    birth_resource = ϵ̄^-params.longevity * ϵ̄real^-params.survival *
                     min(K / E, params.boost)
    death_resource = ϵ̄^-params.longevity * ϵ̄real^params.survival * (E / K)
    return birth_resource, death_resource
end

# As above but for a two-demand environment (e.g. solar resource and water),
# combining the two supplies. The species is limited by whichever resource is
# scarcest: births use the `min` of the two availability ratios (`K1/E1`, `K2/E2`,
# still capped at `params.boost`) and deaths the `max` of the two demand ratios
# (`E1/K1`, `E2/K2`), so both demands must be met for the population to grow.
function _resource_adjustment(eco::AbstractEcosystem,
                              supply::SupplyCollection2,
                              i::Int64,
                              sp::Int64)
    width = getdimension(eco)[1]
    (x, y) = convert_coords(eco, i, width)
    params = eco.spplist.params
    K1 = _getsupply(eco.habitat.supply.one)[x, y] *
         eco.spplist.demand.one.exchange_rate
    K2 = _getsupply(eco.habitat.supply.two)[x, y] *
         eco.spplist.demand.two.exchange_rate
    # Get abundances of square we are interested in
    # Get resource supplies of species in square
    ϵ̄1 = eco.spplist.demand.one.resource[sp] *
         eco.spplist.demand.one.exchange_rate
    ϵ̄2 = eco.spplist.demand.two.resource[sp] *
         eco.spplist.demand.two.exchange_rate
    E1 = eco.cache.totalE[i, 1]
    E2 = eco.cache.totalE[i, 2]
    ϵ̄real1 = 1 / suitability(eco, i, sp)
    ϵ̄real2 = 1 / suitability(eco, i, sp)
    # Alter rates by resource available in current pop & own demands
    birth_resource = (ϵ̄1 * ϵ̄2)^-params.longevity *
                     (ϵ̄real1 * ϵ̄real2)^-params.survival *
                     min(K1 / E1, K2 / E2, params.boost)
    death_resource = (ϵ̄1 * ϵ̄2)^-params.longevity *
                     (ϵ̄real1 * ϵ̄real2)^params.survival *
                     max(E1 / K1, E2 / K2)
    return birth_resource, death_resource
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
                (eco.habitat.active[lookup.x[i] + x, lookup.y[i] + y])

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
                (eco.habitat.active[newx, lookup.y[i] + y])

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
        valid = eco.habitat.active[newx, newy]

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
    moves = lookup.moves
    for j in eachindex(lookup.x)
        newx = mod(lookup.x[j] + x - 1, width) + 1
        newy = mod(lookup.y[j] + y - 1, height) + 1
        loc = convert_coords(eco, (newx, newy), width)
        grd[sp, loc] += moves[j]
    end
    return eco
end

# Return the two ingredients the population routines share when spreading
# individuals across the grid:
#   `grid`     - a vector of the linear indices `1:ncells` of every grid cell,
#                used both to size/flatten the supply and to sample cell locations;
#   `activity` - a flattened copy of `habitat.active`, the boolean mask of which cells
#                are habitable. Callers zero the supply of inactive (`false`) cells so
#                that no individuals are ever placed outside the active region.
function _gridactivity(habitat::AbstractHabitat)
    dim = _getdimension(habitat.regime)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    activity = reshape(copy(habitat.active), len)
    return grid, activity
end

"""
    populate!(ml::GridLandscape, spplist::SpeciesList, habitat::AbstractHabitat,
              nichefit::AbstractNicheFit, rngs::Vector{Random.Xoshiro})
    populate!(ml::GridLandscape, spplist::SpeciesList,
              habitat::GridHabitat{H, SupplyCollection2{B1, B2}}, nichefit, rngs)

Populate the grid landscape `ml` by randomly scattering each species' total
abundance (taken from `spplist.abun`) across the grid cells, choosing each cell
with probability proportional to its available resource supply. Inactive cells are
given zero probability, so no individuals are placed outside the habitable
region. Each species is drawn from its own generator in `rngs`, so the result is
reproducible and independent of the number of threads or MPI processes.

`nichefit` is unused by these resource-based methods; it is accepted only so that they
share a signature with [`tolerancepopulate!`](@ref) and can be passed
interchangeably as the population function when constructing an
[`Ecosystem`](@ref). For a two-supply environment (`SupplyCollection2`) the
sampling weight of a cell is the product of its two separately normalised
supplies.
"""
function populate!(ml::GridLandscape,
                   spplist::SpeciesList,
                   habitat::AB,
                   nichefit::R,
                   rngs::Vector{Random.Xoshiro}) where {AB <: AbstractHabitat,
                                                        R <:
                                                        AbstractNicheFit}
    grid, activity = _gridactivity(habitat)
    # Set up copy of supply
    b = reshape(ustrip.(_getsupply(habitat.supply)), length(grid))
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
                   habitat::GridHabitat{H, SupplyCollection2{B1, B2}},
                   nichefit::R,
                   rngs::Vector{Random.Xoshiro}) where {H <: AbstractRegime,
                                                        B1 <: AbstractSupply,
                                                        B2 <: AbstractSupply,
                                                        R <:
                                                        AbstractNicheFit}
    # Calculate size of regime
    grid, activity = _gridactivity(habitat)
    # Set up copy of supply
    b1 = reshape(copy(_getsupply(habitat.supply, :one)), length(grid))
    b2 = reshape(copy(_getsupply(habitat.supply, :two)), length(grid))
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
    eco.abundances = emptygridlandscape(eco.habitat, eco.spplist)
    eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun),
                                        length(eco.spplist.abun)))
    return populate!(eco.abundances, eco.spplist, eco.habitat, eco.nichefit,
                     eco.rngs)
end

function repopulate!(eco::Ecosystem, abun::Int64)
    grid, activity = _gridactivity(eco.habitat)
    # Set up copy of supply
    b = reshape(copy(_getsupply(eco.habitat.supply)), length(grid))
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
    tolerancepopulate!(ml::GridLandscape, spplist::SpeciesList, habitat::AbstractHabitat,
                   nichefit::AbstractNicheFit, rngs::Vector{Random.Xoshiro})

Populate the grid landscape `ml` by scattering each species' total abundance
(taken from `spplist.abun`) across the grid cells with probability proportional
to how well the species tolerances match each cell's environment, as scored by the
trait nichefit `nichefit` applied to `spplist.tolerance` and `habitat.regime`. Where a
species matches no cell the distribution falls back to uniform. Only native
species (those flagged in `spplist.native`) are placed; non-native species are
left empty.

This is the trait-based counterpart of [`populate!`](@ref), which instead weights
cells by their available resource supply.
"""
function tolerancepopulate!(ml::GridLandscape,
                            spplist::SpeciesList,
                            habitat::AB,
                            nichefit::R,
                            rngs::Vector{Random.Xoshiro}) where {AB <:
                                                                 AbstractHabitat,
                                                                 R <:
                                                                 AbstractNicheFit}
    # Calculate size of regime
    dim = _getdimension(habitat.regime)
    numsquares = dim[1] * dim[2]
    numspp = length(spplist.names)
    regime = reshape(habitat.regime.matrix, numsquares)
    probabilities = [_suitability(habitat.regime, spplist.tolerance, nichefit,
                                  i, sp)
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
    tolerancerepopulate!(eco::Ecosystem)

Repopulate an ecosystem `eco` according to how well species tolerances match their
environment, redistributing the total abundance across species at random.
"""
function tolerancerepopulate!(eco::Ecosystem)
    eco.abundances = emptygridlandscape(eco.habitat, eco.spplist)
    eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun),
                                        length(eco.spplist.abun)))
    return tolerancepopulate!(eco.abundances, eco.spplist, eco.habitat,
                              eco.nichefit, eco.rngs)
end

"""
    emptypopulate!(ml::GridLandscape, spplist::SpeciesList, habitat::AB, nichefit::R,
                   rngs::Vector{Random.Xoshiro}) where {AB <: EcoSISTEM.AbstractHabitat, R <: EcoSISTEM.AbstractNicheFit}

Placeholder population function that leaves the landscape empty and warns.
"""
function emptypopulate!(ml::GridLandscape,
                        spplist::SpeciesList,
                        habitat::AB,
                        nichefit::R,
                        rngs::Vector{Random.Xoshiro}) where {AB <:
                                                             EcoSISTEM.AbstractHabitat,
                                                             R <:
                                                             EcoSISTEM.AbstractNicheFit}
    @warn "Ecosystem not populated!"
end
"""
    resupply!(eco::Ecosystem, supply::Union{Float64, Unitful.Quantity{Float64}}, grid::Tuple{Int64, Int64})
Refill an ecosystem `eco`, with resource from a supply value, `supply` and a grid
size.
"""
function resupply!(eco::Ecosystem,
                   supply::Union{Float64, Unitful.Quantity{Float64}},
                   grid::Tuple{Int64, Int64})
    return fill!(eco.habitat.supply.matrix, supply / (grid[1] * grid[2]))
end
function resupply!(eco::Ecosystem,
                   supply::Tuple{Unitful.Quantity{Float64},
                                 Unitful.Quantity{Float64}},
                   grid::Tuple{Int64, Int64})
    fill!(eco.habitat.supply.one.matrix, supply[1] / (grid[1] * grid[2]))
    return fill!(eco.habitat.supply.two.matrix, supply[2] / (grid[1] * grid[2]))
end
