# SPDX-License-Identifier: LGPL-3.0-or-later
#
# HOW ABUNDANCES CHANGE. Everything that writes to an ecosystem's abundance matrix: `populate!` and
# its relatives, which set it up, and `update!`, which is one timestep of the model — the demographic
# draws, then dispersal, then the layer changes.
#
# `update!` is what `simulate!` calls repeatedly; the order within it is load-bearing and stated at
# the function itself. Nothing here decides *what* a species needs or *where* it can live — that is
# the species requirements and the niche fit; this is only what those imply for the numbers.
#
# The MPI extension **duplicates** `update!`, `move!` and `populate!` rather than calling them, so
# a rule phrased as "the hot loop" has two places to apply and the second is in `ext/`.

using StatsBase

using LinearAlgebra

using Random

"""
    update!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update an ecosystem's abundances and environment for one timestep, with no intervention scheduled.

`eco` is the ecosystem to advance and `timestep` how far to advance it. Equivalent to
`update!(eco, timestep, nothing)`: the three-argument method is where the work happens, and each
concrete ecosystem supplies its own, so this one covers every kind — serial, distributed, and any
later addition.
"""
function update!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return update!(eco, timestep, nothing)
end

"""
    update!(eco::Ecosystem, timestep::Unitful.Time, intervention)

Update an ecosystem for one timestep, applying any scheduled [`Intervention`](@ref).

**The ordering is the point.** Interventions run *after* the population dynamics and *before* the
layer update, with the clock advanced between — so a [`SetChange`](@ref) installed this step takes
effect **this** step rather than one step late.
"""
function update!(eco::Ecosystem, timestep::Unitful.Time, intervention)

    # Calculate dimenions of regime and number of species
    numsc = countsubcommunities(eco.habitat.regime)
    numsp = size(eco.abundances.grid, 1)
    params = eco.spplist.params
    height = getgridshape(eco)[1]

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
    nblocks = cld(numsp, block)
    # :greedy hands the cache-line-sized species blocks to cores as they free up
    # (dynamic load balancing); blocks are independent so results are unchanged.
    Threads.@threads :greedy for b in 1:nblocks
        spstart = (b - 1) * block + 1
        spend = min(b * block, numsp)
        # Loop through grid squares
        for sc in 1:numsc
            # Convert 1D dimension to 2D coordinates
            (y, x) = convert_coords(eco, sc, height)
            # Check if grid cell currently active
            (eco.habitat.active[y, x] && (eco.cache.totaldemand[sc, 1] > 0)) ||
                continue
            for sp in spstart:spend
                rng = getrng(eco, sp)
                # Calculate how much birth and death should be adjusted
                adjusted_birth, adjusted_death = resource_adjustment(eco,
                                                                     eco.habitat.supply,
                                                                     sc, sp)

                # Both are per-individual rates over the timestep. Only the death one becomes
                # a probability: deaths are drawn per individual (Binomial), while births are a
                # count (Poisson) whose mean is the rate itself. `NoUnits` is needed on the birth
                # rate because it reaches `Poisson` as a bare number; the death rate is made
                # dimensionless by `exp`.
                birthrate = params.birth[sp] * timestep * adjusted_birth |>
                            NoUnits
                deathrate = params.death[sp] * timestep * adjusted_death

                deathprob = 1.0 - exp(-deathrate)

                (birthrate >= 0) & (deathprob >= 0) ||
                    error("Birth: $birthrate \n Death: $deathprob \n \n sc: $sc \n sp: $sp")
                # Calculate how many births and deaths
                births = rand(rng,
                              Poisson(eco.abundances.matrix[sp, sc] * birthrate))
                deaths = rand(rng,
                              Binomial(eco.abundances.matrix[sp, sc],
                                       deathprob))

                # Update population
                eco.abundances.matrix[sp, sc] += (births - deaths)

                # Calculate moves and write to cache
                move!(eco, eco.spplist.movement, sc, sp, eco.cache.netmigration,
                      births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.matrix .+= eco.cache.netmigration

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Advance the simulation clock before the layers change, so a layer sees the time it is
    # changing *to*.
    _advanceclock!(eco, timestep)

    # Interventions land here — after the dynamics, after the clock, before the layers — so a
    # `SetChange` installed now is applied by the very next line rather than a step later.
    applyinterventions!(eco, intervention, simulationtime(eco), timestep,
                        _stepnumber(eco, timestep))

    # Update environment - regime and resource supplies
    regimeupdate!(eco, timestep)
    return supplyupdate!(eco, timestep)
end

"""
    populate!(ml::GridLandscape, spplist::SpeciesList, habitat::AbstractHabitat,
              nichefit::AbstractNicheFit, rngs::Vector{Random.Xoshiro})
    populate!(ml::GridLandscape, spplist::SpeciesList,
              habitat::GridHabitat{H, <:LayerCollection{Resource}}, nichefit, rngs)

Populate the grid landscape `ml` by randomly scattering each species' total
abundance (taken from `spplist.abun`) across the grid cells, choosing each cell
with probability proportional to its available resource supply. Inactive cells are
given zero probability, so no individuals are placed outside the habitable
region. Each species is drawn from its own generator in `rngs`, so the result is
reproducible and independent of the number of threads or MPI processes.

`nichefit` is unused by these resource-based methods; it is accepted only so that they
share a signature with [`populate_by_tolerance!`](@ref) and can be passed
interchangeably as the population function when constructing an
[`Ecosystem`](@ref). For a multi-supply environment (a [`LayerCollection`](@ref)) the
sampling weight of a cell is the product of its two separately normalised
supplies.
"""
function populate!(ml::GridLandscape,
                   spplist::SpeciesList,
                   habitat::AB,
                   nichefit::NF,
                   rngs::Vector{Random.Xoshiro}) where {AB <: AbstractHabitat,
                                                        NF <:
                                                        AbstractNicheFit}
    grid, activity = _gridactivity(habitat)
    # Set up copy of supply
    b = reshape(parent(ustrip.(_getsupply(habitat.supply))), length(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    B = b ./ sum(b)
    # Loop through species, drawing from each species' own RNG stream. Three per-species things
    # walked in lockstep, so there is no index to keep: `eachrow` yields exactly the writable
    # `@view ml.matrix[i, :]` that `rand!` needs.
    for (rng, abun, row) in zip(rngs, spplist.abun, eachrow(ml.matrix))
        rand!(rng, Multinomial(abun, B), row)
    end
end

function populate!(ml::GridLandscape,
                   spplist::SpeciesList,
                   habitat::GridHabitat{H, <:LayerCollection{Resource}},
                   nichefit::NF,
                   rngs::Vector{Random.Xoshiro}) where {H <: AbstractRegime,
                                                        NF <:
                                                        AbstractNicheFit}
    # Calculate size of regime
    grid, activity = _gridactivity(habitat)
    # Set up a copy of each supply, zeroed outside the active cells and normalised to its own
    # total, then multiply them together so a cell needs every resource to be worth populating.
    fractions = _zipmap(values(habitat.supply)) do supply
        b = reshape(parent(copy(_getsupply(supply))), length(grid))
        b[.!activity] .= zero(eltype(b))
        return b ./ sum(b)
    end
    B = _fold(fractions) do f1, f2
        return f1 .* f2
    end
    # Loop through species, drawing from each species' own RNG stream
    for sp in eachindex(spplist.abun)
        rand!(rngs[sp], Multinomial(spplist.abun[sp], B ./ sum(B)),
              (@view ml.matrix[sp, :]))
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
    b = reshape(parent(copy(_getsupply(eco.habitat.supply))), length(grid))
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
    populate_by_tolerance!(ml::GridLandscape, spplist::SpeciesList, habitat::AbstractHabitat,
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
function populate_by_tolerance!(ml::GridLandscape,
                                spplist::SpeciesList,
                                habitat::AB,
                                nichefit::NF,
                                rngs::Vector{Random.Xoshiro}) where {AB <:
                                                                     AbstractHabitat,
                                                                     NF <:
                                                                     AbstractNicheFit}
    # Calculate size of regime
    dim = getgridshape(habitat)
    numsquares = dim[1] * dim[2]
    numspp = length(spplist.names)
    probabilities = [_suitability(habitat.regime, spplist.tolerance, nichefit,
                                  sc, sp)
                     for sc in 1:numsquares,
                         sp in 1:numspp]
    # Loop through species, drawing from each species' own RNG stream
    for sp in eachindex(spplist.abun)
        if spplist.native[sp]
            # Get abundance of species
            probs = probabilities[:, sp] ./ sum(probabilities[:, sp])
            probs[isnan.(probs)] .= 1 / numsquares
            abun = rand(rngs[sp], Multinomial(spplist.abun[sp], probs))
            # Add individual to this location
            ml.matrix[sp, :] .+= abun
        end
    end
end

"""
    repopulate_by_tolerance!(eco::Ecosystem)

Repopulate an ecosystem `eco` according to how well species tolerances match their
environment, redistributing the total abundance across species at random.
"""
function repopulate_by_tolerance!(eco::Ecosystem)
    eco.abundances = emptygridlandscape(eco.habitat, eco.spplist)
    eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun),
                                        length(eco.spplist.abun)))
    return populate_by_tolerance!(eco.abundances, eco.spplist, eco.habitat,
                                  eco.nichefit, eco.rngs)
end

"""
    emptypopulate!(ml::GridLandscape, spplist::SpeciesList, habitat::AB, nichefit::NF,
                   rngs::Vector{Random.Xoshiro}) where {AB <: EcoSISTEM.AbstractHabitat, NF <: EcoSISTEM.AbstractNicheFit}

Placeholder population function that leaves the landscape empty and warns.
"""
function emptypopulate!(ml::GridLandscape,
                        spplist::SpeciesList,
                        habitat::AB,
                        nichefit::NF,
                        rngs::Vector{Random.Xoshiro}) where {AB <:
                                                             EcoSISTEM.AbstractHabitat,
                                                             NF <:
                                                             EcoSISTEM.AbstractNicheFit}
    @warn "Ecosystem not populated!"
end

# **The first coordinate is the row**, `y`, matching the `(y, x)` order used throughout: the guard
# reads `dims[1]` from it and the returned pairs index `mat[n[1], n[2]]`. A caller that believes
# otherwise and passes `(x, y)` reads the transposed neighbourhood, which a square grid cannot
# distinguish and any other grid rejects with *"Coordinates outside grid"* — which is why the
# parameter names have to say which is which.
# Get the neighbours of a grid square in a matrix in 4 or 8 directions.
#
# The coordinates are `(y, x)` — row first, then column — as everywhere else in the package, and the
# rows of the returned matrix index `mat` in that same order (`mat[n[1], n[2]]`).
function _getneighbours(mat::Matrix, y_coord::Int64, x_coord::Int64,
                        chess::Int64 = 4)
    # Calculate dimensions
    dims = size(mat)
    y_coord <= dims[1] && x_coord <= dims[2] ||
        error("Coordinates outside grid")
    # Include 4 directions
    if chess == 4
        neighbour_vec = [y_coord x_coord-1
                         y_coord x_coord+1
                         y_coord-1 x_coord
                         y_coord+1 x_coord]
        # Include 8 directions
    elseif chess == 8
        neighbour_vec = [y_coord x_coord-1
                         y_coord x_coord+1
                         y_coord-1 x_coord
                         y_coord+1 x_coord
                         y_coord-1 x_coord-1
                         y_coord-1 x_coord+1
                         y_coord+1 x_coord-1
                         y_coord+1 x_coord+1]
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

# As the scalar method above, but accepting vectors of coordinates and returning the combined
# neighbours for all positions.
function _getneighbours(mat::Matrix,
                        y_coord::Vector{Int64},
                        x_coord::Vector{Int64},
                        chess::Int64 = 4)
    neighbours = map(n -> _getneighbours(mat, y_coord[n], x_coord[n], chess),
                     eachindex(y_coord))
    return vcat(neighbours...)
end

# Which step the clock has just completed. Derived from elapsed time rather than counted, so that
# it is the same on every MPI rank without communicating and survives a `CachedEcosystem` rebuilding
# an ecosystem at an arbitrary elapsed time.
function _stepnumber(eco::AbstractEcosystem, timestep::Unitful.Time)
    return round(Int, uconvert(NoUnits, simulationtime(eco) / timestep))
end

"""
    update_resource_usage!(eco::Ecosystem)
Calculate how much resource has been used up by the current species in each grid
square in the ecosystem, `eco`. This function is parameterised on whether the
species have one type of resource demand or two.
"""
function update_resource_usage!(eco::AbstractEcosystem{Part,
                                                       SpeciesList{TL, DM, MO,
                                                                   T,
                                                                   P}, NF}) where {Part,
                                                                                   MO,
                                                                                   T,
                                                                                   P,
                                                                                   NF,
                                                                                   TL,
                                                                                   DM <:
                                                                                   Demand}
    !eco.cache.valid || return true

    # Get resource supplies of species in square
    ϵ̄ = eco.spplist.demand.resource

    # Loop through grid squares
    Threads.@threads for sc in Base.axes(eco.abundances.matrix, 2)
        eco.cache.totaldemand[sc, 1] = ((@view eco.abundances.matrix[:, sc]) ⋅
                                        ϵ̄) *
                                       eco.spplist.demand.exchange_rate
    end
    return eco.cache.valid = true
end

function update_resource_usage!(eco::AbstractEcosystem{Part,
                                                       SpeciesList{TL, DM, MO,
                                                                   T,
                                                                   P}, NF}) where {Part,
                                                                                   MO,
                                                                                   T,
                                                                                   P,
                                                                                   NF,
                                                                                   TL,
                                                                                   DM <:
                                                                                   SpeciesRequirementCollection{Resource}}
    !eco.cache.valid || return true

    # Get resource supplies of species in square
    ds = values(eco.spplist.demand)

    # Loop through grid squares
    Threads.@threads for sc in Base.axes(eco.abundances.matrix, 2)
        _totaldemand!(eco.cache.totaldemand, sc,
                      (@view eco.abundances.matrix[:, sc]),
                      ds,
                      1)
    end
    return eco.cache.valid = true
end

# The three ways `totaldemand`'s per-resource columns are written, one column per member of a
# a demand collection. Each recurses over the demands rather than looping over them, so the walk is
# unrolled at compile time and allocates nothing whatever the number of resources. `k` is the
# column the head of `ds` writes to, so it tracks the recursion.
#
# `_totaldemand!` is the whole-cell write used by the shared-memory loop; MPI accumulates a cell's total
# across species blocks instead, so it needs the `_zerototaldemand!`/`_addtotaldemand!` pair.
_totaldemand!(totaldemand, sc, abun, ::Tuple{}, k::Int) = nothing

function _totaldemand!(totaldemand, sc, abun, ds::Tuple, k::Int)
    totaldemand[sc, k] = (abun ⋅ first(ds).resource) * first(ds).exchange_rate
    return _totaldemand!(totaldemand, sc, abun, Base.tail(ds), k + 1)
end

# Clear one cell's demand accumulator, one entry per resource. Written as a recursion over the
# demand tuple rather than a loop so the resource count is known at compile time and the whole thing
# unrolls: this runs once per cell per timestep, and must not allocate.
_zerototaldemand!(totaldemand, sc, ::Tuple{}, k::Int) = nothing

function _zerototaldemand!(totaldemand, sc, ds::Tuple, k::Int)
    totaldemand[sc, k] = 0.0
    return _zerototaldemand!(totaldemand, sc, Base.tail(ds), k + 1)
end

# Add a block of species' demand into one cell's accumulator, unrolled over the resources exactly as
# `_zerototaldemand!` is. `from`/`to` are the block this process owns, which is the whole range
# serially and one rank's share under MPI.
function _addtotaldemand!(totaldemand, sc, abun, from::Int, to::Int, ::Tuple{},
                          k::Int)
    return nothing
end

function _addtotaldemand!(totaldemand, sc, abun, from::Int, to::Int, ds::Tuple,
                          k::Int)
    ϵ̄ = @view first(ds).resource[from:to]
    totaldemand[sc, k] += (abun ⋅ ϵ̄) * first(ds).exchange_rate
    return _addtotaldemand!(totaldemand, sc, abun, from, to, Base.tail(ds),
                            k + 1)
end

"""
    resource_adjustment(eco::Ecosystem, supply::AbstractSupply, sc::Int64, sp::Int64)

Calculate how much birth and death rates should be adjusted by, according to how
much resource is available, `supply`, in the grid square, `sc`, and how much resource
the species, `sp`, requires.
"""
function resource_adjustment(eco::AbstractEcosystem, supply::AbstractSupply,
                             sc::Int64, sp::Int64)
    return _resourceadjustmentbytype(eco.spplist.params, eco, supply, sc, sp)
end

# NoGrowth freezes the population; anything else adjusts birth/death rates by the available
# resource — dispatched on `params`'s type rather than an `isa` branch in `resource_adjustment`.
_resourceadjustmentbytype(::NoGrowth, eco, supply, sc, sp) = (0.0, 0.0)

function _resourceadjustmentbytype(::AbstractParams, eco, supply, sc, sp)
    return _resourceadjustment(eco, supply, sc, sp)
end

# Birth and death rate multipliers for a single-demand environment. Weighs
# the species' own resource demand (`ϵ̄`) and how well its tolerances match the cell
# (`ϵ̄real`) against the resource available in the cell (`K`) relative to the total
# demand there (`E`): births are boosted when resource is plentiful (`K/E`, capped at
# `params.boost`) and deaths rise as demand approaches the supply (`E/K`). Called
# only for growing populations — [`resource_adjustment`](@ref) short-circuits NoGrowth.
function _resourceadjustment(eco::AbstractEcosystem, supply::AbstractSupply,
                             sc::Int64, sp::Int64)
    params = eco.spplist.params
    height = getgridshape(eco)[1]
    (y, x) = convert_coords(eco, sc, height)
    K = _getsupply(eco.habitat.supply)[y, x] * eco.spplist.demand.exchange_rate
    # Get resource supplies of species in square
    ϵ̄ = eco.spplist.demand.resource[sp] *
        eco.spplist.demand.exchange_rate
    E = eco.cache.totaldemand[sc, 1]
    # Traits
    ϵ̄real = 1 / suitability(eco, sc, sp)
    # Alter rates by resource available in current pop & own demands
    birth_resource = ϵ̄^-params.longevity * ϵ̄real^-params.survival *
                     min(K / E, params.boost)
    death_resource = ϵ̄^-params.longevity * ϵ̄real^params.survival * (E / K)
    return birth_resource, death_resource
end

# As above but for a multi-resource environment (e.g. solar resource and water), combining the
# supplies. The species is limited by whichever resource is scarcest: births use the `min` of the
# availability ratios (`K/E`, still capped at `params.boost`) and deaths the `max` of the demand
# ratios (`E/K`), so every demand must be met for the population to grow. Per-resource quantities
# are built and combined by compile-time-unrolled folds, so this stays allocation-free at any arity.
#
# The suitability term enters **once**, whatever the arity, which is what makes this reduce to the
# single-supply method above rather than merely resembling it.
#
# Raising it to the number of supplies instead would be wrong three ways over, and each is worth
# knowing because the arithmetic looks plausible:
#
#   * **a resource that cannot limit anything would change the answer.** Add a second resource with
#     infinite supply and identical demand for every species, and birth and death move by
#     `suit^∓survival` — from an input that by construction cannot matter.
#   * **it would index one collection by another's length.** Suitability is a **regime** quantity;
#     the supply count is a **resource** one, and the two have had no arity relationship since layer
#     collections gained members. Two regimes over one supply would give an exponent of 1, and one
#     regime over two supplies an exponent of 2.
#   * **the regime side is already combined.** `nichefitcombine` folds the per-layer suitabilities
#     before this ever sees them, so there is nothing left here to combine.
#
# `demanded` below is still a **product** over resources, and so fails the first of those tests. That
# is deliberate rather than overlooked: unlike suitability it has no single obviously right
# replacement — the limiting resource's demand, a sum, a mean? — so it is a separate question with
# its own answer to find.
function _resourceadjustment(eco::AbstractEcosystem,
                             supply::LayerCollection{Resource},
                             sc::Int64,
                             sp::Int64)
    height = getgridshape(eco)[1]
    (y, x) = convert_coords(eco, sc, height)
    params = eco.spplist.params
    ds = values(eco.spplist.demand)
    bs = values(supply)
    _samearity(ds, bs)
    # Resource available in the cell, the species' own demand, and the cell's total demand —
    # one entry per resource, in the collections' shared order.
    K = _zipmap(bs, ds) do b, d
        return _getsupply(b)[y, x] * d.exchange_rate
    end
    ϵ̄ = _zipmap(ds) do d
        return d.resource[sp] * d.exchange_rate
    end
    E = ntuple(k -> eco.cache.totaldemand[sc, k], Val(length(ds)))
    # Once, not once per resource — see the note above.
    ϵ̄real = 1 / suitability(eco, sc, sp)
    # Alter rates by resource available in current pop & own demands
    demanded = _fold(*, ϵ̄)
    birth_resource = demanded^-params.longevity * ϵ̄real^-params.survival *
                     min(_fold(min, _zipmap(/, K, E)), params.boost)
    death_resource = demanded^-params.longevity * ϵ̄real^params.survival *
                     _fold(max, _zipmap(/, E, K))
    return birth_resource, death_resource
end

"""
    convert_coords(eco, sc::Int64, height::Int64)
    convert_coords(eco, y::Int64, x::Int64, height::Int64)
Convert coordinates from two-dimensional (`y`,`x`) format to one dimension
(`sc`), or vice versa, using the `height` (dimension 1) of the grid. This
function can also be applied to arrays of coordinates.
"""
function convert_coords(eco::AbstractEcosystem,
                        sc::Int64,
                        height::Int64 = getgridshape(eco)[1])
    y = ((sc - 1) % height) + 1
    x = div((sc - 1), height) + 1
    return (y, x)
end

function convert_coords(eco::AbstractEcosystem,
                        pos::Tuple{Int64, Int64},
                        height::Int64 = getgridshape(eco)[1])
    sc = pos[1] + height * (pos[2] - 1)
    return sc
end

function convert_coords(sc::Int64, height::Int64)
    y = ((sc - 1) % height) + 1
    x = div((sc - 1), height) + 1
    return (y, x)
end

function convert_coords(y::Int64, x::Int64, height::Int64)
    sc = y + height * (x - 1)
    return sc
end

# Normalise the per-destination weights `lookup.pnew` (already zeroed for every unreachable or
# inactive destination by the caller) and draw `abun` moves from them. Shared tail of all three
# `calc_lookup_moves!` boundary methods.
#
# If *no* destination is reachable the weights sum to zero, and normalising would give `NaN`s and a
# throw from `Multinomial`. That case is defined as "no move this step": `moves` is zeroed and the
# species' RNG stream is left untouched, so a stranded cell costs no draws and cannot re-phase later
# ones. Defensive only — it is not actually reachable, because staying put is one of the kernel's
# own destinations and a cell only disperses while it is itself active (`update!`'s gate above), so
# the self-offset's weight is always non-zero.
#
# `disperse_safely` is the boundary's flag for what happens to an individual aimed at a dead cell.
# When `true`, the weights are renormalised and it is redistributed among the reachable destinations.
# When `false` it is lost instead: `lost` is the kernel weight that was blocked, so the survivors are
# drawn binomially against `total / (total + lost)` and only they are then distributed. Drawing
# the survivor *count* first is what keeps this allocation-free; extending `pnew` with a "lost"
# category would allocate a vector per call, in the hot loop.
#
# The `iszero(lost)` guard is **not** an optimisation — it is what makes the flag a no-op when
# there is nothing to lose. Without it the extra `Binomial` draw would consume from the species' RNG
# stream and re-phase every later draw, so a `Torus` with no inactive cells would give different
# results with the flag on and off despite nothing being lost either way. Measured before the guard
# existed: 5883 against 5805 on an identical landscape and seed.
function _drawmoves!(lookup::Lookup, sp::Int64, eco::AbstractEcosystem,
                     abun::Int64, disperse_safely::Bool = true,
                     lost::Float64 = 0.0)
    total = sum(lookup.pnew)
    if iszero(total)
        fill!(lookup.moves, 0)
        return lookup.moves
    end
    surviving = (disperse_safely || iszero(lost)) ? abun :
                rand(getrng(eco, sp), Binomial(abun, total / (total + lost)))
    lookup.pnew ./= total
    dist = Multinomial(surviving, lookup.pnew)
    return rand!(getrng(eco, sp), dist, lookup.moves)
end

# Resolve one coordinate against its axis's boundary condition. Two tiny methods, dispatched on the
# condition, so the rule is stated once per axis rather than branched on — which is what lets the
# three near-identical `calc_lookup_moves!` methods this replaces become one.
#
# `Bounded` returns `nothing` for a step off the edge, and the caller reads that as "dead cell";
# `Periodic` always lands somewhere, wrapping with `mod`.
function _stepto(::Type{Bounded}, here::Int64, delta::Int64, extent::Int64)
    return (0 < here + delta <= extent) ? here + delta : nothing
end

function _stepto(::Type{Periodic}, here::Int64, delta::Int64, extent::Int64)
    return mod(here + delta - 1, extent) + 1
end

"""
    calc_lookup_moves!(topology::EdgeTopology, y, x, sp, eco, abun)

Calculate the number of moves taken by a species, `sp`, from a specific grid
square location (`y`, `x`). `topology` says how the grid's edges join — one
[`AbstractBoundaryCondition`](@ref) per axis — and the total abundance of
individuals to move is `abun`, which may be the number of births in the timestep
or the total individuals.

**One method for all four topologies.** They differ only in whether each axis wraps, which is the
two coordinate rules in `_stepto`, so the fourth combination — Y wrapping, X bounded — comes free
rather than as a fourth copy.
"""
function calc_lookup_moves!(::EdgeTopology{BCY, BCX},
                            y::Int64,
                            x::Int64,
                            sp::Int64,
                            eco::AbstractEcosystem,
                            abun::Int64) where {BCY, BCX}
    lookup = getlookup(eco, sp)
    maxY, maxX = getgridshape(eco)
    # `lost` accumulates the kernel weight aimed at dead cells, which `_drawmoves!` needs exactly
    # (not inferred from a float comparison) to tell "nothing was blocked" from "a rounding step
    # short of nothing".
    lost = 0.0
    for i in eachindex(lookup.y)
        newy = _stepto(BCY, y, lookup.y[i], maxY)
        newx = _stepto(BCX, x, lookup.x[i], maxX)
        valid = !isnothing(newy) && !isnothing(newx) &&
                eco.habitat.active[newy, newx]

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
        lost += valid ? 0.0 : lookup.p[i]
    end
    return _drawmoves!(lookup, sp, eco, abun,
                       dispersesafely(eco.spplist.movement, sp), lost)
end

"""
    move!(eco::Ecosystem, ::AbstractMovement, sc::Int64, sp::Int64, grd::Matrix{Int64}, abun::Int64)

Calculate the movement of species `sp` from a given position in the landscape
`sc`, using the lookup table found in the [`Ecosystem`](@ref) and updating the
movement patterns on a cached grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process,
instead of the entire population
"""
function move!(eco::AbstractEcosystem,
               ::AlwaysMovement,
               sc::Int64,
               sp::Int64,
               grd::Matrix{Int64},
               ::Int64)
    # "Animal-like": the whole current population disperses
    return _move!(eco, sc, sp, grd, eco.abundances.matrix[sp, sc])
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
               sc::Int64,
               sp::Int64,
               grd::Matrix{Int64},
               births::Int64)
    # "Plant-like": only the newly born individuals disperse
    return _move!(eco, sc, sp, grd, births)
end

# Common dispersal code: move `amount` individuals of species `sp` from
# position `sc`, scattering them across the landscape according to the
# species' lookup table.
function _move!(eco::AbstractEcosystem,
                sc::Int64,
                sp::Int64,
                grd::Matrix{Int64},
                amount::Int64)
    height, width = getgridshape(eco)
    (y, x) = convert_coords(eco, sc, height)
    lookup = getlookup(eco, sp)
    calc_lookup_moves!(eco.habitat.topology, y, x, sp, eco, amount)
    # Lose moves from current grid square
    grd[sp, sc] -= amount
    # Map moves to location in grid
    moves = lookup.moves
    for i in eachindex(lookup.y)
        newy = mod(lookup.y[i] + y - 1, height) + 1
        newx = mod(lookup.x[i] + x - 1, width) + 1
        loc = convert_coords(eco, (newy, newx), height)
        grd[sp, loc] += moves[i]
    end
    return eco
end

# Return the two ingredients the population routines share when spreading
# individuals across the grid:
# `grid`     - a vector of the linear indices `1:ncells` of every grid cell,
# used both to size/flatten the supply and to sample cell locations;
# `activity` - a flattened copy of `habitat.active`, the boolean mask of which cells
# are habitable. Callers zero the supply of inactive (`false`) cells so
# that no individuals are ever placed outside the active region.
function _gridactivity(habitat::AbstractHabitat)
    dim = getgridshape(habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # `reshape` on a `DimArray` isn't a clean zero-copy operation the way a plain `Array`
    # reshape is — unwrap via `parent()` first.
    activity = reshape(copy(parent(habitat.active)), len)
    return grid, activity
end
