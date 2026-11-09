# SPDX-License-Identifier: LGPL-3.0-or-later
#
# HOW ABUNDANCES CHANGE. Everything that writes to an ecosystem's abundance matrix: placing its
# individuals, at build and by `repopulate!`, and `update!`, which is one timestep of the model - the
# demographic draws, then dispersal, then the layer changes.
#
# `update!` is what `simulate!` calls repeatedly; the order within it is load-bearing and stated at
# the function itself. Nothing here decides *what* a species needs or *where* it can live - that is
# the species requirements and the niche fit; this is only what those imply for the numbers.
#
# The MPI extension **duplicates** `update!` and `update_resource_usage!` rather than calling them, so
# a rule phrased as "the hot loop" has two places to apply and the second is in `ext/`.

using StatsBase

using LinearAlgebra

using KahanSummation: sum_kbn

using Random

# == Functions ==================================================================================

"""
    repopulate!(eco::AbstractEcosystem, placement = eco.placement; seed = eco.seed,
                verbosity = :normal)

Place an ecosystem's individuals again from scratch: its abundances are cleared, and `placement`
places them over its active cells and layers as they now are, as [`build_ecosystem`](@ref) placed
them. The ecosystem keeps `placement` as its own, so a later `repopulate!(eco)` repeats it. The
placement draws from streams of its own made from `seed`, so on an unchanged habitat, with the
ecosystem's own placement and seed, it places what its build placed; the clock, and the streams the
dynamics draw from, are left as they are. Under MPI with more than one rank this is a collective
call.

# Arguments

  - `eco`: the ecosystem.
  - `placement`: any placement [`build_ecosystem`](@ref) takes; one reading a table that can be read
    only once is refused.
  - `seed`: the seed of the placement's streams.
  - `verbosity`: `:silent` reports only what the placement could not place; anything else reports
    its notices too.
"""
function repopulate!(eco::AbstractEcosystem, placement = eco.placement;
                     seed::Integer = eco.seed, verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    isnothing(placement) &&
        error("this ecosystem records no placement to place by again; give one, as " *
              "`repopulate!(eco, BySupply())`.")
    _checkplacement(placement)
    _reusable(placement) ||
        error("$(_placementphrase(placement)) reads a table that can be read only once, so it " *
              "cannot place again; give `FromTable` a path, a record or a table of columns.")
    context = _placementcontext(eco, seed)
    _placeinto!(eco, placement, context, verbosity, clear = true)
    invalidatecaches!(eco)
    return eco
end

# The context of a placement in `eco`: its active cells, a placement stream per species and the
# community stream, both made from `seed`, and its nichefit, epoch and calendar, at its start.
function _placementcontext(eco::AbstractEcosystem, seed::Integer)
    streams = _placementrngs(seed, eco.spplist)
    return PlacementContext(cells = findall(parent(eco.habitat.active)),
                            streams = Dict(zip(eco.spplist.names, streams)),
                            community = _abundancerng(seed),
                            nichefit = eco.nichefit, epoch = eco.epoch,
                            calendar = eco.calendar)
end

# Place `placement`'s individuals in `eco`: called and checked first, its problems reported, then -
# its abundances cleared, where `clear` - written, and its records kept, the placement too where
# `keep` (as at build, not at an intervention's firing). A dated table is refused where `undated`.
function _placeinto!(eco::AbstractEcosystem, placement,
                     context::PlacementContext, verbosity::Symbol; clear::Bool,
                     keep::Bool = true, undated::Bool = false)
    outcome = _resolveplacement(eco, placement, context, verbosity,
                                undated = undated)
    _reportproblems(vcat(outcome.placed.problems, outcome.resolved.problems),
                    verbosity)
    clear && (_ownedabundances(eco).rows .= 0)
    _writeplaced!(eco, outcome.resolved)
    keep && (eco.placement = placement)
    eco.inputs = _uniqueinputs(vcat(eco.inputs, outcome.placed.records))
    return eco
end

# What `placement` places in `eco`, called and checked but not yet written: its result, its table in
# one form with its species matched to the species list, and its rows resolved against the grid.
# Under MPI the root alone calls it and matches its names - so only the root reads a name authority -
# once every rank has built the layers it weights by together, and every rank takes the root's result.
function _resolveplacement(eco::AbstractEcosystem, placement,
                           context::PlacementContext, verbosity::Symbol;
                           undated::Bool = false)
    prepared = _prepared(placement, eco.habitat)
    placed = _sharedresult(eco, context) do
        raw = _callplacement(prepared, eco.spplist, eco.habitat, context)
        undated && :time in Tables.columnnames(Tables.columns(raw.table)) &&
            error("`AddAbundance` adds this placement's individuals at each firing, and its table is " *
                  "dated: a table served by its dates is a `FromTable` with a `time` column or a " *
                  "`FromOccurrences` given `date` or `time`, alone in an `Intervention` with " *
                  "`EveryStep()`.")
        matched = _matchedtable(_canonicaltable(raw.table, eco.habitat,
                                                context),
                                eco.spplist, context, verbosity)
        return PlacedAbundances(matched.table, records = raw.records,
                                problems = vcat(raw.problems, matched.problems))
    end
    return (placed = placed,
            resolved = _resolveplaced(placed, eco.habitat, context))
end

# A placement's result, serially; under MPI the extension's method has the root build it and every
# rank take it.
_sharedresult(build, ::AbstractEcosystem, ::PlacementContext) = build()

# Placed individuals added to an ecosystem's abundances, each rank its own species'.
function _writeplaced!(eco::AbstractEcosystem, resolved::NamedTuple)
    owned = _ownedabundances(eco)
    held = axes(owned.rows, 1)
    for (sp, column, n) in zip(resolved.species, resolved.cells,
                               resolved.counts)
        row = sp - owned.firstspecies + 1
        row in held || continue
        owned.rows[row, column] += n
    end
    return _rowswritten!(eco)
end

# What follows writing the rows of an ecosystem's abundances: nothing, serially; under MPI the
# extension's method rebuilds the columns from them.
_rowswritten!(::AbstractEcosystem) = nothing

"""
    update!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update an ecosystem's abundances and environment for one timestep, with no intervention scheduled.

`eco` is the ecosystem to advance and `timestep` how far to advance it. Equivalent to
`update!(eco, timestep, nothing)`.
"""
function update!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return update!(eco, timestep, nothing)
end

"""
    update!(eco::AbstractEcosystem, timestep::Unitful.Time, intervention)

Update an ecosystem for one timestep, applying any scheduled [`Intervention`](@ref).

**The ordering is the point.** Interventions run *after* the population dynamics and *before* the
layer update, with the clock advanced between - so a [`SetChange`](@ref) installed this step takes
effect **this** step rather than one step late. On the first step, with the clock at zero, those
due at the start run once more *before* the dynamics, so they act on the starting state.

`eco` is the ecosystem to advance, serial or distributed, `timestep` how far to advance it, and
`intervention` an [`Intervention`](@ref), an [`InterventionSet`](@ref) or `nothing`.
"""
function update!(eco::AbstractEcosystem, timestep::Unitful.Time, intervention)
    iszero(simulationtime(eco)) && _startrun!(eco, intervention, timestep)
    _resumerun!(eco, intervention, timestep)
    return _step!(eco, timestep, intervention)
end

# The start of a run: what is due at elapsed zero acts on the starting state, before the first
# step's dynamics, which then read the abundances it wrote. A driver that observes the starting state
# calls this itself and then `_step!`, so the start is applied once.
function _startrun!(eco::Ecosystem, intervention, timestep::Unitful.Time)
    applyinterventions!(eco, intervention, simulationtime(eco), timestep, 0)
    invalidatecaches!(eco)
    return eco
end

# The first step of a run started from a saved state, after anything due at the start: a state keeps
# a layer's values but not its change rule, so a run continuing one installs again the rules its
# interventions would have installed by then, and is warned of any layer whose rule still differs from
# the one it was saved with. A run restarted at zero has its rules installed as it reaches them. Under
# MPI every rank installs the same rules, as every rank holds the layers.
function _resumerun!(eco::AbstractEcosystem, intervention,
                     timestep::Unitful.Time)
    saved = eco.savedchanges
    isnothing(saved) && return eco
    iszero(simulationtime(eco)) ||
        _replaychanges!(eco, intervention, timestep)
    _isrootrank(eco) && _warnchanges(saved, eco.habitat)
    eco.savedchanges = nothing
    return eco
end

# One timestep of an ecosystem: the dynamics, then the clock, the interventions due and the layers.
function _step!(eco::Ecosystem, timestep::Unitful.Time, intervention)

    # Calculate dimenions of regime and number of species
    numsc = countsubcommunities(eco.habitat.regime)
    numsp = size(eco.abundances.grid, 1)
    params = eco.spplist.params
    height = getgridshape(eco)[1]

    # Set the overall resource supply of that square
    update_resource_usage!(eco)

    # Loop through species in cache-line-sized contiguous blocks (see
    # `species_blocksize`): each thread owns whole blocks, and the cell loop sits
    # outside the inner species loop so a block's species - adjacent rows of the
    # column-major (species, cells) matrix - are touched as one cache line. The
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
                adjusted_birth,
                adjusted_death = resource_adjustment(eco,
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
                move!(eco, eco.spplist.movement[sp], sc, sp,
                      eco.cache.netmigration, births)
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

    # Interventions land here - after the dynamics, after the clock, before the layers - so a
    # `SetChange` installed now is applied by the very next line rather than a step later.
    applyinterventions!(eco, intervention, simulationtime(eco), timestep,
                        _stepnumber(eco, timestep))

    # Update environment - regime and resource supplies
    regimeupdate!(eco, timestep)
    return supplyupdate!(eco, timestep)
end

# **The first coordinate is the row**, `y`, matching the `(y, x)` order used throughout: the guard
# reads `dims[1]` from it and the returned pairs index `mat[n[1], n[2]]`. A caller that believes
# otherwise and passes `(x, y)` reads the transposed neighbourhood, which a square grid cannot
# distinguish and any other grid rejects with *"Coordinates outside grid"* - which is why the
# parameter names have to say which is which.
# Get the neighbours of a grid square in a matrix in 4 or 8 directions.
#
# The coordinates are `(y, x)` - row first, then column - as everywhere else in the package, and the
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
                                                                                   AbstractDemand}
    !eco.cache.valid || return true

    # Get resource supplies of species in square
    ds = _demandtuple(eco.spplist.demand)

    # Loop through grid squares
    Threads.@threads for sc in Base.axes(eco.abundances.matrix, 2)
        _totaldemand!(eco.cache.totaldemand, sc,
                      (@view eco.abundances.matrix[:, sc]),
                      ds,
                      1)
    end
    return eco.cache.valid = true
end

# One demand or several, as a tuple, so everything downstream is written once against the
# many-resource case and the single-resource case is simply its arity-1 instance. Measured rather
# than assumed: routing a lone demand through `_totaldemand!` is *cheaper* than the direct dot
# product it replaced -- 18.1 ns against 21.2 ns per cell, and fewer bytes -- because the recursion
# unrolls at compile time and a one-member tuple costs nothing to build.
_demandtuple(demand::Demand) = (demand,)

_demandtuple(demands::SpeciesRequirementCollection{Resource}) = values(demands)

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
# resource - dispatched on `params`'s type rather than an `isa` branch in `resource_adjustment`.
_resourceadjustmentbytype(::NoGrowth, eco, supply, sc, sp) = (0.0, 0.0)

function _resourceadjustmentbytype(::AbstractParams, eco, supply, sc, sp)
    return _resourceadjustment(eco, supply, sc, sp)
end

# Birth and death rate multipliers for a single-demand environment. Weighs
# the species' own resource demand (`ϵ̄`) and how well its tolerances match the cell
# (`ϵ̄real`) against the resource available in the cell (`K`) relative to the total
# demand there (`E`): births are boosted when resource is plentiful (`K/E`, capped at 1, so a
# species reproduces no faster than its baseline rate however plentiful the resource) and deaths
# rise as demand approaches the supply (`E/K`). Called
# only for growing populations - [`resource_adjustment`](@ref) short-circuits NoGrowth.
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
    birth_resource = ϵ̄^-params.longevity[sp] * ϵ̄real^-params.survival *
                     min(K / E, 1.0)
    death_resource = ϵ̄^-params.longevity[sp] * ϵ̄real^params.survival * (E / K)
    return birth_resource, death_resource
end

# As above but for a multi-resource environment (e.g. solar resource and water), combining the
# supplies. The species is limited by whichever resource is scarcest: births use the `min` of the
# availability ratios (`K/E`, still capped at 1) and deaths the `max` of the demand
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
#     `suit^∓survival` - from an input that by construction cannot matter.
#   * **it would index one collection by another's length.** Suitability is a **regime** quantity;
#     the supply count is a **resource** one, and the two have had no arity relationship since layer
#     collections gained members. Two regimes over one supply would give an exponent of 1, and one
#     regime over two supplies an exponent of 2.
#   * **the regime side is already combined.** `nichefitcombine` folds the per-layer suitabilities
#     before this ever sees them, so there is nothing left here to combine.
#
# `demanded` below is still a **product** over resources, and so fails the first of those tests. That
# is deliberate rather than overlooked: unlike suitability it has no single obviously right
# replacement - the limiting resource's demand, a sum, a mean? - so it is a separate question with
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
    # Resource available in the cell, the species' own demand, and the cell's total demand -
    # one entry per resource, in the collections' shared order.
    K = _zipmap(bs, ds) do b, d
        return _getsupply(b)[y, x] * d.exchange_rate
    end
    ϵ̄ = _zipmap(ds) do d
        return d.resource[sp] * d.exchange_rate
    end
    E = ntuple(k -> eco.cache.totaldemand[sc, k], Val(length(ds)))
    # Once, not once per resource - see the note above.
    ϵ̄real = 1 / suitability(eco, sc, sp)
    # Alter rates by resource available in current pop & own demands
    demanded = _fold(*, ϵ̄)
    birth_resource = demanded^-params.longevity[sp] * ϵ̄real^-params.survival *
                     min(_fold(min, _zipmap(/, K, E)), 1.0)
    death_resource = demanded^-params.longevity[sp] * ϵ̄real^params.survival *
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
# ones. Defensive only - it is not actually reachable, because staying put is one of the kernel's
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
# The `iszero(lost)` guard is **not** an optimisation - it is what makes the flag a no-op when
# there is nothing to lose. Without it the extra `Binomial` draw would consume from the species' RNG
# stream and re-phase every later draw, so a `Torus` with no inactive cells would give different
# results with the flag on and off despite nothing being lost either way. Measured before the guard
# existed: 5883 against 5805 on an identical landscape and seed.
function _drawmoves!(lookup::Lookup, sp::Int64, eco::AbstractEcosystem,
                     abun::Int64, disperse_safely::Bool = true,
                     lost::Float64 = 0.0)
    # A compensated sum, for the correctly rounded total: `sum` may reassociate once vectorised,
    # and this sum affects reproducibility. For speed, we should revert to sum().
    total = sum_kbn(lookup.pnew)
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
# condition, so the rule is stated once per axis rather than branched on - which is what lets the
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
square location (`y`, `x`). `topology` says how the grid's edges join - one
[`AbstractBoundaryCondition`](@ref) per axis - and the total abundance of
individuals to move is `abun`, which may be the number of births in the timestep
or the total individuals.

**One method for all four topologies.** They differ only in whether each axis wraps, which is the
two coordinate rules in `_stepto`, so the fourth combination - Y wrapping, X bounded - comes free
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
    move!(eco::AbstractEcosystem, movement::SpeciesDispersal, sc::Int64, sp::Int64,
          grd::Matrix{Int64}, births::Int64)

Calculate the movement of species `sp` from a given position in the landscape
`sc`, using the lookup table found in the [`Ecosystem`](@ref) and updating the
movement patterns on a cached grid, `grd`. The number of births are passed in
as `births` so that the movement can be restricted to only the newly born
individuals. `movement` is the species' own, and its mode says which
individuals disperse: only the `births` under [`SeedDispersal`](@ref) and the
other modes dispersing at birth, the whole standing population under
[`AlwaysDispersal`](@ref), and none under [`NoDispersal`](@ref).
"""
function move!(eco::AbstractEcosystem,
               ::SpeciesDispersal{<:AlwaysDispersal},
               sc::Int64,
               sp::Int64,
               grd::Matrix{Int64},
               ::Int64)
    # "Animal-like": the whole current population disperses
    return _move!(eco, sc, sp, grd, _standingpopulation(eco, sc, sp))
end

function move!(eco::AbstractEcosystem,
               ::SpeciesDispersal{NoDispersal},
               ::Int64,
               ::Int64,
               ::Matrix{Int64},
               ::Int64)
    return eco
end

function move!(eco::AbstractEcosystem,
               ::SpeciesDispersal{<:_BirthDispersal},
               sc::Int64,
               sp::Int64,
               grd::Matrix{Int64},
               births::Int64)
    # "Plant-like": only the newly born individuals disperse
    return _move!(eco, sc, sp, grd, births)
end

# Which row of the landscape's abundance and migration matrices belongs to species `sp`. Serially
# that is `sp` itself; a distributed run holds only its own block of species, so the MPI extension
# maps the global index onto this rank's local row. Written as a hook rather than a branch so that
# `_move!` below is ONE body serving both, instead of two that have to be kept in step -- which is
# the failure this whole file's header warns about.
_landscaperow(::AbstractEcosystem, sp::Int64) = sp

# How many individuals of species `sp` are standing in cell `sc` - what `AlwaysDispersal` disperses
# in full, as against the newborns a mode dispersing at birth is handed directly. Serially the whole
# abundance matrix is to hand; a distributed run holds only its own block of species, so the MPI
# extension reads the same count from its rank-local rows. A hook rather than a branch, for the same
# reason `_landscaperow` above is one: `move!` stays a single body serving both sides.
function _standingpopulation(eco::AbstractEcosystem, sc::Int64, sp::Int64)
    return eco.abundances.matrix[sp, sc]
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
    # `grd` is indexed by the row this process owns for `sp`, which is `sp` serially and a
    # rank-local row under MPI; `lookup` and `calc_lookup_moves!` take the global index either way.
    row = _landscaperow(eco, sp)
    # Lose moves from current grid square
    grd[row, sc] -= amount
    # Map moves to location in grid
    moves = lookup.moves
    for i in eachindex(lookup.y)
        newy = mod(lookup.y[i] + y - 1, height) + 1
        newx = mod(lookup.x[i] + x - 1, width) + 1
        loc = convert_coords(eco, (newy, newx), height)
        grd[row, loc] += moves[i]
    end
    return eco
end
