# SPDX-License-Identifier: LGPL-3.0-or-later
#
# HOW TO RUN AN EXPERIMENT. This file is the workflow and nothing else: the verbs a user calls, in the
# order they call them. It is the last of the model's own files - only `deprecations.jl` loads
# after it - and, like `Ecology.jl`, which is the first, it should be
# readable on its own by someone who never opens another.
#
# `Ecology.jl` says what the model IS. This says what you DO with it:
#
#     investigate_study_area   decide the grid, and argue with it, before building anything on it
#     build_habitat            the environment: what each place is like, and what it provides
#     build_species            the species: what they tolerate, what they need, how they move
#     build_ecosystem          put the species in the environment and check the two agree
#     simulate!                run it
#     simulate!(recorder, ...) run it and keep abundances or a diversity measure on a schedule, with
#                              a recorder from `Recorder.jl`
#     simulate!(f, ...)        run it and call back on a schedule, for anything else
#     generate_storage         allocate the array a recorder writes into
#
# WHAT IS NOT HERE, and the test that decides it: these are the verbs a USER calls. Everything the
# package calls on their behalf stays with its machinery. So `update!` is absent - `simulate!` calls
# it, and it is `public` rather than exported for that reason - and so is every private helper,
# including the defaults behind the builders
# (`DefaultEcosystem.jl`) and the grid-decision machinery behind `investigate_study_area`
# (`StudyArea.jl`).
#
# The `populate!` family is the exception, and it is absent for the opposite reason: those five
# share one signature so that any of them can be passed to `Ecosystem(popfun, ...)`, so they are the
# shipped values of a user-selectable strategy rather than verbs in their own right. They stay
# exported, and they live with the abundance machinery they are values for.
#
# THE STRICT COUNTERPARTS. Each builder here is the forgiving face of a constructor that refuses to
# guess, and those live with their own types:
#
#     StudyArea(...)      `StudyArea.jl`      decides a grid from the layers it is given
#     GridHabitat(...)    `GridHabitat.jl`    one construction route, and no defaulted inputs
#     SpeciesList(...)    `Species.jl`        every per-species field stated explicitly
#     Ecosystem(...)      `Ecosystem.jl`      species and environment, already checked
#
# The contrast is deliberate rather than duplication. A constructor cannot supply a default, because
# no default is right for every model; a builder can, and announces every value it chose. Reach for
# the constructor when you know what you want, and the builder when you want to be told.

using Unitful
using DimensionalData: Ti, At
using Unitful.DefaultSymbols
using Diversity
using EcoSISTEM.Units

# == Deciding the grid ==========================================================================

# `Missing` and `Nothing` are both in the unions and mean **different** things: `missing` is "not
# given - derive it, or inherit it from `base`", `nothing` is "clear what would have been inherited".
# Neither is a layer, which is why each has to be admitted explicitly rather than falling out of
# [`LayerInput`](@ref).
#
# **An already-built layer must NOT be accepted here, in either role** - unlike `GridHabitat`,
# which takes a built `AbstractSupply`. The difference is what the function does with it: this pair
# *decides the grid from its layers*, so a layer must be able to answer `_shapesgrid`/`_layerfacts`
# (extent, CRS, cell size), and a built layer cannot. `GridHabitat` only has to *match* a grid that
# has already been decided, which a built layer can. Widening the union would give
# `MethodError: no method matching _asraster(::Supply{...})` - a signature promising what the body
# cannot do, which is worse than not naming the type at all.
"""
    investigate_study_area(base = missing; regime::Union{LayerInput, Missing, Nothing} = missing,
        supply::Union{LayerInput, Missing, Nothing} = missing,
                           within = missing, crs = missing, cellsize = missing,
                           extent = missing, align = missing, simulate_safely = missing)

Report what a [`StudyArea`](@ref) built from these arguments would be, without building it.

Takes exactly the arguments [`StudyArea`](@ref) does - see there for what each one decides - so you
can investigate, adjust, and then commit to the result. In brief: `regime`/`supply` are the layers
allowed to shape the grid, `within` **positions** it, `crs` and `cellsize` fix the projection and
resolution, `align` names the layer kept exactly, `extent` sizes a *synthetic* area only, and
`simulate_safely` says whether a partly covered cell may be simulated.

Returns a [`StudyAreaReport`](@ref), which displays as a readable summary and can equally be
inspected programmatically (`report.problems`, `report.layers`, `report.footprint`).

**Under MPI with more than one rank this is a collective call**, made on every rank in the same
order; see [Building is collective](@ref).

The report can itself be handed back as the `base` of a `StudyArea`, which reuses its cache, so
committing to an investigated grid re-reads nothing:

```julia
report = investigate_study_area(; regime, supply, within = ShapeSpec("scotland.shp"))
area = StudyArea(report)            # same grid, no re-reading
```
"""
function investigate_study_area(base = missing;
                                regime::Union{LayerInput, Missing,
                                              Nothing} = missing,
                                supply::Union{LayerInput, Missing,
                                              Nothing} = missing,
                                within = missing, crs = missing,
                                cellsize = missing, extent = missing,
                                align = missing,
                                simulate_safely::Union{Bool, Missing,
                                                       Nothing} = missing)
    inputs = _resolveinputs(base; regime, supply, within, crs, cellsize, extent,
                            align, simulate_safely)
    return _analyse(inputs.layers; inputs.constraints..., cache = inputs.cache)
end

# == Building the pieces ========================================================================

"""
    build_habitat(source = DefaultEcosystem(); regime, supply, area, topology,
                  verbosity = :normal)

Create a [`GridHabitat`](@ref), filling anything not named from `source`.

[`GridHabitat`](@ref) itself is strict: omitting a required input is an error, because no default
could be right for every model. This fills those omissions from somewhere explicit instead, and
announces every value it chose.

  - `source`: where an unnamed input comes from. [`DefaultEcosystem`](@ref) gives a small synthetic
    toy grid; an existing [`GridHabitat`](@ref) gives its own layers, grid and topology, so
    `build_habitat(habitat, supply = ...)` is a rebuild with one thing changed.
  - `regime`: the **Condition** layer(s) - see [`GridHabitat`](@ref) for what is accepted.
  - `supply`: the **Resource** layer(s).
  - `area`: the [`StudyArea`](@ref) that decides the grid. From a habitat this is the grid it was
    *actually built on*, narrowing included - not the one originally investigated.
  - `topology`: how the grid's edges join - see [`EdgeTopology`](@ref).
  - `verbosity`: `:silent` announces nothing; anything else announces each filled-in value.

Under MPI it is a collective call, as [`GridHabitat`](@ref) is.

```julia
toy   = build_habitat()                         # every input defaulted, and said so
wetter = build_habitat(toy, supply = rainier)   # same grid, same regime, new supply
```
"""
function build_habitat(source = DefaultEcosystem();
                       # **`verbosity` must be declared FIRST, because the four below reference
                       # it**: a keyword default may only name keywords declared before it. Defaults
                       # evaluate in declaration order, and one passed explicitly is never evaluated
                       # - so naming a keyword suppresses exactly its own message.
                       verbosity::Symbol = :normal,
                       regime = _getdefaultvalue(Val(:regime), source,
                                                 verbosity, _HABITAT),
                       supply = _getdefaultvalue(Val(:supply), source,
                                                 verbosity, _HABITAT),
                       area = _getdefaultvalue(Val(:area), source, verbosity,
                                               _HABITAT),
                       topology = _getdefaultvalue(Val(:topology), source,
                                                   verbosity, _HABITAT))
    return GridHabitat(; regime, supply, area, topology)
end

"""
    build_species(numspecies::Integer; tolerance, toleranceaxis, demand, demandaxis, dispersal = 10.0km,
        pthresh = 1.0e-9, movement = InsectDispersal(),
        disperse_safely = nothing, birth = 0.6/year, death = 0.6/year, longevity = 1.0, survival = 0.2,
        abundance = 1000 * numspecies, native = true, names = nothing, seed = nothing,
        provenance = ProvenanceRecord[])

Build a `SpeciesList` of `numspecies` species. `tolerance` (the environmental **Condition** a
species is matched to) and `demand` (the **Resource** it consumes) are **required**, and so is the
[`NicheAxis`](@ref) each is measured on - `toleranceaxis` and `demandaxis`. Omitting any of the four
errors (use `build_species(DefaultEcosystem(); ...)` to fill omissions with announced defaults). Every per-species keyword accepts either a scalar (applied to all species) or a
length-`numspecies` vector (validated, with a clear error naming the argument on a length
mismatch):

  - `tolerance = (mean, width)` - Gaussian niche tolerance (a [`NicheTolerance`](@ref) with a `Normal`
    response). A tuple of `(mean, width)` pairs, e.g. `((298.0K, 2.0K), (50.0mm, 10.0mm))`, gives one
    tolerance per regime (a [`SpeciesRequirementCollection`](@ref)), to match a
    multi-variable environment.
  - `toleranceaxis` - the [`NicheAxis`](@ref) the tolerance is on. A `NicheTolerance`'s unit comes
    from its axis, so the axis is what says which environmental variable the species responds to and
    in what unit its `(mean, width)` are read. For a per-regime tolerance tuple, pass a matching
    tuple of axes.
  - `demand` - per-species resource demand rate, in the axis's own unit
    (`Demand{SolarRadiation}` is `kJ/day`, `Demand{Precipitation}` `L/day`,
    `Demand{CarbonFlux}` `g/day`); any scale of the right dimension is converted. A tuple of two
    demands (e.g. `(10.0kJ/day, 5.0L/day)`) gives a [`SpeciesRequirementCollection`](@ref) to match a
    two-resource supply.
  - `demandaxis` - the [`NicheAxis`](@ref) the demand is on, matching `demand`'s arity as
    `toleranceaxis` matches `tolerance`'s.

!!! note "The axis is named, never inferred"
    Both axes used to be guessable: `axis` defaulted to `Temperature`, and a demand's type was
    chosen from its **unit** (`kJ/day` -> solar, and so on). Neither is true now. A unit cannot say
    what a value *means* - `m/s` and `mm/day` share a dimension, which is how a wind-speed layer
    once became a water supply - so the axis is the single declaration, on both sides, and it is
    required.
  - `dispersal` - each species' mean dispersal distance, from which its kernel is built, cut off
    below `pthresh`.
  - `movement` - the dispersal mode, an [`AbstractMovement`](@ref), by default
    [`InsectDispersal`](@ref); its family must be built from a mean distance alone. Or the species' movements themselves, one [`SpeciesDispersal`](@ref) each,
    where it needs more, in which case `dispersal`, `pthresh` and `disperse_safely` are not read.
  - `disperse_safely` - per species, what becomes of an individual whose dispersal is aimed at a
    **dead cell** (off the grid, or an inactive one): `true` redistributes it among the reachable
    destinations, `false` loses it. By default the mode's own. A property of the *disperser*, so
    it may differ between species - a wind-dispersed seed blown out to sea is gone, an
    animal-dispersed one is not. The **grid's** topology (whether the world wraps) is a different
    question, answered by `GridHabitat(..., topology = ...)`.
  - `birth`, `death`, `longevity`, `survival` - demographic rates
    (scalar rates give an [`EqualPop`](@ref), vectors a [`PopGrowth`](@ref)).

`abundance` is either a total number of individuals split at random across
species or an explicit per-species vector; `seed` fixes the split, from a stream of its own, so
the same seed gives the same starting abundances whatever else has been drawn. `native` marks
species as native (default all `true`). `names` gives the species their names, one unique string
each - binomials, say - which also label their Diversity types; left unset they are `"1"` to
`"numspecies"`. `provenance` is an [`ProvenanceRecord`](@ref), or a vector of them, for the published
data the species were built from - an occurrence download, a trait table - which
[`provenance`](@ref) of the ecosystem then lists.
"""
function build_species(numspecies::Integer;
                       tolerance = _require(:tolerance),
                       # Required *unless* `tolerance` is already built - a pre-built one
                       # carries its own axis, so demanding a second statement of it would be
                       # ceremony that could also disagree.
                       toleranceaxis = nothing,
                       demand = _require(:demand),
                       demandaxis = _require(:demandaxis),
                       dispersal = _default(Val(:dispersal)).value,
                       pthresh::Real = 1.0e-9,
                       movement = _default(Val(:movement)).value.mode,
                       disperse_safely = nothing,
                       birth = _default(Val(:birth)).value,
                       death = _default(Val(:death)).value,
                       longevity::Real = _default(Val(:longevity)).value,
                       survival::Real = _default(Val(:survival)).value,
                       abundance = 1000 * numspecies,
                       native = true,
                       names = nothing,
                       seed = nothing,
                       provenance = ProvenanceRecord[])
    n = Int64(numspecies)
    _prebuilt = tolerance isa AbstractTolerance ||
                (tolerance isa Union{Tuple, NamedTuple} &&
                 all(t -> t isa AbstractTolerance, values(tolerance)))
    (_prebuilt || !isnothing(toleranceaxis)) || _require(:toleranceaxis)

    # **A pre-built tolerance must be checked FIRST**, because the branches below reach inside
    # `tolerance` with `first(...)`, and on a `NicheTolerance` object that gives
    # `MethodError: no method matching iterate(::NicheTolerance{...})` - an error naming `iterate`
    # rather than the mistake.
    # **And accepting one is better than merely reporting it well**, because it is the only way to
    # build a species against a layer whose frame `build_species` cannot otherwise express: a
    # tolerance on `NicheAxis` carrying a real unit, say, which `NicheTolerance(..., support = ...)`
    # can construct and no combination of keywords here could.
    traits = if tolerance isa AbstractTolerance
        _checktolerancespecies(tolerance, n)
        tolerance
    elseif tolerance isa Union{Tuple, NamedTuple} &&
           all(t -> t isa AbstractTolerance, values(tolerance))
        foreach(t -> _checktolerancespecies(t, n), values(tolerance))
        SpeciesRequirementCollection(tolerance)
    elseif first(tolerance) isa Tuple
        axes = toleranceaxis isa Union{Tuple, NamedTuple} ?
               values(toleranceaxis) :
               ntuple(_ -> toleranceaxis, length(tolerance))
        # Built positionally, then put back into whichever backing `tolerance` had, so a named
        # `tolerance = (temperature = ..., rainfall = ...)` gives a named collection
        # however `axis` was written.
        built = _zipmap(axes, values(tolerance)) do ax, tol
            return _nichetolerance(ax, tol, n)
        end
        SpeciesRequirementCollection(_asbacking(tolerance, built))
    else
        _nichetolerance(toleranceaxis, tolerance, n)
    end

    move = _speciesmovement(movement, dispersal, disperse_safely, pthresh, n)

    # A tuple `demand` gives one demand per resource -> a collection; the backing is put back
    # so a named `demand` gives a named collection. `demandaxis` is zipped against it exactly as
    # `toleranceaxis` is against `tolerance` - including `_zipmap`'s arity check, so a mismatched
    # pair is refused rather than silently truncated.
    demands = if demand isa Union{Tuple, NamedTuple}
        daxes = demandaxis isa Union{Tuple, NamedTuple} ? values(demandaxis) :
                ntuple(_ -> demandaxis, length(demand))
        built = _zipmap(daxes, values(demand)) do ax, d
            return _demand(ax, _tofield(d, n, "demand"))
        end
        SpeciesRequirementCollection(_asbacking(demand, built))
    else
        _demand(demandaxis, _tofield(demand, n, "demand"))
    end
    param = _params(birth, death, longevity, survival, n)
    nat = _tofield(native, n, "native")
    abun = _abundances(abundance, n, seed)

    return SpeciesList(n, traits, abun, demands, move, param, nat,
                       names = names, provenance = provenance)
end

"""
    build_species(pool::SpeciesPool; tolerance, demand,
                  movement = MovementSpec(InsectDispersal(), meandist = 10.0km),
                  birth = 0.6 / year, death = 0.6 / year, longevity = 1.0, survival = 0.2,
                  seed = nothing, verbosity = :normal, provenance = ProvenanceRecord[])

Build a `SpeciesList` of every species of `pool`, one spec per trait. `tolerance` and `demand` are
required: they must match the habitat the species will live in, axis for axis. Every other trait
has a default, announced when it is taken. Each parameter is a trait value - one value for every
species, a [`TraitColumn`](@ref) of the pool, a [`DrawFrom`](@ref) a distribution, or a function of
the species' row of the pool - and the pool's missing values are filled first, as
[`fill_traits`](@ref) fills them; the filled values, and which were filled, are kept on the species
list. Nothing is taken by position.

A deactivated species is built with the rest and starts with no individuals; the active ones share
1000 individuals each, split at random, from a stream of `seed`'s own. The pool's records, and its
authority's, are the species list's inputs, which [`provenance`](@ref) of an ecosystem lists.

# Arguments

  - `pool`: the species, a [`SpeciesPool`](@ref).
  - `tolerance`: each species' niche, a [`ToleranceSpec`](@ref), a tuple or named tuple of them for
    a habitat with several regimes, or an [`AbstractTolerance`](@ref) built already.
  - `demand`: what each species takes from the shared supply, a [`DemandSpec`](@ref), a tuple or
    named tuple of them for several supplies, or an [`AbstractDemand`](@ref) built already.
  - `movement`: how each species disperses, a [`MovementSpec`](@ref), or a [`ChosenBy`](@ref)
    choosing one per species.
  - `birth`, `death`: each species' rates, per unit time.
  - `longevity`: the exponent on each species' own demand, which sets the tempo of its turnover.
  - `survival`: the exponent on niche suitability, one number for every species.
  - `seed`: an integer fixing every value drawn and the split of the starting individuals, each
    from a stream of its own; left unset, both come from the global generator.
  - `verbosity`: `:silent` says nothing, `:normal` announces each default taken and notes each column
    read and filled, `:verbose` names every name a note concerns.
  - `provenance`: a [`ProvenanceRecord`](@ref), or a vector of them, for what nothing else records -
    a paper parameters were taken from.
"""
function build_species(pool::SpeciesPool;
                       # Declared first, because the defaults below announce themselves through it.
                       verbosity::Symbol = :normal,
                       tolerance = _require(:tolerance),
                       demand = _require(:demand),
                       movement = _getdefaultvalue(Val(:movement),
                                                   DefaultEcosystem(),
                                                   verbosity, _POOLSPECIES),
                       birth = _getdefaultvalue(Val(:birth), DefaultEcosystem(),
                                                verbosity, _POOLSPECIES),
                       death = _getdefaultvalue(Val(:death), DefaultEcosystem(),
                                                verbosity, _POOLSPECIES),
                       longevity = _getdefaultvalue(Val(:longevity),
                                                    DefaultEcosystem(),
                                                    verbosity,
                                                    _POOLSPECIES),
                       survival = _getdefaultvalue(Val(:survival),
                                                   DefaultEcosystem(),
                                                   verbosity, _POOLSPECIES),
                       seed::Union{Nothing, Integer} = nothing,
                       provenance = ProvenanceRecord[])
    _checkverbosity(verbosity)
    n = length(pool)
    n > 0 || error("the pool holds no species to build.")
    slots = _TraitSlot[]
    _slots!(slots, pool, "tolerance", tolerance, Colon())
    _slots!(slots, pool, "demand", demand, Colon())
    _slots!(slots, pool, "movement", movement, Colon())
    for (name, value) in ((:birth, birth), (:death, death),
        (:longevity, longevity))
        push!(slots, _TraitSlot(name, value, Colon(), nothing))
    end
    filled = _fillpool(pool, slots, seed, verbosity)
    table = filled.values
    rows = 1:n
    problems = Problem[]
    tolerances = _buildtrait(tolerance, pool, table, "tolerance", rows,
                             problems,
                             verbosity)
    demands = _buildtrait(demand, pool, table, "demand", rows, problems,
                          verbosity)
    movements = _buildmovement(movement, pool, table, "movement", rows,
                               problems,
                               verbosity)
    params = _buildparams(birth, death, longevity, survival, pool, table, rows,
                          problems, verbosity)
    _reportpool(problems, verbosity)
    active = table.active
    abun = zeros(Int64, n)
    abun[active] = _abundances(1000 * count(active), count(active), seed)
    return SpeciesList(n, tolerances, abun, demands, movements, params,
                       fill(true, n), names = table.name,
                       provenance = vcat(pool.inputs,
                                         _recordvector(provenance)),
                       pool = (traits = table, filled = filled.filled,
                               authority = pool.authority,
                               namemap = _speciesnamemap(pool)))
end

# As `build_habitat` fills its inputs: `numspecies`, `tolerance` and `demand` and their axes from
# `DefaultEcosystem()` where omitted, announced, then delegated. Any other keyword (dispersal,
# movement, rates, abundance, ...) passes straight through.
function build_species(source::DefaultEcosystem;
                       verbosity::Symbol = :normal,
                       numspecies = _getdefaultvalue(Val(:numspecies), source,
                                                     verbosity,
                                                     _DEFAULTSPECIES),
                       tolerance = _getdefaultvalue(Val(:tolerance), source,
                                                    verbosity, _DEFAULTSPECIES),
                       toleranceaxis = _getdefaultvalue(Val(:toleranceaxis),
                                                        source,
                                                        verbosity,
                                                        _DEFAULTSPECIES),
                       demand = _getdefaultvalue(Val(:demand), source,
                                                 verbosity,
                                                 _DEFAULTSPECIES),
                       demandaxis = _getdefaultvalue(Val(:demandaxis), source,
                                                     verbosity,
                                                     _DEFAULTSPECIES),
                       kw...)
    return build_species(numspecies; tolerance = tolerance,
                         toleranceaxis = toleranceaxis, demand = demand,
                         demandaxis = demandaxis, kw...)
end

"""
    build_ecosystem(species::SpeciesList, environment::GridHabitat;
        nichefit = nothing, seed = nothing, distributed = :auto, epoch = nothing,
        calendar = ExactDates(), provenance = ProvenanceRecord[])

Assemble an ecosystem from a `species` list and an `environment`. When
`nichefit` is not given it is inferred from the trait type (`NicheTolerance` ->
`NicheSuitability`, any `AbstractCategoricalTolerance` -> `CategoricalSuitability`). Checks that the
species' resource demand matches the environment's supply before building,
and passes `seed` through for reproducibility.

`distributed` selects the ecosystem type: `:auto` (the default) builds a
distributed `MPIEcosystem` when the MPI extension is loaded and the process is
running with more than one rank (`MPI.Init()` already called), and a serial
`Ecosystem` otherwise; `true` forces an `MPIEcosystem` (erroring if the MPI
extension is not loaded); `false` forces a serial `Ecosystem`. Note that
auto-selection only covers construction and `simulate!` - under MPI the
abundances are partitioned across ranks, so collecting results still needs
[`gatherabundance`](@ref)/[`gatherdiversity`](@ref).

`epoch` is the real date the run begins - the date elapsed time zero means. It fixes the run's
calendar for output (see [`simulationdate`](@ref)) and, more importantly, *phases* every series the
environment carries: a [`DatedSeries`](@ref) starts at the slice covering the epoch, and a
[`MonthOfYearSeries`](@ref) climatology at the slice for the epoch's own month, so a run beginning in
July starts in July rather than in January.

It is resolved the way [`StudyArea`](@ref) resolves a CRS - adopt if unambiguous, ask if not. An
explicit `epoch` always wins; otherwise the environment's own series supply it, and if exactly one
real start date is found it is used and everything else is phased to it. Series that disagree are an
error naming the candidates, and an environment with no dated series has no epoch at all, which is
the behaviour of a run that never mentions dates.

`calendar` is how dates become elapsed time: [`ExactDates`](@ref), the default, places every dated
slice and the epoch by the real time between them, and [`MeanMonths`](@ref) counts each calendar
month as `month_mean_duration`, so a dated monthly series stepped by `month_mean_duration` makes
every month current once. Under `ExactDates` such a run is refused by [`simulate!`](@ref), since
real months are 28 to 31 days long and the step skips some. The calendar also decides the date
[`simulationdate`](@ref) reports.

`provenance` is an [`ProvenanceRecord`](@ref), or a vector of them, for published data that belongs to
the run as a whole rather than to its species or its environment - the records a starting
population was seeded from, say, or the study's own DOI. [`provenance`](@ref) of the ecosystem lists
them beside everything else it was built from.
"""
function build_ecosystem(species::SpeciesList, environment::GridHabitat;
                         nichefit = nothing, seed = nothing,
                         distributed = :auto,
                         epoch::Union{Nothing, Dates.TimeType} = nothing,
                         calendar::AbstractRunCalendar = ExactDates(),
                         provenance = ProvenanceRecord[])
    records = _recordvector(provenance)
    _checksimulatable(environment)
    # Checked here as well as in the `Ecosystem` constructor, and *before* the nichefit is inferred:
    # `_defaultsuitability` pairs tolerances with regimes member by member, so a mismatch reaches
    # it first and would be reported from inside the inference rather than about the two inputs.
    _checkaligned("species demand and environment supply",
                  _demandside(species.demand),
                  _supplyside(environment.supply))
    _checkaligned("species tolerances and environment regime",
                  _toleranceside(species.tolerance),
                  _regimeside(environment.regime))
    nichefit = isnothing(nichefit) ?
               _defaultsuitability(species.tolerance, environment.regime) :
               nichefit
    # Resolved and applied *before* the ecosystem is built, because re-pointing rewrites the
    # layers' changes and the ecosystem holds the very same habitat object. Mutating in place is the
    # existing semantics here - an `Ecosystem` shares its habitat with the caller rather than copying
    # it - and re-pointing is idempotent, since an origin is computed from the calendar and the epoch
    # rather than accumulated.
    resolved = _resolveepoch(environment, epoch)
    _repointseries!(environment.regime, resolved, calendar)
    _repointseries!(environment.supply, resolved, calendar)
    # Then write those values in, so the run *starts* in the state its series and epoch describe
    # rather than reaching it one timestep late. Ordered after re-pointing for the obvious reason,
    # and safe against data gaps because `GridHabitat` has already cleaned both a supply's matrix
    # and the slices its change is holding to write later.
    _primeseries!(environment.regime, zero(1.0s))
    _primeseries!(environment.supply, zero(1.0s))
    eco = if _resolvedistributed(distributed)
        isnothing(seed) ? MPIEcosystem(species, environment, nichefit) :
        MPIEcosystem(species, environment, nichefit, seed = seed)
    else
        isnothing(seed) ? Ecosystem(species, environment, nichefit) :
        Ecosystem(species, environment, nichefit, seed = seed)
    end
    eco.epoch = resolved
    eco.calendar = calendar
    eco.inputs = _uniqueinputs(records)
    return eco
end

# Complete the trio: build the default environment and default species (each announcing its own
# defaults) and assemble them into a runnable `Ecosystem`. `seed`/`distributed` pass to the assembly.
function build_ecosystem(::DefaultEcosystem; seed = nothing,
                         distributed = :auto,
                         epoch::Union{Nothing, Dates.TimeType} = nothing,
                         calendar::AbstractRunCalendar = ExactDates(),
                         provenance = ProvenanceRecord[])
    environment = build_habitat()
    species = build_species(DefaultEcosystem())
    return build_ecosystem(species, environment, seed = seed,
                           distributed = distributed, epoch = epoch,
                           calendar = calendar, provenance = provenance)
end

# == Running it =================================================================================

"""
    simulate!(eco::AbstractEcosystem, times::Unitful.Time, timestep::Unitful.Time)

Run an ecosystem, `eco` for a specified length of time, `times`, for a
particular timestep, `timestep`.
"""
function simulate!(eco::AbstractEcosystem, duration::Unitful.Time,
                   timestep::Unitful.Time; intervention = nothing)
    checkcoverage(eco, duration, timestep)
    check_bounds(eco, duration, timestep)
    _checkschedules(intervention, eco, duration, timestep)
    for _ in 1:_stepcount(duration, timestep)
        update!(eco, timestep, intervention)
    end
end

"""
    simulate!(f, eco::AbstractEcosystem, duration::Unitful.Time, timestep::Unitful.Time;
              every = timestep, intervention = nothing)

Run `eco` for `duration` in steps of `timestep`, as the three-argument `simulate!` does, calling `f`
whenever the schedule `every` fires - usually written as a `do` block:

```julia
totals = Int[]
simulate!(eco, 10year, 1month_mean_duration, every = EveryInterval(1year)) do occurrence
    push!(totals, sum(eco.abundances.matrix))
end
```

`f` is handed one named tuple per occurrence: `count`, which occurrence of this call it is, from 1;
`elapsed`, the simulation time it fires at; and `date`, that time as a date under the run's epoch and
calendar, or `nothing` for a run with no epoch. The ecosystem is reached by closure. Each occurrence
sees the state at exactly its time: on a run starting from elapsed zero the first is the starting
state, after any intervention due at the start, and each later one comes after its step's dynamics,
interventions and layer update.

`f` may also be a recorder - an [`AbstractRecorder`](@ref) such as [`RecordAbundance`](@ref) -
which is called with the ecosystem and the occurrence; to keep several things from one run, call
each from a callback of your own.

**The callback is for observing the run** - recording and logging. It must not draw from the global
random number generator or write a layer's values, either of which would make the result depend on
how the work was divided. Under MPI every rank runs the loop and the callback, so a collective the
callback calls, such as [`gatherabundance`](@ref) or [`gatherdiversity`](@ref), must be reached on
every rank - a callback returning before it on the ranks other than the root leaves every rank waiting
for ever. A callback writing abundances reaches only its own rank's species there, and the demand the
next step computes does not include what it wrote. To change the ecosystem, declare an
[`Intervention`](@ref), which is applied once and identically everywhere.

# Arguments

  - `f`: called at each occurrence with the named tuple above.
  - `eco`: the ecosystem to run.
  - `duration`, `timestep`: the length of the run and of each step; the run takes
    `duration / timestep` steps - the nearest whole number for a step of a day or less, and
    otherwise exactly a whole number.
  - `every`: when `f` is called - an [`AbstractSchedule`](@ref), or a duration meaning
    [`EveryInterval`](@ref) of it. Every timestep by default, the starting state included.
  - `intervention`: an [`Intervention`](@ref) or [`InterventionSet`](@ref) applied as the run
    proceeds.
"""
function simulate!(f, eco::AbstractEcosystem, duration::Unitful.Time,
                   timestep::Unitful.Time; every = timestep,
                   intervention = nothing)
    schedule = _everyschedule(every)
    checkcoverage(eco, duration, timestep)
    check_bounds(eco, duration, timestep)
    _checkschedules(intervention, eco, duration, timestep)
    _checkschedule(schedule, eco, duration, timestep)
    count = 0
    if iszero(simulationtime(eco))
        _startrun!(eco, intervention, timestep)
        count = _occur(f, schedule, eco, timestep, count)
    end
    for _ in 1:_stepcount(duration, timestep)
        _step!(eco, timestep, intervention)
        count = _occur(f, schedule, eco, timestep, count)
    end
    return nothing
end

function simulate!(recorder::AbstractRecorder, eco::AbstractEcosystem,
                   duration::Unitful.Time, timestep::Unitful.Time;
                   every = timestep, intervention = nothing)
    simulate!(eco, duration, timestep, every = every,
              intervention = intervention) do occurrence
        return recorder(eco, occurrence)
    end
    return nothing
end

"""
    simulate!(cache::CachedEcosystem, srt::Unitful.Time, timestep::Unitful.Time)

Run a cached ecosystem, `cache` at a specified timepoint, `srt`, for a
particular timestep, `timestep`.
"""
function simulate!(cache::CachedEcosystem, srt::Unitful.Time,
                   timestep::Unitful.Time)
    eco = Ecosystem{typeof(cache.habitat), typeof(cache.spplist),
                    typeof(cache.nichefit)}(GridLandscape(copy(cache.abundances.matrix[Ti(At(srt))].matrix),
                                                          cache.spplist.names,
                                                          getcoords(cache.habitat)),
                                            cache.spplist,
                                            cache.habitat,
                                            cache.nichefit,
                                            cache.lookup,
                                            cache.cache,
                                            cache.rngs,
                                            uconvert(s, float(srt)),
                                            cache.seed,
                                            cache.epoch)
    eco.calendar = cache.calendar
    update!(eco, timestep)
    return cache.abundances.matrix[Ti(At(srt + timestep))] = eco.abundances
end

"""
    generate_storage(eco::Ecosystem, times::Int64, reps::Int64;
                     maxspecies = length(eco.spplist.abun))

Allocate an integer array of shape `(maxspecies, gridSize, times, reps)` for
recording species abundances across the ecosystem `eco` over multiple timesteps
and replicate runs.

`maxspecies` defaults to the number of species `eco` starts with, which is what
a run that never gains one needs. **Raise it to make room for species that
arrive during the run** - an [`AddSpecies`](@ref) intervention, an invasion - and
the recorder then has somewhere to put them:

```julia
storage = generate_storage(eco, times, reps, maxspecies = 12)
```

The array is sized **once**, before the run, and cannot grow afterwards: it is
a dense array indexed by species number, so growing it mid-run would mean
reallocating and copying every recorded step. A run that gains more species than
`maxspecies` allows for says so when it tries to record them, rather than
silently dropping them.
"""
function generate_storage(eco::Ecosystem, times::Int64, reps::Int64;
                          maxspecies::Int64 = length(eco.spplist.abun))
    maxspecies >= length(eco.spplist.abun) ||
        error("`maxspecies` is $maxspecies but the ecosystem already has " *
              "$(length(eco.spplist.abun)) species - the recording has to hold at least what is " *
              "there before the run starts.")
    gridSize = countsubcommunities(eco.habitat.regime)
    return abun = Array{Int64, 4}(undef, maxspecies, gridSize, times, reps)
end

"""
    generate_storage(eco::Ecosystem, qs::Int64, times::Int64, reps::Int64)

Allocate a float array of shape `(gridSize, qs, times, reps)` for recording
diversity values across the ecosystem `eco` for `qs` diversity orders over
multiple timesteps and replicate runs.
"""
function generate_storage(eco::Ecosystem, qs::Int64, times::Int64, reps::Int64)
    gridSize = countsubcommunities(eco.habitat.regime)
    return abun = Array{Float64, 4}(undef, gridSize, qs, times, reps)
end

# The schedule a callback's `every` names: a duration means every whole multiple of it.
_everyschedule(every::Unitful.Time) = EveryInterval(every)

_everyschedule(every::AbstractSchedule) = every

# Call `f` for an occurrence if `schedule` fires at the clock as it now stands, and return how many
# occurrences there have been.
function _occur(f, schedule::AbstractSchedule, eco::AbstractEcosystem,
                timestep::Unitful.Time, count::Integer)
    _fires(schedule, eco, simulationtime(eco), timestep) || return count
    count += 1
    f((count = count, elapsed = simulationtime(eco),
       date = simulationdate(eco)))
    return count
end
