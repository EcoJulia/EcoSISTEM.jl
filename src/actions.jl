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
# The placements are absent for the opposite reason: they are values `build_ecosystem` takes rather
# than verbs in their own right, so they live with the placement machinery (`Placement.jl`).
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
    build_species(pool::SpeciesPool; tolerance, demand,
                  movement = MovementSpec(InsectDispersal(), meandist = 10.0km),
                  birth = 0.6 / year, death = 0.6 / year, longevity = 1.0, survival = 0.2,
                  seed = nothing, verbosity = :normal, provenance = ProvenanceRecord[])
    build_species(DefaultEcosystem(); numspecies = 10, tolerance, demand, ...)

Build a `SpeciesList` of every species of `pool`, one spec per trait. `tolerance` and `demand` are
required: they must match the habitat the species will live in, axis for axis. Every other trait
has a default, announced when it is taken. Each parameter is a trait value - one value for every
species, a [`TraitColumn`](@ref) of the pool, a [`DrawFrom`](@ref) a distribution, or a function of
the species' row of the pool - and the pool's missing values are filled first, as
[`fill_traits`](@ref) fills them; the filled values, and which were filled, are kept on the species
list. Nothing is taken by position.

A deactivated species is built with the rest, and a placement places none of it; where the
individuals start is the placement [`build_ecosystem`](@ref) is given. The pool's records, and its
authority's, are the species list's inputs, which [`provenance`](@ref) of an ecosystem lists.

# Arguments

  - `pool`: the species, a [`SpeciesPool`](@ref); or `DefaultEcosystem()` for a pool of `numspecies`
    species named `"1"` onwards, whose `tolerance` and `demand`, like every other trait, default
    and are announced.
  - `numspecies`: with `DefaultEcosystem()`, how many species.
  - `tolerance`: each species' niche, a [`ToleranceSpec`](@ref), a tuple or named tuple of them for
    a habitat with several regimes, or an [`AbstractTolerance`](@ref) built already.
  - `demand`: what each species takes from the shared supply, a [`DemandSpec`](@ref), a tuple or
    named tuple of them for several supplies, or an [`AbstractDemand`](@ref) built already.
  - `movement`: how each species disperses, a [`MovementSpec`](@ref), or a [`ChosenBy`](@ref)
    choosing one per species.
  - `birth`, `death`: each species' rates, per unit time.
  - `longevity`: the exponent on each species' own demand, which sets the tempo of its turnover.
  - `survival`: the exponent on niche suitability, one number for every species.
  - `seed`: an integer fixing every value drawn, each from a stream of its own; left unset, they
    come from the global generator.
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
    _reportproblems(problems, verbosity)
    return SpeciesList(n, tolerances, demands, movements, params,
                       names = table.name,
                       provenance = vcat(pool.inputs,
                                         _recordvector(provenance)),
                       pool = (traits = table, filled = filled.filled,
                               authority = pool.authority,
                               namemap = _speciesnamemap(pool),
                               phylogeny = pool.phylogeny))
end

# As `build_habitat` fills its inputs: `numspecies`, `tolerance` and `demand` from
# `DefaultEcosystem()` where omitted, announced, then built from a pool of that many anonymous
# species. Any other keyword (movement, rates, seed, ...) passes straight through; one only the
# removed `build_species(numspecies; ...)` took is refused with its rewrite.
function build_species(source::DefaultEcosystem;
                       verbosity::Symbol = :normal,
                       numspecies = _getdefaultvalue(Val(:numspecies), source,
                                                     verbosity,
                                                     _DEFAULTSPECIES),
                       tolerance = _getdefaultvalue(Val(:tolerance), source,
                                                    verbosity, _DEFAULTSPECIES),
                       demand = _getdefaultvalue(Val(:demand), source,
                                                 verbosity,
                                                 _DEFAULTSPECIES),
                       kw...)
    _refuseremovedkeywords(numspecies, tolerance, demand, kw)
    return build_species(SpeciesPool(numspecies, verbosity = :silent);
                         tolerance = tolerance, demand = demand,
                         verbosity = verbosity, kw...)
end

"""
    build_ecosystem(placement, species::SpeciesList, environment::GridHabitat;
        nichefit = nothing, seed = nothing, distributed = :auto, epoch = nothing,
        calendar = ExactDates(), provenance = ProvenanceRecord[], verbosity = :normal)
    build_ecosystem(species::SpeciesList, environment::GridHabitat; ...)
    build_ecosystem(species::SpeciesList, environment::GridHabitat; from, elapsed = nothing, ...)
    build_ecosystem(species, environment; ...) do species, habitat, context
        ...
    end

Assemble an ecosystem from a `species` list and an `environment`, its individuals placed by
`placement`. When `nichefit` is not given it is inferred from the trait type (`NicheTolerance` ->
`NicheSuitability`, any `AbstractCategoricalTolerance` -> `CategoricalSuitability`). Checks that the
species' resource demand matches the environment's supply before building,
and passes `seed` through for reproducibility.

`placement` says how many individuals of each species start in each cell: a shipped placement -
[`BySupply`](@ref), [`BySuitability`](@ref), [`ByLayer`](@ref), [`Weighted`](@ref),
[`FromTable`](@ref), [`FromOccurrences`](@ref), [`NoAbundance`](@ref), any filled to capacity by
[`ToCapacity`](@ref), or several added with `+` -
or any function of `(species, habitat, context)`, as a `do` block, returning a
[`PlacedAbundances`](@ref) or a table of the individuals it places. It is called with the species list, the environment and a
[`PlacementContext`](@ref) - the cells it may place in (the active ones), a random stream per
species, the run's dates - and its table is checked before anything is written: a row naming a
species the list does not hold or that is not active, or a place off the grid or in an inactive
cell, is dropped and reported, and a count that is not a whole number stops the build. Under MPI the
first rank alone calls it. Left out, it is `DefaultEcosystem()`, 1000 individuals of each active
species, split at random and placed by supply, and is announced. The ecosystem keeps it, for
[`repopulate!`](@ref) and the run's provenance.

`verbosity` is `:silent`, `:normal` or `:verbose`; `:silent` announces nothing, and reports only what
could not be placed.

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

`from` starts the ecosystem from a state [`write_state`](@ref) saved, a path or its record, in place of
a placement. The species and habitat are built from the same pool and study area as the run that
saved it, on its grid: every species it saved is in the list, matched by name and key in the pool,
and the list may hold more, which start with none and as their pool says. The state gives every
saved species its abundances and whether it is in the run, the habitat its active cells and every
layer's values, and the run its epoch and calendar, which a run given others refuses. On the
state's own seed, the default, each saved species continues its own random stream, so the run goes
on as the saved one would have; another `seed` gives fresh streams, for replicates from one state.
The clock starts at `elapsed`, the state's own unless given (`elapsed = 0.0s` starts a run afresh
from it), and any series a layer follows takes its values at that time. A layer's change rule is not
saved: at its first step a run continuing a state installs again the rules its interventions would
have installed by then, so given the saved run's interventions it goes on with the saved run's rules,
and a run whose rules still differ from those saved is warned. The state's record is among the
ecosystem's inputs.
"""
function build_ecosystem(placement, species::SpeciesList,
                         environment::GridHabitat;
                         nichefit = nothing, seed = nothing,
                         distributed = :auto,
                         epoch::Union{Nothing, Dates.TimeType} = nothing,
                         calendar::AbstractRunCalendar = ExactDates(),
                         provenance = ProvenanceRecord[],
                         verbosity::Symbol = :normal, from = nothing,
                         elapsed = nothing)
    isnothing(from) ||
        error("`from` starts an ecosystem from a saved state, which says where its individuals " *
              "are, so it takes no placement: give `build_ecosystem(species, habitat; from = " *
              "...)`.")
    isnothing(elapsed) ||
        error("`elapsed` sets the clock of an ecosystem started from a saved state, so it goes " *
              "with `from`.")
    _checkverbosity(verbosity)
    _checkplacement(placement)
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
        _emptympiecosystem(species, environment, nichefit, seed = seed)
    else
        isnothing(seed) ? _emptyecosystem(species, environment, nichefit) :
        _emptyecosystem(species, environment, nichefit, seed = seed)
    end
    eco.epoch = resolved
    eco.calendar = calendar
    eco.inputs = _uniqueinputs(records)
    context = _placementcontext(eco, eco.seed)
    return _placeinto!(eco, placement, context, verbosity, clear = false)
end

# The placement left out: a saved state, or the default, announced.
function build_ecosystem(species::SpeciesList, environment::GridHabitat;
                         verbosity::Symbol = :normal, from = nothing, kw...)
    isnothing(from) ||
        return _startfromstate(species, environment, from;
                               verbosity = verbosity, kw...)
    placement = _getdefaultvalue(Val(:placement), DefaultEcosystem(), verbosity,
                                 "build_ecosystem")
    return build_ecosystem(placement, species, environment;
                           verbosity = verbosity, kw...)
end

# Complete the trio: build the default environment and default species (each announcing its own
# defaults) and assemble them into a runnable `Ecosystem`. `seed`/`distributed` pass to the assembly.
function build_ecosystem(::DefaultEcosystem; seed = nothing,
                         distributed = :auto,
                         epoch::Union{Nothing, Dates.TimeType} = nothing,
                         calendar::AbstractRunCalendar = ExactDates(),
                         provenance = ProvenanceRecord[],
                         verbosity::Symbol = :normal)
    environment = build_habitat(verbosity = verbosity)
    species = build_species(DefaultEcosystem(), verbosity = verbosity)
    return build_ecosystem(species, environment, seed = seed,
                           distributed = distributed, epoch = epoch,
                           calendar = calendar, provenance = provenance,
                           verbosity = verbosity)
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
    _prepareinterventions!(intervention, eco)
    for _ in 1:_stepcount(duration, timestep)
        update!(eco, timestep, intervention)
    end
    _finishinterventions!(eco, intervention)
    return nothing
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
    _prepareinterventions!(intervention, eco)
    _checkschedule(schedule, eco, duration, timestep)
    count = 0
    if iszero(simulationtime(eco))
        _startrun!(eco, intervention, timestep)
        count = _occur(f, schedule, eco, timestep, count)
    end
    _resumerun!(eco, intervention, timestep)
    for _ in 1:_stepcount(duration, timestep)
        _step!(eco, timestep, intervention)
        count = _occur(f, schedule, eco, timestep, count)
    end
    _finishinterventions!(eco, intervention)
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
    eco.placement = cache.placement
    update!(eco, timestep)
    return cache.abundances.matrix[Ti(At(srt + timestep))] = eco.abundances
end

"""
    generate_storage(eco::Ecosystem, times::Int64, reps::Int64)

Allocate an integer array of shape `(numspecies, gridSize, times, reps)` for recording species
abundances across the ecosystem `eco` over multiple timesteps and replicate runs: a row for every
species of its pool, those an [`AddSpecies`](@ref) brings in during the run included.

# Arguments

  - `eco`: the ecosystem whose species and cells size the storage.
  - `times`: how many occurrences to leave room for.
  - `reps`: how many replicate runs.
  - `maxspecies`: deprecated, and ignored: every species a run may hold has its row already.
"""
function generate_storage(eco::Ecosystem, times::Int64, reps::Int64;
                          maxspecies = nothing)
    isnothing(maxspecies) || _maxspeciesdeprecated(:generate_storage)
    gridSize = countsubcommunities(eco.habitat.regime)
    return abun = Array{Int64, 4}(undef, length(eco.spplist.names), gridSize,
                                  times, reps)
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
