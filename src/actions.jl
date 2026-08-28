# SPDX-License-Identifier: LGPL-3.0-or-later
#
# HOW TO RUN AN EXPERIMENT. This file is the workflow and nothing else: the verbs a user calls, in the
# order they call them. It is the last of the model's own files — only `deprecations.jl` and the
# `ClimatePref` shims load after it — and, like `Ecology.jl`, which is the first, it should be
# readable on its own by someone who never opens another.
#
# `Ecology.jl` says what the model IS. This says what you DO with it:
#
#     investigate_study_area   decide the grid, and argue with it, before building anything on it
#     build_habitat            the environment: what each place is like, and what it provides
#     build_species            the species: what they tolerate, what they need, how they move
#     build_ecosystem          put the species in the environment and check the two agree
#     simulate!                run it
#     simulate_record!         run it and record abundances, into `generate_storage`'s array
#     simulate_record_diversity!   the same, recording a diversity measure instead
#     simulate_action!         run it and call back at intervals, for anything else
#     generate_storage         allocate the array the two recorders write into
#
# WHAT IS NOT HERE, and the test that decides it: these are the verbs a USER calls. Everything the
# package calls on their behalf stays with its machinery. So `update!` is absent — `simulate!` calls
# it, and it is `public` rather than exported for that reason — and so is every private helper,
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
using Printf
using JLD2
import Diversity.Gamma

# ══ Deciding the grid ══════════════════════════════════════════════════════════════════════════

# `Missing` and `Nothing` are both in the unions and mean **different** things: `missing` is "not
# given — derive it, or inherit it from `base`", `nothing` is "clear what would have been inherited".
# Neither is a layer, which is why each has to be admitted explicitly rather than falling out of
# [`LayerInput`](@ref).
#
# **An already-built layer must NOT be accepted here, in either role** — unlike `GridHabitat`,
# which takes a built `AbstractSupply`. The difference is what the function does with it: this pair
# *decides the grid from its layers*, so a layer must be able to answer `_shapesgrid`/`_layerfacts`
# (extent, CRS, cell size), and a built layer cannot. `GridHabitat` only has to *match* a grid that
# has already been decided, which a built layer can. Widening the union would give
# `MethodError: no method matching _asraster(::Supply{…})` — a signature promising what the body
# cannot do, which is worse than not naming the type at all.
"""
    investigate_study_area(base = missing; regime::Union{LayerInput, Missing, Nothing} = missing,
        supply::Union{LayerInput, Missing, Nothing} = missing,
                           within = missing, crs = missing, cellsize = missing,
                           extent = missing, align = missing, simulate_safely = missing)

Report what a [`StudyArea`](@ref) built from these arguments would be, without building it.

Takes exactly the arguments [`StudyArea`](@ref) does — see there for what each one decides — so you
can investigate, adjust, and then commit to the result. In brief: `regime`/`supply` are the layers
allowed to shape the grid, `within` **positions** it, `crs` and `cellsize` fix the projection and
resolution, `align` names the layer kept exactly, `extent` sizes a *synthetic* area only, and
`simulate_safely` says whether a partly covered cell may be simulated.

Returns a [`StudyAreaReport`](@ref), which displays as a readable summary and can equally be
inspected programmatically (`report.problems`, `report.layers`, `report.footprint`).

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

# ══ Building the pieces ════════════════════════════════════════════════════════════════════════

"""
    build_habitat(source = DefaultEcosystem(); regime, supply, area, topology,
                  verbosity = :normal)

Create a [`GridHabitat`](@ref), filling anything not named from `source`.

[`GridHabitat`](@ref) itself is strict: omitting a required input is an error, because no default
could be right for every model. This fills those omissions from somewhere explicit instead, and
announces every value it chose.

  - `source`: where an unnamed input comes from. [`DefaultEcosystem`](@ref) gives a small synthetic
    toy grid; an existing [`GridHabitat`](@ref) gives its own layers, grid and topology, so
    `build_habitat(habitat, supply = …)` is a rebuild with one thing changed.
  - `regime`: the **Condition** layer(s) — see [`GridHabitat`](@ref) for what is accepted.
  - `supply`: the **Resource** layer(s).
  - `area`: the [`StudyArea`](@ref) that decides the grid. From a habitat this is the grid it was
    *actually built on*, narrowing included — not the one originally investigated.
  - `topology`: how the grid's edges join — see [`EdgeTopology`](@ref).
  - `verbosity`: `:silent` announces nothing; anything else announces each filled-in value.

```julia
toy   = build_habitat()                         # every input defaulted, and said so
wetter = build_habitat(toy, supply = rainier)   # same grid, same regime, new supply
```
"""
function build_habitat(source = DefaultEcosystem();
                       # **`verbosity` must be declared FIRST, because the four below reference
                       # it**: a keyword default may only name keywords declared before it. Defaults
                       # evaluate in declaration order, and one passed explicitly is never evaluated
                       # — so naming a keyword suppresses exactly its own message.
                       verbosity::Symbol = :normal,
                       regime = _getdefaultvalue(AbstractRegime, source,
                                                 verbosity),
                       supply = _getdefaultvalue(AbstractSupply, source,
                                                 verbosity),
                       area = _getdefaultvalue(StudyArea, source, verbosity),
                       topology = _getdefaultvalue(AbstractTopology, source,
                                                   verbosity))
    return GridHabitat(; regime, supply, area, topology)
end

"""
    build_species(numspecies::Integer; tolerance, toleranceaxis, demand, demandaxis, dispersal = 10.0km,
        pthresh = 1.0e-9, movement = BirthOnlyMovement, disperse_safely = true, birth = 0.6/year,
        death = 0.6/year, longevity = 1.0, survival = 0.2, boost = 1.0,
        abundance = 1000 * numspecies, native = true, seed = nothing)

Build a `SpeciesList` of `numspecies` species. `tolerance` (the environmental **Condition** a
species is matched to) and `demand` (the **Resource** it consumes) are **required**, and so is the
[`NicheAxis`](@ref) each is measured on — `toleranceaxis` and `demandaxis`. Omitting any of the four
errors (use `build_species(DefaultEcosystem(); …)` to fill omissions with announced defaults). Every per-species keyword accepts either a scalar (applied to all species) or a
length-`numspecies` vector (validated, with a clear error naming the argument on a length
mismatch):

  - `tolerance = (mean, width)` — Gaussian niche tolerance (a [`NicheTolerance`](@ref) with a `Normal`
    response). A tuple of `(mean, width)` pairs, e.g. `((298.0K, 2.0K), (50.0mm, 10.0mm))`, gives one
    tolerance per regime (a [`SpeciesRequirementCollection`](@ref)), to match a
    multi-variable environment.
  - `toleranceaxis` — the [`NicheAxis`](@ref) the tolerance is on. A `NicheTolerance`'s unit comes
    from its axis, so the axis is what says which environmental variable the species responds to and
    in what unit its `(mean, width)` are read. For a per-regime tolerance tuple, pass a matching
    tuple of axes.
  - `demand` — per-species resource demand rate, in the axis's own unit
    (`Demand{SolarRadiation}` is `kJ/day`, `Demand{Precipitation}` `L/day`,
    `Demand{CarbonFlux}` `g/day`); any scale of the right dimension is converted. A tuple of two
    demands (e.g. `(10.0kJ/day, 5.0L/day)`) gives a [`SpeciesRequirementCollection`](@ref) to match a
    two-resource supply.
  - `demandaxis` — the [`NicheAxis`](@ref) the demand is on, matching `demand`'s arity as
    `toleranceaxis` matches `tolerance`'s.

!!! note "The axis is named, never inferred"
    Both axes used to be guessable: `axis` defaulted to `Temperature`, and a demand's type was
    chosen from its **unit** (`kJ/day` → solar, and so on). Neither is true now. A unit cannot say
    what a value *means* — `m/s` and `mm/day` share a dimension, which is how a wind-speed layer
    once became a water supply — so the axis is the single declaration, on both sides, and it is
    required.
  - `dispersal` — mean Gaussian dispersal distance, cut off below `pthresh`.
  - `movement` — the movement type: [`BirthOnlyMovement`](@ref) (default), [`AlwaysMovement`](@ref)
    or [`NoMovement`](@ref).
  - `disperse_safely` — per species, what becomes of an individual whose dispersal is aimed at a
    **dead cell** (off the grid, or an inactive one): `true` (the default) redistributes it among
    the reachable destinations, `false` loses it. A property of the *disperser*, so it may differ
    between species — a wind-dispersed seed blown out to sea is gone, an animal-dispersed one is
    not. The **grid's** topology (whether the world wraps) is a different question, answered by
    `GridHabitat(…, topology = …)`.
  - `birth`, `death`, `longevity`, `survival`, `boost` — demographic rates
    (scalar rates give an [`EqualPop`](@ref), vectors a [`PopGrowth`](@ref)).

`abundance` is either a total number of individuals split at random across
species (seedable via `seed`) or an explicit per-species vector. `native` marks
species as native (default all `true`).
"""
function build_species(numspecies::Integer;
                       tolerance = _require(:tolerance),
                       # Required *unless* `tolerance` is already built — a pre-built one
                       # carries its own axis, so demanding a second statement of it would be
                       # ceremony that could also disagree.
                       toleranceaxis = nothing,
                       demand = _require(:demand),
                       demandaxis = _require(:demandaxis),
                       dispersal = 10.0km,
                       pthresh::Real = 1.0e-9,
                       movement::Type{<:AbstractMovement} = BirthOnlyMovement,
                       disperse_safely = true,
                       birth = 0.6 / year,
                       death = 0.6 / year,
                       longevity::Real = 1.0,
                       survival::Real = 0.2,
                       boost::Real = 1.0,
                       abundance = 1000 * numspecies,
                       native = true,
                       seed = nothing)
    n = Int64(numspecies)
    _prebuilt = tolerance isa AbstractTolerance ||
                (tolerance isa Union{Tuple, NamedTuple} &&
                 all(t -> t isa AbstractTolerance, values(tolerance)))
    (_prebuilt || !isnothing(toleranceaxis)) || _require(:toleranceaxis)

    # **A pre-built tolerance must be checked FIRST**, because the branches below reach inside
    # `tolerance` with `first(…)`, and on a `NicheTolerance` object that gives
    # `MethodError: no method matching iterate(::NicheTolerance{…})` — an error naming `iterate`
    # rather than the mistake.
    # **And accepting one is better than merely reporting it well**, because it is the only way to
    # build a species against a layer whose frame `build_species` cannot otherwise express: a
    # tolerance on `NicheAxis` carrying a real unit, say, which `NicheTolerance(…, support = …)`
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
        # `tolerance = (temperature = …, rainfall = …)` gives a named collection
        # however `axis` was written.
        built = _zipmap(axes, values(tolerance)) do ax, tol
            return _nichetolerance(ax, tol, n)
        end
        SpeciesRequirementCollection(_asbacking(tolerance, built))
    else
        _nichetolerance(toleranceaxis, tolerance, n)
    end

    kernels = GaussianKernel.(_tofield(dispersal, n, "dispersal"),
                              Float64(pthresh))
    move = movement(kernels, _tofield(disperse_safely, n, "disperse_safely"))

    # A tuple `demand` gives one demand per resource → a collection; the backing is put back
    # so a named `demand` gives a named collection. `demandaxis` is zipped against it exactly as
    # `toleranceaxis` is against `tolerance` — including `_zipmap`'s arity check, so a mismatched
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
    param = _params(birth, death, longevity, survival, boost, n)
    nat = _tofield(native, n, "native")
    abun = _abundances(abundance, n, seed)

    return SpeciesList(n, traits, abun, demands, move, param, nat)
end

# As above for `build_species`: default `numspecies`, `tolerance` and `demand` if omitted, announce,
# then delegate. Any other keyword (dispersal, movement, rates, abundance, …) passes straight through.
function build_species(::DefaultEcosystem; numspecies = nothing,
                       tolerance = nothing, toleranceaxis = nothing,
                       demand = nothing, demandaxis = nothing,
                       verbosity::Symbol = :normal, kw...)
    if isnothing(numspecies)
        numspecies = 10
        _announce(verbosity,
                  "build_species(DefaultEcosystem()): `numspecies` defaulted to 10.")
    end
    if isnothing(tolerance)
        tolerance = (298.0K, 2.0K)
        _announce(verbosity,
                  "build_species(DefaultEcosystem()): `tolerance` defaulted to a (298.0 K, 2.0 K) temperature niche.")
    end
    if isnothing(toleranceaxis)
        toleranceaxis = Temperature
        _announce(verbosity,
                  "build_species(DefaultEcosystem()): `toleranceaxis` defaulted to `Temperature`.")
    end
    if isnothing(demand)
        demand = 10.0 * canonicalunit(Resource, SolarRadiation)
        _announce(verbosity,
                  "build_species(DefaultEcosystem()): `demand` defaulted to 10× the canonical solar rate.")
    end
    if isnothing(demandaxis)
        demandaxis = SolarRadiation
        _announce(verbosity,
                  "build_species(DefaultEcosystem()): `demandaxis` defaulted to `SolarRadiation`.")
    end
    return build_species(numspecies; tolerance = tolerance,
                         toleranceaxis = toleranceaxis, demand = demand,
                         demandaxis = demandaxis, kw...)
end

"""
    build_ecosystem(species::SpeciesList, environment::GridHabitat;
        nichefit = nothing, seed = nothing, distributed = :auto, epoch = nothing)

Assemble an ecosystem from a `species` list and an `environment`. When
`nichefit` is not given it is inferred from the trait type (`NicheTolerance` →
`NicheSuitability`, any `AbstractCategoricalTolerance` → `CategoricalSuitability`). Checks that the
species' resource demand matches the environment's supply before building,
and passes `seed` through for reproducibility.

`distributed` selects the ecosystem type: `:auto` (the default) builds a
distributed `MPIEcosystem` when the MPI extension is loaded and the process is
running with more than one rank (`MPI.Init()` already called), and a serial
`Ecosystem` otherwise; `true` forces an `MPIEcosystem` (erroring if the MPI
extension is not loaded); `false` forces a serial `Ecosystem`. Note that
auto-selection only covers construction and `simulate!` — under MPI the
abundances are partitioned across ranks, so collecting results still needs
[`gatherabundance`](@ref)/[`gatherdiversity`](@ref).

`epoch` is the real date the run begins — the date elapsed time zero means. It fixes the run's
calendar for output (see [`simulationdate`](@ref)) and, more importantly, *phases* every series the
environment carries: a [`DatedSeries`](@ref) starts at the slice covering the epoch, and a
[`MonthOfYearSeries`](@ref) climatology at the slice for the epoch's own month, so a run beginning in
July starts in July rather than in January.

It is resolved the way [`StudyArea`](@ref) resolves a CRS — adopt if unambiguous, ask if not. An
explicit `epoch` always wins; otherwise the environment's own series supply it, and if exactly one
real start date is found it is used and everything else is phased to it. Series that disagree are an
error naming the candidates, and an environment with no dated series has no epoch at all, which is
the behaviour of a run that never mentions dates.
"""
function build_ecosystem(species::SpeciesList, environment::GridHabitat;
                         nichefit = nothing, seed = nothing,
                         distributed = :auto,
                         epoch::Union{Nothing, Dates.TimeType} = nothing)
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
    # existing semantics here — an `Ecosystem` shares its habitat with the caller rather than copying
    # it — and re-pointing is idempotent, since an origin is computed from the calendar and the epoch
    # rather than accumulated.
    resolved = _resolveepoch(environment, epoch)
    _repointseries!(environment.regime, resolved)
    _repointseries!(environment.supply, resolved)
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
    return eco
end

# Complete the trio: build the default environment and default species (each announcing its own
# defaults) and assemble them into a runnable `Ecosystem`. `seed`/`distributed` pass to the assembly.
function build_ecosystem(::DefaultEcosystem; seed = nothing,
                         distributed = :auto,
                         epoch::Union{Nothing, Dates.TimeType} = nothing)
    environment = build_habitat()
    species = build_species(DefaultEcosystem())
    return build_ecosystem(species, environment, seed = seed,
                           distributed = distributed, epoch = epoch)
end

# ══ Running it ═════════════════════════════════════════════════════════════════════════════════

"""
    simulate!(eco::AbstractEcosystem, times::Unitful.Time, timestep::Unitful.Time)

Run an ecosystem, `eco` for a specified length of time, `times`, for a
particular timestep, `timestep`.
"""
function simulate!(eco::AbstractEcosystem, duration::Unitful.Time,
                   timestep::Unitful.Time; intervention = nothing)
    checkcoverage(eco, duration, timestep)
    check_bounds(eco, duration, timestep)
    times = length((0s):timestep:duration)
    for i in 1:times
        update!(eco, timestep, intervention)
    end
end

"""
    simulate!(cache::CachedEcosystem, srt::Unitful.Time, timestep::Unitful.Time)

Run a cached ecosystem, `cache` at a specified timepoint, `srt`, for a
particular timestep, `timestep`.
"""
function simulate!(cache::CachedEcosystem, srt::Unitful.Time,
                   timestep::Unitful.Time)
    eco = Ecosystem{typeof(cache.habitat), typeof(cache.spplist),
                    typeof(cache.nichefit)}(copy(cache.abundances.matrix[Ti(At(srt))]),
                                            cache.spplist,
                                            cache.habitat,
                                            cache.nichefit,
                                            cache.lookup,
                                            cache.cache,
                                            cache.rngs,
                                            uconvert(s, float(srt)),
                                            cache.seed,
                                            cache.epoch)
    update!(eco, timestep)
    return cache.abundances.matrix[Ti(At(srt + timestep))] = eco.abundances
end

"""
    simulate!(eco::Ecosystem, times::Unitful.Time, timestep::Unitful.Time,
              cacheInterval::Unitful.Time, cacheFolder::String,
              scenario_name::String)

Run an ecosystem, `eco` for specified length of times, `duration`, for a
particular timestep, 'timestep'. A cache interval and folder/file name are
specified for saving output.
"""
function simulate!(eco::Ecosystem,
                   times::Unitful.Time,
                   timestep::Unitful.Time,
                   cacheInterval::Unitful.Time,
                   cacheFolder::String,
                   scenario_name::String)
    checkcoverage(eco, times, timestep)
    check_bounds(eco, times, timestep)
    time_seq = zero(times):timestep:times
    for i in eachindex(time_seq)
        update!(eco, timestep)
        # Save cache of abundances
        if mod(time_seq[i], cacheInterval) == zero(time_seq[i])
            @save joinpath(cacheFolder,
                           scenario_name *
                           (@sprintf "%02d.jld2" uconvert(NoUnits,
                                                          time_seq[i] /
                                                          cacheInterval))) abun=eco.abundances.matrix
        end
    end
end

"""
    generate_storage(eco::Ecosystem, times::Int64, reps::Int64;
                     maxspecies = length(eco.spplist.abun))

Allocate an integer array of shape `(maxspecies, gridSize, times, reps)` for
recording species abundances across the ecosystem `eco` over multiple timesteps
and replicate runs.

`maxspecies` defaults to the number of species `eco` starts with, which is what
a run that never gains one needs. **Raise it to make room for species that
arrive during the run** — an [`AddSpecies`](@ref) intervention, an invasion — and
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
              "$(length(eco.spplist.abun)) species — the recording has to hold at least what is " *
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

"""
    simulate_action!(action!::Function, eco::AbstractEcosystem, times::Unitful.Time,
                     interval::Unitful.Time, timestep::Unitful.Time;
                     intervention = nothing, offset = false)

Run an ecosystem `eco` up to time `times` in steps of `timestep`, calling the
user-supplied `action!` at regular intervals so that any periodic task can be
performed as the simulation proceeds — recording a quantity, logging progress,
applying a management intervention, checking a stopping condition, and so on.
This is the general engine behind the [`simulate_record!`](@ref) and
[`simulate_record_diversity!`](@ref) recorders; use it directly when you want to
do something they do not.

At each step the ecosystem is advanced with [`update!`](@ref), which applies any
[`Intervention`](@ref) passed as `intervention`. Whenever the elapsed time falls on a multiple of `interval`, `action!(counting)`
is called, where `counting` is the 1-based index of that occurrence (handy as a
storage slot when the action is recording); `interval` must be a whole multiple of
`timestep`. Anything the action needs to read or update — the ecosystem, an output
array, an external counter — is captured by the closure, typically written as a
`do` block:

```julia
totals = zeros(Int, length((0s):interval:times))
simulate_action!(eco, times, interval, timestep) do counting
    totals[counting] = sum(eco.abundances.matrix)
end
```

`offset` shifts the action grid to start at `timestep` rather than `0`, which drops
the first occurrence (one fewer action in total). Use it to make the number of
actions match a pre-allocated array (such as one from [`generate_storage`](@ref));
the built-in diversity recorders pass `offset = iseven(size(storage, 3))`.

Returns the ecosystem `eco`, now advanced to `times`.
"""
function simulate_action!(action!::F,
                          eco::AbstractEcosystem,
                          times::Unitful.Time,
                          interval::Unitful.Time,
                          timestep::Unitful.Time;
                          intervention = nothing,
                          offset = false) where {F <: Function}
    iszero(mod(interval, timestep)) ||
        error("Interval must be a multiple of timestep")
    checkcoverage(eco, times, timestep)
    check_bounds(eco, times, timestep)
    action_seq = offset ? (timestep:interval:times) : ((0s):interval:times)
    time_seq = offset ? (timestep:timestep:times) : ((0s):timestep:times)
    counting = 0
    for i in eachindex(time_seq)
        update!(eco, timestep, intervention)
        if time_seq[i] in action_seq
            counting += 1
            action!(counting)
        end
    end
    return eco
end

"""
    simulate_record!(storage::AbstractArray, eco::Ecosystem, times::Unitful.Time,
         interval::Unitful.Time, timestep::Unitful.Time)

Run an ecosystem, `eco` for a specified length of time, `times`, for a
particular timestep, `timestep`, recording abundances into `storage` at each
time interval `interval`.

Pre-allocate `storage` with [`generate_storage`](@ref)`(eco, ntimes, reps)`,
where `ntimes = length((0s):interval:times)` is the number of recordings.

An `intervention` keyword takes an [`Intervention`](@ref) or
[`InterventionSet`](@ref). If it
can add species ([`AddSpecies`](@ref)), size `storage` for them with
`generate_storage(eco, ntimes, reps, maxspecies = …)`: the array is allocated
before the run and cannot grow.

To record diversity rather than raw abundances, see
[`simulate_record_diversity!`](@ref); to perform an arbitrary action at regular
intervals via a callback, see [`simulate_action!`](@ref).
"""
function simulate_record!(storage::AbstractArray,
                          eco::Ecosystem,
                          times::Unitful.Time,
                          interval::Unitful.Time,
                          timestep::Unitful.Time;
                          intervention = nothing)
    iszero(mod(interval, timestep)) ||
        error("Interval must be a multiple of timestep")
    checkcoverage(eco, times, timestep)
    check_bounds(eco, times, timestep)
    record_seq = (0s):interval:times
    time_seq = (0s):timestep:times
    _record!(storage, eco, 1)
    counting = 1
    for i in 2:length(time_seq)
        update!(eco, timestep, intervention)
        if time_seq[i] in record_seq
            counting = counting + 1
            _record!(storage, eco, counting)
        end
    end
    return storage
end

"""
    simulate_record_diversity!(storage, eco, times, interval, timestep,
                               divfun, qs::Vector{Float64})
    simulate_record_diversity!(substorage, metastorage, eco, times, interval, timestep,
                               qs::Vector{Float64})
    simulate_record_diversity!(storage, eco, times, interval, timestep,
                               divfuns::Array{Function}, q::Float64)

Run an ecosystem `eco` up to `times` in steps of `timestep`, recording diversity
into `storage` (and, for the alpha/beta/gamma form, `substorage`/`metastorage`) every `interval`,
which must be a whole multiple of `timestep`. These are all thin wrappers over
[`simulate_action!`](@ref) — see it for the recording mechanics — and differ only
in what diversity they record:

  - `divfun, qs` — a single diversity function `divfun` (which returns a
    `DataFrame` with a `:diversity` column) evaluated over the diversity orders
    `qs`, reshaped into `storage`;
  - `substorage, metastorage, …, qs` — normalised alpha, normalised beta and gamma
    diversity over `qs`; subcommunity-level values are written to `substorage`
    (gridSize × 3 × timepoints × qs) and metacommunity-level values to `metastorage`
    (3 × timepoints × qs). This form returns **both**, as the named tuple
    `(subcommunity = substorage, metacommunity = metastorage)`, so a caller need not
    remember which of the two came first;
  - `divfuns, q` — several diversity functions at a single diversity order `q`,
    one per column of `storage`.

For the `divfun`/`divfuns` forms, pre-allocate `storage` with
[`generate_storage`](@ref)`(eco, ncols, ntimes, reps)`, where `ncols` is
`length(qs)` (or `length(divfuns)`) and `ntimes = length((0s):interval:times)`.
"""
function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    divfun::F,
                                    qs::Vector{Float64}) where {F <: Function}
    simulate_action!(eco, times, interval, timestep,
                     offset = iseven(size(storage, 3))) do counting
        diversity = divfun(eco, qs)[!, :diversity]
        return storage[:, :, counting] = reshape(diversity,
                                                 Int(length(diversity) /
                                                     length(qs)),
                                                 length(qs))
    end
    return storage
end

function simulate_record_diversity!(substorage::AbstractArray,
                                    metastorage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    qs::Vector{Float64})
    simulate_action!(eco, times, interval, timestep,
                     offset = iseven(size(substorage, 3))) do counting
        measures = [NormalisedAlpha, NormalisedBeta, Gamma]
        for (i, msr) in enumerate(measures)
            dm = msr(eco)
            diversity = subdiv(dm, qs)[!, :diversity]
            diversity2 = metadiv(dm, qs)[!, :diversity]
            substorage[:, :, i, counting] = reshape(diversity,
                                                    Int(length(diversity) /
                                                        length(qs)),
                                                    length(qs))
            metastorage[:, i, counting] = diversity2
        end
    end
    return (subcommunity = substorage, metacommunity = metastorage)
end

function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    divfuns::Array{Function},
                                    q::Float64)
    simulate_action!(eco, times, interval, timestep) do counting
        # `j` is a position: it addresses `storage`, allocated by `generate_storage`, as well as
        # picking the measure.
        for (j, divfun) in enumerate(divfuns)
            storage[:, j, counting] .= divfun(eco, q)[!, :diversity][1]
        end
    end
    return storage
end
