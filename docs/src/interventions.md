```@meta
CurrentModule = EcoSISTEM
```

# Interventions

An **intervention** changes the ecosystem: it destroys habitat, moves individuals, converts land,
or alters what a layer will do next. That is a different thing from a
[layer change](@ref "Layers that change over time"), which changes one layer's *values* as a pure
function of elapsed time - and the two are deliberately separate mechanisms.

| | changes | is a function of | applied |
|---|---|---|---|
| [`AbstractLayerChange`](@ref) | one layer's values | elapsed time only | redundantly on every MPI rank |
| [`Intervention`](@ref) | the ecosystem | a schedule | once, identically everywhere |

An intervention answers three questions, and each is a type rather than a callback - so what one
does is visible without running it:

```@setup iv
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

area = StudyArea(extent = (10km, 10km), cellsize = 1km, verbosity = :silent)
makeeco() = begin
    env = GridHabitat(regime = UniformSpec(10.0°C, axis = Temperature),
                            supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                                 axis = SolarRadiation),
                            area = area)
    species = build_species(SpeciesPool(5),
                            tolerance = ToleranceSpec(Temperature, μ = 10.0°C, σ = 5.0K),
                            demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                            seed = 1, verbosity = :silent)
    build_ecosystem(BySupply(total = 10_000), species, env, seed = 1)
end
```

```@example iv
eco = makeeco()
simulate!(eco, 5year, 1month_mean_duration,
          intervention = Intervention(AtTime(2year), RandomCells(20), Deactivate()))
count(parent(eco.habitat.active))          # 20 cells destroyed
```

## When - the schedule

[`EveryStep`](@ref), [`AtTime`](@ref), [`AtTimes`](@ref), [`BetweenTimes`](@ref),
[`EveryInterval`](@ref) (every whole multiple of an interval, counted from the start of the run),
[`AtDates`](@ref) and [`EveryYear`](@ref) (real dates, placed through the run's epoch and
`calendar`, so a run needs an epoch to use them) and [`NeverScheduled`](@ref) (for disabling one
without removing it).

A one-off schedule fires on the step that **reaches** its instant, not on one that equals it:
elapsed time accumulates as a float and a run's steps need not land on the instant exactly, so an
equality test would silently never fire. An instant at or before the start of the run acts on the
starting state, before the first step's births and deaths. `EveryStep` and `BetweenTimes` act over
the steps taken, so they first act after the first step.

A date schedule is held to the steps more strictly, because a date is meant to be a date. Each date
the run reaches must fall at the end of a step, or the step must be no longer than a day; otherwise
`simulate!` refuses before the first step, naming the date and the one it would have acted on. So
`EveryYear()` at steps of `month_mean_duration` needs `calendar = MeanMonths()` - under exact dates
1 January drifts against the steps and would be acted on in late January in three years out of four
- and `AtDates([Date(2000, 3, 15)])` needs daily steps, since no monthly step ends on the fifteenth.

## Where - the region

[`AllCells`](@ref), [`ActiveCells`](@ref), [`CellMask`](@ref) (your own boolean grid),
[`CellsAt`](@ref) (the cells at named places), [`RandomCells`](@ref) and [`SpreadingCells`](@ref)
(a contiguous cluster) - or a shape.

A site is usually known by where it is, not by its row and column. `CellsAt` takes places - a
[`LatLong`](@ref), a position in the grid's own frame, or a grid key - and refuses one off the grid
before the run starts. It takes the cells whether active or not, so a site can be reactivated by
naming it; with `snap`, a place that falls in an inactive cell (a release point on the shore of a
coarse grid, say) moves to the nearest active cell within that distance:

```julia
Intervention(AtTime(10year), CellsAt(LatLong(56.0°, -3.2°), snap = 5km),
             AddAbundance("Rhododendron", 500))
```

A shape is a region as it is a study area's `within`: a [`NaturalEarthSpec`](@ref), a
[`ShapeSpec`](@ref) of your own boundary file or a combination of them means the cells it covers at
least half of, and a [`ShapeMaskSpec`](@ref) says which cells under another rule. It is placed on the
habitat's own grid once, before the run's first step, so the grid needs a coordinate reference
system:

```julia
Intervention(AtTime(20year), ShapeMaskSpec(NaturalEarthSpec("Orkney Is."), AnyOverlap()),
             Deactivate())
```

A count may be **exact or a rate**. A rate means each candidate is taken independently over the
step - a binomial draw - which is how habitat is actually lost:

```@example iv
eco = makeeco()
simulate!(eco, 10year, 1month_mean_duration,
          intervention = Intervention(EveryStep(), RandomCells(0.05 / year), Deactivate()))
100 - count(parent(eco.habitat.active))    # cells lost, drawn not dictated
```

## What - the operations

A **closed set**: [`Deactivate`](@ref), [`Reactivate`](@ref), [`SetLandCover`](@ref),
[`SetChange`](@ref), [`AddAbundance`](@ref), [`RemoveAbundance`](@ref) and [`AddSpecies`](@ref).
Closed on purpose - a callback could do anything, including the things
that break reproducibility and MPI, whereas named operations can each be checked once and then
trusted.

**`Deactivate` kills what lives there.** Destroying a cell is not a pause: a deactivated cell is
skipped by the simulation entirely, so anything left in it would neither breed nor die - a frozen
population with no ecological meaning. [`Reactivate`](@ref) does not bring them back either; it makes
the cell habitable again so that dispersal can recolonise it, as vegetation returns to a slag heap
once it stops being used. To restock deliberately, add [`AddAbundance`](@ref) alongside it.

### Several operations, one region

Operations after the first share the **same resolved region**, applied in order - clear ground and
plant a crop on exactly the ground you cleared:

```@example iv
eco = makeeco()
simulate!(eco, 5year, 1month_mean_duration,
          intervention = Intervention(EveryStep(), RandomCells(0.02 / year),
                                      Deactivate(), AddAbundance(5, 500)))
converted = findall(.!vec(parent(eco.habitat.active)))
all(eco.abundances.matrix[5, converted] .>= 500)
```

Two interventions could **not** do this: each resolves its own region, and a random one draws its
own cells, so the crop would be planted somewhere other than the cleared ground.

### Adding individuals by a placement

`AddAbundance` also takes a placement - anything [`build_ecosystem`](@ref) takes first - and adds
the individuals it places over the region's active cells each time the intervention fires. Wrapped in
[`ToCapacity`](@ref) it tops each cell up to what its resources support, counting what already lives
there, and a species' name places that species alone, so sowing a crop among the weeds and hedgerows
of some fields is

```julia
Intervention(AtDates([Date(2030, 3, 1)]), CellMask(fields), AddAbundance(ToCapacity("wheat")))
```

A dated placement - a [`FromTable`](@ref) with a `time` column, or a [`FromOccurrences`](@ref) given
its records' `date` or `time` - is served by its dates instead: alone in an `Intervention` with
`EveryStep()`, its rows are added as the clock reaches them, those at or before
the start onto the starting state. It is read once for an ecosystem, over the cells the ecosystem was
built with; each step adds the rows falling due in the region's cells still active, and the run's end
notes what it did not add - rows dated after the end, rows a continuing run started past, and rows due
in cells deactivated since:

```julia
arrivals = Intervention(EveryStep(), ActiveCells(),
                        AddAbundance(FromOccurrences(records, date = (:year, :month))))
```

A placement may be a function of the user's, which the operations otherwise refuse, because a placement
only says where individuals go: it returns a table, which is checked before anything is written, and
under MPI the first rank alone calls it and every rank adds its own species' rows. Each firing gives it
random streams of its own, apart from the region's, so drawing only from its context's streams keeps
a run replayable.

### Adding individuals from a table

A table listing which species, which cell, how many and when seeds a population from dated records by
one operation rather than an intervention per species and month: [`FromTable`](@ref) reads it, a row's
place by its grid key (`CartesianIndex(y, x)`) or by `latitude` and `longitude`, and `AddAbundance`
serves it by its dates, the region filtering its rows - `ActiveCells()` leaves out any in an inactive
cell:

```@example iv
eco = makeeco()
records = (species = fill(eco.spplist.names[1], 3),
           cell = [CartesianIndex(1, 1), CartesianIndex(1, 1), CartesianIndex(2, 3)],
           time = [1.0, 1.0, 6.0] .* month_mean_duration)
simulate!(eco, 1year, 1month_mean_duration,
          intervention = Intervention(EveryStep(), ActiveCells(),
                                      AddAbundance(FromTable(records))))
```

`AddAbundanceTable` and `build_abundance_table` are deprecated in its favour: a table goes to
`FromTable` as it is, naming each species by name and each place by its grid key or location.

### Bringing a species into a run

Every species a run may hold is declared in its [`SpeciesPool`](@ref) before the run is built, so
one arriving later - an invader, a reintroduction, a crop sown for the first time - is declared with
traits of its own and left out of the run with [`deactivate_species!`](@ref). [`AddSpecies`](@ref)
brings it in, its first individuals placed over the region by a placement, as `AddAbundance` of a
placement places them, the arrivals the only species placed. Nothing about the run changes size, so
it works under MPI as everywhere else:

```@example iv
pool = SpeciesPool(6)
deactivate_species!(pool, "6", verbosity = :silent)     # declared, but not in the run
species = build_species(pool, tolerance = ToleranceSpec(Temperature, μ = 10.0°C, σ = 5.0K),
                        demand = DemandSpec(SolarRadiation, 1.0e9kJ / day), seed = 1,
                        verbosity = :silent)
habitat = GridHabitat(regime = UniformSpec(10.0°C, axis = Temperature),
                      supply = UniformSpec(1.0e5kJ / (m^2 * day), axis = SolarRadiation),
                      area = area)
eco = build_ecosystem(BySupply(total = 10_000), species, habitat, seed = 1)
simulate!(eco, 2year, 1month_mean_duration,
          intervention = Intervention(AtTime(1year), CellsAt(CartesianIndex(5, 5)),
                                      AddSpecies("6", BySupply(each = 200))))
sum(eco.abundances.matrix[6, :])         # the arrival, a year after it was introduced
```

[`RandomSpecies`](@ref) draws the arrivals from the species not yet in the run instead - a number,
or a rate at which each arrives - from the intervention's own stream. Names are looked up before the
run's first step, and one reaching no species of the run stops it unless a lower `level` is given; a
species already in the run gets its arrivals all the same, and the run's end says so.

## Several interventions

[`InterventionSet`](@ref) applies them in the order written, every step.

## Reproducibility

Selections come from a **counter-based** stream, seeded from `(seed, :intervention, k, step)`, which
generalises the per-species scheme. So a run replays exactly, every MPI rank and thread computes the
same selection without communicating, and species streams stay reserved for birth, death and
dispersal, meaning adding an intervention cannot re-phase the demography.

## Ordering

Interventions run **inside** [`update!`](@ref): after the population dynamics, after the clock
advances, and **before** the layer update. So a [`SetChange`](@ref) installed on a step takes effect
that same step rather than one step late.

## Worked examples

[`examples/interventions/`](https://github.com/EcoJulia/EcoSISTEM.jl/tree/main/examples/interventions) recreates a published set of scenarios - steady warming and drying,
a seasonal cycle, random and clustered habitat loss, land conversion, and generalist and specialist
invasions - and [`examples/interventions.jl`](https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/interventions.jl) runs them as part of the test suite.
[`examples/models.jl`](https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/models.jl) reuses the same scenarios verbatim, following them through time with a
full suite of diversity measures.

**Contrast with a layer change**, which is the *other* mechanism: [Time in EcoSISTEM](@ref) and
[`examples/VaryingClimate.jl`](https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/VaryingClimate.jl) cover it. The distinction is not stylistic. A layer change is
a pure function of elapsed time, so it can be applied redundantly on every MPI rank and still agree;
an intervention *mutates* the ecosystem, so it must be applied once and identically everywhere. That
is also why the operation set here is closed rather than a callback - a callback could draw from the
global RNG or write a layer's matrix, and either would silently desynchronise ranks.
