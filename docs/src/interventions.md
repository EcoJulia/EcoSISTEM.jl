```@meta
CurrentModule = EcoSISTEM
```

# Interventions

An **intervention** changes the ecosystem: it destroys habitat, moves individuals, converts land,
or alters what a layer will do next. That is a different thing from a
[layer change](@ref "Layers that change over time"), which changes one layer's *values* as a pure
function of elapsed time — and the two are deliberately separate mechanisms.

| | changes | is a function of | applied |
|---|---|---|---|
| [`AbstractLayerChange`](@ref) | one layer's values | elapsed time only | redundantly on every MPI rank |
| [`Intervention`](@ref) | the ecosystem | a schedule | once, identically everywhere |

An intervention answers three questions, and each is a type rather than a callback — so what one
does is visible without running it:

```@setup iv
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

area = StudyArea(extent = (10km, 10km), cellsize = 1km, verbosity = :silent)
makeeco() = begin
    env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                            supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                                 axis = SolarRadiation),
                            area = area)
    species = build_species(5, tolerance = (285.0K, 5.0K), toleranceaxis = Temperature, demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                            abundance = 10_000, seed = 1)
    build_ecosystem(species, env, seed = 1)
end
```

```@example iv
eco = makeeco()
simulate!(eco, 5year, 1month_mean_duration,
          intervention = Intervention(AtTime(2year), RandomCells(20), Deactivate()))
count(parent(eco.habitat.active))          # 20 cells destroyed
```

## When — the schedule

[`EveryStep`](@ref), [`AtTime`](@ref), [`AtTimes`](@ref), [`BetweenTimes`](@ref) and
[`NeverScheduled`](@ref) (for disabling one without removing it).

A one-off schedule fires on the step that **reaches** its instant, not on one that equals it:
elapsed time accumulates as a float and a run's steps need not land on the instant exactly, so an
equality test would silently never fire.

## Where — the region

[`AllCells`](@ref), [`ActiveCells`](@ref), [`CellMask`](@ref) (your own boolean grid),
[`RandomCells`](@ref) and [`SpreadingCells`](@ref) (a contiguous cluster).

A count may be **exact or a rate**. A rate means each candidate is taken independently over the
step — a binomial draw — which is how habitat is actually lost:

```@example iv
eco = makeeco()
simulate!(eco, 10year, 1month_mean_duration,
          intervention = Intervention(EveryStep(), RandomCells(0.05 / year), Deactivate()))
100 - count(parent(eco.habitat.active))    # cells lost, drawn not dictated
```

## What — the operations

A **closed set**: [`Deactivate`](@ref), [`Reactivate`](@ref), [`SetLandCover`](@ref),
[`SetChange`](@ref), [`AddAbundance`](@ref) and [`RemoveAbundance`](@ref). Closed on purpose — a
callback could do anything, including the things that break reproducibility and MPI, whereas six
named operations can each be checked once and then trusted.

**`Deactivate` kills what lives there.** Destroying a cell is not a pause: a deactivated cell is
skipped by the simulation entirely, so anything left in it would neither breed nor die — a frozen
population with no ecological meaning. [`Reactivate`](@ref) does not bring them back either; it makes
the cell habitable again so that dispersal can recolonise it, as vegetation returns to a slag heap
once it stops being used. To restock deliberately, add [`AddAbundance`](@ref) alongside it.

### Several operations, one region

Operations after the first share the **same resolved region**, applied in order — clear ground and
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

## Several interventions

[`InterventionSet`](@ref) applies them in the order written, every step.

## Reproducibility

Selections come from a **counter-based** stream — `hash((seed, :intervention, k, step))` — which
generalises the per-species scheme. So a run replays exactly, every MPI rank and thread computes the
same selection without communicating, and species streams stay reserved for birth, death and
dispersal, meaning adding an intervention cannot re-phase the demography.

## Ordering

Interventions run **inside** [`update!`](@ref): after the population dynamics, after the clock
advances, and **before** the layer update. So a [`SetChange`](@ref) installed on a step takes effect
that same step rather than one step late.

## Worked examples

[`examples/interventions/`](https://github.com/EcoJulia/EcoSISTEM.jl/tree/main/examples/interventions) recreates a published set of scenarios — steady warming and drying,
a seasonal cycle, random and clustered habitat loss, land conversion, and generalist and specialist
invasions — and [`examples/interventions.jl`](https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/interventions.jl) runs them as part of the test suite.
[`examples/models.jl`](https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/models.jl) reuses the same scenarios verbatim, following them through time with a
full suite of diversity measures.

**Contrast with a layer change**, which is the *other* mechanism: [Time in EcoSISTEM](@ref) and
[`examples/VaryingClimate.jl`](https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/VaryingClimate.jl) cover it. The distinction is not stylistic. A layer change is
a pure function of elapsed time, so it can be applied redundantly on every MPI rank and still agree;
an intervention *mutates* the ecosystem, so it must be applied once and identically everywhere. That
is also why the operation set here is closed rather than a callback — a callback could draw from the
global RNG or write a layer's matrix, and either would silently desynchronise ranks.
