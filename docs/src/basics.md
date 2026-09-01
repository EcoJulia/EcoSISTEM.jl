```@meta
CurrentModule = EcoSISTEM
```

# The basics of EcoSISTEM.jl

## Install

```julia
using Pkg
Pkg.add("EcoSISTEM")
```

## An ecosystem in three parts

An EcoSISTEM simulation runs on an [`Ecosystem`](@ref), assembled from three things:

| part | built by | holds |
| --- | --- | --- |
| the **environment** | [`GridHabitat`](@ref) | what each grid cell offers - conditions species are matched to, and resources they compete for |
| the **species** | [`build_species`](@ref) | what each species needs, how it moves, and how fast it breeds and dies |
| the **fit** between them | [`build_ecosystem`](@ref) | how well a species' tolerance suits a cell's conditions |

The grid itself is decided *first*, by a [`StudyArea`](@ref), and everything else is built
onto it. That ordering is deliberate: a layer introduced later can mark cells unusable, but it
can never move or resize the grid under a simulation that is already set up.

## A first simulation

Everything below runs as written. Start with the grid - here a synthetic one, 50 km square in
cells of 10 km, with no geography attached:

```@example basics
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

area = StudyArea(extent = (50km, 50km), cellsize = 10km, verbosity = :silent)
```

Now the environment. Every cell has the same temperature - a **condition**, which species are
matched to - and the same solar energy supply, a **resource**, which they compete for:

```@example basics
env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                        supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                             axis = SolarRadiation),
                        area = area)
```

Then the species: ten of them, each suited to temperatures around 285 K with a spread of 3 K,
each individual drawing energy at a fixed rate, and 10,000 individuals shared between them.

```@example basics
species = build_species(10,
                        tolerance = (285.0K, 3.0K), toleranceaxis = Temperature,
                        demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                        dispersal = 5.0km,
                        abundance = 10_000,
                        seed = 1)
```

Assemble, and run for ten years in monthly steps:

```@example basics
eco = build_ecosystem(species, env)
simulate!(eco, 10year, 1month_mean_duration)
sum(eco.abundances.matrix)          # individuals still alive
```

Passing a `seed` makes the run reproducible: it builds one random number stream per species,
so the result does not depend on how many threads or processes the simulation happens to use.

`eco.abundances.matrix` is species × cells. `eco.abundances.grid` is the same memory seen as
species × y × x, which is usually the more useful view:

```@example basics
size(eco.abundances.grid)
```

Both are plain arrays, so they index and reduce as you would expect. Alongside them,
`eco.abundances.dimmatrix` and `eco.abundances.dimgrid` are labelled views of that same memory,
carrying the species names and, for `dimgrid`, the grid's real coordinates - so you can ask for a
species by name or a cell by position rather than by index. They cost nothing to keep: the labels
are computed on demand, and no data is copied.

## Recording a run

[`simulate!`](@ref) leaves only the final state. To keep the whole time series, allocate
storage and use [`simulate_record!`](@ref), which records at an interval you choose:

```@example basics
times = 10year
interval = 1year
timestep = 1month_mean_duration

storage = generate_storage(eco, length((0year):interval:times), 1)
simulate_record!(storage, eco, times, interval, timestep)
size(storage)                       # species × cells × recordings × replicates
```

The `interval` must be a whole multiple of the `timestep` - see
[Time in EcoSISTEM](@ref) for why a year and a day do not qualify.

## Varying the model

Every per-species keyword takes either one value for all species or a vector with one entry
each, so a community of differing species is a matter of passing vectors:

```@example basics
varied = build_species(10,
                       tolerance = (collect(280.0:1.0:289.0)K, 3.0K), toleranceaxis = Temperature,
                       demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                       dispersal = collect(1.0:1.0:10.0)km,
                       abundance = 10_000,
                       seed = 1)
length(varied.tolerance.dists)
```

A multi-variable environment is a tuple - several conditions, several resources, or both -
and the species then carry a tolerance and a demand for each. Name the members and the two
sides are matched up by name; leave both sides positional and they are matched by position.
Naming one side and not the other is rejected rather than guessed at:

```@example basics
env2 = GridHabitat(regime = (temperature = UniformSpec(285.0K,
                                                            axis = Temperature),
                                   rain = UniformSpec(2.0mm / day,
                                                      axis = Precipitation)),
                         supply = (sunlight = UniformSpec(1.0e5kJ / (m^2 * day),
                                                          axis = SolarRadiation),
                                   water = UniformSpec(20.0mm / day,
                                                       axis = Precipitation)),
                         area = area)
species2 = build_species(10,
                         tolerance = (temperature = (285.0K, 3.0K),
                                      rain = (2.0mm / day, 1.0mm / day)),
                         toleranceaxis = (temperature = Temperature,
                                 rain = Precipitation),
                         demand = (sunlight = 1.0e9kJ / day,
                                   water = 1.0e6Unitful.L / day), demandaxis = (sunlight = SolarRadiation, water = Precipitation),
                         abundance = 10_000, seed = 1)
eco2 = build_ecosystem(species2, env2)
simulate!(eco2, 1year, 1month_mean_duration)
sum(eco2.abundances.matrix)
```

Real climate data goes in the same place, as a [`SourceSpec`](@ref) rather than a synthetic
spec - and the study area should be told about it, so that the grid comes from the data
instead of being invented:

```julia
area = StudyArea(regime = SourceSpec(WorldClim{BioClim}, :bio1),
                 within = LatLong(-10.0°, 3.0°, 50.0°, 60.0°),
                 cellsize = 10km, crs = EPSG(27700))
env = GridHabitat(regime = SourceSpec(WorldClim{BioClim}, :bio1),
                        supply = SourceSpec(WorldClim{Climate}, :srad),
                        area = area)
```

## Where to go next

- [How the model works](@ref) - the population model underneath: what regulates abundance, where
  carrying capacity comes from, and what `longevity` and `survival` each change. **Worth reading
  before interpreting any output.**
- [Layers, conditions and resources](@ref) - what a layer means, why temperature is a
  condition and light a resource, and how to find out what a dataset can be used for.
- [Time in EcoSISTEM](@ref) - the simulation clock, choosing a timestep, and environments
  that change as a simulation runs.
- [Integration with Diversity.jl](@ref) - measuring diversity on a built ecosystem.
- [Running at scale](@ref "Running a simulation at scale") - sizing a large run before building it,
  and spreading one across MPI ranks.
