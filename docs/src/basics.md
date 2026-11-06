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

Then the species: a pool of ten, each suited to temperatures around 285 K with a spread of 3 K,
each individual drawing energy at a fixed rate and dispersing 5 km on average. Each trait is a spec,
and every trait not given - here the birth and death rates, `longevity` and `survival` - takes a
default, announced unless told `verbosity = :silent`, as the examples here are:

```@example basics
species = build_species(SpeciesPool(10),
                        tolerance = ToleranceSpec(Temperature, μ = 285.0K, σ = 3.0K),
                        demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                        movement = MovementSpec(InsectDispersal(), meandist = 5.0km),
                        seed = 1, verbosity = :silent)
```

Assemble, with 10,000 individuals shared between the species and placed by supply, and run for ten
years in monthly steps:

```@example basics
eco = build_ecosystem(BySupply(total = 10_000), species, env, seed = 1)
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

[`simulate!`](@ref) leaves only the final state. To keep the whole time series, hand it a
[`RecordAbundance`](@ref) and say how often to record with `every`:

```@example basics
times = 10year
interval = 1year
timestep = 1month_mean_duration

recording = RecordAbundance(eco, length((0year):interval:times))
simulate!(recording, eco, times, timestep, every = EveryInterval(interval))
size(recording.storage)             # species × cells × recordings × replicates
```

The first recording is the starting state and each later one the state at exactly that multiple of
the interval. A run of `times` in steps of `timestep` takes `times / timestep` steps and ends at
`times`, so it needs a slot for each multiple of the interval from zero to `times`.
[`RecordDiversity`](@ref) records a diversity measure the same way, and [`SaveAbundance`](@ref)
writes each recording to a file; each keeps the run's [`provenance`](@ref) as of its last write.

## Where the individuals start

[`build_ecosystem`](@ref) places the starting individuals by a **placement**, its first argument:
above, [`BySupply`](@ref) shares 10,000 between the species and spreads each one's share over the
cells by supply. Left out, the placement is `DefaultEcosystem()` - 1,000 individuals per species,
shared between them at random and placed by supply - and is announced unless told
`verbosity = :silent`. Here 5,000 are shared between the species, each spread over the cells by how
well each cell suits it:

```@example basics
eco3 = build_ecosystem(BySuitability(total = 5_000), species, env)
sum(eco3.abundances.matrix)
```

[`BySupply`](@ref), [`BySuitability`](@ref) and [`ByLayer`](@ref) (any layer of your own, such as a
vegetation index) each weight the cells, and [`Weighted`](@ref) multiplies several. Each takes a
community `total` or `each` species' own, as a count or as a density (`2.0/km^2`).
[`FromTable`](@ref) reads rows of species, place and count; [`FromOccurrences`](@ref) reads
occurrence records, one individual each by default; [`NoAbundance`](@ref) places nothing;
[`ToCapacity`](@ref) fills each cell to what its resources support, in the shares another placement
gives;
and placements add with `+`. Any function of the species, the habitat and a
[`PlacementContext`](@ref) is a placement too, written as a `do` block returning a table:

```@example basics
eco4 = build_ecosystem(species, env) do species, habitat, context
    return (species = ["1", "2"],
            cell = [CartesianIndex(1, 1), CartesianIndex(3, 3)],
            count = [100, 50])
end
sum(eco4.abundances.matrix)
```

Every row is checked before anything is written: a species the list does not hold, or a place off
the grid or in an inactive cell, is dropped and reported. The ecosystem keeps its placement, and
[`repopulate!`](@ref) places by it again.

## Varying the model

Every trait value is either one value for all species or read per species from the pool, never
taken by position. A pool made from a table holds a column per value, and a [`TraitColumn`](@ref)
reads one, so a community of differing species is a table of how they differ (see
[Species](@ref) for pools and the other trait values):

```@example basics
using DataFrames

pool = SpeciesPool(DataFrame(species = string.(1:10),
                             optimum = collect(280.0:1.0:289.0)K,
                             dispersal = collect(1.0:1.0:10.0)km))
varied = build_species(pool,
                       tolerance = ToleranceSpec(Temperature, μ = TraitColumn(:optimum),
                                                 σ = 3.0K),
                       demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                       movement = MovementSpec(InsectDispersal(),
                                               meandist = TraitColumn(:dispersal)),
                       seed = 1, verbosity = :silent)
length(varied.tolerance.dists)
```

A multi-variable environment is a tuple - several conditions, several resources, or both -
and the species then carry a tolerance and a demand for each, a tuple of specs. Name the members
and the two sides are matched up by name; leave both sides positional and they are matched by
position. Naming one side and not the other is rejected rather than guessed at:

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
species2 = build_species(SpeciesPool(10),
                         tolerance = (temperature = ToleranceSpec(Temperature,
                                                                  μ = 285.0K, σ = 3.0K),
                                      rain = ToleranceSpec(Precipitation,
                                                           μ = 2.0mm / day, σ = 1.0mm / day)),
                         demand = (sunlight = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                                   water = DemandSpec(Precipitation, 1.0e6Unitful.L / day)),
                         seed = 1, verbosity = :silent)
eco2 = build_ecosystem(BySupply(total = 10_000), species2, env2, seed = 1)
simulate!(eco2, 1year, 1month_mean_duration)
sum(eco2.abundances.matrix)
```

Real climate data goes in the same place, as a [`SourceSpec`](@ref) rather than a synthetic
spec - and the study area should be told about it, so that the grid comes from the data
instead of being invented:

```julia
area = StudyArea(regime = SourceSpec(WorldClim{BioClim}, :bio1),
                 within = EcoSISTEM.boundingbox("Scotland"),
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
