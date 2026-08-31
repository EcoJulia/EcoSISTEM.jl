```@meta
CurrentModule = EcoSISTEM
```

# Integration with Diversity.jl

EcoSISTEM is integrated with the [Diversity](https://github.com/EcoJulia/Diversity.jl)
package, so diversity measures can be calculated directly on an ecosystem - an
[`Ecosystem`](@ref) *is* a Diversity metacommunity, with one subcommunity per grid cell.

See [The basics of EcoSISTEM.jl](@ref) for how to set one up; the example below builds a small
one to measure.

```@setup diversity
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

area = StudyArea(extent = (50km, 50km), cellsize = 10km, verbosity = :silent)
env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                        supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                             axis = SolarRadiation),
                        area = area)
species = build_species(10, tolerance = (285.0K, 3.0K), toleranceaxis = Temperature, demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                        abundance = 10_000, seed = 1)
eco = build_ecosystem(species, env)
simulate!(eco, 1year, 1month_mean_duration)
```

**Every** grid cell is a subcommunity, including inactive ones - sea off a coastline, or a cell a
[`Deactivate`](@ref) intervention has destroyed. That costs nothing for the metacommunity measures,
which weight each subcommunity by its share of the total abundance and so give an empty cell zero
weight; the partitioned results are identical whether or not the empty cells are there. But a
**subcommunity** measure returns one row per cell, so on a small island in a large grid most of those
rows describe empty sea, and an index of an empty community is not defined. Select the occupied cells
before asking a per-cell question.

```@example diversity
using Diversity

norm_sub_alpha(eco, 1.0)      # subcommunity measures - one row per grid cell
```

```@example diversity
norm_meta_alpha(eco, 1.0)     # or metacommunity measures - one row for the landscape
```

Any measure takes several values of the viewpoint parameter `q` at once - `0` counts rare and
common species alike, and larger values weight towards the commonest:

```@example diversity
norm_sub_beta(eco, 0.0:3.0)
```
