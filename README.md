# EcoSISTEM

| **Documentation** | **Build Status** | **DOI** |
|:-----------------:|:----------------:|:-------:|
| [![stable docs][docs-stable-img]][docs-stable-url] | [![build tests][actions-img]][actions-url] [![JuliaNightly][nightly-img]][nightly-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![dev docs][docs-dev-img]][docs-dev-url] | [![codecov][codecov-img]][codecov-url] | |

## Package for running dynamic ecosystem simulations

**EcoSISTEM** (Ecosystem Simulation through Integrated Species-Trait Environment Modelling) is a
[Julia](http://www.julialang.org) package for simulating species undergoing birth, death,
competition and dispersal, under climate and habitat that change over time.

The package was primarily developed for global scale simulations of plant biodiversity. The
underlying model is described in the arXiv paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale responses of biodiversity to
environmental and land-use change*.

### What the model does

A simulation is a **spatial metacommunity**: abundances per species per cell on a grid, rather
than individuals. Each cell offers **conditions** a species is more or less suited to, and
**resources** every species present competes for from a shared pool. Births and deaths respond to
how the community's total demand compares with what a cell supplies, so a cell's carrying capacity
is an outcome of who is living there rather than a number set in advance — conditions decide where
a species can persist, and resources how many of it can.

The environment is built from layers, which may be uniform, drawn from real climate and land-cover
data, or made to change over the course of a run. See
[How the model works](https://docs.ecojulia.org/EcoSISTEM.jl/stable/model/) for the mechanism in
full.

### Getting started

```julia
using Pkg
Pkg.add("EcoSISTEM")
```

A first simulation needs a grid, an environment on it, a set of species, and a run:

```julia
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

# A synthetic grid, 50 km square, in cells of 10 km
area = StudyArea(extent = (50km, 50km), cellsize = 10km)

# One condition species are matched to, and one resource they compete for
env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                  supply = UniformSpec(1.0e5kJ / (m^2 * day), axis = SolarRadiation),
                  area = area)

# Ten species suited to 285 K, sharing 10,000 individuals
species = build_species(10,
                        tolerance = (285.0K, 3.0K), toleranceaxis = Temperature,
                        demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                        dispersal = 5.0km,
                        abundance = 10_000,
                        seed = 1)

eco = build_ecosystem(species, env)
simulate!(eco, 10year, 1month_mean_duration)

sum(eco.abundances.matrix)      # individuals still alive
```

The [documentation][docs-stable-url] walks through this in more detail, and covers what an
environmental layer means, where climate data comes from, how time works in a simulation, how to
intervene in a running ecosystem, and how to measure the diversity of the result.

### Interactive introduction

Two Pluto notebooks come with the source, in their own environment. With the repository cloned:

```julia
julia --project=notebooks -e 'import Pluto; Pluto.run()'
```

Then open `notebooks/Introduction.jl` for a guided tour, or `notebooks/InteractiveAfrica.jl` to
watch an invasive species colonise Africa. Both may be slow to start, as the second downloads
climate data.

### Contributing

This package is in beta now, so please raise an issue if you find any problems. For more
information on how to contribute, please read [our contributing guidelines](CONTRIBUTING.md). We
are supported by NERC's Landscape Decisions [small][NERC-small] and [large][NERC-big] maths grants
and an [EPSRC][EPSRC-stu] studentship.

[paper-url]: https://arxiv.org/abs/1911.12257

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://docs.ecojulia.org/EcoSISTEM.jl/stable/

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://docs.ecojulia.org/EcoSISTEM.jl/dev/

[actions-img]: https://github.com/EcoJulia/EcoSISTEM.jl/actions/workflows/testing.yaml/badge.svg?branch=main
[actions-url]: https://github.com/EcoJulia/EcoSISTEM.jl/actions/workflows/testing.yaml?branch=main

[nightly-img]: https://github.com/EcoJulia/EcoSISTEM.jl/actions/workflows/nightly.yaml/badge.svg?branch=main
[nightly-url]: https://github.com/EcoJulia/EcoSISTEM.jl/actions/workflows/nightly.yaml?branch=main

[codecov-img]: https://codecov.io/gh/EcoJulia/EcoSISTEM.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/EcoJulia/EcoSISTEM.jl?branch=main

[zenodo-img]: https://zenodo.org/badge/251665824.svg
[zenodo-url]: https://zenodo.org/badge/latestdoi/251665824

[NERC-small]: https://gtr.ukri.org/projects?ref=NE%2FT004193%2F1
[NERC-big]: https://gtr.ukri.org/projects?ref=NE%2FT010355%2F1
[EPSRC-stu]: https://gtr.ukri.org/projects?ref=EP%2FM506539%2F1
