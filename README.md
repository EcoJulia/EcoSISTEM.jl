# Simulation
[![][docs-dev-img]][docs-dev-url] [![][actions-img]][actions-url] [![][codecov-img]][codecov-url]

*Package for running dynamic ecosystem and epidemiological simulations*

## Summary

**Simulation** is a [Julia](http://www.julialang.org) package that
provides functionality for simulating species undergoing dynamic
biological processes such as birth, death, competition and dispersal, as well as environmental changes in climate and habitat. This branch has now been adapted to include basic functionality for spatially explicit, dynamic, stochastic epidemiological models.

Examples for running small test compartmental simulations and scaled up versions (UK and Scotland sized) are available in `examples/Epidemiology`. Documentation of model structure can be found [here][model-struct-url] and ideas for model development [here][model-dev-url].

The package was primarily developed for global scale simulations of
plant biodiversity. The underlying model for this is described in the arXiv
paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale
responses of biodiversity to environmental and land-use change*.
Future updates to the package functionality involve incorporating
age-structure and epidemiological models (see the [SCRC fork](https://github.com/ScottishCovidResponse/Simulation.jl) for more details).

This package is in alpha now, so please raise an issue if you find any problems.

[paper-url]: https://arxiv.org/abs/1911.12257
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://boydorr.github.io/Simulation.jl/dev/
[actions-img]: https://github.com/boydorr/Simulation.jl/workflows/Simulation%20testing/badge.svg
[actions-url]: https://github.com/boydorr/Simulation.jl/actions
[codecov-img]: https://codecov.io/gh/boydorr/Simulation.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/boydorr/Simulation.jl?branch=master
[model-struct-url]: https://boydorr.github.io/Simulation.jl/dev/model_structure/
[model-dev-url]: https://boydorr.github.io/Simulation.jl/dev/model_structure/
