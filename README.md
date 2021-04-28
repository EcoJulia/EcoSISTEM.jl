# EcoSISTEM

| **Documentation** | **Build Status** | **DOI** |
|:-----------------:|:----------------:|:-------:|
| [![stable docs][docs-stable-img]][docs-stable-url] | [![build tests][actions-img]][actions-url] [![JuliaNightly][nightly-img]][nightly-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![dev docs][docs-dev-img]][docs-dev-url] | [![codecov][codecov-img]][codecov-url] | |

## Package for running dynamic ecosystem and epidemiological simulations

**EcoSISTEM** (Ecosystem Simulation through Integrated Species-Trait Environment Modelling) is a [Julia](http://www.julialang.org) package that
provides functionality for simulating species undergoing dynamic
biological processes such as birth, death, competition and dispersal, as well as environmental changes in climate and habitat. This branch has now been adapted to include basic functionality for spatially explicit, dynamic, stochastic epidemiological models.

Examples for running small test compartmental simulations and scaled up versions (UK and Scotland sized) are available in `examples/Epidemiology`. Documentation of model structure can be found [here][model-struct-url] and ideas for model development [here][model-dev-url].

The package was primarily developed for global scale simulations of plant biodiversity. The underlying model for this is described in the arXiv paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale
responses of biodiversity to environmental and land-use change*.
Future updates to the package functionality involve incorporating
age-structure and epidemiological models (see the [SCRC fork](https://github.com/ScottishCovidResponse/Simulation.jl) for more details).

This package is in beta now, so please raise an issue if you find any problems. For more information on how to contribute, please read [our contributing guidelines](CONTRIBUTING.md).

[paper-url]: https://arxiv.org/abs/1911.12257

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://boydorr.github.io/EcoSISTEM.jl/stable/

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://boydorr.github.io/EcoSISTEM.jl/dev/

[actions-img]: https://github.com/boydorr/EcoSISTEM.jl/workflows/EcoSISTEM%20testing/badge.svg?branch=dev
[actions-url]: https://github.com/boydorr/EcoSISTEM.jl/actions

[nightly-img]: https://github.com/boydorr/EcoSISTEM.jl/actions/workflows/nightly.yaml/badge.svg?branch=dev
[nightly-url]: https://github.com/boydorr/EcoSISTEM.jl/actions/workflows/nightly.yaml?branch=dev

[codecov-img]: https://codecov.io/gh/boydorr/EcoSISTEM.jl/branch/dev/graph/badge.svg
[codecov-url]: https://codecov.io/gh/boydorr/EcoSISTEM.jl?branch=dev

[model-struct-url]: https://boydorr.github.io/EcoSISTEM.jl/dev/model_structure/
[model-dev-url]: https://boydorr.github.io/EcoSISTEM.jl/dev/model_structure/

[zenodo-img]: https://zenodo.org/badge/251665824.svg
[zenodo-url]: https://zenodo.org/badge/latestdoi/251665824
