# EcoSISTEM
[![][docs-main-img]][docs-main-url] [![][actions-img]][actions-url] [![][codecov-img]][codecov-url]
[![][docs-dev-img]][docs-dev-url]

*Package for running dynamic ecosystem simulations*

## Summary

**EcoSISTEM** (Ecosystem Simulation through Integrated Species-Trait Environment Modelling) is a [Julia](http://www.julialang.org) package that
provides functionality for simulating species undergoing dynamic
biological processes such as birth, death, competition and dispersal, as well as environmental changes in climate and habitat

The package was primarily developed for global scale simulations of
plant biodiversity. The underlying model for this is described in the arXiv
paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale
responses of biodiversity to environmental and land-use change*.

This package is in alpha now, so please raise an issue if you find any problems.

[paper-url]: https://arxiv.org/abs/1911.12257
[docs-main-img]: https://img.shields.io/badge/docs-main-blue.svg
[docs-main-url]: https://boydorr.github.io/EcoSISTEM.jl/main/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://boydorr.github.io/EcoSISTEM.jl/dev/
[actions-img]: https://github.com/boydorr/EcoSISTEM.jl/workflows/EcoSISTEM%20testing/badge.svg?branch=main
[actions-url]: https://github.com/boydorr/EcoSISTEM.jl/actions
[codecov-img]: https://codecov.io/gh/boydorr/EcoSISTEM.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/boydorr/EcoSISTEM.jl?branch=main
