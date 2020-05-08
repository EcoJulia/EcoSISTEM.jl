# Simulation

*Package for running dynamic ecosystem and epidemiological simulations*

## Summary

**Simulation** is a [Julia](http://www.julialang.org) package that
provides functionality for simulating species undergoing dynamic
biological processes such as birth, death, competition and dispersal, as well as environmental changes in climate and habitat. This branch has now been adapted to include basic functionality for spatially explicit, dynamic, stochastic SIR models.

Examples for running small test SIR simulations and scaled up versions (UK and Scotland sized) are available in `examples/Epidemiology`. Documentation of model structure can be found [here](docs/Structure_gh.md) and ideas for model development [here](docs/Design.md).

The package was primarily developed for global scale simulations of
plant biodiversity. The underlying model for this is described in the arXiv paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale responses of biodiversity to environmental and land-use change*. Future updates to the package functionality involve incorporating
age-structure and more complex epidemiological models.

This package is in alpha now, so please raise an issue if you find any
problems.

[paper-url]: https://arxiv.org/abs/1911.12257

<img src="test/examples/ScotlandSIRSim.gif" width="300" height="400" />
