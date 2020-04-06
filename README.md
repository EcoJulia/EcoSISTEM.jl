# Simulation

*Package for running dynamic ecosystem simulations*

![](test/examples/Simulation.gif)

## Summary

**Simulation** is a [Julia](http://www.julialang.org) package that
provides functionality for simulating species undergoing dynamic
biological processes such as birth, death, competition and dispersal, as well as
environmental changes in climate and habitat.

The package was primarily developed for global scale simulations of
plant biodiversity. The underlying model for this is described in the arXiv
paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale
responses of biodiversity to environmental and land-use change*. This has now been adapted to include basic functionality for spatially explicit, dynamic, stochastic SIR models. Future updates to the package functionality involve incorporating
age-structure and more complex epidemiological models.

This package is in beta now, so please raise an issue if you find any
problems. The code has been substantially optimised for speed and can be fully
run in parallel through use of multi-threading. Over 30 cores, a century worth of simulation of
South America for 80km grid cells and ~60,000 species takes roughly 12 hours. This has been extended to run
on HPC using MPI. Scaling up to 36 nodes using MPI worked with 95% efficiency.

[paper-url]: https://arxiv.org/abs/1911.12257
