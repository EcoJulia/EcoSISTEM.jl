# Simulation

*Package for running dynamic ecosystem simulations*

## Summary

**Simulation** is a [Julia](http://www.julialang.org) package that
provides functionality for simulating species undergoing dynamic 
biological processes such as birth, death, competition and dispersal, as well as 
environmental changes in climate and habitat.

The package was primarily developed for global scale simulations of 
plant biodiversity. The underlying model for this is described in the arXiv
paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale 
responses of biodiversity to environmental and land-use change*.
Future updates to the package functionality involve incorporating
age-structure and epidemiological models.

This package is in beta now, so please raise an issue if you find any
problems. The code has been substantially optimised for speed and can be fully
run in parallel through use of multi-threading. This can be extended to run 
on HPCC using MPI. Over 30 nodes, a century worth of simulation of 
South America for 80km grid cells and ~60,000 species takes roughly 12 hours.

[paper-url]: https://arxiv.org/abs/1911.12257