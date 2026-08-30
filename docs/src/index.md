# EcoSISTEM.jl

**EcoSISTEM** is a [Julia](http://www.julialang.org) package that provides functionality for simulating species undergoing dynamic biological processes such as birth, death, competition and dispersal, as well as environmental changes in climate and habitat.

The package was primarily developed for global scale simulations of plant biodiversity. The underlying model for this is described in the arXiv paper [arXiv:1911.12257 (q-bio.QM)](https://arxiv.org/abs/1911.12257)
*Dynamic virtual ecosystems as a tool for detecting large-scale responses of biodiversity to environmental and land-use change*. The paper is the reference for the model; [How the model works](@ref) is
the reference for what the package currently implements.

Species are regulated by **competition for a shared pool of resources**: in each cell the demands of
everything present are summed against what the cell supplies, so births fall and deaths rise as that
total approaches the supply. Carrying capacity is not a parameter anywhere in the package — it is
what those two pressures leave behind. Conditions such as temperature decide *where* a species can
persist; resources decide *how many*.

A simulation is assembled from three things — an environment, a list of species, and the fit
between them — on a grid decided in advance. [The basics](@ref "The basics of EcoSISTEM.jl")
walks through a first run and [How the model works](@ref) explains the population model underneath;
[Layers, conditions and resources](@ref) covers what an environmental layer means and where
climate data comes from; [Time in EcoSISTEM](@ref) covers the simulation clock and
environments that change as a run proceeds;
[Interventions](@ref) covers changing the ecosystem itself part way through a run;
[Integration with Diversity.jl](@ref) covers what to measure once it has run; and
[Running at scale](@ref "Running a simulation at scale") covers sizing a large run and
spreading one across MPI ranks. [Examples](@ref) is the shortest route to something running.

This package is in beta now, so please raise an issue if you find any problems. For more information on how to contribute, please read [our contributing guidelines](https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/CONTRIBUTING.md).
