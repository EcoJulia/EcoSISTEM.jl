# EcoSISTEM.jl

**EcoSISTEM** is a [Julia](http://www.julialang.org) package that provides functionality for simulating species undergoing dynamic biological processes such as birth, death, competition and dispersal, as well as environmental changes in climate and habitat.

The package was primarily developed for global scale simulations of plant biodiversity. The underlying model for this is described in the arXiv paper [arXiv:1911.12257 (q-bio.QM)][paper-url]
*Dynamic virtual ecosystems as a tool for detecting large-scale responses of biodiversity to environmental and land-use change*.

There are substantial changes to the package introduced through the [`dev`][dev-url] branch ([docs][docs-dev-url]), including epidemiological simulations and refactoring of the code base for further flexibility.

This package is in beta now, so please raise an issue if you find any problems. For more information on how to contribute, please read [our contributing guidelines][contrib-url].

![](Simulation.gif)

[paper-url]: https://arxiv.org/abs/1911.12257

[dev-url]: https://github.com/boydorr/EcoSISTEM.jl/tree/dev
[docs-dev-url]: https://boydorr.github.io/EcoSISTEM.jl/dev/

[contrib-url]: https://github.com/boydorr/EcoSISTEM.jl/blob/main/CONTRIBUTING.md
