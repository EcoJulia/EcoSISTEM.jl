# EcoSISTEM

| **Documentation** | **Build Status** | **DOI** |
|:-----------------:|:----------------:|:-------:|
| [![stable docs][docs-stable-img]][docs-stable-url] | [![build tests][actions-img]][actions-url] [![JuliaNightly][nightly-img]][nightly-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![dev docs][docs-dev-img]][docs-dev-url] | [![codecov][codecov-img]][codecov-url] | |

## Package for running dynamic ecosystem simulations

### Summary

**EcoSISTEM** (Ecosystem Simulation through Integrated Species-Trait Environment Modelling) is a [Julia](http://www.julialang.org) package that provides functionality for simulating species undergoing dynamic biological processes such as birth, death, competition and dispersal, as well as environmental changes in climate and habitat.

The package was primarily developed for global scale simulations of plant biodiversity. The underlying model for this is described in the arXiv paper [arXiv:1911.12257 (q-bio.QM)][paper-url],
*Dynamic virtual ecosystems as a tool for detecting large-scale
responses of biodiversity to environmental and land-use change*.

There are substantial changes to the package introduced through the [`dev`][dev-url] branch ([docs][docs-dev-url]), including epidemiological simulations and refactoring of the code base for further flexibility.

This package is in beta now, so please raise an issue if you find any problems. For more information on how to contribute, please read [our contributing guidelines](CONTRIBUTING.md). We are supported by NERC's Landscape Decisions [small][NERC-small] and [large][NERC-big] maths grants and an [EPSRC][EPSRC-stu] studentship.

## Introduction to EcoSISTEM

You can now run through a full introduction to EcoSISTEM with Pluto.jl if you have the source of the package cloned. To get started (if you are in the root of the package):

```julia
(EcoSISTEM) pkg> activate --temp
  Activating new project at `/var/folders/fv/1rqrvwq14ssggm_gc_x4v1cw0000gq/T/jl_IFxVXO`

(jl_IFxVXO) pkg> dev .
   Resolving package versions...
    Updating `/private/var/folders/fv/1rqrvwq14ssggm_gc_x4v1cw0000gq/T/jl_IFxVXO/Project.toml`
  [ed2dc23b] + EcoSISTEM
    Updating `/private/var/folders/fv/1rqrvwq14ssggm_gc_x4v1cw0000gq/T/jl_IFxVXO/Manifest.toml`

(jl_IFxVXO) pkg> add Pluto
   Resolving package versions...
    Updating `/private/var/folders/fv/1rqrvwq14ssggm_gc_x4v1cw0000gq/T/jl_IFxVXO/Project.toml`
  [c3e4b0f8] + Pluto
  No Changes to `/private/var/folders/fv/1rqrvwq14ssggm_gc_x4v1cw0000gq/T/jl_IFxVXO/Manifest.toml`

julia> import Pluto

julia> Pluto.run()
┌ Info: 
└ Opening http://localhost:1235/?secret=xxxxxxxx in your default browser... ~ have fun!
┌ Info: 
│ Press Ctrl+C in this terminal to stop Pluto
└ 
```

This should open a Pluto window in your browser - from there you can find `notebooks/Introduction.jl` in the `Open a notebook` box. You can also test `notebooks/InteractiveAfrica.jl` to see an invasive species colonising Africa. Note that it may be slow on first launch as it must install packages and in the latter case download climate data from the internet. If you are using a different Julia version you may need to add a block at the start of the examples to update the manifest:

```julia
begin
    using Pkg
    Pkg.update()
end
```

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

[dev-url]: https://github.com/EcoJulia/EcoSISTEM.jl/tree/dev
[NERC-small]: https://gtr.ukri.org/projects?ref=NE%2FT004193%2F1
[NERC-big]: https://gtr.ukri.org/projects?ref=NE%2FT010355%2F1
[EPSRC-stu]: https://gtr.ukri.org/projects?ref=EP%2FM506539%2F1
