# SPDX-License-Identifier: LGPL-3.0-or-later

using Pkg
Pkg.resolve()

using Documenter
using EcoSISTEM
using RasterDataSources

DocMeta.setdocmeta!(EcoSISTEM, :DocTestSetup,
                    :(using EcoSISTEM), recursive = true)

makedocs(modules = [EcoSISTEM],
         sitename = "EcoSISTEM.jl",
         format = Documenter.HTML(canonical = "https://docs.ecojulia.org/EcoSISTEM.jl/stable",
                                  edit_link = "main",
                                  size_threshold_ignore = ["api.md"]),
         pages = [
             "Home" => "index.md",
             # Grouped by what a reader is trying to do. The single "Biodiversity" heading this
             # replaces was a leftover from when the repository had two arms and this was one of
             # them; everything had since accumulated under the survivor.
             "Getting started" => [
                 "Basics" => "basics.md",
                 "Examples" => "examples.md"
             ],
             "Building a model" => [
                 "Layers" => "layers.md",
                 "Axes, units and roles" => "units.md",
                 "Time" => "time.md",
                 "Interventions" => "interventions.md"
             ],
             "How the model works" => "model.md",
             "Analysing results" => ["Diversity" => "diversity.md"],
             "Going further" => [
                 "Running at scale" => "scale.md",
                 "Africa" => "africa.md",
                 "Data Pipeline" => "pipeline.md"
             ],
             "API" => "api.md"
         ])

deploydocs(repo = "github.com/EcoJulia/EcoSISTEM.jl.git",
           push_preview = true,
           devbranch = "main",
           # `dev` is the ecosystem convention, and what other people's links assume. This was
           # `main`, which left every link written the usual way pointing at a stale path.
           devurl = "dev")
