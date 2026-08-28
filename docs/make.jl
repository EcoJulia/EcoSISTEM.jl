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
             "Biodiversity" => [
                 "Basics" => "basics.md",
                 "How the model works" => "model.md",
                 "Layers" => "layers.md",
                 "Axes, units and roles" => "units.md",
                 "Time" => "time.md",
                 "Interventions" => "interventions.md",
                 "Diversity" => "diversity.md",
                 "Running at scale" => "scale.md",
                 "Examples" => "examples.md",
                 "Africa" => "africa.md",
                 "Data Pipeline" => "pipeline.md"
             ],
             "API" => "api.md"
         ])

deploydocs(repo = "github.com/EcoJulia/EcoSISTEM.jl.git",
           push_preview = true,
           devbranch = "main",
           devurl = "main")
