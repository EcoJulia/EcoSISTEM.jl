# SPDX-License-Identifier: LGPL-3.0-or-later

using Pkg
"EcoSISTEM" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("EcoSISTEM")
Pkg.develop(path = joinpath(@__DIR__, ".."))

using Documenter
using EcoSISTEM

DocMeta.setdocmeta!(EcoSISTEM, :DocTestSetup,
                    :(using EcoSISTEM); recursive = true)

makedocs(modules = [EcoSISTEM],
         sitename = "EcoSISTEM.jl",
         format = Documenter.HTML(canonical = "https://docs.ecojulia.org/EcoSISTEM.jl/stable",
                                  edit_link = "main",
                                  size_threshold_ignore = ["api.md"]),
         pages = [
             "Home" => "index.md",
             "Biodiversity" => [
                 "Basics" => "basics.md",
                 "Diversity" => "diversity.md",
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
