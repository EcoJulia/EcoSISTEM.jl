using Documenter
using EcoSISTEM

makedocs(modules = [EcoSISTEM],
         sitename = "EcoSISTEM.jl",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) ==
                                               "true",
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
