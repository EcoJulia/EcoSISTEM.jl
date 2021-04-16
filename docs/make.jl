using Documenter
using EcoSISTEM

makedocs(
    modules = [EcoSISTEM],
    sitename = "EcoSISTEM.jl",
    format=Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Biodiversity" => [
        "Basics" => "basics.md",
        "Diversity" => "diversity.md",
        "Examples" => "examples.md",
        "Africa" => "africa.md"
        ]
    ],
    strict=true,
    checkdocs=:none,
)

deploydocs(
    repo = "github.com/boydorr/EcoSISTEM.jl.git",
    push_preview = true,
    devbranch = "main",
    devurl = "main",
)
