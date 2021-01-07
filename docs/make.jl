using Documenter
using Simulation

makedocs(
    modules = [Simulation],
    sitename = "Simulation.jl",
    format=Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
    ],
    strict=true,
    checkdocs=:none,
)

deploydocs(
    repo = "github.com/boydorr/Simulation.jl.git",
    push_preview = true,
)
