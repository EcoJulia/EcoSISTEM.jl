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
        "Biodiversity" => [
        "Basics" => "basics.md",
        "Diversity" => "diversity.md"
        ]
        "Epidemiology" => [
        "SCRC" => "epi.md",
        "Model Structure" => "model_structure.md",
        "Model Development" => "model_development.md",
        "API" => "api.md",
        "HPC" => "HPC.md",
        "Data" => "data.md"
        ]
    ],
    strict=true,
    checkdocs=:none,
)

deploydocs(
    repo = "github.com/boydorr/Simulation.jl.git",
    push_preview = true
)
