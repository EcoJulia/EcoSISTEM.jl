# sort out dependencies
using Pkg
tempdir = mktempdir()
path_to_repo = joinpath(@__DIR__, "..")
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=path_to_repo))
Pkg.add(["BenchmarkTools", "PkgBenchmark"])  # for benchmarking
Pkg.add(["StatsBase", "Unitful", "Plots", "PlotlyJS", "ORCA"])  # for running experiments. NOTE: keep plotting pkgs otherwise some code in examples can't complie (e.g. @layout)
Pkg.resolve()

using BenchmarkTools

const SUITE = BenchmarkGroup()
const PATH_TO_REPO = "../"
const PATH_TO_EXAMPLES = "examples/Epidemiology/" # from repo

function run_example(file)
    include(file)
end

# set to repo path as some examples use the relative path from repo to load the data
cd(joinpath(@__DIR__, PATH_TO_REPO))
for file in readdir(joinpath(@__DIR__, PATH_TO_REPO, PATH_TO_EXAMPLES))
    # temporarily restrict to a few files only as a prototype
    if file in [
        "Benchmarking.jl",
        "Small_SIR.jl",
        "SEI2HRD_example.jl",
        # "UK_SIR.jl"
        ]
        SUITE[file[1:end - length(".jl")]] =
            @benchmarkable $(run_example)(joinpath(@__DIR__, PATH_TO_REPO, PATH_TO_EXAMPLES, $(file)))
    end
end
cd(@__DIR__)
