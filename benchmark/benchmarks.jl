# benchmark on files in "examples/Epidemiology/"

using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "PkgBenchmark"])  # for benchmarking
# NOTE: keep plotting pkgs otherwise some code in examples can't complie (e.g. @layout)
Pkg.add(["StatsBase", "Unitful", "Plots", "PlotlyJS", "ORCA"])  # for running experiments.
Pkg.resolve()

using BenchmarkTools

const SUITE = BenchmarkGroup()
const PATH_TO_REPO = "../"
const PATH_TO_EXAMPLES = "examples/Epidemiology/" # from repo

function run_example(file)
    # need to change the dir to repo as assumed by files in `/examples/Epidemiology/`,
    #   as `@benchmarkable` will automatically change the working dir to `repo/benchmark/`
    cd(joinpath(@__DIR__, PATH_TO_REPO))
    include(file)
end

# set to repo path as some examples use the relative path from repo to load the data
for file in readdir(joinpath(@__DIR__, PATH_TO_REPO, PATH_TO_EXAMPLES))
    SUITE[file[1:end - length(".jl")]] =
        @benchmarkable $(run_example)(joinpath(@__DIR__, PATH_TO_REPO, PATH_TO_EXAMPLES, $(file)))
end
