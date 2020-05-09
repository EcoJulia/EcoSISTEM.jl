# benchmark on files in "../test/canonical/"
using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "PkgBenchmark"])  # for benchmarking
Pkg.add(["JLSO", "StatsBase", "Unitful", "Test"])  # for running this file and examples
Pkg.resolve()

using BenchmarkTools
using JLSO
using Simulation

# whether or not to regenerate simulate inputs
gen_epi = false

const SUITE = BenchmarkGroup()
# - absolute path to various folders
const PATH_TO_REPO = joinpath(@__DIR__, "../")
const PATH_TO_EXAMPLES = joinpath(@__DIR__, "../test/canonical/")
# folder that stores precomputed inputs to simulate
const PATH_TO_EPI = joinpath(@__DIR__, "episystem", "numthread_$(Threads.nthreads())/")


##### functions #####
"""
    save_epi(file, case_name)

Save the EpiSystem and other inputs needed for simulate as JLSO files
"""
function save_epi(file, case_name)
    println("save episystem for case: ", case_name)
    # prepare saving arguments that will be used in canonical tests through `simulate!`
    #   or `simulate_record!`. These variables are passed into files being included later
    # NOTE: declare `global` variables as `include` will always evaluate files in the global
    #   scope
    global do_save = true
    global save_path = normpath(joinpath(
        PATH_TO_EPI,
        case_name,
    ))
    # run file to save epi
    include(joinpath(PATH_TO_EXAMPLES, file))
end

"""
    bench_example(case_name)

Load precomputed inputs for simulate and benchmark simulate
"""
function bench_example(case_name)
    # load precomputed episystem and run benchmark
    for (root, dirs, files) in walkdir(joinpath(PATH_TO_EPI, case_name))
        # get the subfolders that contain jlso files to run simulate
        if isempty(dirs) && !isempty(files)
            reldir = relpath(root, PATH_TO_EPI)
            println("\n benchmark episystem: ", reldir)
            # load simulate inputs
            @time epi = Simulation.load(joinpath(root, "initial_system.jlso"), EpiSystem)
            @time configs = JLSO.load(joinpath(root, "configuration.jlso"))
            # run simulate
            SUITE[reldir] = if :storage in keys(configs)
                # run simulate_record!
                @benchmarkable $(simulate_record!)(
                    $(configs[:storage]),
                    $epi,
                    $(configs[:times]),
                    $(configs[:interval]),
                    $(configs[:timestep]),
                )
            else
                # run simulate!
                @benchmarkable $(simulate!)(
                    $epi,
                    $(configs[:duration]),
                    $(configs[:timestep]),
                )
            end
        end
    end
end
##### END of functions #####

for file in readdir(PATH_TO_EXAMPLES)
    case_name = file[length("test_")+1 : end-length(".jl")]
    # generate episystem and other simulate inputs if necessary
    if gen_epi || !isdir(PATH_TO_EPI, case_name)
        save_epi(file, case_name)
    end
    # run benchmark
    bench_example(case_name)
end
