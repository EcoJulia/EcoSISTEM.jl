# profile files in "../test/canonical/"
using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["StatProfilerHTML"])  # for profiling
Pkg.add(["JLSO", "StatsBase", "Unitful", "Test"])  # for running this file and examples
Pkg.resolve()

using StatProfilerHTML
using JLSO
using Simulation

include(joinpath(@__DIR__, "utils.jl"))

# whether or not to regenerate simulate inputs
gen_epi = false

# - absolute path to various folders
const PATH_TO_REPO = joinpath(@__DIR__, "..")
const PATH_TO_EXAMPLES = joinpath(@__DIR__, "..", "test", "canonical")
# folder that stores precomputed inputs to simulate
const PATH_TO_EPI = joinpath(@__DIR__, "episystem", "numthread_$(Threads.nthreads())")

canonical_file = ARGS[1]

function profile_example(case_name)
    for (root, dirs, files) in walkdir(joinpath(PATH_TO_EPI, case_name))
        # get the subfolders that contain jlso files to run simulate
        if isempty(dirs) && !isempty(files)
            reldir = relpath(root, PATH_TO_EPI)
            println("\n profile episystem: ", reldir)
            # load simulate inputs
            @time epi = Simulation.load(joinpath(root, "initial_system.jlso"), EpiSystem)
            @time configs = JLSO.load(joinpath(root, "configuration.jlso"))
            # run simulate
            if :storage in keys(configs)
                # run simulate_record!
                @profilehtml simulate_record!(
                    configs[:storage],
                    epi,
                    configs[:times],
                    configs[:interval],
                    configs[:timestep],
                )
            else
                # run simulate!
                @profilehtml simulate!(
                    epi,
                    configs[:duration],
                    configs[:timestep],
                )
            end
        end
    end
end


case_name = canonical_file[length("test_")+1 : end-length(".jl")]
# generate episystem and other simulate inputs if necessary
if gen_epi || !isdir(PATH_TO_EPI, case_name)
    save_epi(
        joinpath(PATH_TO_EXAMPLES, canonical_file),
        normpath(joinpath(PATH_TO_EPI, case_name))
     )
end
# run benchmark
profile_example(case_name)
