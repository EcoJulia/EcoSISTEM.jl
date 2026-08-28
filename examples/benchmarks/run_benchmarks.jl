# SPDX-License-Identifier: LGPL-3.0-or-later

# Benchmark orchestrator: compares EcoSISTEM's simulation speed across a small,
# fixed set of thread/process configurations chosen to isolate two questions on a
# machine with `N` logical cores (N = Sys.CPU_THREADS, override with
# ECOSISTEM_BENCH_MAXCORES):
#
#   1. threads vs processes at full core count:
#        * 1 process  x N threads         (pure threading, no MPI)
#        * N processes x 1 thread          (pure MPI, no threading)
#
#   2. the speedup from doubling one axis while holding the other, always going
#      from N/2 to N cores (so the ideal speedup is 2):
#        * processes 1 -> 2  at N/2 threads each   (1p N/2t  vs  2p N/2t)
#        * processes N/2 -> N at 1 thread each      (N/2p 1t  vs  Np 1t)
#        * threads   1 -> 2  at N/2 processes each  (N/2p 1t  vs  N/2p 2t)
#        * threads   N/2 -> N at 1 process          (1p N/2t  vs  1p Nt)
#
# The six distinct configurations that cover all of the above are:
#
#        1p x Nt   (threaded)   1p x (N/2)t (threaded)
#        Np x 1t   (mpi)        (N/2)p x 1t (mpi)
#        2p x (N/2)t (mpi)      (N/2)p x 2t (mpi)
#
# Because Julia fixes its thread count (`-t`) and MPI rank count (`-n`) at process
# launch, each configuration is timed in its own fresh subprocess
# (`benchmark_worker.jl`).
#
# The model is deliberately started with ~1 billion individuals (see
# benchmark_worker.jl) so it is compute-bound from the first timestep and a short
# run gives a stable estimate; the grid and species count are as large as before.
#
# The integer type storing species counts can be varied with
# ECOSISTEM_BENCH_COUNT_TYPE (default Int64); the orchestrator applies it via the
# EcoSISTEM `count_type` preference so the workers compile with that
# `EcoSISTEM.Count` (e.g. Int32 to halve the abundance memory), restoring the
# preference afterwards. Outputs are suffixed with the type so runs coexist.
#
# Outputs (in ECOSISTEM_BENCH_OUTDIR, default this directory): the summary table,
# the comparison report, benchmark_results_<Count>.csv, a per-configuration time
# bar chart and a doubling-speedup bar chart.
#
# Usage:
#     julia --project=examples examples/benchmarks/run_benchmarks.jl
#     ECOSISTEM_BENCH_COUNT_TYPE=Int32 julia --project=examples \
#         examples/benchmarks/run_benchmarks.jl

using Printf
using Plots
using TOML

const PROJECT = normpath(joinpath(@__DIR__, ".."))
const WORKER = joinpath(@__DIR__, "benchmark_worker.jl")
const OUTDIR = get(ENV, "ECOSISTEM_BENCH_OUTDIR", @__DIR__)
const JULIA = first(Base.julia_cmd())
const MPI_PKG = Base.PkgId(Base.UUID("da04e1cc-30fd-572f-bb4f-1f8673147195"),
                           "MPI")

# Integer type the workers should store species counts as (see EcoSISTEM.Count).
# Selected through the EcoSISTEM `count_type` preference below; the workers
# inherit this env var and assert they actually compiled with it.
const COUNT_TYPE = get(ENV, "ECOSISTEM_BENCH_COUNT_TYPE", "Int64")
const PREFS_FILE = joinpath(PROJECT, "LocalPreferences.toml")

# Load MPI lazily so a broken/unbuilt MPI does not stop the threaded runs.
get_mpi() = Base.require(MPI_PKG)

# Run `f` with the EcoSISTEM `count_type` preference set to `count_type` in the
# worker project's LocalPreferences.toml, restoring the file afterwards. `Count`
# is a compile-time type, so the preference (not an env var) is what makes the
# worker subprocesses compile with it; the first worker of a new type pays a
# recompile, the rest reuse the cache. Other packages' preferences are preserved.
function with_count_pref(f, count_type)
    backup = isfile(PREFS_FILE) ? read(PREFS_FILE, String) : nothing
    prefs = isnothing(backup) ? Dict{String, Any}() : TOML.parse(backup)
    get!(prefs, "EcoSISTEM", Dict{String, Any}())["count_type"] = count_type
    open(io -> TOML.print(io, prefs), PREFS_FILE, "w")
    try
        return f()
    finally
        isnothing(backup) ? (isfile(PREFS_FILE) && rm(PREFS_FILE)) :
        write(PREFS_FILE, backup)
    end
end

# Extract the "BENCH_RESULT,mode,procs,threads,seconds" line from worker output.
function parse_result(output::AbstractString)
    for line in eachsplit(output, '\n')
        if startswith(line, "BENCH_RESULT,")
            fields = split(line, ',')
            return parse(Float64, fields[5])
        end
    end
    return error("no BENCH_RESULT line found in worker output:\n$output")
end

# Run one threaded configuration in a fresh subprocess and return seconds taken.
function run_threaded(threads::Int)
    cmd = `$JULIA -t $threads --startup-file=no --project=$PROJECT $WORKER`
    cmd = addenv(cmd, "ECOSISTEM_BENCH_MODE" => "threaded")
    return parse_result(read(cmd, String))
end

# Run one MPI configuration (procs ranks x threads threads each) and return
# seconds taken, as reported by rank 0.
function run_mpi(procs::Int, threads::Int)
    # Use the no-argument `mpiexec()` form, which just returns the launcher path.
    # The do-block form runs the command inside a `withenv` that rewrites
    # DYLD_FALLBACK_LIBRARY_PATH to MPI's own libraries, which breaks the worker's
    # loading of other Julia artifacts (e.g. PROJ/SQLite) on macOS. invokelatest
    # is needed because MPI is loaded lazily via `get_mpi`.
    launcher = Base.invokelatest(get_mpi().mpiexec)
    cmd = `$launcher -n $procs $JULIA -t $threads --startup-file=no --project=$PROJECT $WORKER`
    cmd = addenv(cmd, "ECOSISTEM_BENCH_MODE" => "mpi")
    return parse_result(read(cmd, String))
end

# The fixed set of configurations to time, given `N` cores and `half = N ÷ 2`.
# Returned as (mode, procs, threads) tuples, de-duplicated (small N can collapse
# some of them, e.g. N == 2 gives half == 1).
function configurations(N::Int)
    half = max(1, N ÷ 2)
    wanted = [("threaded", 1, N),
        ("threaded", 1, half),
        ("mpi", N, 1),
        ("mpi", half, 1),
        ("mpi", 2, half),
        ("mpi", half, 2)]
    seen = Set{Tuple{String, Int, Int}}()
    return [c for c in wanted if !(c in seen) && (push!(seen, c); true)]
end

function main()
    maxcores = parse(Int,
                     get(ENV, "ECOSISTEM_BENCH_MAXCORES",
                         string(Sys.CPU_THREADS)))

    @info "Benchmarking on $maxcores cores" grid=get(ENV,
                                                     "ECOSISTEM_BENCH_GRID",
                                                     "100") species=get(ENV,
                                                                        "ECOSISTEM_BENCH_SPECIES",
                                                                        "10000") individuals=get(ENV,
                                                                                                 "ECOSISTEM_BENCH_INDIVIDUALS",
                                                                                                 "1000000000") years=get(ENV,
                                                                                                                         "ECOSISTEM_BENCH_YEARS",
                                                                                                                         "3")

    @info "Storing counts as $COUNT_TYPE (EcoSISTEM.Count)"

    # Each result: (mode, procs, threads, cores, seconds). seconds is NaN if the
    # configuration failed (e.g. MPI is not built). All runs share one count-type
    # preference so the workers compile with the requested `EcoSISTEM.Count`.
    results = with_count_pref(COUNT_TYPE) do
        rs = NamedTuple[]
        for (mode, procs, threads) in configurations(maxcores)
            @info "$mode run: $procs process(es) x $threads thread(s)"
            secs = try
                mode == "mpi" ? run_mpi(procs, threads) : run_threaded(threads)
            catch err
                @warn "$mode run with $procs process(es) x $threads " *
                      "thread(s) failed" exception=err
                NaN
            end
            push!(rs,
                  (mode = mode, procs = procs, threads = threads,
                   cores = procs * threads, seconds = secs))
        end
        return rs
    end

    report(results, maxcores)
    return nothing
end

# Seconds for a specific (mode, procs, threads) configuration, or NaN if it was
# not run or failed.
function timeof(results, mode, procs, threads)
    idx = findfirst(r -> r.mode == mode && r.procs == procs &&
                         r.threads == threads, results)
    return isnothing(idx) ? NaN : results[idx].seconds
end

function report(results, N::Int)
    half = max(1, N ÷ 2)

    println("\n", "="^68)
    println("EcoSISTEM thread-vs-process benchmark  ",
            "(N = $N cores, Count = $COUNT_TYPE)")
    println("="^68)

    # Per-configuration timing table.
    println("\nConfigurations:")
    @printf("%-10s %-6s %-8s %-6s %12s\n",
            "mode", "procs", "threads", "cores", "seconds")
    for r in results
        @printf("%-10s %-6d %-8d %-6d %12.3f\n",
                r.mode, r.procs, r.threads, r.cores, r.seconds)
    end

    # 1. threads vs processes at the full core count.
    t_threaded = timeof(results, "threaded", 1, N)
    t_mpi = timeof(results, "mpi", N, 1)
    println("\nThreads vs processes at $N cores:")
    @printf("  1 process  x %d threads (threaded): %8.3f s\n", N, t_threaded)
    @printf("  %d processes x 1 thread  (MPI)     : %8.3f s\n", N, t_mpi)
    if !isnan(t_threaded) && !isnan(t_mpi)
        @printf("  -> MPI is %.2fx the speed of threading (threaded/MPI time)\n",
                t_threaded / t_mpi)
    end

    # 2. doubling one axis from N/2 to N cores; ideal speedup is 2.
    # Each entry: (label, before (mode,procs,threads), after (mode,procs,threads)).
    doublings = [
        ("processes 1->2 (N/2 threads each)",
         ("threaded", 1, half), ("mpi", 2, half)),
        ("processes N/2->N (1 thread each)",
         ("mpi", half, 1), ("mpi", N, 1)),
        ("threads 1->2 (N/2 processes each)",
         ("mpi", half, 1), ("mpi", half, 2)),
        ("threads N/2->N (1 process)",
         ("threaded", 1, half), ("threaded", 1, N))]

    println("\nSpeedup from doubling cores N/2 -> N (ideal = 2.00):")
    @printf("%-36s %10s %10s %8s\n",
            "doubling", "before(s)", "after(s)", "speedup")
    speedups = Tuple{String, Float64}[]
    for (label, before, after) in doublings
        tb = timeof(results, before...)
        ta = timeof(results, after...)
        sp = tb / ta
        push!(speedups, (label, sp))
        @printf("%-36s %10.3f %10.3f %8.2f\n", label, tb, ta, sp)
    end

    # Machine-readable CSV of the raw timings (suffixed by count type so runs
    # with different `Count` types do not overwrite each other).
    csv_path = joinpath(OUTDIR, "benchmark_results_$COUNT_TYPE.csv")
    open(csv_path, "w") do io
        println(io, "mode,procs,threads,cores,seconds")
        for r in results
            @printf(io, "%s,%d,%d,%d,%.6f\n",
                    r.mode, r.procs, r.threads, r.cores, r.seconds)
        end
    end
    println("\nWrote ", csv_path)

    make_plots(results, speedups)
    return nothing
end

# A short label like "1p x 12t\n(threaded)" for a configuration.
conf_label(r) = @sprintf("%dp x %dt\n(%s)", r.procs, r.threads, r.mode)

function make_plots(results, speedups)
    # Per-configuration wall-clock time.
    kept = filter(r -> !isnan(r.seconds), results)
    if !isempty(kept)
        p1 = bar([conf_label(r) for r in kept], [r.seconds for r in kept],
                 legend = false, ylabel = "seconds",
                 title = "Simulation time by configuration ($COUNT_TYPE)",
                 rotation = 0)
        time_path = joinpath(OUTDIR, "config_times_$COUNT_TYPE.png")
        savefig(p1, time_path)
        println("Wrote ", time_path)
    end

    # Doubling speedups with the ideal (2x) reference line.
    kept_sp = filter(s -> !isnan(s[2]), speedups)
    if !isempty(kept_sp)
        p2 = bar([s[1] for s in kept_sp], [s[2] for s in kept_sp],
                 legend = false, ylabel = "speedup (before / after)",
                 title = "Speedup from doubling cores N/2 -> N ($COUNT_TYPE)",
                 rotation = 20)
        hline!(p2, [2.0], linestyle = :dash, color = :black)
        speedup_path = joinpath(OUTDIR, "doubling_speedup_$COUNT_TYPE.png")
        savefig(p2, speedup_path)
        println("Wrote ", speedup_path)
    end
    return nothing
end

main()
