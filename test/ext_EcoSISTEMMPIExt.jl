# SPDX-License-Identifier: LGPL-3.0-or-later

module ExtEcoSISTEMMPIExt

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using Diversity
using MPI
using JLD2
using Test

if !MPI.Initialized()
    MPI.Init()
end

@testset "MPI" begin
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println(Threads.nthreads())
    numSpecies = 100
    grid = (10, 10)
    demand = 10.0kJ / day
    individuals = 1_000
    area = 100.0 * km^2
    totalK = 100.0kJ / km^2 / day
    # Set up initial parameters for ecosystem

    # Set up how much resource each species consumes
    resource_vec = Demand{SolarRadiation}(fill(demand, numSpecies))

    # Set probabilities
    birth = 0.6 / year
    death = 0.6 / year
    longevity = 1.0
    survival = 0.2

    # Collect model parameters together
    param = EqualPop(birth, death, longevity, survival)

    # Create kernel for movement
    kernel = fill(GaussianKernel(1.0km, 10e-10), numSpecies)
    movement = BirthOnlyMovement(kernel)

    # Create species list, including their temperature preferences, seed abundance and native status
    opts = fill(274.0K, numSpecies)
    vars = fill(0.5K, numSpecies)
    tolerance = NicheTolerance(Temperature, Normal, opts, vars)
    native = fill(true, numSpecies)
    # abun = rand(Multinomial(individuals, numSpecies))
    abun = fill(div(individuals, numSpecies), numSpecies)
    sppl = SpeciesList(numSpecies, tolerance, abun, resource_vec,
                       movement, param, native)

    # Create abiotic environment - an even grid of one temperature. The study area decides the
    # grid; `GridHabitat` only samples onto it.
    studyarea = StudyArea(extent = (sqrt(area), sqrt(area)),
                          cellsize = sqrt(area) / grid[1], verbosity = :silent)
    habitat = GridHabitat(regime = UniformSpec(274.0K,
                                               axis = Temperature),
                          supply = UniformSpec(totalK,
                                               axis = SolarRadiation),
                          area = studyarea)

    # Set nichefit between species and environment (gaussian)
    nichefit = NicheSuitability{Temperature, typeof(1.0K)}()

    # build_ecosystem auto-selection: with MPI initialised but a single rank, `:auto` stays serial
    # (Comm_size == 1) while `distributed = true` forces the distributed type - from the same
    # SpeciesList + GridHabitat that MPIEcosystem takes directly.
    @test build_ecosystem(sppl, habitat, nichefit = nichefit) isa Ecosystem
    @test build_ecosystem(sppl, habitat, nichefit = nichefit,
                          distributed = true) isa
          MPIEcosystem

    # Create ecosystem
    @test_nowarn eco = MPIEcosystem(sppl, habitat, nichefit)
    eco = MPIEcosystem(sppl, habitat, nichefit)
    @test sum(eco.sppcounts) == length(eco.spplist.names)
    @test eco.firstsp == 1
    @test sum(eco.sccounts) == prod(size(eco.habitat.regime.matrix))
    @test eco.firstsc == 1

    # Simulation Parameters
    burnin = 1month_mean_duration
    times = 3month_mean_duration
    timestep = 1month_mean_duration
    record_interval = 1month_mean_duration
    repeats = 1
    lensim = length((0year):record_interval:times)
    # Burnin
    MPI.Barrier(comm)
    @time simulate!(eco, burnin, timestep)

    # Set columns vector to zero and check synchronise from rows
    eco.abundances.cols_vector .= 0
    EcoSISTEM.synchronise_from_rows!(eco.abundances)
    @test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

    # Set rows matrix to zero and check synchronise from cols
    eco.abundances.rows_matrix .= 0
    EcoSISTEM.synchronise_from_cols!(eco.abundances)
    @test sum(eco.abundances.cols_vector) == sum(eco.abundances.rows_matrix)

    # `getabundance` refuses on a distributed ecosystem whatever the rank count: it asks for the
    # whole metacommunity's abundances, which no rank holds, and Diversity's consumers each reduce
    # that matrix over a different axis. The refusal stands even here, where a single rank does
    # happen to hold everything - an answer that appeared at one rank and vanished at two would be
    # exactly the kind of rank-dependence the distributed code exists to avoid.
    @test_throws "no single rank holds them" getabundance(eco)

    # The genuinely global quantities, asserted on their length as well as their total: a sum of one
    # is true of any block normalised by its own total, which is how the per-species totals came to
    # cover only the calling rank's species without any test noticing.
    @test length(getmetaabundance(eco)) == numSpecies
    @test sum(getmetaabundance(eco)) ≈ 1.0
    @test sum(getordinariness!(eco)) ≈ 1.0
    @test sum(getweight(eco)) ≈ 1.0
    @test length(getweight(eco)) == countsubcommunities(eco)

    # Gather abundances and check against rows matrix - should be same for 1 process
    abuns = gatherabundance(eco)
    @test abuns == eco.abundances.rows_matrix
end

# The output of `cmd`, or a line saying why there is none. It only describes a launch that has
# already failed, so it must neither hang nor throw: it waits at most `limit` seconds.
function capture(cmd::Cmd; limit = 30)
    path, io = mktemp()
    try
        process = run(pipeline(cmd, stdout = io, stderr = io), wait = false)
        if timedwait(() -> process_exited(process), limit, pollint = 0.5) !==
           :ok
            kill(process)
            return "($(first(cmd.exec)) gave no answer within $limit s)"
        end
    catch err
        return "($(first(cmd.exec)) could not run: $(sprint(showerror, err)))"
    finally
        close(io)
    end
    return read(path, String)
end

# Every process descended from `root`, by pid, with its command line. `mpiexec` starts the ranks
# through a proxy, so they are its grandchildren rather than its children.
function descendants(root::Integer)
    table = Dict{Int, Tuple{Int, String}}()
    for line in eachline(IOBuffer(capture(`ps -ww -A -o pid= -o ppid= -o command=`)))
        fields = split(strip(line), limit = 3)
        length(fields) == 3 || continue
        pid, ppid = tryparse(Int, fields[1]), tryparse(Int, fields[2])
        (isnothing(pid) || isnothing(ppid)) && continue
        table[pid] = (ppid, fields[3])
    end
    found = Dict{Int, String}()
    frontier = [Int(root)]
    while !isempty(frontier)
        parent = pop!(frontier)
        for (pid, (ppid, command)) in table
            if ppid == parent && !haskey(found, pid)
                found[pid] = command
                push!(frontier, pid)
            end
        end
    end
    return found
end

# The parts of a macOS `sample` report that say where a process was: the first thread's call graph,
# which is the main thread's, and the frames it was most often found in.
function mainstack(report::AbstractString)
    lines = split(report, '\n')
    graph = findfirst(==("Call graph:"), lines)
    isnothing(graph) && return report
    nextthread = findnext(l -> occursin(r"^\s+\d+ Thread_", l), lines,
                          graph + 2)
    keep = lines[graph:min(something(nextthread, graph + 60) - 1,
                           lastindex(lines))]
    top = findfirst(l -> startswith(l, "Sort by top of stack"), lines)
    if !isnothing(top)
        images = findnext(l -> startswith(l, "Binary Images"), lines, top)
        append!(keep,
                lines[top:min(something(images, top + 20) - 1,
                              lastindex(lines))])
    end
    return join(keep, '\n')
end

# A Linux process's kernel wait channel, and the lines of its status saying whether it is running,
# sleeping or swapped out.
function procstate(pid::Integer)
    try
        wchan = read("/proc/$pid/wchan", String)
        status = filter(readlines("/proc/$pid/status")) do line
            return occursin(r"^(State|VmRSS|VmSwap|voluntary_ctxt|nonvoluntary_ctxt)",
                            line)
        end
        return join(["wchan: $wchan"; status], '\n')
    catch err
        return "(could not read /proc/$pid: $(sprint(showerror, err)))"
    end
end

# What a launch was doing when it reached its limit: each process under `mpiexec`, listed twice five
# seconds apart to tell a rank still using CPU from a blocked one; the machine's memory and swap; and
# where each Julia rank was, as a native stack sample on macOS or the kernel's view on Linux, where
# attaching a debugger needs privileges a CI job lacks. Lines are cut to 200 characters, since a
# rank's command runs to several hundred. Nothing on Windows, which has no `ps`.
function launchstate(root::Integer)
    Sys.iswindows() && return "(not described on Windows)"
    procs = descendants(root)
    isempty(procs) && return "(nothing is running under mpiexec)"
    report = IOBuffer()
    listing = `ps -ww -o pid,ppid,stat,%cpu,time,rss,etime,command -p $(join(keys(procs), ","))`
    println(report, capture(listing))
    sleep(5)
    println(report, capture(listing))
    if Sys.isapple()
        println(report, capture(`sysctl vm.swapusage`), capture(`vm_stat`))
    else
        println(report, capture(`free -m`))
    end
    for (pid, command) in procs
        basename(first(split(command))) == "julia" || continue
        location = Sys.isapple() ? mainstack(capture(`sample $pid 1`)) :
                   procstate(pid)
        println(report, "process $pid:\n", location)
    end
    lines = split(String(take!(report)), '\n')
    return join((first(line, 200) for line in lines), '\n')
end

# Run one `mpiexec` launch, and fail rather than hang if its ranks never finish: a rank left waiting
# in a collective waits forever. A launch takes a few minutes, far inside the limit. Killing
# `mpiexec` stops its ranks too. The ranks' output goes to `logfile`, and the end of it is reported
# only if the launch fails, after what its processes were doing if it was stopped.
function runmpi(cmd, logfile; limit = 1800)
    state = nothing
    finished = open(logfile, "w") do io
        process = run(pipeline(cmd, stdout = io, stderr = io), wait = false)
        waited = timedwait(() -> process_exited(process), limit, pollint = 1.0)
        waited === :ok && return success(process)
        state = launchstate(getpid(process))
        kill(process)
        wait(process)
        return false
    end
    finished && return true
    # In the message, which is printed whole: a logged value is printed with `show`, which escapes
    # the newlines and elides all but the ends of a long string.
    report = isnothing(state) ? "`mpiexec` failed: $cmd" :
             "`mpiexec` was stopped after $limit s without finishing: $cmd\n\n" *
             "What its processes were doing:\n$state"
    output = join(last(readlines(logfile), 200), '\n')
    @error "$report\n\nThe last of its output:\n$output"
    return false
end

# Hydra's flag labelling each line of output with the rank that wrote it, on MPICH, whose launcher
# has it; Microsoft MPI's has no equivalent.
rankprefix() = MPI.MPI_LIBRARY == "MPICH" ? `-prepend-rank` : ``

# The Julia command a launch runs, with coverage instrumentation only where it is wanted. It roughly
# doubles a launch - measured on a runner, this file took 861 s on the instrumented job against 402 s
# uninstrumented - and each launch runs the same package code, so one carrying it covers the same
# lines as three would. The parent process, which runs this file, keeps whatever it was given.
function launchcmd(; coverage::Bool)
    cmd = Base.julia_cmd()
    coverage && return cmd
    return Cmd(filter(a -> !startswith(a, "--code-coverage"), cmd.exec))
end

@testset "mpirun" begin
    # Keep the MPI outputs in a temp dir the OS cleans up (no manual `rm` needed for hygiene).
    # The child `mpiexec` processes read its path as their first command-line argument.
    datadir = mktempdir()
    # Compare 1 thread 4 processes vs. 4 threads 1 process vs. 2 threads 2 processes
    withenv("JULIA_NUM_THREADS" => "4") do
        nprocs = 1
        function cmd(n = nprocs)
            return `$(mpiexec()) $(rankprefix()) -n $nprocs $(launchcmd(coverage = false)) --startup-file=no $(pkgdir(EcoSISTEM, "test", "SmallMPItest.jl")) $datadir`
        end
        @test runmpi(cmd(), joinpath(datadir, "mpiexec1.log"))
    end
    withenv("JULIA_NUM_THREADS" => "2") do
        nprocs = 2
        function cmd(n = nprocs)
            return `$(mpiexec()) $(rankprefix()) -n $nprocs $(launchcmd(coverage = true)) --startup-file=no $(pkgdir(EcoSISTEM, "test", "SmallMPItest.jl")) $datadir`
        end
        @test runmpi(cmd(), joinpath(datadir, "mpiexec2.log"))
    end
    withenv("JULIA_NUM_THREADS" => "1") do
        nprocs = 4
        function cmd(n = nprocs)
            return `$(mpiexec()) $(rankprefix()) -n $nprocs $(launchcmd(coverage = false)) --startup-file=no $(pkgdir(EcoSISTEM, "test", "SmallMPItest.jl")) $datadir`
        end
        @test runmpi(cmd(), joinpath(datadir, "mpiexec4.log"))
    end

    ## All answers should be the same across process/thread splits. Load the
    ## saved values (note: `@load file var` binds `var` but returns the symbol
    ## list, so use `load(file, "abuns")` to get the data itself).
    abuns1thread = load(joinpath(datadir, "Test_abuns1.jld2"), "abuns")
    abuns2thread = load(joinpath(datadir, "Test_abuns2.jld2"), "abuns")
    abuns4thread = load(joinpath(datadir, "Test_abuns4.jld2"), "abuns")

    # Each launch also checks its own run against the blessed values, so a launch that does not
    # reproduce itself fails there, apart from any disagreement between the three here.
    @test abuns1thread == abuns2thread == abuns4thread
end

if !MPI.Finalized()
    MPI.Finalize()
end

end
