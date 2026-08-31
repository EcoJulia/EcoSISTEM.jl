# HPC scripts

Everything here runs from the root of a checked-out [EcoSISTEM][ecosistem-git], with the project
environment in `examples/`:

```sh
julia --project=examples examples/HPC/<script>.jl
```

> Install the packages on a login node (or anywhere with internet access) first - a compute node
> will otherwise try to do it at startup and may not be able to.

**None of these is run by the test suite.** `test/extras_examples.jl` executes only the top level
of `examples/`, so these are exercised by hand. Each takes `ECOSISTEM_SCALE=small` (except
`SmallMPIRun.jl`, which is already small) to give a fast smoke test of the machinery.

## What is here

| | |
|---|---|
| `memory.jl` | works out what a run can allocate, and what a configuration would cost. Included by the others; run it directly to see what this machine or allocation offers |
| `Africa.jl` | 50,000 plant species over Africa on real WorldClim temperature. **Picks its grid** from the memory available, serial or MPI |
| `MPIRun.jl` | the scaling benchmark - 65,536 species on 256 × 256, timed over two burn-ins, one timing file per rank |
| `SmallMPIRun.jl` | the small smoke test - 100 species on 50 × 50. Try this first on a new machine or MPI build |

**Two ways of using `memory.jl`, and the difference is whether the size is negotiable.**
`Africa.jl` has a *choice* of resolutions and takes `choose_cellsize` to find the finest that fits.
The two MPI scripts are benchmarks whose size is the point, so they take `check_memory` and **error**
rather than shrink - before anything is allocated, because on a batch system an out-of-memory kill an
hour in reports only that the job died.

Neither writes down how much anything costs: `investigate_study_area` resolves a grid without
building it, and `EcoSISTEM.getspeciesstorage` says what one species' abundances occupy on it.

## Memory

`memory.jl` sums the physical memory of the distinct nodes an MPI job actually landed on
(`MPI.Comm_split_type` with `COMM_TYPE_SHARED`), preferring Slurm's allocation to the machine's RAM
where `SLURM_MEM_PER_NODE` says so - a job cgroup-limited to a fraction of a large node would
otherwise size itself for the whole thing and be killed.

A run holds several full `species × cells` `Int64` arrays: **two** serially (abundances and net
migration) and **three** distributed (the row- and column-partitioned landscapes plus net migration),
and one more is allowed as slack.

## Multithreaded runs

Any run is multithreaded and uses every available thread. On a node with two 64-core processors:

```sh
sbatch examples/HPC/demo-threads.bash
```

## MPI runs

Julia's own MPI libraries are used here. On a real cluster its native MPI is usually much faster,
being built for that machine's interconnect - see [MPIPreferences][mpiprefs] for pointing Julia at
it.

```sh
sbatch examples/HPC/demo-MPI-threads.bash      # 2 tasks x 32 threads
sbatch examples/HPC/demo-MPI-processes.bash    # 64 tasks x 1 thread
sbatch examples/HPC/demo-MPI-nodes.bash        # 4 nodes, 32 tasks x 8 threads
```

The `#SBATCH` lines in those files name an account, partition and memory that are specific to one
cluster - edit them for yours. `MPIRun.jl` writes its timing files to `./MPIRun-<ranks>x<threads>`,
overridable with `ECOSISTEM_MPIRUN_SAVEDIR`.

[ecosistem-git]: https://github.com/EcoJulia/EcoSISTEM.jl.git
[mpiprefs]: https://juliaparallel.org/MPI.jl/stable/configuration/
