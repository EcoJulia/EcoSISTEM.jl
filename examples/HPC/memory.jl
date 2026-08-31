# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Memory guidance for the HPC scripts: how much memory this run can plan around, what a given
# configuration would cost, and therefore the finest grid it can afford.
#
# **One helper for both cases.** A laptop run is not started by `mpirun`, so `MPI.Initialized()`
# is false (or there is a single rank) and the budget is simply this machine's. Under `mpirun`/`srun`
# with several ranks it becomes the sum across the nodes the job actually landed on. Nothing above
# has to know which it got - the same call answers both, and reports which it was.
#
# **The assumption, stated because it decides the arithmetic**: ranks do *not* share the landscape
# (each holds its own private slice), and each node is exclusively this job's. The per-rank figure is
# therefore its node's memory divided by the ranks on that node, and the global figure takes the
# **minimum** of those - the weakest rank sets the chunk size for a uniformly distributed array, so
# a heterogeneous allocation is costed at its worst node rather than its average.
#
# Include it and ask:
#
#     include(joinpath(@__DIR__, "memory.jl"))
#     budget = MemoryGuidance.memory_budget()
#     chosen = MemoryGuidance.choose_cellsize(50_000, (5.0km, 10.0km, 20.0km, 50.0km)) do cellsize
#         prod(size(StudyArea(; regime, supply, within, crs, cellsize).report.active))
#     end
#
# Run the demonstration below directly:
#
#     julia --project=examples examples/HPC/memory.jl
#     mpiexecjl --project=examples -n 128 julia examples/HPC/memory.jl
#     srun -N 4 --ntasks-per-node=32 --exclusive julia --project=examples examples/HPC/memory.jl

module MemoryGuidance

using MPI
using Printf

# --- what a run costs -------------------------------------------------------------------

# The number of full `species × cells` `Int64` arrays a run keeps resident, which is what the
# estimate is built on. Serial holds the abundance matrix and `cache.netmigration`
# (`src/Ecosystem.jl:307`). MPI holds **three**, not two: `rows_matrix` (species-partitioned),
# `cols_vector` (location-partitioned) and its own `netmigration`, sized like `rows_matrix`
# (`MPIGridLandscape`'s fields in `ext/EcoSISTEMMPIExt/Landscape.jl`, and the `nm` allocation in its
# `MPIEcosystem` constructor). Each is partitioned across ranks,
# so each totals `species × cells` over the whole job.
const SERIAL_ARRAYS = 2
const MPI_ARRAYS = 3

# Fraction of nominal RAM that can actually be allocated at all. The rest goes to the kernel, page
# tables, MPI's internal and communication buffers, `/dev/shm`, and slack for Julia's GC not
# returning freed pages promptly. Lower it if the run does a lot of large transient allocation.
const HEADROOM = 0.85

# One array's worth of slack on top of what the run actually holds, for transients, for the
# dispersal lookups and per-species vectors, and for Julia's GC not returning freed pages promptly.
#
# **This is why there is no separate "use less than a fraction of memory" constant.** Requiring
# room for `arrays + 1` *is* that rule: at two arrays the real ones may occupy no more than a third
# of the budget, and at three no more than a quarter. Writing the fraction out as well would restate
# the array count and then have to be kept in step with it.
#
# In particular MPI needs no extra allowance for communication. `Alltoallv!` is handed
# `rows_matrix` and `cols_vector` directly as its send and receive buffers
# (`synchronise_from_rows!`/`synchronise_from_cols!` in `ext/EcoSISTEMMPIExt/Landscape.jl`), and
# `MPI.VBuffer` wraps a pointer rather than copying,
# so nothing is allocated at the Julia level; the library's own per-message chunking is bounded and
# belongs to `HEADROOM` above.
const SPARE_ARRAYS = 1

# --- what we have -----------------------------------------------------------------------

# The memory this process may use on its own node, in bytes. Prefers the batch scheduler's
# allocation to the machine's physical RAM: a job on a 1 TB node may be held to a fraction of it,
# and sizing from `Sys.total_memory()` there gets the run OOM-killed rather than slowed down.
# `Sys.total_memory()` rather than `Sys.total_physical_memory()` for the fallback - it is the one
# documented to respect a Linux control group, which is how the limit usually arrives.
function _nodememory()
    permb = get(ENV, "SLURM_MEM_PER_NODE", "")
    isempty(permb) || return parse(Int64, permb) * 1024^2
    percpu = get(ENV, "SLURM_MEM_PER_CPU", "")
    cpus = get(ENV, "SLURM_CPUS_ON_NODE", "")
    (isempty(percpu) || isempty(cpus)) && return Int64(Sys.total_memory())
    return parse(Int64, percpu) * parse(Int64, cpus) * 1024^2
end

# Whether this process is part of a live multi-rank MPI session - the same test the package's own
# `_should_mpi` makes, so `build_ecosystem(distributed = :auto)` and this agree on which run it is.
_isdistributed() = MPI.Initialized() && MPI.Comm_size(MPI.COMM_WORLD) > 1

# Measured once and kept. Under MPI `memory_budget` is **collective** - it reduces across
# `COMM_WORLD` - so calling it again from a subset of ranks (say, inside a `rank == 0` branch) would
# hang every rank with no error at all. Caching makes every call after the first free and safe to
# make from anywhere, leaving only the first call needing to be reached by all ranks.
const _BUDGET = Ref{Union{Nothing, NamedTuple}}(nothing)

"""
    memory_budget()

Return what this run can plan around, as a named tuple of `total` and `per_rank` bytes, the number
of `nodes` and `ranks` it spans, and whether it is `distributed`.

Serial runs report this machine alone. Under MPI, ranks are grouped by the memory they share - one
group per node - and each rank's budget is its node's memory divided by the ranks on it; `per_rank`
is the **minimum** across the job and `total` is that minimum times the rank count, so an uneven
allocation is costed at its weakest node.

 **The first call is collective and every rank must reach it** - call it once at startup, before
any `rank == 0` branch. The result is cached, so later calls cost nothing and may be made from a
single rank.
"""
function memory_budget()
    isnothing(_BUDGET[]) || return _BUDGET[]

    if !_isdistributed()
        node = _nodememory()
        return _BUDGET[] = (total = node, per_rank = node, nodes = 1, ranks = 1,
                            distributed = false)
    end

    comm = MPI.COMM_WORLD
    ranks = MPI.Comm_size(comm)

    # Ranks that can share memory are ranks on the same node, for every mainstream MPI
    # implementation - and unlike counting distinct `MPI.Get_processor_name()` strings, this stays
    # right where one node runs several containers.
    nodecomm = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, 0)
    hereranks = MPI.Comm_size(nodecomm)
    isleader = MPI.Comm_rank(nodecomm) == 0
    MPI.free(nodecomm)

    perrank = _nodememory() ÷ Int64(hereranks)
    minperrank = MPI.Allreduce(perrank, min, comm)
    nodes = MPI.Allreduce(Int64(isleader), +, comm)

    return (total = minperrank * Int64(ranks), per_rank = minperrank,
            nodes = nodes, ranks = ranks, distributed = true)
end

"""
    usable_memory([budget])

Return the bytes this run may allocate, after allowing for what the operating system, the MPI
library and Julia's GC take off the top. Takes a [`memory_budget`](@ref) result, measuring one if
not given.
"""
function usable_memory(budget = memory_budget())
    return floor(Int64, budget.total * HEADROOM)
end

"""
    required_memory(numspecies, perspecies; distributed)

Return the bytes `numspecies` species would need: every full `species × cells` array the run holds -
two serially, three under MPI - plus one more as slack.

`perspecies` is the bytes **one** species' abundances occupy over the whole grid.  Ask the package
rather than working it out: `investigate_study_area(...).footprint.perspecies` is exactly this
figure, decided by the same analysis that `StudyArea` would use, so it can never describe a grid
other than the one that would be built.

 It covers the **whole** grid, not the active cells: the abundance matrix is allocated over all of
it, and an inactive cell costs the same as an active one.
"""
function required_memory(numspecies::Integer, perspecies::Integer;
                         distributed::Bool = _isdistributed())
    arrays = (distributed ? MPI_ARRAYS : SERIAL_ARRAYS) + SPARE_ARRAYS
    return Int64(arrays) * Int64(numspecies) * Int64(perspecies)
end

"""
    check_memory(numspecies, perspecies; budget = memory_budget(), what = "this run")

Return the bytes required, or throw if they will not fit in `budget`.

For a run whose grid is **fixed** - a benchmark that must be the size it says, not whatever happens
to fit. Where the size is negotiable, [`choose_cellsize`](@ref) picks one instead.

- `numspecies`, `perspecies` - as for [`required_memory`](@ref); ask the package for the second with
  `investigate_study_area(...)` and `EcoSISTEM.getspeciesstorage`.
- `budget` - a [`memory_budget`](@ref) result, measured if not given.
- `what` - how to name the run in the error, so a script with several says which one failed.

 Checked **before** anything is allocated, so an impossible configuration is a message rather than
an OOM kill an hour into a batch job - which is the whole point on a shared machine, where the
scheduler's report says only that the job died.
"""
function check_memory(numspecies::Integer, perspecies::Integer;
                      budget = memory_budget(), what = "this run")
    required = required_memory(numspecies, perspecies,
                               distributed = budget.distributed)
    usable = usable_memory(budget)
    required ≤ usable && return required
    return error("$what needs $(_gib(required)) but only $(_gib(usable)) is usable " *
                 "($(budget.distributed ? "$(budget.ranks) MPI ranks over $(budget.nodes) node(s)" : "serial"), " *
                 "$(_gib(budget.total)) total). Reduce the species count or the grid, or ask for " *
                 "more nodes.")
end

"""
    choose_cellsize(footprint, numspecies, candidates; budget = memory_budget())

Return the finest of `candidates` that `numspecies` species will fit in, as a named tuple of the
chosen `cellsize`, its `perspecies` bytes, the `required` and `usable` bytes, and whether it `fits`.

`footprint` is called with a cell size and must return the bytes one species' abundances would
occupy on that grid - `investigate_study_area(...).footprint.perspecies`. It is taken first so it can
be written as a do-block. Candidates are tried finest first and the first that fits is returned, so
it is called once per candidate at most. If none fits, the coarsest is returned with `fits = false`
rather than an error: the caller decides whether to shrink the species pool, ask for more nodes, or
run it anyway.
"""
function choose_cellsize(footprint, numspecies::Integer, candidates;
                         budget = memory_budget())
    usable = usable_memory(budget)
    for cellsize in sort(collect(candidates))
        perspecies = footprint(cellsize)
        required = required_memory(numspecies, perspecies,
                                   distributed = budget.distributed)
        required ≤ usable &&
            return (cellsize = cellsize, perspecies = perspecies,
                    required = required, usable = usable, fits = true)
    end

    coarsest = maximum(candidates)
    perspecies = footprint(coarsest)
    return (cellsize = coarsest, perspecies = perspecies,
            required = required_memory(numspecies, perspecies,
                                       distributed = budget.distributed),
            usable = usable, fits = false)
end

# --- saying so --------------------------------------------------------------------------

# Bytes as GiB to one decimal place - every figure here is a planning estimate, so more precision
# than that would be claiming an accuracy the headroom fraction does not have.
_gib(bytes) = @sprintf("%.1f GiB", bytes / 2^30)

"""
    describe([budget]; io = stdout)

Print what was detected: serial or distributed, over how many ranks and nodes, the memory it may
allocate, and how many `species × cells` arrays a run there has to fit into it.

 Every rank returns the same figures, so under MPI this should be called from one rank only -
otherwise the report appears once per rank.
"""
function describe(budget = memory_budget(); io = stdout)
    place = budget.distributed ?
            "MPI, $(budget.ranks) ranks over $(budget.nodes) node(s)" : "serial"
    arrays = budget.distributed ? MPI_ARRAYS : SERIAL_ARRAYS
    println(io, "Memory budget ($place):")
    println(io, "  per rank      ", _gib(budget.per_rank))
    println(io, "  total         ", _gib(budget.total))
    println(io, "  allocatable   ", _gib(usable_memory(budget)),
            "  (fits $arrays arrays + $SPARE_ARRAYS spare)")
    return nothing
end

"""
    describe(chosen, numspecies; io = stdout)

Print the outcome of a [`choose_cellsize`](@ref), naming what it picked and what that costs - or, if
nothing fitted, saying so rather than letting the run discover it by being killed.
"""
function describe(chosen::NamedTuple, numspecies::Integer; io = stdout)
    if chosen.fits
        println(io, "Chose $(chosen.cellsize) cells: ",
                Base.format_bytes(chosen.perspecies), " per species × ",
                "$(numspecies) species needs ", _gib(chosen.required),
                " of ", _gib(chosen.usable), ".")
    else
        println(io,
                "⚠️ Nothing fits. The coarsest candidate ($(chosen.cellsize)) ",
                "still needs ", _gib(chosen.required), " of ",
                _gib(chosen.usable),
                " - reduce the species pool or ask for more nodes.")
    end
    return nothing
end

end # module MemoryGuidance

# --- run directly to see what this machine or allocation offers -------------------------
# Only when executed as a script, so an `include` from another example stays silent.
if abspath(PROGRAM_FILE) == @__FILE__
    using MPI
    using Unitful, Unitful.DefaultSymbols

    MPI.Initialized() || MPI.Init()

    # Every rank, and before the `rank == 0` branch below: this is the collective call.
    budget = MemoryGuidance.memory_budget()

    if !budget.distributed || MPI.Comm_rank(MPI.COMM_WORLD) == 0
        MemoryGuidance.describe(budget)
        println()

        # **A canned illustration, not the real mechanism.** These are the Africa grids as
        # measured, hardcoded so this file stays standalone - it depends on nothing but MPI and can
        # be run anywhere. A real script asks the package instead, and
        # `examples/HPC/Africa.jl` is the worked example: `investigate_study_area` resolves the grid
        # and `EcoSISTEM.getspeciesstorage` prices it, so nothing is written down.
        const AFRICA_CELLS = Dict(5.0km => 2_138_808, 10.0km => 535_108,
                                  20.0km => 133_650, 50.0km => 21_252,
                                  100.0km => 5_280)

        # `budget` passed rather than defaulted - inside a single-rank branch that is the habit
        # worth having, even though the cache above now makes the default safe too.
        numspecies = 50_000
        chosen = MemoryGuidance.choose_cellsize(numspecies, keys(AFRICA_CELLS),
                                                budget = budget) do cellsize
            return AFRICA_CELLS[cellsize] * sizeof(Int64)
        end
        MemoryGuidance.describe(chosen, numspecies)
    end
end
