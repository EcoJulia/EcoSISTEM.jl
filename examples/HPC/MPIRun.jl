# SPDX-License-Identifier: LGPL-3.0-or-later
#
# **The MPI scaling benchmark**: 65,536 species on a 256 × 256 grid, timed over two identical
# burn-ins so the first pays for compilation and the second measures the run. Each rank writes its
# own timing file, so a set of jobs at different process/thread splits can be compared afterwards —
# which is what `demo-MPI-threads.bash`, `demo-MPI-processes.bash` and `demo-MPI-nodes.bash` launch.
#
# **It goes through `StudyArea` → `GridHabitat` → `build_species` → `build_ecosystem`.** Building an
# environment with a deprecated builder, a hand-assembled `SpeciesList` and a direct `MPIEcosystem`
# call decides no grid, and such a run cannot start.
#
# **The grid is fixed, and checked rather than chosen.** This is a benchmark: the point is to time
# *this* configuration, so a run that will not fit is an error, not an invitation to shrink. That is
# `MemoryGuidance.check_memory`, and it is asked **before** anything is allocated — on a batch system
# an OOM kill an hour in reports only that the job died.
#
# At the full size this needs about **128 GiB across the job**, so it wants either a large node or
# several. Smoke-test the machinery in seconds with `ECOSISTEM_SCALE=small`.
#
#     ECOSISTEM_SCALE=small mpiexecjl --project=examples -n 2 julia --project=examples examples/HPC/MPIRun.jl
#     mpiexecjl --project=examples -n 32 julia -t 8 --project=examples examples/HPC/MPIRun.jl
#     sbatch examples/HPC/demo-MPI-nodes.bash
#
# | variable | default | |
# |---|---|---|
# | `ECOSISTEM_SCALE` | `large` | `small` for a fast smoke test |
# | `ECOSISTEM_MPIRUN_SAVEDIR` | `./MPIRun-<ranks>x<threads>` | where the timing files go |

start = time()

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using MPI
using Statistics

include(joinpath(@__DIR__, "memory.jl"))

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const RANKS = MPI.Comm_size(COMM)
const ISROOT = RANK == 0

@info "Total memory: $(Sys.total_memory() / 2^30) GiB, threads: $(Threads.nthreads())"

# --- configuration ----------------------------------------------------------------------

const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"

# Species and cells are both powers of two at full size (2^16 each), which divides evenly across
# any power-of-two rank count — deliberate for a *benchmark*, where an uneven split would add a load
# imbalance to whatever is being measured. Do not copy that into a test: `test/SmallMPItest.jl`
# uses 7 species on 77 cells precisely because even division hides partitioning bugs.
const NUMSPECIES = SMALL ? 64 : 2^16
const GRID = SMALL ? 16 : 256
const INDIVIDUALS = SMALL ? 2^16 : 2^26

# The grid the original expressed as `(256, 256)` cells over `1_000_000 km²` — the same thing said
# the way a `StudyArea` needs it, as an extent and a cell size.
const SIDE = SMALL ? 100.0km : 1000.0km
const CELLSIZE = SIDE / GRID

const BURNIN = SMALL ? 1month_mean_duration : 1year
const TIMESTEP = 1month_mean_duration

const SAVEDIR = get(ENV, "ECOSISTEM_MPIRUN_SAVEDIR",
                    joinpath(pwd(),
                             "MPIRun-$(RANKS)x$(Threads.nthreads())"))

ISROOT && println("using: $((time() - start) * s)")

# --- will it fit? -----------------------------------------------------------------------

# Collective, and every rank must reach it — see `memory.jl`.
const BUDGET = MemoryGuidance.memory_budget()

# Costed from the *report*, so nothing is built to find out whether it can be: the analysis behind
# `investigate_study_area` is the one `StudyArea` would run, and `getspeciesstorage` prices its
# grid. A synthetic area needs no data at all, so this is free.
const REPORT = investigate_study_area(extent = (SIDE, SIDE),
                                      cellsize = CELLSIZE)

if ISROOT
    MemoryGuidance.describe(BUDGET)
end

MemoryGuidance.check_memory(NUMSPECIES,
                            EcoSISTEM.getspeciesstorage(REPORT),
                            budget = BUDGET,
                            what = "MPIRun ($(NUMSPECIES) species on $(GRID) × $(GRID) cells)")

# --- build ------------------------------------------------------------------------------

const AREA = StudyArea(extent = (SIDE, SIDE), cellsize = CELLSIZE,
                       verbosity = ISROOT ? :normal : :silent)

# A flat environment: one temperature everywhere, one solar supply everywhere. That is the point
# of a scaling benchmark — nothing about the landscape should vary the work per cell.
const ENVIRONMENT = GridHabitat(regime = UniformSpec(274.0K,
                                                     axis = Temperature),
                                supply = UniformSpec(1000.0kJ / km^2 /
                                                     day,
                                                     axis = SolarRadiation),
                                area = AREA)

# `BirthOnlyMovement` is not a stylistic choice here: `AlwaysMovement` is unimplemented under MPI
# and refuses (see `ext/EcoSISTEMMPIExt/dynamics.jl`).
const SPECIES = build_species(NUMSPECIES,
                              tolerance = (274.0K, 0.5K),
                              toleranceaxis = Temperature,
                              demand = 1.0kJ / day,
                              demandaxis = SolarRadiation,
                              dispersal = 1.0km,
                              movement = BirthOnlyMovement,
                              survival = 0.2,
                              abundance = INDIVIDUALS,
                              seed = 1)

ISROOT && println("Startup: $((time() - start) * s)")

const ECO = build_ecosystem(SPECIES, ENVIRONMENT, seed = 1,
                            distributed = :auto)

# --- timing -----------------------------------------------------------------------------

# One line per rank per event, appended. Each rank owns its own file, so no rank waits on another to
# write and the files can be collated afterwards without interleaving.
function _record(message)
    mkpath(SAVEDIR)
    cores = RANKS * Threads.nthreads()
    open(joinpath(SAVEDIR, "output-cores$(cores)-np$(RANKS)-$(RANK).txt"),
         append = true) do io
        return write(io, "$message @ $((time() - start) * s)\n")
    end
    return nothing
end

# How many species this rank holds and how their abundances are spread — the load-balance check the
# original made, and worth keeping: an uneven partition shows up here before it shows up in a timing.
#
# `rows_matrix` exists only on the distributed landscape, so this must ask which it has: a
# single-rank job (`:auto` makes that serial) otherwise dies here rather than in anything it measures.
function _describelocal()
    rows = ECO isa EcoSISTEM.MPIEcosystem ? ECO.abundances.rows_matrix :
           ECO.abundances.matrix
    counts = sum(rows, dims = 2)
    return "numspp = $(size(rows, 1)) mean = $(mean(counts)) std = $(std(counts))"
end

_record("rank $RANK / $RANKS: $(Threads.nthreads()) threads")
_record(_describelocal())

# **Two identical burn-ins, and only the second is the measurement.** The first compiles the whole
# hot loop for these concrete types; timing it would measure the compiler.
for pass in (:warmup, :measured)
    ISROOT && println("Start $pass burnin: $((time() - start) * s)")
    MPI.Barrier(COMM)
    elapsed = @elapsed simulate!(ECO, BURNIN, TIMESTEP)
    _record("$pass burnin: time = $(elapsed * s)")
    _record(_describelocal())
end

ISROOT && println("End: $((time() - start) * s)")
ISROOT && println("Timings written to $(SAVEDIR)")

MPI.Barrier(COMM)
