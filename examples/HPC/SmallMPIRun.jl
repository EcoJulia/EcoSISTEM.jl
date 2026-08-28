# SPDX-License-Identifier: LGPL-3.0-or-later
#
# **The small MPI smoke test**: 100 species on a 50 × 50 grid, one burn-in, timed. Small enough to
# run anywhere, so it is what to try first on a new machine or a new MPI build — if this does not
# work, `MPIRun.jl` will not either, and it fails in seconds rather than after a queue wait.
#
# Like its larger sibling, it goes through `StudyArea` → `GridHabitat` → `build_species` →
# `build_ecosystem`: a deprecated builder with a hand-assembled `SpeciesList` and a direct
# `MPIEcosystem` call decides no grid, and such a run cannot start.
#
# **No `Pkg.activate("examples")` here, deliberately.** It assumes the process was started from the
# repository root and silently activates the wrong environment from anywhere else. Pass
# `--project=examples` instead, which is what every other script here does and what the `sbatch`
# demos pass.
#
# **The grid is fixed and checked, not chosen.** `MemoryGuidance.check_memory` refuses a
# configuration that cannot fit *before* anything is allocated. At this size it never will refuse —
# which is the point of having it here: the check costs nothing and the habit is what matters when
# the same lines appear in `MPIRun.jl`, where it does bite.
#
#     mpiexecjl --project=examples -n 4 julia --project=examples examples/HPC/SmallMPIRun.jl
#     srun --ntasks=4 julia --project=examples examples/HPC/SmallMPIRun.jl

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using MPI

include(joinpath(@__DIR__, "memory.jl"))

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const RANKS = MPI.Comm_size(COMM)
const ISROOT = RANK == 0

ISROOT && println("$(RANKS) rank(s), $(Threads.nthreads()) thread(s) each.")

# --- configuration ----------------------------------------------------------------------

const NUMSPECIES = 100
const GRID = 50
const INDIVIDUALS = 1_000_000

# `(50, 50)` cells over `100_000 km²`, said the way a `StudyArea` needs it. `sqrt` because the
# original gave a total area and a cell count and left the cell size implicit.
const SIDE = sqrt(100_000.0km^2)
const CELLSIZE = SIDE / GRID

const BURNIN = 20year
const TIMESTEP = 1month_mean_duration

# --- will it fit? -----------------------------------------------------------------------

# Collective, and every rank must reach it — see `memory.jl`.
const BUDGET = MemoryGuidance.memory_budget()

# Priced from the report rather than the built area, so an impossible configuration costs nothing
# to discover. A synthetic grid needs no data, so this reads nothing.
const REPORT = investigate_study_area(extent = (SIDE, SIDE),
                                      cellsize = CELLSIZE)

ISROOT && MemoryGuidance.describe(BUDGET)

MemoryGuidance.check_memory(NUMSPECIES,
                            EcoSISTEM.getspeciesstorage(REPORT),
                            budget = BUDGET,
                            what = "SmallMPIRun ($(NUMSPECIES) species on $(GRID) × $(GRID) cells)")

# --- build ------------------------------------------------------------------------------

const AREA = StudyArea(extent = (SIDE, SIDE), cellsize = CELLSIZE,
                       verbosity = ISROOT ? :normal : :silent)

const ENVIRONMENT = GridHabitat(regime = UniformSpec(274.0K,
                                                     axis = Temperature),
                                supply = UniformSpec(100.0kJ / km^2 /
                                                     day,
                                                     axis = SolarRadiation),
                                area = AREA)

# `BirthOnlyMovement` because `AlwaysMovement` is unimplemented under MPI and refuses — see
# `ext/EcoSISTEMMPIExt/generate.jl`.
const SPECIES = build_species(NUMSPECIES,
                              tolerance = (274.0K, 0.5K),
                              toleranceaxis = Temperature,
                              demand = 10.0kJ / day,
                              demandaxis = SolarRadiation,
                              dispersal = 1.0km,
                              movement = BirthOnlyMovement,
                              survival = 0.2,
                              abundance = INDIVIDUALS,
                              seed = 1)

const ECO = build_ecosystem(SPECIES, ENVIRONMENT, seed = 1,
                            distributed = :auto)

# --- run --------------------------------------------------------------------------------

# How many species this rank holds and how many individuals they have between them — the load-balance
# check the original printed. `rows_matrix` only exists on the distributed type, so a single-rank
# run (which `:auto` makes serial) reads the ordinary landscape instead.
function _localshare(eco)
    rows = eco isa EcoSISTEM.MPIEcosystem ? eco.abundances.rows_matrix :
           eco.abundances.matrix
    return (species = size(rows, 1), individuals = sum(rows))
end

# Reported in rank order through a barrier rather than by sleeping for `rank` seconds, which is
# what the original did — that cost `RANKS` seconds and still raced whenever a rank was slow.
function _reportinturn(message)
    for r in 0:(RANKS - 1)
        r == RANK && println(message)
        flush(stdout)
        RANKS > 1 && MPI.Barrier(COMM)
    end
    return nothing
end

before = _localshare(ECO)
_reportinturn("rank $RANK: $(before.species) species, " *
              "$(before.individuals) individuals before burn-in")

MPI.Barrier(COMM)
elapsed = @elapsed simulate!(ECO, BURNIN, TIMESTEP)

after = _localshare(ECO)
_reportinturn("rank $RANK: $(after.individuals) individuals after $(BURNIN), " *
              "$(round(elapsed, digits = 2)) s")

MPI.Barrier(COMM)
