# SPDX-License-Identifier: LGPL-3.0-or-later
#
# **Africa at whatever resolution this machine can afford.** 50,000 plant species - roughly the real
# number - over the whole continent, on real WorldClim temperature, run serially or across MPI ranks
# without changing a line.
#
# **Nothing here asserts a memory figure: it *measures* what is available and picks a resolution to
# match.** A hardcoded `@assert Sys.total_memory() >= 100GB` guards a run that actually needs about
# **8 GiB** - measured - so it is out by an order of magnitude and refuses most machines that could
# run it perfectly well.
#
# **The grid is chosen, not fixed, and nothing about its size is written down here.**
# `examples/HPC/memory.jl` works out what this run can allocate - across every node when launched
# under MPI - and each candidate cell size is costed by asking the package:
# `investigate_study_area` resolves the grid without building it, and `getspeciesstorage` says
# what one species' abundances would occupy on it. The finest that fits wins. So the same file gives
# a few-GiB run on a laptop and a multi-TiB one across HPC nodes, and there is no table to fall out
# of date.
#
# **5 km is multi-node only.** At 50,000 species it needs over 3 TiB, so no single node reaches it
# however generous - it is there for an MPI job spanning several.
#
# **Not part of the test suite**, and deliberately: `test/extras_examples.jl` runs only top-level
# `examples/*.jl`. Smoke-test it in seconds with `ECOSISTEM_SCALE=small`, which drops to 200 species
# on the coarsest grid.
#
#     ECOSISTEM_SCALE=small julia --project=examples examples/HPC/Africa.jl
#     julia -t 8 --project=examples examples/HPC/Africa.jl
#     mpiexecjl --project=examples -n 32 julia -t 8 --project=examples examples/HPC/Africa.jl
#     srun -N 4 --ntasks-per-node=8 julia -t 8 --project=examples examples/HPC/Africa.jl
#
# Everything is configured by environment variable, so a batch script needs no edits:
#
# | variable | default | |
# |---|---|---|
# | `ECOSISTEM_SCALE` | `large` | `small` for a fast smoke test |
# | `ECOSISTEM_AFRICA_SPECIES` | 50000 | size of the species pool |
# | `ECOSISTEM_AFRICA_CELLSIZE` | chosen | `20km` to override the memory-based choice |
# | `ECOSISTEM_AFRICA_YEARS` | 100 | simulated years after burn-in |
# | `ECOSISTEM_AFRICA_RECORD` | `none` | `none`, `interval` or `dates` |
# | `ECOSISTEM_AFRICA_SAVEDIR` | `./Africa_run` | where recordings go |

using EcoSISTEM
# `simulationtime`/`simulationdate` are `public` but not exported, so `using EcoSISTEM` alone
# does not bring them into scope - they have to be named.
using EcoSISTEM: simulationtime, simulationdate
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Unitful
using Unitful.DefaultSymbols
using Dates: Dates
using JLD2
using MPI

include(joinpath(@__DIR__, "memory.jl"))

# --- configuration ----------------------------------------------------------------------

# MPI is initialised even on a laptop, where it simply reports one rank. That keeps a single code
# path: `MemoryGuidance` and `build_ecosystem(distributed = :auto)` both ask the same question and
# both answer "serial" when there is only one rank.
MPI.Initialized() || MPI.Init()
const RANK = MPI.Comm_rank(MPI.COMM_WORLD)
const ISROOT = RANK == 0

const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"

# The cell sizes this run may choose between. Only the *sizes* are listed - what each costs is
# asked of the package below, never written down, because the reprojected extent does not scale
# linearly with cell size (each grid is rounded outwards to whole cells independently) and a
# hand-kept table would drift the moment anything about the area changed.
const CANDIDATES = (5.0km, 10.0km, 20.0km, 50.0km, 100.0km)

const NUMSPECIES = SMALL ? 200 :
                   parse(Int, get(ENV, "ECOSISTEM_AFRICA_SPECIES", "50000"))
const YEARS = (SMALL ? 2 :
               parse(Int, get(ENV, "ECOSISTEM_AFRICA_YEARS", "100"))) * year
const BURNIN = SMALL ? 1year : 10year
const TIMESTEP = 1month_mean_duration
const RECORD = Symbol(get(ENV, "ECOSISTEM_AFRICA_RECORD", "none"))
const SAVEDIR = get(ENV, "ECOSISTEM_AFRICA_SAVEDIR",
                    joinpath(pwd(), "Africa_run"))

# The epoch only matters for `RECORD = :dates`, but is set unconditionally so that
# `simulationdate` answers at all - without one a run has elapsed time and nothing else.
const EPOCH = Dates.Date(2000, 1, 1)
# Scaled with the run, or a smoke test would finish decades before the first one and `:dates`
# would appear to work while having saved nothing.
const SAVE_DATES = [Dates.Date(2000 + y, 1, 1)
                    for y in (SMALL ? (1, 2) : (10, 25, 50, 100))]

# --- the layers, named before the grid is decided ----------------------------------------

# Real temperature rather than the original's uniform 274 K: bio1 is WorldClim's annual mean, and
# its unit and axis come from the shipped catalogue rather than being asserted here. The supply is
# left uniform, so what varies across the continent is the *condition* species are matched to.
#
# These come **before** the grid is chosen, because choosing it means costing each candidate, and
# costing one means resolving the area these very layers would land on.
const TEMPERATURE = SourceSpec(WorldClim{BioClim}, :bio1)
const SUNLIGHT = UniformSpec(1000.0kJ / km^2 / day, axis = SolarRadiation)
# `within` positions the area and is not optional: WorldClim is global, so without it the grid
# would be the world.
const WITHIN = EcoSISTEM.boundingbox("Africa")
# A projected CRS is required to simulate - dispersal assumes one uniform cell size, and a degree
# cell's real extent shrinks towards the poles. EPSG:10592 (WGS 84 / GLANCE Africa) is the package's
# own suggestion for this extent.
const CRS = EPSG(10592)

# --- decide the grid --------------------------------------------------------------------

# Collective, and every rank must reach it - see `memory.jl`. Called before any `ISROOT` branch
# for exactly that reason.
const BUDGET = MemoryGuidance.memory_budget()

# What one species' abundances would occupy at `cellsize`, asked of the package rather than
# calculated here. `investigate_study_area` runs the *same* analysis `StudyArea` would - the report
# it returns "can never describe a grid other than the one that would be built" - so this cannot drift
# from what the run actually allocates, and it reads no raster data to answer.
function _perspecies(cellsize)
    report = investigate_study_area(regime = TEMPERATURE, supply = SUNLIGHT,
                                    within = WITHIN, crs = CRS,
                                    cellsize = cellsize)
    return EcoSISTEM.getspeciesstorage(report)
end

# Honour an explicit cell size if one was asked for, otherwise take the finest that fits. Parsed
# rather than `eval`ed: the value arrives as a string from the environment, and only these five are
# offered anyway.
function _requestedcellsize()
    want = get(ENV, "ECOSISTEM_AFRICA_CELLSIZE", "")
    isempty(want) && return nothing
    for cs in CANDIDATES
        string(cs) == want || "$(Int(ustrip(cs)))km" == want || continue
        return cs
    end
    return error("ECOSISTEM_AFRICA_CELLSIZE=$want is not one of " *
                 join(CANDIDATES, ", "))
end

# A smoke run takes the coarsest grid outright rather than the finest that fits - the whole point
# is to be quick, and on a large machine "what fits" is a grid that takes hours.
const CHOSEN = let want = _requestedcellsize()
    pool = SMALL ? (maximum(CANDIDATES),) :
           isnothing(want) ? CANDIDATES : (want,)
    MemoryGuidance.choose_cellsize(_perspecies, NUMSPECIES, pool,
                                   budget = BUDGET)
end

if ISROOT
    MemoryGuidance.describe(BUDGET)
    MemoryGuidance.describe(CHOSEN, NUMSPECIES)
    println()
end

# Stop rather than let the batch system discover this by killing the job an hour in. An explicit
# `ECOSISTEM_AFRICA_CELLSIZE` that does not fit lands here too, which is the point of checking it
# through the same path as the automatic choice.
CHOSEN.fits ||
    error("$(NUMSPECIES) species at $(CHOSEN.cellsize) does not fit in the memory available; " *
          "ask for more nodes, a coarser ECOSISTEM_AFRICA_CELLSIZE, or fewer species.")

# --- build ------------------------------------------------------------------------------

# The same four inputs the costing above used, now at the chosen resolution - so the grid built
# here is the one that was priced, not a second opinion about it.
const AREA = StudyArea(regime = TEMPERATURE, supply = SUNLIGHT, within = WITHIN,
                       crs = CRS, cellsize = CHOSEN.cellsize,
                       verbosity = ISROOT ? :normal : :silent)

const ENVIRONMENT = GridHabitat(regime = TEMPERATURE, supply = SUNLIGHT,
                                area = AREA)

# Generalists: a wide 50 K tolerance about 274 K, as in the original, so the pool coexists rather
# than partitioning the continent between narrow specialists.
const SPECIES = build_species(NUMSPECIES,
                              tolerance = (274.0K, 50.0K),
                              toleranceaxis = Temperature,
                              demand = 10.0kJ / day,
                              demandaxis = SolarRadiation,
                              dispersal = 15.0km,
                              # **`BirthOnlyMovement`, not the original's `AlwaysMovement`, and
                              # that is a package bug rather than a preference.** `move!`'s
                              # `AlwaysMovement` method (`src/Generate.jl:511-518`) is typed on
                              # `AbstractEcosystem`, so it catches an `MPIEcosystem` too, and reads
                              # `eco.abundances.matrix` - which an `MPIGridLandscape` has not got
                              # (`rows_matrix`/`cols_vector`). It fails on the first timestep of any
                              # distributed run. `test/SmallMPItest.jl:63` uses
                              # `BirthOnlyMovement`, so nothing exercised it.
                              movement = BirthOnlyMovement,
                              survival = 0.1,
                              abundance = SMALL ? 10^6 : 3 * 10^8,
                              seed = 1)

# `distributed = :auto` is the whole point of one file for both: an `MPIEcosystem` when the
# process is a live multi-rank MPI session, a serial `Ecosystem` otherwise.
const ECO = build_ecosystem(SPECIES, ENVIRONMENT, seed = 1, epoch = EPOCH,
                            distributed = :auto)

if ISROOT
    ny, nx = size(ENVIRONMENT.active)
    println("Africa: $(ny) × $(nx) cells of $(CHOSEN.cellsize), ",
            "$(count(ENVIRONMENT.active)) active, $(NUMSPECIES) species, ",
            BUDGET.distributed ? "$(BUDGET.ranks) MPI ranks." : "serial.")
end

# --- recording --------------------------------------------------------------------------

# The abundance matrix as a plain `species × cells` array on the root rank. Collective under MPI,
# and it allocates a **whole extra copy on rank 0** - at 50,000 species over a 20 km grid that is
# 53 GiB on one rank, against a per-rank share of far less. That is why `RECORD = :none` is the
# default and why the large configurations are not meant to be recorded: the run itself fits when
# the recording does not.
function _snapshot(eco)
    eco isa EcoSISTEM.MPIEcosystem && return gatherabundance(eco)
    return Array(eco.abundances.matrix)
end

# Write one recording. Only the root holds meaningful data after a gather, so only it writes.
function _save(eco, label)
    abun = _snapshot(eco)
    ISROOT || return nothing
    mkpath(SAVEDIR)
    path = joinpath(SAVEDIR, "Africa_$(label).jld2")
    jldsave(path, abun = abun, elapsed = simulationtime(eco),
            date = simulationdate(eco))
    println("  saved $(basename(path))")
    return nothing
end

# Run with no recording at all - the only mode the largest configurations can use, and the reason
# `simulate!` without a storage argument exists.
function _runsilently(eco)
    ISROOT && println("Running $(YEARS) with no recording.")
    simulate!(eco, YEARS, TIMESTEP)
    return nothing
end

# Run recording at a regular elapsed interval. `simulate_action!` takes its callback first, so
# this reads as a do-block, and it accepts an `MPIEcosystem` where `simulate_record!` does not.
function _runatintervals(eco, interval)
    ISROOT && println("Running $(YEARS), recording every $(interval).")
    simulate_action!(eco, YEARS, interval, TIMESTEP) do count
        return _save(eco, "step$(count)")
    end
    return nothing
end

# Run recording on specific real dates. The callback fires every timestep and saves when the clock
# has reached the next target date, which is how an irregular calendar schedule is expressed with
# only a regular one available. `simulationdate` returns `nothing` without an epoch, so this mode
# needs one - `EPOCH` above.
function _runatdates(eco, dates)
    isnothing(simulationdate(eco)) &&
        error("RECORD = :dates needs an epoch; build the ecosystem with `epoch = ...`.")
    pending = sort(filter(>(simulationdate(eco)), dates))
    ISROOT &&
        println("Running $(YEARS), recording at ", join(pending, ", "), ".")
    simulate_action!(eco, YEARS, TIMESTEP, TIMESTEP) do _
        (isempty(pending) || simulationdate(eco) < first(pending)) &&
            return nothing
        return _save(eco, Dates.format(popfirst!(pending), "yyyy-mm-dd"))
    end
    return nothing
end

# --- run --------------------------------------------------------------------------------

if ISROOT
    println("Burn-in $(BURNIN)...")
end
@time simulate!(ECO, BURNIN, TIMESTEP)

if RECORD === :none
    @time _runsilently(ECO)
elseif RECORD === :interval
    @time _runatintervals(ECO, SMALL ? 1year : 12month_mean_duration)
elseif RECORD === :dates
    @time _runatdates(ECO, SAVE_DATES)
else
    error("ECOSISTEM_AFRICA_RECORD=$(RECORD) is not one of none, interval, dates.")
end

# --- what happened ----------------------------------------------------------------------

# Collected once at the end whatever the recording mode, so even a silent run reports something.
# Collective, so every rank calls it and only the root reads the result.
const FINAL = _snapshot(ECO)

if ISROOT
    occupied = count(>(0), sum(FINAL, dims = 1))
    surviving = count(>(0), sum(FINAL, dims = 2))
    println("\nAfter $(BURNIN + YEARS) (", simulationdate(ECO), "):")
    println("  total abundance   ", sum(FINAL))
    println("  species surviving ", surviving, " of ", NUMSPECIES)
    println("  cells occupied    ", occupied, " of ",
            length(ENVIRONMENT.active))
end

BUDGET.distributed && MPI.Barrier(MPI.COMM_WORLD)
