```@meta
CurrentModule = EcoSISTEM
```

# Running a simulation at scale

A continental run at fine resolution is bound by memory, not by processor time, and it may not fit
on one machine. This page covers the two things that follow from that: **sizing a run before you
build it**, and **spreading one across MPI ranks**.

## Size the run before building it

The abundance matrix is allocated over the **whole** grid, one entry per species per cell, and an
inactive cell costs exactly as much as an active one. So the memory a run needs is decided the
moment the grid is, and it can be asked before anything is constructed.

[`investigate_study_area`](@ref) reports what a [`StudyArea`](@ref) *would* be without building it,
and [`getspeciesstorage`](@ref) can be asked of that report:

```julia
using EcoSISTEM, EcoSISTEM.Units, Unitful, Unitful.DefaultSymbols

report = investigate_study_area(regime = SourceSpec(WorldClim{BioClim}, :bio1),
                                cellsize = 10.0km)
perspecies = EcoSISTEM.getspeciesstorage(report)   # bytes, for ONE species
```

Two things to know before trusting the number:

  - it is **one array**, and a run holds several — abundances and net migration serially, more when
    distributed — so multiply by the species count *and* by the number of arrays;
  - it counts every cell, active or not, which is why a coarser `cellsize` is the effective lever.

That is how `examples/HPC/Africa.jl` picks its own resolution: it measures the memory available,
including across the nodes of an MPI job, and chooses the finest `cellsize` that fits rather than
asserting one.

## Running across MPI ranks

Load `MPI` and the extension activates. [`build_ecosystem`](@ref) then chooses for you:

```julia
using MPI, EcoSISTEM
MPI.Init()

eco = build_ecosystem(species, env, seed = 1)      # distributed = :auto, the default
```

`distributed = :auto` builds a distributed ecosystem when MPI is initialised **and** there is more
than one rank; `mpirun -n 1` or a bare `using MPI` still gives an ordinary [`Ecosystem`](@ref).
Pass `distributed = true` or `false` to force the choice.

Each rank owns a block of species and a block of grid cells, because the work needs both views: the
demographics run per species and dispersal runs per cell. Bringing a result back to one rank for
recording is explicit — `gatherabundance` for the abundances, `gatherdiversity` for a diversity
measure.

### The same seed gives the same answer

Reproducibility across ranks is a design requirement, not a convenience. Each species has its own
deterministic random stream addressed by its **global** index, so a result does not depend on how
the work was divided: one rank or sixteen, one thread or many, the same seed gives the same numbers.

That requirement is also why a layer's change over time must be a pure function of elapsed time.
Layers are updated redundantly on every rank, so anything drawn from a shared random stream, or
depending on the ecosystem's own state, would let the ranks drift apart silently. Change that
*does* depend on the ecosystem goes through [Interventions](@ref) instead, which are applied once
and identically everywhere.

### Movement

All three movement types work distributed. Each rank owns a block of species across the *whole*
grid while the dynamics run, so an individual dispersing between cells stays rank-local whether it
is a newborn under [`BirthOnlyMovement`](@ref) or an established individual under
[`AlwaysMovement`](@ref); [`NoMovement`](@ref) disperses nothing at all.

### What is not supported

The [`AddSpecies`](@ref) operation is refused in a distributed run. Species are what the ranks are
partitioned by, so adding one changes the partition itself, where every other operation acts within
a rank. Add the species before the run starts, or run that scenario serially.

## Threads

Threading needs nothing declared: start Julia with `-t` and the per-species loop is shared out.
The same reproducibility guarantee applies, for the same reason — the random stream is chosen by
species, never by which thread happened to take the work.

```bash
julia -t 8 --project myrun.jl
```

## Practical notes

  - **Measure before you scale.** A run that does not fit is an error worth getting at startup
    rather than after an hour, which is what sizing from the report is for.
  - **Reading real rasters is the other memory peak**, and it is separate from the simulation's:
    a continental read can transiently need more than the run itself.
  - **Recording costs memory too.** [`generate_storage`](@ref) allocates
    `species × cells × recordings × replicates` up front and cannot grow afterwards, so an
    intervention that adds species needs `maxspecies` set for them in advance.

## Where to go next

  - [Virtual plant simulations of Africa](@ref) — a worked continental run, including the
    MPI invocation.
  - [Interventions](@ref) — ecosystem-level change, and why it is a separate mechanism from
    layer change.
