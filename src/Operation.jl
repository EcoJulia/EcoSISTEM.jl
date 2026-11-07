# SPDX-License-Identifier: LGPL-3.0-or-later
#
# *What* an intervention does - a deliberately closed set. Closed because a callback could draw from
# the global RNG or write a continuous layer's matrix, either of which desynchronises MPI ranks;
# every operation here is reproducible from the run's seed alone.

using Unitful

using Tables: Tables

"""
    AbstractOperation

**What** an [`Intervention`](@ref) does - a **closed set of seven**: [`Deactivate`](@ref),
[`Reactivate`](@ref), [`SetLandCover`](@ref), [`SetChange`](@ref), [`AddAbundance`](@ref),
[`RemoveAbundance`](@ref) and [`AddSpecies`](@ref).

**Closed on purpose.** A callback would let an intervention do anything, including the things that
break reproducibility and MPI - writing a continuous layer's matrix directly, drawing from the global
RNG, resizing the landscape mid-run. Seven named operations can each be checked once and then
trusted.

**A bound is not among them**, and does not need to be. It is a property of the *quantity* rather
than of an intervention - a supply cannot be negative because of what a resource is - so it is
enforced on the layer automatically; see `bounds` and `_enforcebounds!`.
"""
abstract type AbstractOperation end

"""
    Deactivate() <: AbstractOperation

Destroy the region's cells: **everything living there dies**, and nothing grows there again until the
cells are reactivated.

**Killing is the point, not a side effect.** A deactivated cell is skipped entirely by the hot loop,
so anything left in it would neither breed nor die - a frozen population with no ecological meaning.
Left alive, stranded individuals come to outnumber the living, because they hold their numbers while
the real population turns over around them.

[`Reactivate`](@ref) does **not** bring them back, and should not: it makes the cell habitable again
so that dispersal can recolonise it, as vegetation returns to a slag heap once it stops being used.
"""
struct Deactivate <: AbstractOperation end

"""
    Reactivate() <: AbstractOperation

Mark the region's cells active again.

A deactivated cell **keeps its supply** - `_zerogaps` zeroes data gaps only - so a reactivated cell
has resource waiting. That matters: a cell whose supply had been zeroed would be lethal rather than
merely barren, since the death term is `E/K` and `K = 0` gives a death probability of exactly one.

**It does not restock the cell**, and that asymmetry with [`Deactivate`](@ref) is deliberate: the
habitat becomes available again and dispersal recolonises it from the surviving neighbours, which is
how a site actually recovers. Actively replanting is a different act - combine [`AddAbundance`](@ref)
with this over the same region - and is a restoration scenario rather than a release of constraint.
"""
struct Reactivate <: AbstractOperation end

"""
    SetLandCover(class)

Set the region's cells to a land-cover class.

**The only operation permitted to write a layer's matrix directly, and only a categorical one.** A
continuous layer's values are owned by its change rule, which is a pure function of elapsed time;
writing them from outside would break that and desynchronise MPI ranks. A categorical layer has no
such rule, so land cover is the one thing an intervention may set.

# Arguments

  - `class`: the class to set the cells to, by **name** (`:open_water`) or by its code.
"""
struct SetLandCover{C} <: AbstractOperation
    class::C
end

"""
    SetChange(target, spec)

Install a new change rule on one layer, part-way through a run:

```julia
SetChange(:rainfall, IncrementBy(-0.1mm / day / year))
```

Region-independent by nature - a layer's change applies to the whole layer - so an intervention
carrying it should use [`AllCells`](@ref). Any other region is accepted but ignored, and says so.

# Arguments

  - `target`: which layer, by name. This addresses a **member of a collection**, so a multi-variable
    environment can say which variable is the one that changes.
  - `spec`: the new change, as an [`AbstractChangeSpec`](@ref).
"""
struct SetChange{T, S} <: AbstractOperation
    target::T
    spec::S
end

"""
    AddAbundance(species, count)
    AddAbundance(placement)

Add individuals to **each** cell of the region - an invasion, or a reintroduction - or, given a
placement, the individuals it places over the region's active cells at each firing, as
[`build_ecosystem`](@ref) places them at the start: `AddAbundance(ToCapacity("wheat"))` sows a crop,
filling each cell to what its resources support given what already grows there. A dated placement - a
[`FromTable`](@ref) with a `time` column, or a [`FromOccurrences`](@ref) given `date` or `time` - is
served instead: alone in an `Intervention` with `EveryStep()`, its rows are added as the clock
reaches them, those at or before the start onto the starting state.

# Arguments

  - `species`: which species, by name or by index.
  - `count`: how many per cell. An `Integer` is an exact number; a quantity per unit time is a mean
    rate, drawn afresh each step.
  - `placement`: a placement - a shipped one, or a function of `(species, habitat, context)` - called
    at each firing with the region's active cells as its context's `cells` and random streams of that
    firing's own, so it draws nothing a region's selection or a species' dynamics draw. Under MPI the
    first rank calls it, and every rank adds its own species' rows. A dated placement is called once
    for an ecosystem, at its first firing, over the cells the ecosystem was built with; each step adds
    its rows falling due in the region's cells still active, and a run's end notes what it did not
    add - rows dated after it, rows a continuing run started past, and rows due in cells deactivated
    since.
"""
struct AddAbundance{S, C <: CELLCOUNT} <: AbstractOperation
    species::S
    count::C

    # `CELLCOUNT`, shared with the regions, which face the same "an exact number or a rate" choice.
    # It admits any `Integer` - a perfectly ordinary `typemax(Int32)` would be a `MethodError`
    # against `Int` - and pins a rate's dimension to `𝐓^-1`, so a length is refused at the call site
    # rather than where the rate is multiplied by the timestep.
    function AddAbundance(species::S, count::CELLCOUNT) where {S}
        return new{S, typeof(count)}(species, count)
    end
end

# The placement form of `AddAbundance`: a placement's individuals added over the region at each firing,
# or, where its table is `dated`, its rows as they fall due, resolved once for the ecosystem `served`
# holds them for.
struct _PlacementAddition{P} <: AbstractOperation
    placement::P
    dated::Bool
    served::Base.RefValue{Any}
end

function AddAbundance(placement)
    _checkplacement(placement)
    return _PlacementAddition(placement, _isdated(placement), Ref{Any}(nothing))
end

"""
    RemoveAbundance(species, count)

Remove individuals from each cell of the region, never taking a cell below zero - a cull, or a
harvest.

# Arguments

  - `species`: which species, by name or by index.
  - `count`: at most how many per cell. An `Integer` is an exact number; a quantity per unit time is
    a mean rate, drawn afresh each step.
"""
struct RemoveAbundance{S, C <: CELLCOUNT} <: AbstractOperation
    species::S
    count::C

    # `CELLCOUNT`, for the reason given on `AddAbundance`.
    function RemoveAbundance(species::S, count::CELLCOUNT) where {S}
        return new{S, typeof(count)}(species, count)
    end
end

"""
    AddSpecies(tolerance = nothing, demand = nothing, dispersal = nothing,
               birth = nothing, death = nothing, abundance, name = nothing)

Introduce a **new species** - an invasion, a reintroduction, a crop sown for the first time -
scattered across the intervention's region.

**It carries the arrival's own traits**, which is what it is for: an invader can be given a niche of
its own rather than inheriting one. Any trait left as `nothing` clones the last species instead,
which is what a plain reintroduction of something already present wants.

An existing species' dispersal can be handed straight over, since [`speciesdispersal`](@ref) returns
exactly what `dispersal` takes:

```julia
AddSpecies(dispersal = EcoSISTEM.speciesdispersal(eco, 3), abundance = 500)
```

**The recording has to have room.** A [`RecordAbundance`](@ref) writes into an array sized before the
run starts, so a run that gains species needs `RecordAbundance(eco, ntimes, maxspecies = ...)`, or
storage from `generate_storage(eco, times, reps, maxspecies = ...)`; see [`generate_storage`](@ref). **Not supported under MPI**: species are partitioned across ranks, so
adding one changes the partition itself, and unlike the abundance operations this is not rank-local.

# Arguments

  - `abundance`: how many individuals arrive in total, scattered over the region. The one required
    argument.
  - `tolerance`, `demand`: the arrival's niche and its resource requirements. Cloned from the last
    species if omitted.
  - `dispersal`: a dispersal **kernel**, an [`AbstractKernel`](@ref), matching
    `build_species`' keyword of the same name. The arrival disperses in the last species' mode,
    with this kernel.
  - `birth`, `death`: the arrival's demographic rates.
  - `name`: what to call it.
"""
struct AddSpecies{T, D, M, B, X} <: AbstractOperation
    tolerance::T
    demand::D
    dispersal::M
    birth::B
    death::X
    abundance::Int
    name::Union{Nothing, String}

    function AddSpecies(; tolerance = nothing, demand = nothing,
                        dispersal = nothing, birth = nothing, death = nothing,
                        abundance::Integer, name = nothing)
        return new{typeof(tolerance), typeof(demand), typeof(dispersal),
                   typeof(birth), typeof(death)}(tolerance, demand, dispersal,
                                                 birth, death, Int(abundance),
                                                 name)
    end
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# As for the schedules and regions: the call the caller wrote, not the type signature. `AddSpecies`
# carries five optional fields, so it reports only the ones that were set - which is what makes it
# readable inside an `Intervention` alongside the others.
Base.show(io::IO, o::SetLandCover) = print(io, "SetLandCover($(o.class))")
function Base.show(io::IO, o::SetChange)
    return print(io, "SetChange($(repr(o.target)), ", sprint(show, o.spec), ")")
end
function Base.show(io::IO, o::AddAbundance)
    return print(io, "AddAbundance($(repr(o.species)), $(o.count))")
end
function Base.show(io::IO, o::_PlacementAddition)
    return print(io, "AddAbundance(", _placementphrase(o.placement), ")")
end
function Base.show(io::IO, o::RemoveAbundance)
    return print(io, "RemoveAbundance($(repr(o.species)), $(o.count))")
end

function Base.show(io::IO, o::AddSpecies)
    parts = String[]
    for f in fieldnames(typeof(o))
        v = getfield(o, f)
        isnothing(v) || push!(parts, "$(f) = $(v)")
    end
    return print(io, "AddSpecies(", join(parts, ", "), ")")
end
