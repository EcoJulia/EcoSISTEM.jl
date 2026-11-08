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
    RandomSpecies(count)

`count` species drawn at random from those not yet in the run, each time an [`AddSpecies`](@ref)
acts: exactly that many (`RandomSpecies(2)`), or each species not yet in the run arriving
independently at a rate over the step (`RandomSpecies(0.1/year)`), a binomial draw. They are drawn
from the intervention's own stream, as [`RandomCells`](@ref) are, so every rank and thread draws the
same species. Where fewer are left than it asks for, those there are arrive, and the run's end says
how many it could not take.

# Arguments

  - `count`: how many species, a whole number or a rate per unit time, neither negative.
"""
struct RandomSpecies{C <: CELLCOUNT}
    count::C

    # `CELLCOUNT`, for the reason given on `AddAbundance`.
    function RandomSpecies(count::CELLCOUNT)
        count >= zero(count) ||
            error("`RandomSpecies` takes how many species arrive, a number or a rate that is not " *
                  "negative; got $count.")
        return new{typeof(count)}(count)
    end
end

"""
    AddSpecies(species, placement; level = ProblemError())

Bring species into a run - an invasion, a reintroduction, a crop sown for the first time: each of
`species` joins the run, and its first individuals are placed by `placement` over the intervention's
region.

```julia
Intervention(AtTime(5year), CellsAt(LatLong(56.0°, -3.2°)),
             AddSpecies("Rhododendron ponticum", BySupply(each = 500)))
Intervention(EveryStep(), ActiveCells(), AddSpecies(RandomSpecies(0.1 / year), BySupply(each = 50)))
```

**Every species a run may hold is declared before it is built**, in its [`SpeciesPool`](@ref), and one
not in the run at the start is deactivated there ([`deactivate_species!`](@ref)): so an arrival brings
traits of its own, and nothing about the run changes size, under MPI or otherwise.

The placement places the arrivals alone, as though they were the only species in the run, worked out
over the region's active cells with random streams of the firing's own, as [`AddAbundance`](@ref) of
a placement is; its rows naming any other species are dropped, and reported. A species already in the
run gets its arrivals all the same, and the run's end says which did.

Its names are looked up before the run's first step, as a placement's are - a species' own name, a
name of its pool's data, an identifier, or through its pool's name authority - and a name reaching no
species of the run is a problem at `level`.

# Arguments

  - `species`: which species - a name, a vector of names, or [`RandomSpecies`](@ref) drawn from those
    not yet in the run.
  - `placement`: where and how many arrive, as any undated placement [`build_ecosystem`](@ref) takes -
    `BySupply(each = 500)`, a [`FromTable`](@ref), a [`ToCapacity`](@ref), or a function of
    `(species, habitat, context)`. A species with dated records is declared in the run with none
    placed at the start, and the records added by [`AddAbundance`](@ref) of the dated placement.
  - `level`: what a name reaching no species of the run is - a [`ProblemError`](@ref), stopping the
    run before it starts, by default; `ProblemWarning()` or `ProblemNotice()` to report it and go on
    without it.
"""
struct AddSpecies{S, P, L <: AbstractProblemSeverity} <: AbstractOperation
    species::S
    placement::P
    level::L
    # What it found for the one ecosystem it last acted on (`_ArrivalState`), `nothing` before.
    arrivals::Base.RefValue{Any}

    function AddSpecies(species::Union{AbstractString, Symbol, AbstractVector,
                                       RandomSpecies}, placement;
                        level::AbstractProblemSeverity = ProblemError())
        _checkplacement(placement)
        _isdated(placement) &&
            error("`AddSpecies` adds its arrivals when it acts, so it takes an undated placement, " *
                  "and $(_placementphrase(placement)) is dated. Declare the species in the run with " *
                  "none placed at the start, and add its records with `AddAbundance` of the dated " *
                  "placement, in an `Intervention` of its own with `EveryStep()`.")
        names = _arrivingnames(species)
        return new{typeof(names), typeof(placement), typeof(level)}(names,
                                                                    placement,
                                                                    level,
                                                                    Ref{Any}(nothing))
    end
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# As for the schedules and regions: the call the caller wrote, not the type signature. An
# `AddSpecies` names its `level` only where it is not the default.
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

Base.show(io::IO, r::RandomSpecies) = print(io, "RandomSpecies($(r.count))")
function Base.show(io::IO, o::AddSpecies)
    print(io, "AddSpecies(", _arrivalsphrase(o.species), ", ",
          _placementphrase(o.placement))
    o.level isa ProblemError ||
        print(io, ", level = ", nameof(typeof(o.level)), "()")
    return print(io, ")")
end

# The species an `AddSpecies` names, as their names: one, or several, each a name; or a draw.
_arrivingnames(name::Union{AbstractString, Symbol}) = [String(name)]

_arrivingnames(r::RandomSpecies) = r

function _arrivingnames(names::AbstractVector)
    isempty(names) &&
        error("`AddSpecies` needs at least one species: a name, a vector of names, or " *
              "`RandomSpecies(count)`.")
    all(n -> n isa Union{AbstractString, Symbol}, names) ||
        error("`AddSpecies` names its species by name; got $(repr(names)).")
    return unique(String.(names))
end

# An `AddSpecies`' species as written: a name, a vector of names, or the draw.
function _arrivalsphrase(names::Vector{String})
    return length(names) == 1 ? repr(only(names)) : repr(names)
end

_arrivalsphrase(r::RandomSpecies) = sprint(show, r)
