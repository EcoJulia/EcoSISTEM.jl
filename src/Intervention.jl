# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The three answers bound together - when, where and what - and a set of them.

using Unitful

using Unitful.DefaultSymbols

using EcoSISTEM.Units

using DimensionalData

using Random

using StatsBase: sample

using Distributions: Binomial, Poisson

"""
    Intervention(schedule, region, operations...)

A declared change to the **ecosystem** - as against an [`AbstractLayerChange`](@ref), which changes a
single layer. The parts answer *when* ([`AbstractSchedule`](@ref)), *where*
([`AbstractRegion`](@ref), or a shape) and *what* ([`AbstractOperation`](@ref)):

```julia
Intervention(AtTime(50year), RandomCells(20), Deactivate())
Intervention(EveryStep(), AllCells(), SetChange(:rainfall, IncrementBy(-0.1mm/day/year)))
Intervention(AtTime(10year), NaturalEarthSpec("Fife"), Deactivate())
```

A shape is a region as a study area's `within` is: the cells it covers at least half of, or those a
[`ShapeMaskSpec`](@ref)'s rule picks. It is placed on the grid the habitat was built on, which must
have a coordinate reference system, once for each ecosystem before the run's first step, and the
records of its files join the ecosystem's when the intervention first acts.

**Several operations share one resolved region**, applied in the order written - clear an area and
plant a crop on it, or reactivate a site and restock it:

```julia
Intervention(AtTime(10year), RandomCells(20), Deactivate(), AddAbundance(:wheat, 500))
```

That sharing is the whole reason for the varargs, and it cannot be had from two interventions in a
set: each resolves its **own** region, and a random one draws its own cells from its own stream (`k`
is its index in the set), so two `RandomCells(20)` would pick different cells. A region that must be
acted on twice has to be resolved once.

Pass one, or an [`InterventionSet`](@ref), as `simulate!`'s `intervention` keyword.

An intervention built from published data - a land-use scenario, the occurrence records a
population is seeded from - says so with `provenance`, an [`ProvenanceRecord`](@ref) or a vector of them.
The records join the ecosystem's own the first time the intervention acts, so [`provenance`](@ref)
of the ecosystem lists an intervention that has acted and not one that never fired.
"""
struct Intervention{S <: AbstractSchedule, RG <: AbstractRegion, O <: Tuple}
    schedule::S
    region::RG
    operations::O
    inputs::Vector{ProvenanceRecord}

    # The sole constructor. Varargs, so the one-operation form is unchanged and needs no second
    # spelling; at least one is required, since an intervention that does nothing is not one.
    function Intervention(schedule::AbstractSchedule,
                          region::Union{AbstractRegion, AbstractShapeSpec,
                                        ShapeMaskSpec},
                          operations::AbstractOperation...;
                          provenance = ProvenanceRecord[])
        isempty(operations) &&
            error("an `Intervention` needs at least one operation: a schedule and a region with " *
                  "nothing to do are not an intervention. Use `NeverScheduled()` to disable one.")
        _checktableintervention(schedule, operations)
        cells = _asregion(region)
        return new{typeof(schedule), typeof(cells),
                   typeof(operations)}(schedule, cells, operations,
                                       _recordvector(provenance))
    end
end

"""
    InterventionSet(interventions...)

Several [`Intervention`](@ref)s applied in the order written, every step. Any number of them, of any
kinds, including several of the same kind.
"""
struct InterventionSet{T <: Tuple}
    interventions::T

    # The sole constructor: a bare tuple field would otherwise be a second, unlabelled spelling.
    function InterventionSet(interventions::Intervention...)
        return new{typeof(interventions)}(interventions)
    end
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# The expression that builds it, as for the specs. The schedule, region and operations are all
# small fieldless or few-field types that already print as they are written - `AtTime(5.0 yr)`,
# `RandomCells(20)`, `Deactivate()` - so nesting them costs nothing and keeps the line honest.
function Base.show(io::IO, i::Intervention)
    ops = join(map(o -> sprint(show, o), i.operations), ", ")
    return print(io,
                 "Intervention(", sprint(show, i.schedule), ", ",
                 sprint(show, i.region), ", ", ops, ")")
end

function Base.show(io::IO, s::InterventionSet)
    n = length(s.interventions)
    return print(io,
                 "InterventionSet($(n) intervention$(n == 1 ? "" : "s"))")
end

function Base.show(io::IO, ::MIME"text/plain", s::InterventionSet)
    println(io, sprint(show, s))
    for i in s.interventions
        println(io, "  ", sprint(show, i))
    end
    return nothing
end

# A dated placement's rows for one ecosystem, this rank's species only, in time order: each one's row
# in this rank's abundances, storage column, count and elapsed time; the elapsed time served through,
# `nothing` before the first firing; and the individuals passed over at a continuing run's start, due
# in cells deactivated since, and added at the start to a state the run started from, until a run's
# end notes them.
mutable struct _ServedRows
    rows::Vector{Int}
    cells::Vector{Int}
    counts::Vector{Int}
    times::Vector{typeof(1.0s)}
    servedthrough::Union{Nothing, typeof(1.0s)}
    skipped::Int
    deactivated::Int
    ontostate::Int
end

# What an `AddSpecies` has found for the one ecosystem it last acted on: the rows of the species it
# names, found before the run's first step, and what the run's end notes - each species it found in
# the run already, with how often, and how many species a `RandomSpecies` asked for and found none
# left to take.
mutable struct _ArrivalState
    owner::WeakRef
    rows::Vector{Int}
    already::Dict{String, Int}
    short::Int
end

using Unitful

using Unitful.DefaultSymbols

using EcoSISTEM.Units

using DimensionalData

using Random

using StatsBase: sample

using Distributions: Binomial, Poisson

# == Functions ==================================================================================

"""
    applyinterventions!(eco::AbstractEcosystem, intervention, elapsed, timestep, step)

Apply every scheduled [`Intervention`](@ref) for the step that has just reached `elapsed`. Called by
[`update!`](@ref) - after the population dynamics and **before** the layer update, so that a
[`SetChange`](@ref) takes effect in the same step rather than one step late - and, at elapsed zero,
once before the first step's dynamics, so an intervention due at the start acts on the starting
state.

**This is also the supported imperative route**: it is the way to act on a running ecosystem
without having declared anything in advance, and it takes any operation, not just some privileged
subset. To bring a species into a run mid-way, having passed no `intervention` to
[`simulate!`](@ref):

```julia
simulate!(eco, 1year, 1month_mean_duration)
applyinterventions!(eco,
                     Intervention(EveryStep(), AllCells(),
                                  AddSpecies("invader", BySupply(total = 500))),
                     simulationtime(eco), 1month_mean_duration, 12)
simulate!(eco, 1year, 1month_mean_duration)     # carries on with the new species
```

**Prefer declaring an intervention and passing it to `simulate!`.** Calling this yourself gives up
the two guarantees the declarative form provides, and `step` is where both bite:

  - **Reproducibility.** A selection's random stream is seeded from `(seed, :intervention, k,
    step)`, so a `step` you have already used **reuses that stream** - a [`RandomCells`](@ref) region would draw
    the very same cells again - and a wrong `step` means the run no longer follows from its seed.
    Pass the step number the simulation has actually reached.
  - **Determinism across MPI ranks.** An intervention mutates the ecosystem, so it must be applied
    once and identically everywhere. `update!` guarantees that; a hand call does not.

For a *conditional* intervention - one whose firing depends on the state of the run, which no
schedule can express - give [`simulate!`](@ref) a callback: it closes over the ecosystem and is
handed the elapsed time it fires at, so it can decide and then act.
"""
function applyinterventions!(eco::AbstractEcosystem, intervention,
                             elapsed::Unitful.Time, timestep::Unitful.Time,
                             step::Integer)
    for (k, iv) in enumerate(_interventions(intervention))
        _fires(iv, eco, elapsed, timestep) || continue
        rng = _interventionrng(eco, k, step)
        # Resolved **once**, then every operation acts on the same cells - see `Intervention`.
        cells = _regioncells(iv.region, eco, rng, timestep)
        _applyoperations!(iv, eco, cells, rng, elapsed, timestep,
                          (k = k, step = step))
        _recordintervention!(eco, iv)
    end
    return eco
end

# The change rules a run's interventions would have installed by the time its clock stands at, for a
# run continuing a saved state: every firing of a `SetChange` applied again in the order
# `applyinterventions!` applies it - firing by firing, the interventions of one step in order - over
# steps of `timestep` from zero, the last ending on the clock wherever this run's steps would fall.
# Only rules are installed, since every other operation's effect is in the state, so nothing is drawn,
# no cells are resolved and every rank installs the same.
function _replaychanges!(eco::AbstractEcosystem, intervention,
                         timestep::Unitful.Time)
    changing = filter(iv -> any(_installsrule, iv.operations),
                      collect(_interventions(intervention)))
    isempty(changing) && return eco
    final = simulationtime(eco)
    step = uconvert(s, float(timestep))
    elapsed = zero(final)
    window = step
    while true
        for iv in changing
            _fires(iv, eco, elapsed, window) || continue
            foreach(op -> _installrule!(op, eco), iv.operations)
            _recordintervention!(eco, iv)
        end
        elapsed < final || return eco
        reached = elapsed + step
        next = final - reached <= _DRIFT * step ? final : reached
        window = next - elapsed
        elapsed = next
    end
end

# Whether an operation installs a change rule, which a saved state does not keep.
_installsrule(::AbstractOperation) = false

_installsrule(::SetChange) = true

# Install an operation's change rule, and do nothing for an operation installing none.
_installrule!(::AbstractOperation, eco::AbstractEcosystem) = eco

function _installrule!(op::SetChange, eco::AbstractEcosystem)
    return setchange!(_targetlayer(eco, op.target), op.spec)
end

# ---------------------------------------------------------------------------
# Ecosystem-level intervention
# ---------------------------------------------------------------------------
# **An intervention is a declaration, not a callback.** The three questions it answers - **when**,
# **where** and **what** - are each a type, so an intervention can be dispatched on, validated,
# reported and reproduced. A function reference could do none of those: two interventions differing
# only in their function would be the same type, and what either did would be invisible until it ran.
#
# This is the same move `AbstractLayerChange` made for layer-level change (see `Layer.jl`), and for
# the same reason. The two mechanisms are deliberately separate: a *layer* change is a pure function
# of elapsed time and may be applied redundantly on every MPI rank, whereas an intervention mutates
# the ecosystem - abundances, the active mask - and so must be applied once, identically, everywhere.

# Whether a schedule fires on the step that has just advanced the clock to `elapsed`, having covered
# `timestep`. The step *covers* `(elapsed - timestep, elapsed]`, and a one-off schedule fires when
# its instant falls in that half-open window - so it fires exactly once however the steps are sized,
# and never falls between two steps.
#
# At elapsed zero no step has been taken, so the start is not one: a schedule naming an instant at or
# before the start fires there, on the starting state, while one that acts over a step - every step,
# or every step within a span - waits for the first.
_fires(::EveryStep, elapsed, _) = !iszero(elapsed)

_fires(::NeverScheduled, _, _) = false

function _fires(s::AtTime, elapsed::Unitful.Time, timestep::Unitful.Time)
    return _reached(s.time, elapsed, timestep)
end

function _fires(s::AtTimes, elapsed::Unitful.Time, timestep::Unitful.Time)
    return any(t -> _reached(t, elapsed, timestep), s.times)
end

function _fires(s::BetweenTimes, elapsed::Unitful.Time, ::Unitful.Time)
    return !iszero(elapsed) && s.from <= elapsed <= s.to
end

# Only the multiples that could lie in this step's window are tried, each through `_reached` exactly
# as `AtTimes` tries its list, so the two agree step for step.
function _fires(s::EveryInterval, elapsed::Unitful.Time,
                timestep::Unitful.Time)
    lower = max(0,
                floor(Int, ustrip(NoUnits, (elapsed - timestep) / s.interval)))
    upper = ceil(Int, ustrip(NoUnits, elapsed / s.interval))
    return any(k -> _reached(k * s.interval, elapsed, timestep), lower:upper)
end

# A schedule asked with the ecosystem to hand. One reading elapsed time alone needs nothing from it;
# a date schedule reads the epoch and calendar that turn its dates into elapsed time.
function _fires(schedule::AbstractSchedule, ::AbstractEcosystem, elapsed,
                timestep)
    return _fires(schedule, elapsed, timestep)
end

function _fires(s::AtDates, eco::AbstractEcosystem, elapsed::Unitful.Time,
                timestep::Unitful.Time)
    epoch = _scheduleepoch(s, eco)
    return any(s.dates) do date
        return _reached(_elapsedat(epoch, date, eco.calendar), elapsed,
                        timestep)
    end
end

# Only the anniversaries in the years this step's window spans are tried, and none before the epoch.
function _fires(s::EveryYear, eco::AbstractEcosystem, elapsed::Unitful.Time,
                timestep::Unitful.Time)
    epoch = _scheduleepoch(s, eco)
    window = max(elapsed - timestep, zero(elapsed))
    fromyear = Dates.year(_shiftdate(epoch, window, eco.calendar))
    toyear = Dates.year(_shiftdate(epoch, elapsed, eco.calendar))
    return any(fromyear:toyear) do year
        t = _elapsedat(epoch, Dates.DateTime(year, s.month, s.day),
                       eco.calendar)
        return t >= zero(t) && _reached(t, elapsed, timestep)
    end
end

# The epoch a date schedule is placed from, refused where the run has none.
function _scheduleepoch(schedule::AbstractSchedule, eco::AbstractEcosystem)
    isnothing(eco.epoch) || return eco.epoch
    return error("`$(sprint(show, schedule))` names real dates, but this run has no epoch to " *
                 "place them from. Pass `epoch = ` to `build_ecosystem` to say which date the " *
                 "run starts on.")
end

# Find, before a run's first step, the cells of every region that names places or a shape, and the
# species every `AddSpecies` names: a place off the grid or a name reaching no species is refused
# here rather than when it fires, and a shape is placed on the grid once for the ecosystem, by every
# rank at the same point, since under MPI that is a collective build, as looking up a name is.
function _prepareinterventions!(intervention, eco::AbstractEcosystem)
    for iv in _interventions(intervention)
        _prepareregion!(iv.region, eco)
        foreach(op -> _prepareoperation!(op, eco), iv.operations)
    end
    return nothing
end

_prepareoperation!(::AbstractOperation, ::AbstractEcosystem) = nothing

function _prepareoperation!(op::AddSpecies, eco::AbstractEcosystem)
    _arrivalstate!(op, eco)
    return nothing
end

_prepareregion!(::AbstractRegion, ::AbstractEcosystem) = nothing

function _prepareregion!(r::CellsAt, eco::AbstractEcosystem)
    _placecells(r, eco)
    return nothing
end

function _prepareregion!(r::_ShapeCells, eco::AbstractEcosystem)
    _shapecells!(r, eco)
    return nothing
end

# Refuse, before a run's first step, a date schedule the run cannot honour: one on a run with no
# epoch, and one naming a date that falls inside a step longer than a day rather than at the end of
# one, which would be acted on late - and, as the steps drift against the calendar, by a different
# amount each time.
function _checkschedules(intervention, eco::AbstractEcosystem,
                         duration::Unitful.Time, timestep::Unitful.Time)
    foreach(iv -> _checkschedule(iv.schedule, eco, duration, timestep),
            _interventions(intervention))
    return nothing
end

_checkschedule(::AbstractSchedule, ::AbstractEcosystem, _, _) = nothing

# A step of a day or less shows every date on its own day; a longer one must end on each date the run
# reaches, to the same tolerance a series' slice is found to.
function _checkschedule(schedule::Union{AtDates, EveryYear},
                        eco::AbstractEcosystem, duration::Unitful.Time,
                        timestep::Unitful.Time)
    epoch = _scheduleepoch(schedule, eco)
    timestep <= 1.0u"d" && return nothing
    start = simulationtime(eco)
    final = _finalelapsed(eco, duration, timestep)
    step = uconvert(s, float(timestep))
    for (date, t) in _scheduleddates(schedule, epoch, eco.calendar, start,
                                     final)
        start < t <= final || continue
        k = ustrip(NoUnits, (t - start) / step)
        abs(k - round(k)) <= _DRIFT && continue
        acted = _shiftdate(epoch, start + ceil(Int, k) * step, eco.calendar)
        error("`$(sprint(show, schedule))` names $(Dates.Date(date)), which falls inside a " *
              "step of $timestep rather than at the end of one, so it would be acted on " *
              "$(Dates.Date(acted)). " * _scheduleremedy(eco.calendar))
    end
    return nothing
end

# The dates a date schedule names that a run from `start` to `final` could reach, each with its
# elapsed time under the run calendar.
function _scheduleddates(schedule::AtDates, epoch, calendar, _, _)
    return [(date = d, elapsed = _elapsedat(epoch, d, calendar))
            for d in schedule.dates]
end

function _scheduleddates(schedule::EveryYear, epoch, calendar, start, final)
    fromyear = Dates.year(_shiftdate(epoch, start, calendar))
    toyear = Dates.year(_shiftdate(epoch, final, calendar))
    dates = [Dates.DateTime(y, schedule.month, schedule.day)
             for y in fromyear:toyear]
    return [(date = d, elapsed = _elapsedat(epoch, d, calendar)) for d in dates]
end

# What to do about a date that falls inside a step, under each run calendar.
function _scheduleremedy(::ExactDates)
    return "Step by no more than a day, or - for dates on the first of a month - build the " *
           "ecosystem with `calendar = MeanMonths()` and step by a whole number of " *
           "`month_mean_duration`s, which then end on the first of every month."
end

function _scheduleremedy(::MeanMonths)
    return "Step by no more than a day, or schedule dates on the first of a month and step by a " *
           "whole number of `month_mean_duration`s."
end

# `(elapsed - timestep, elapsed]` - half-open below so consecutive steps cannot both claim the same
# instant, closed above so an instant landing exactly on a step boundary fires on that step.
#
# **At elapsed zero the window is everything up to zero**, which `update!` asks before the first
# step's dynamics: an instant at or before the start acts on the starting state, and the first step's
# window, `(0, timestep]`, then excludes it, so it fires once. Without that start an `AtTime(0)` would
# never fire, which is the failure this function exists to avoid.
function _reached(t::Unitful.Time, elapsed::Unitful.Time,
                  timestep::Unitful.Time)
    t <= elapsed || return false
    return iszero(elapsed) || t > (elapsed - timestep)
end

# One intervention answers as a set of one, so the apply path needs no separate leaf case - the same
# shape a bare layer takes against `LayerCollection`.
_interventions(set::InterventionSet) = set.interventions

_interventions(intervention::Intervention) = (intervention,)

_interventions(::Nothing) = ()

# Add the records of an intervention that has just acted to the ecosystem's own, so what a run says
# it was built from names the interventions that changed it, the files of a shape it acted over among
# them. Every rank fires the same interventions on the same step, so every rank's records agree.
function _recordintervention!(eco::AbstractEcosystem, iv::Intervention)
    inputs = vcat(iv.inputs, _regioninputs(iv.region, eco))
    isempty(inputs) && return eco
    eco.inputs = _uniqueinputs(vcat(eco.inputs, inputs))
    return eco
end

_regioninputs(::AbstractRegion, ::AbstractEcosystem) = ProvenanceRecord[]

function _regioninputs(r::_ShapeCells, eco::AbstractEcosystem)
    return _shapecells!(r, eco).inputs
end

# ---------------------------------------------------------------------------
# Applying an intervention
# ---------------------------------------------------------------------------
# **Reproducibility is the constraint that shapes all of this.** Counter-based per step: a stream
# seeded from `(seed, :intervention, k, step)` for intervention `k` on step `step`. That
# generalises the per-species `(seed, j)` scheme, and buys three things at once - every MPI rank and
# every thread computes bit-identical selections without communicating; a run replays exactly from
# any step; and species streams stay reserved for birth/death/dispersal, so adding an intervention
# cannot re-phase the demography.
#

# The word standing for `:intervention` in a stream's seed - "inter" in ASCII - so that an
# intervention's seed can never coincide with a species' `(seed, j)`. A literal rather than
# `hash(:intervention)`, because `Base.hash` changes between Julia versions and the seed must not.
const _INTERVENTION_STREAM = 0x696e746572

# The stream for intervention `k` on the step reaching `elapsed`. Keyed on the *step number*, not
# on elapsed time: a float is a poor seed, and the step index is what "replay from step n" means.
# Seeded through the generator's own seeding, never `Base.hash`, so the stream is the same on every
# Julia version.
function _interventionrng(eco::AbstractEcosystem, k::Integer, step::Integer)
    return Random.Xoshiro(UInt64[eco.seed, _INTERVENTION_STREAM, k, step])
end

# **A count may be exact or drawn.** An integer means exactly that many; a **rate** means each
# candidate is taken independently at that rate over the step, which is a binomial draw and is what a
# fixed count cannot express. The rate is per unit time, so it is the timestep that turns it into a
# probability, and a rate high enough to exceed 1 over a step is clamped - "more certain than
# certain" has no meaning.
_drawcount(count::Integer, _, _, _) = Int(count)

function _drawcount(rate::Unitful.Quantity, candidates::Integer,
                    timestep::Unitful.Time, rng)
    p = min(one(Float64), ustrip(NoUnits, rate * timestep))
    return rand(rng, Binomial(candidates, p))
end

# The linear cell indices a region resolves to.
function _regioncells(::AllCells, eco::AbstractEcosystem, _, _)
    return eachindex(parent(eco.habitat.active))
end

function _regioncells(::ActiveCells, eco::AbstractEcosystem, _, _)
    return findall(vec(parent(eco.habitat.active)))
end

function _regioncells(r::CellMask, eco::AbstractEcosystem, _, _)
    size(r.mask) == size(parent(eco.habitat.active)) ||
        error("a `CellMask` must be on the habitat's own grid: the mask is $(size(r.mask)) but " *
              "the grid is $(size(parent(eco.habitat.active))).")
    return findall(vec(r.mask))
end

# Each place's cell, in order, and a place in an inactive cell moved to the nearest active one within
# `snap`, where there is one; each cell once, in grid order, so the places' order changes nothing.
function _regioncells(r::CellsAt, eco::AbstractEcosystem, _, _)
    active = parent(eco.habitat.active)
    cells = _placecells(r, eco)
    if !iszero(r.snap)
        yx = _gridyx(eco)
        cells = map(cells) do cell
            active[cell] && return cell
            target = _nearestactive(yx, active, cell, r.snap)
            return isnothing(target) ? cell : target
        end
    end
    linear = LinearIndices(active)
    return sort!(unique(linear[cell] for cell in cells))
end

function _regioncells(r::_ShapeCells, eco::AbstractEcosystem, _, _)
    return _shapecells!(r, eco).cells
end

# The cell of each of a `CellsAt`'s places in `eco`'s grid, a place off the grid refused by name.
function _placecells(r::CellsAt, eco::AbstractEcosystem)
    shape = size(parent(eco.habitat.active))
    return map(place -> _placecell(place, eco, shape), r.places)
end

function _placecell(key::CartesianIndex{2}, _, shape::Tuple)
    checkbounds(Bool, CartesianIndices(shape), key) ||
        error("`CellsAt` names the cell $(Tuple(key)), which is off the $(shape[1]) by " *
              "$(shape[2]) grid.")
    return key
end

function _placecell(place::SpatialLocation, eco::AbstractEcosystem, ::Tuple)
    cell = only(getcellat(eco, [place]))
    isnothing(cell) &&
        error("`CellsAt` names $(sprint(show, place)), which is off the grid.")
    return cell
end

# A region as an `Intervention` keeps it: a shape, bare or under a rule, as the cells it covers.
_asregion(region::AbstractRegion) = region

function _asregion(shape::AbstractShapeSpec)
    return _ShapeCells(shape, ShapeMaskSpec(shape), Ref{Any}(nothing))
end

_asregion(spec::ShapeMaskSpec) = _ShapeCells(spec, spec, Ref{Any}(nothing))

function _asregion(spec::ShapeMaskSpec{<:AbstractShapeSpec, WholeBoundingBox})
    return error("`WholeBoundingBox()` says where a study area's grid is cut, which an " *
                 "intervention's grid already is; give the shape with `AnyOverlap()` for every " *
                 "cell it touches, or a `FractionWithin` share.")
end

# The cells a shape region covers in `eco`'s grid, and the records of the shape's files, found the
# first time it is asked for that ecosystem and kept until it is asked for another: the shape's
# covered share of each cell of the grid the habitat was built on, under the region's rule.
function _shapecells!(r::_ShapeCells, eco::AbstractEcosystem)
    found = r.found[]
    !isnothing(found) && found.owner.value === eco && return found
    isnothing(_gridcrs(eco)) &&
        error("`$(sprint(show, r))` is a shape, so placing it needs a grid with a coordinate " *
              "reference system, and this ecosystem's is synthetic. Name its cells with a " *
              "`CellMask` or `CellsAt` instead.")
    area = StudyArea(eco.habitat)
    covered = materialise(ShapeCoverage(r.spec.shape, axis = NicheAxis), area)
    mask = map(f -> _rulepasses(r.spec.rule, f), parent(covered.matrix))
    # The copied area's cache is its own, so it holds the shape's records alone; taken whole, since
    # under MPI they are the root's, keyed by the root's own spec.
    inputs = reduce(vcat, values(area.report.cache.inputs),
                    init = ProvenanceRecord[])
    r.found[] = (owner = WeakRef(eco), cells = findall(vec(mask)),
                 inputs = _uniqueinputs(inputs))
    return r.found[]
end

function _regioncells(r::RandomCells, eco::AbstractEcosystem, rng, timestep)
    active = findall(vec(parent(eco.habitat.active)))
    wanted = _drawcount(r.count, length(active), timestep, rng)
    return sample(rng, active, min(wanted, length(active)), replace = false)
end

# A contiguous cluster, grown one neighbour at a time. Re-seeds if the cluster is boxed in by
# inactive cells rather than returning short, so `SpreadingCells(n)` gives `n` cells whenever that
# many are available at all - otherwise a run's outcome would depend on the shape of the mask in a
# way nobody asked for.
function _regioncells(r::SpreadingCells, eco::AbstractEcosystem, rng, timestep)
    height = getgridshape(eco)[1]
    active = Set(findall(vec(parent(eco.habitat.active))))
    wanted = min(_drawcount(r.count, length(active), timestep, rng),
                 length(active))
    chosen = Int[]
    frontier = Int[]
    while length(chosen) < wanted
        if isempty(frontier)
            push!(frontier, rand(rng, collect(active)))
        end
        cell = popfirst!(frontier)
        cell in active || continue
        delete!(active, cell)
        push!(chosen, cell)
        for n in _neighbourcells(eco, cell, height)
            n in active && push!(frontier, n)
        end
    end
    return chosen
end

# The four-connected neighbours of a linear cell index that lie on the grid.
function _neighbourcells(eco::AbstractEcosystem, cell::Int, height::Int)
    (y, x) = convert_coords(eco, cell, height)
    width = getgridshape(eco)[2]
    out = Int[]
    for (dy, dx) in ((-1, 0), (1, 0), (0, -1), (0, 1))
        ny, nx = y + dy, x + dx
        (1 <= ny <= height && 1 <= nx <= width) || continue
        push!(out, convert_coords(eco, (ny, nx), height))
    end
    return out
end

# `active` is indexed through its `parent`: a `DimArray` has no dimensions to match a bare linear
# index against.
function _applyoperation!(::Deactivate, eco::AbstractEcosystem, cells, _, _)
    parent(eco.habitat.active)[cells] .= false
    # ...and everything there dies. See `Deactivate` for why this is not optional.
    _ownedabundances(eco).rows[:, cells] .= 0
    return eco
end

function _applyoperation!(::Reactivate, eco::AbstractEcosystem, cells, _, _)
    parent(eco.habitat.active)[cells] .= true
    return eco
end

# The one permitted direct matrix write, and only on a categorical layer - see `SetLandCover`.
function _applyoperation!(op::SetLandCover, eco::AbstractEcosystem, cells, _, _)
    layer = eco.habitat.regime
    layer isa CategoricalLayer ||
        error("`SetLandCover` writes a categorical layer's values directly, which is the only " *
              "direct write an intervention may make - but this environment's regime is a " *
              "$(nameof(typeof(layer))). A continuous layer's values belong to its change rule; " *
              "use `SetChange` to alter them.")
    parent(layer.matrix)[cells] .= _landcovercode(op.class)
    return eco
end

# A layer's change applies to the whole layer, so the region is meaningless here - said rather than
# silently ignored, because a caller who wrote one expected it to matter.
function _applyoperation!(op::SetChange, eco::AbstractEcosystem, cells, _, _)
    length(cells) == length(parent(eco.habitat.active)) ||
        @warn "`SetChange` alters a whole layer's change rule, so the intervention's region is " *
              "ignored; use `AllCells()` to say so." maxlog=1
    return _installrule!(op, eco)
end

function _applyoperation!(op::AddAbundance, eco::AbstractEcosystem, cells, rng,
                          timestep)
    row = _localrow(eco, op.species)
    isnothing(row) && return eco          # another rank owns this species
    owned = _ownedabundances(eco)
    # Drawn **per cell** when the count is a rate, so arrivals vary across the region as they
    # would in reality, rather than every cell receiving the same number.
    for cell in cells
        owned.rows[row, cell] += _drawnumber(op.count, timestep, rng)
    end
    return eco
end

# Clamped at zero per cell: removing more than are there cannot leave a negative abundance.
function _applyoperation!(op::RemoveAbundance, eco::AbstractEcosystem, cells,
                          rng, timestep)
    row = _localrow(eco, op.species)
    isnothing(row) && return eco
    owned = _ownedabundances(eco)
    for cell in cells
        owned.rows[row, cell] = max(0,
                                    owned.rows[row, cell] -
                                    _drawnumber(op.count, timestep, rng))
    end
    return eco
end

# The abundances this process owns, and the global species index its first row is.
# **At the point interventions run the landscape is row-partitioned**, so a rank holds whole
# *species* across **all** cells (`synchronise_from_rows!` has already run - `ext/EcoSISTEMMPIExt/dynamics.jl`). That
# is what makes every abundance operation rank-local: a rank clears or adds to its own species at the
# named cells, every rank does the same, and between them they cover the whole assemblage. No global
# draw and no partition arithmetic is needed.
# The MPI extension supplies the sole method for its own landscape.
function _ownedabundances(eco::AbstractEcosystem)
    return (rows = eco.abundances.matrix, firstspecies = 1)
end

# A land-cover class as the code the layer stores. A `Symbol` is looked up in the shipped catalogue;
# anything else is already a code and passes through, so an intervention can name a class either way.
_landcovercode(class) = class

_landcovercode(class::Symbol) = landcoverclass(class)

# Address a sub-layer of a collection by name, or the layer itself when unnamed.
_targetlayer(eco::AbstractEcosystem, ::Nothing) = eco.habitat.regime

function _targetlayer(eco::AbstractEcosystem, target::Symbol)
    for side in (eco.habitat.regime, eco.habitat.supply)
        side isa LayerCollection || continue
        target in propertynames(side) && return getproperty(side, target)
    end
    return error("no layer named `$target` in this environment. The regime holds " *
                 "$(keys(eco.habitat.regime)) and the supply " *
                 "$(keys(eco.habitat.supply)).")
end

# An exact count, or a Poisson draw about a mean rate - an arrival process, which a fixed number
# cannot express.
_drawnumber(count::Integer, _, _) = Int(count)

function _drawnumber(rate::Unitful.Quantity, timestep::Unitful.Time, rng)
    return rand(rng, Poisson(ustrip(NoUnits, rate * timestep)))
end

# Which local row a species is, or `nothing` when this process does not own it. Serially every
# species is owned, so this is the identity; under MPI it selects the one rank that acts.
function _localrow(eco::AbstractEcosystem, species)
    owned = _ownedabundances(eco)
    row = _speciesindex(eco, species) - owned.firstspecies + 1
    return (1 <= row <= Base.size(owned.rows, 1)) ? row : nothing
end

# By index or by name, matching how `getdispersaldist` already accepts either. Defined against the
# `SpeciesList` as well as the ecosystem, because the public accessors take either - see
# `Ecosystem.jl` - and a name lookup needs only the names.
_speciesindex(::SpeciesList, sp::Integer) = sp

function _speciesindex(spplist::SpeciesList, sp::AbstractString)
    idx = findfirst(==(sp), spplist.names)
    isnothing(idx) &&
        error("no species named `$sp` in this species list.")
    return idx
end

function _speciesindex(spplist::SpeciesList, sp::Symbol)
    return _speciesindex(spplist, String(sp))
end

_speciesindex(eco::AbstractEcosystem, sp) = _speciesindex(eco.spplist, sp)

# ---------------------------------------------------------------------------
# Adding abundances by their dates
# ---------------------------------------------------------------------------
# An `AddAbundance` of a dated placement carries its own times: its intervention fires every step and
# at the start, and it is handed the elapsed time each step reaches, which says which rows are due.

# A dated placement stands alone in an intervention scheduled every step: a schedule firing less often
# would skip rows, and an operation sharing it would act at the start as well.
function _checktableintervention(schedule::AbstractSchedule, operations::Tuple)
    any(_istimed, operations) || return nothing
    schedule isa EveryStep && length(operations) == 1 && return nothing
    return error("an `AddAbundance` of a dated placement carries its own times, so it goes in an " *
                 "`Intervention` of its own with `EveryStep()`, whose region filters its rows; this " *
                 "one has `$(sprint(show, schedule))` and $(length(operations)) operations.")
end

# Whether an operation carries its own times, its rows served as they fall due.
_istimed(::AbstractOperation) = false

_istimed(op::_PlacementAddition) = op.dated

# Whether an intervention fires on the step reaching `elapsed`: when its schedule does, except that a
# dated placement's fires at the start too, to add the rows due at or before it.
function _fires(iv::Intervention, eco::AbstractEcosystem, elapsed, timestep)
    return _fires(iv.schedule, eco, elapsed, timestep)
end

function _fires(iv::Intervention{EveryStep, <:AbstractRegion,
                                 <:Tuple{_PlacementAddition}},
                eco::AbstractEcosystem, elapsed, timestep)
    return only(iv.operations).dated ||
           _fires(iv.schedule, eco, elapsed, timestep)
end

# Apply an intervention's operations in order to its resolved cells.
function _applyoperations!(iv::Intervention, eco::AbstractEcosystem, cells,
                           rng, elapsed, timestep, firing::NamedTuple)
    at = merge(firing, (elapsed = elapsed,))
    foreach(op -> _applyoperation!(op, eco, cells, rng, timestep, at),
            iv.operations)
    return eco
end

# An operation applied at a firing - which intervention, on which step, at what time. Every operation
# but those placing individuals needs only its cells, the firing's stream and the timestep.
function _applyoperation!(op::AbstractOperation, eco::AbstractEcosystem, cells,
                          rng, timestep, ::NamedTuple)
    return _applyoperation!(op, eco, cells, rng, timestep)
end

# A placement's individuals added over the region's active cells, as `build_ecosystem` places them,
# its table undated: under MPI the root calls it and every rank adds its own species' rows.
function _applyoperation!(op::_PlacementAddition, eco::AbstractEcosystem,
                          cells, _, timestep, firing::NamedTuple)
    op.dated && return _servedated!(op, eco, cells, timestep, firing)
    context = _firingcontext(eco, op.placement, cells, firing)
    return _placeinto!(eco, op.placement, context, :normal, clear = false,
                       keep = false, undated = true)
end

# Species joining the run: those named or drawn, each counted where it was in the run already, made
# active in the species list, then their arrivals placed over the region's active cells as
# `AddAbundance` of a placement places, the arrivals the only species placed. Under MPI every rank
# draws and activates the same species in its copy of the list.
function _applyoperation!(op::AddSpecies, eco::AbstractEcosystem, cells, rng,
                          timestep, firing::NamedTuple)
    state = _arrivalstate!(op, eco)
    rows = _arrivingrows!(op.species, state, eco, rng, timestep)
    isempty(rows) && return eco
    active = eco.spplist.traits.active
    for row in rows
        active[row] || continue
        name = eco.spplist.names[row]
        state.already[name] = get(state.already, name, 0) + 1
    end
    arriving = falses(length(active))
    arriving[rows] .= true
    active[rows] .= true
    context = _firingcontext(eco, op.placement, cells, firing,
                             active = arriving)
    return _placeinto!(eco, op.placement, context, :normal, clear = false,
                       keep = false, undated = true)
end

# The rows of the species arriving: those named, as found before the run; or a draw from the species
# not in the run, from the firing's stream, its shortfall counted.
_arrivingrows!(::Vector{String}, state, ::AbstractEcosystem, _, _) = state.rows

function _arrivingrows!(r::RandomSpecies, state, eco::AbstractEcosystem, rng,
                        timestep)
    candidates = findall(!, eco.spplist.traits.active)
    wanted = _drawcount(r.count, length(candidates), timestep, rng)
    taken = min(wanted, length(candidates))
    state.short += wanted - taken
    return sort!(sample(rng, candidates, taken, replace = false))
end

# What an `AddSpecies` has found for `eco`: kept while it acts on that ecosystem, made afresh for
# another, its named species' rows found then.
function _arrivalstate!(op::AddSpecies, eco::AbstractEcosystem)
    held = op.arrivals[]
    !isnothing(held) && held.owner.value === eco && return held
    state = _ArrivalState(WeakRef(eco), _namedrows(op, eco),
                          Dict{String, Int}(),
                          0)
    op.arrivals[] = state
    return state
end

# The rows of the species an `AddSpecies` names, looked up as a placement's names are - on the first
# rank under MPI, which alone reads a name authority - a name reaching no species a problem at its
# level; none for `RandomSpecies`, which draws its species when it acts.
_namedrows(::AddSpecies{<:RandomSpecies}, ::AbstractEcosystem) = Int[]

function _namedrows(op::AddSpecies{Vector{String}}, eco::AbstractEcosystem)
    index = _matchindex(eco.spplist)
    answers = _madeonroot() do
        answers = _nameanswers(index, op.species)
        _reportproblems(_unnamedproblems(op, answers, index), :normal)
        return answers
    end
    return unique(answers[n].row
                  for n in op.species if !isnothing(answers[n].row))
end

# What an `AddSpecies` reports of the names reaching no species of the run, at its level: those no
# species answers to, and those its pool's name authority recognises but no species of the pool is.
function _unnamedproblems(op::AddSpecies, answers, index)
    lost = [n for n in op.species if isnothing(answers[n].row)]
    problems = Problem[]
    recognised = filter(n -> answers[n].recognised, lost)
    unknown = setdiff(lost, recognised)
    leftout = op.level isa ProblemError ? "" : " It adds no such species."
    isempty(unknown) ||
        push!(problems,
              Problem(op.level, :not_recognised,
                      "`AddSpecies` names " *
                      _namelist(["`$n`" for n in unknown], :normal) *
                      ", " *
                      (isnothing(index.authority) ?
                       "which no species of the run is" :
                       "which $(_authorityphrase(index.authority)) does not recognise") *
                      ": every species a run may hold is declared in its pool before it is " *
                      "built." * leftout))
    isempty(recognised) ||
        push!(problems,
              Problem(op.level, :not_in_pool,
                      "`AddSpecies` names " *
                      _namelist(["`$n`" for n in recognised], :normal) *
                      ", recognised by $(_authorityphrase(index.authority)) but no species of " *
                      "the run's pool: every species a run may hold is declared in its pool " *
                      "before it is built." * leftout))
    return problems
end

# The context of a placement firing in intervention `firing.k` on step `firing.step` over `cells`, the
# region's storage columns: the active ones among them, the species it places where they are not
# those in the run, a stream per species, keyed on its key in its pool as at build, and a community
# stream, all made for that firing and apart from the stream the region was chosen from; and, where
# the placement fills to capacity, what is already in those cells.
function _firingcontext(eco::AbstractEcosystem, placement, cells,
                        firing::NamedTuple; active = nothing)
    live = parent(eco.habitat.active)
    grid = CartesianIndices(size(live))
    incells = CartesianIndex{2}[grid[c] for c in cells if live[c]]
    words = UInt64[eco.seed % UInt64, _INTERVENTION_STREAM, firing.k,
                   firing.step,
                   _PLACEMENT_STREAM]
    streams = Dict(name => Random.Xoshiro(vcat(words, key % UInt64))
                   for (name, key) in zip(eco.spplist.names,
                                          eco.spplist.traits.key))
    present = _needsdemand(placement) ? _presentdemand(eco, incells) : nothing
    return PlacementContext(cells = incells, streams = streams,
                            community = Random.Xoshiro(words),
                            nichefit = eco.nichefit, epoch = eco.epoch,
                            calendar = eco.calendar, time = firing.elapsed,
                            demand = present, active = active)
end

# What the individuals already in each of `cells` demand of each resource, a grid-shaped array for
# each in the order of the supplies, every species counted: each rank sums its own species, and the
# ranks' sums are added.
function _presentdemand(eco::AbstractEcosystem,
                        cells::AbstractVector{CartesianIndex{2}})
    owned = _ownedabundances(eco)
    linear = LinearIndices(getgridshape(eco))
    shape = getgridshape(eco)
    return map(values(eco.spplist.demand)) do d
        u = unit(eltype(d.resource))
        present = zeros(shape)
        for cell in cells, row in axes(owned.rows, 1)
            sp = owned.firstspecies + row - 1
            present[cell] += owned.rows[row, linear[cell]] *
                             ustrip(u, d.resource[sp])
        end
        return _sumacrossranks!(eco, present) .* u
    end
end

# Each rank's partial sums added across the ranks; serially, as they are.
_sumacrossranks!(::AbstractEcosystem, x::AbstractArray) = x

# Whether this is the rank that reports what every rank found; serially, the one there is.
_isrootrank(::AbstractEcosystem) = true

# A dated placement's rows falling due on the step reaching `firing.elapsed` added where they fall in
# the region's still-active cells: at the start every row at or before it, otherwise those in
# `(elapsed - timestep, elapsed]`. At a continuing run's first firing the rows it has passed are
# counted, unless the run started from a saved state, which holds them already; so are those due in
# cells deactivated since the ecosystem was built, and those added at the start to a saved state.
function _servedated!(op::_PlacementAddition, eco::AbstractEcosystem, cells,
                      timestep::Unitful.Time, firing::NamedTuple)
    state = _servedrows!(op, eco, firing)
    elapsed = firing.elapsed
    from = iszero(elapsed) ? 1 :
           searchsortedlast(state.times, elapsed - timestep) + 1
    fromstate = _fromstate(eco)
    isnothing(state.servedthrough) && !fromstate &&
        (state.skipped += sum(view(state.counts, 1:(from - 1)), init = 0))
    ontostate = fromstate && iszero(elapsed)
    to = searchsortedlast(state.times, elapsed)
    owned = _ownedabundances(eco)
    active = vec(parent(eco.habitat.active))
    inregion = falses(length(active))
    inregion[cells] .= true
    for i in from:to
        cell = state.cells[i]
        inregion[cell] || continue
        if active[cell]
            owned.rows[state.rows[i], cell] += state.counts[i]
            ontostate && (state.ontostate += state.counts[i])
        else
            state.deactivated += state.counts[i]
        end
    end
    state.servedthrough = elapsed
    return eco
end

# A dated placement's rows for `eco`, resolved at its first firing on it and kept: another ecosystem
# resolves them afresh.
function _servedrows!(op::_PlacementAddition, eco::AbstractEcosystem,
                      firing::NamedTuple)
    held = op.served[]
    !isnothing(held) && held.owner.value === eco && return held.rows
    rows = _resolvedated(op.placement, eco, firing)
    op.served[] = (owner = WeakRef(eco), rows = rows)
    return rows
end

# A dated placement resolved as `build_ecosystem` resolves a placement - under MPI the root calling it
# and every rank taking the table - but over the cells the ecosystem was built with, keeping every row
# whatever its time, with streams of the intervention's own; this rank's species' rows in time order.
function _resolvedated(placement, eco::AbstractEcosystem, firing::NamedTuple)
    built = parent(eco.habitat.area.report.active)
    words = UInt64[eco.seed % UInt64, _INTERVENTION_STREAM, firing.k, 0,
                   _PLACEMENT_STREAM]
    streams = Dict(name => Random.Xoshiro(vcat(words, key % UInt64))
                   for (name, key) in zip(eco.spplist.names,
                                          eco.spplist.traits.key))
    context = PlacementContext(cells = findall(built), streams = streams,
                               community = Random.Xoshiro(words),
                               nichefit = eco.nichefit, epoch = eco.epoch,
                               calendar = eco.calendar, time = Inf * s)
    outcome = _resolveplacement(eco, placement, context, :normal)
    _reportproblems(vcat(outcome.placed.problems, outcome.resolved.problems),
                    :normal)
    eco.inputs = _uniqueinputs(vcat(eco.inputs, outcome.placed.records))
    resolved = outcome.resolved
    isnothing(resolved.times) &&
        error("$(_placementphrase(placement)) is dated, but its table has no `time` column, so " *
              "no row says when it falls due.")
    owned = _ownedabundances(eco)
    rows = resolved.species .- owned.firstspecies .+ 1
    mine = findall(in(axes(owned.rows, 1)), rows)
    order = mine[sortperm(resolved.times[mine])]
    return _ServedRows(rows[order], resolved.cells[order],
                       resolved.counts[order], resolved.times[order], nothing,
                       0, 0, 0)
end

# At a run's end, each operation's notice of what the run did that it was not told: the root
# reporting, once every rank's counts are in.
function _finishinterventions!(eco::AbstractEcosystem, intervention)
    for iv in _interventions(intervention), op in iv.operations
        _finishoperation!(op, eco)
    end
    return eco
end

_finishoperation!(::AbstractOperation, ::AbstractEcosystem) = nothing

# A dated placement's notice of what it did not add - individuals dated after the run's end, passed
# over at a continuing run's start, or due in cells deactivated since - and of what it added at the
# start to the state the run started from: every rank's counted.
function _finishoperation!(op::_PlacementAddition, eco::AbstractEcosystem)
    op.dated || return nothing
    held = op.served[]
    (!isnothing(held) && held.owner.value === eco) || return nothing
    state = held.rows
    late = sum(state.counts[i]
               for i in eachindex(state.times)
               if state.times[i] > simulationtime(eco); init = 0)
    counts = _sumacrossranks!(eco,
                              [late, state.skipped, state.deactivated,
                                  state.ontostate])
    state.skipped = state.deactivated = state.ontostate = 0
    _isrootrank(eco) || return nothing
    problems = Problem[]
    for (n, what) in ((counts[1],
        "dated after the run's end $(_were(counts[1])) not added"),
        (counts[2],
        "dated at or before the run's start $(_were(counts[2])) not added: a run " *
        "continuing another adds only what falls due as it runs, so give the " *
        "intervention to the run before it too"),
        (counts[3],
        "due in cells deactivated since the ecosystem was built " *
        "$(_were(counts[3])) not added"),
        (counts[4],
        "dated at or before the run's start $(_were(counts[4])) added to the saved state it " *
        "started from"))
        iszero(n) ||
            push!(problems,
                  Problem(ProblemNotice(), :not_added,
                          "$(_individualsword(n)) $what."))
    end
    _reportproblems(problems, :normal)
    return nothing
end

# An `AddSpecies`' notice of the species it found in the run already, which had their arrivals all the
# same, and of the species a `RandomSpecies` asked for and found none left to take. Every rank drew
# and found the same.
function _finishoperation!(op::AddSpecies, eco::AbstractEcosystem)
    held = op.arrivals[]
    (!isnothing(held) && held.owner.value === eco) || return nothing
    already = sort!(collect(held.already))
    short = held.short
    empty!(held.already)
    held.short = 0
    _isrootrank(eco) || return nothing
    problems = Problem[]
    isempty(already) ||
        push!(problems,
              Problem(ProblemNotice(), :already_in_run,
                      "`AddSpecies` found " *
                      _namelist(["`$n` ($(_countword(k, "time")))"
                                 for (n, k) in already],
                                :normal) *
                      " in the run already, and added their arrivals all the same."))
    iszero(short) ||
        push!(problems,
              Problem(ProblemNotice(), :none_left,
                      "`$(sprint(show, op.species))` asked for $short more species than were " *
                      "left out of the run, so fewer arrived."))
    _reportproblems(problems, :normal)
    return nothing
end
