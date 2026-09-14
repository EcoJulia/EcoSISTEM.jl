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
([`AbstractRegion`](@ref)) and *what* ([`AbstractOperation`](@ref)):

```julia
Intervention(AtTime(50year), RandomCells(20), Deactivate())
Intervention(EveryStep(), AllCells(), SetChange(:rainfall, IncrementBy(-0.1mm/day/year)))
```

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
population is seeded from - says so with `provenance`, an [`InputRecord`](@ref) or a vector of them.
The records join the ecosystem's own the first time the intervention acts, so [`provenance`](@ref)
of the ecosystem lists an intervention that has acted and not one that never fired.
"""
struct Intervention{S <: AbstractSchedule, RG <: AbstractRegion, O <: Tuple}
    schedule::S
    region::RG
    operations::O
    inputs::Vector{InputRecord}

    # The sole constructor. Varargs, so the one-operation form is unchanged and needs no second
    # spelling; at least one is required, since an intervention that does nothing is not one.
    function Intervention(schedule::AbstractSchedule, region::AbstractRegion,
                          operations::AbstractOperation...;
                          provenance = InputRecord[])
        isempty(operations) &&
            error("an `Intervention` needs at least one operation: a schedule and a region with " *
                  "nothing to do are not an intervention. Use `NeverScheduled()` to disable one.")
        return new{typeof(schedule), typeof(region),
                   typeof(operations)}(schedule, region, operations,
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
subset. To add a species mid-run, having passed no `intervention` to [`simulate!`](@ref):

```julia
simulate!(eco, 1year, 1month_mean_duration)
applyinterventions!(eco,
                     Intervention(EveryStep(), AllCells(),
                                  AddSpecies(tolerance = Normal(298.0, 2.0),
                                             abundance = 500)),
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
        _fires(iv.schedule, eco, elapsed, timestep) || continue
        rng = _interventionrng(eco, k, step)
        # Resolved **once**, then every operation acts on the same cells - see `Intervention`.
        cells = _regioncells(iv.region, eco, rng, timestep)
        foreach(op -> _applyoperation!(op, eco, cells, rng, timestep),
                iv.operations)
        _recordintervention!(eco, iv)
    end
    return eco
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
# it was built from names the interventions that changed it. Every rank fires the same interventions
# on the same step, so every rank's records agree.
function _recordintervention!(eco::AbstractEcosystem, iv::Intervention)
    isempty(iv.inputs) && return eco
    eco.inputs = _uniqueinputs(vcat(eco.inputs, iv.inputs))
    return eco
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
    return setchange!(_targetlayer(eco, op.target), op.spec)
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

# The region is deliberately unused: the arrival is scattered by the same resource-weighted
# `repopulate!` every species is placed by, so where it lands is the environment's business rather
# than the intervention's. Saying so beats silently ignoring a region the caller wrote meaning
# something.
function _applyoperation!(op::AddSpecies, eco::AbstractEcosystem, cells, _, _)
    length(cells) == length(parent(eco.habitat.active)) ||
        @warn "`AddSpecies` scatters the arrival across the whole grid by resource availability, " *
              "so the intervention's region is ignored; use `AllCells()` to say so." maxlog=1
    return _addspecies!(eco, op.abundance, tolerance = op.tolerance,
                        demand = op.demand, dispersal = op.dispersal,
                        birth = op.birth, death = op.death, name = op.name)
end

# The abundances this process owns, and the global species index its first row is.
# **At the point interventions run the landscape is row-partitioned**, so a rank holds whole
# *species* across **all** cells (`synchronise_from_rows!` has already run - `ext/EcoSISTEMMPIExt/dynamics.jl`). That
# is what makes every abundance operation rank-local: a rank clears or adds to its own species at the
# named cells, every rank does the same, and between them they cover the whole assemblage. No global
# draw and no partition arithmetic is needed, which is what an earlier reading of the plan assumed.
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

_speciesindex(eco::AbstractEcosystem, sp) = _speciesindex(eco.spplist, sp)
