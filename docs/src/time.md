```@meta
CurrentModule = EcoSISTEM
```

# Time in EcoSISTEM

Two clocks meet in an EcoSISTEM simulation. One is the **run's** own: elapsed simulation
time, advancing one timestep at a time. The other belongs to the **data**: a climate layer
knows which month each of its slices describes, and over what interval each value
accumulated. This page covers both, and how they are lined up against each other — which is
what a run's [epoch](@ref "The epoch: giving the run a date") does.

For what a layer *means* — condition or resource, and to which species — see
[Layers, conditions and resources](@ref).

## The simulation clock

An [`Ecosystem`](@ref) carries an elapsed time that starts at zero when it is built and
advances by one timestep on each [`update!`](@ref). Read it with `simulationtime`, and
return it to zero with [`resettime!`](@ref):

```@setup clock
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

area = StudyArea(extent = (100km, 100km), cellsize = 10km, verbosity = :silent)
env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                        supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                             axis = SolarRadiation),
                        area = area)
species = build_species(5, tolerance = (285.0K, 3.0K), toleranceaxis = Temperature, demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                        abundance = 1000, seed = 1)
eco = build_ecosystem(species, env)
```

```@example clock
simulate!(eco, 1year, 1day)           # run for a year in daily steps
EcoSISTEM.simulationtime(eco)        # a year, as seconds
```

```@example clock
resettime!(eco)                       # back to the start
EcoSISTEM.simulationtime(eco)
```

**Every demographic parameter is a rate**, so the timestep is a free choice rather than a
property of the model. A birth rate of `0.6/year` contributes `0.6/year * timestep` to each
step's draw, and halving the timestep halves the per-step probability rather than changing
what a year means. Choose the timestep for numerical reasons — how finely you want events
resolved, and how much detail your environmental data can actually support — not because the
model expects a particular one.

To do something periodically while a simulation runs, use [`simulate_action!`](@ref), which
calls a function of yours at a regular `interval`:

```@example clock
totals = zeros(Int, length((0year):(1year):(10year)))
simulate_action!(eco, 10year, 1year, 1month_mean_duration) do counting
    totals[counting] = sum(eco.abundances.matrix)
end
totals
```

The `interval` must be a whole multiple of the `timestep`, so that the action always lands on
a step boundary — and this is the first place the next section's distinction bites. A year is
**not** a whole number of days: `year` is Unitful's Julian year of 365.25 days, so an interval
of `1year` with a timestep of `1day` is rejected. It divides exactly by
`month_mean_duration`, which is a twelfth of that same year, so monthly steps are used above.

[`simulate_record!`](@ref) is this same engine with recording already written for you, and is
what most runs want.

## Naming time

A `Unitful` unit is by definition a fixed quantity. A calendar month is not: February is 28
days, or 29, and July is 31 whatever else is true. Rather than pretend otherwise,
`EcoSISTEM.Units` gives each named period its own duration, so that a piece of code says
which month it means:

| unit(s) | duration |
| --- | --- |
| `january_duration`, `march_duration`, `may_duration`, `july_duration`, `august_duration`, `october_duration`, `december_duration` | 31 days |
| `april_duration`, `june_duration`, `september_duration`, `november_duration` | 30 days |
| `february_common_year_duration` / `february_leap_year_duration` / `february_mean_duration` | 28 / 29 / 28.25 days |
| `year` (Unitful's Julian year, re-exported for reference) | 365.25 days |
| `month_mean_duration` | 30.4375 days = `year` ÷ 12 |
| `quarter_mean_duration` | 91.3125 days = `year` ÷ 4 |

The eleven exact months plus February's mean sum to exactly one `year`, so the scheme is
self-consistent by construction. `EcoSISTEM.Units.month_duration(n)` returns the duration of
month `n` — take a month by number rather than by assembling its name.

The last two rows are **declared approximations**, and are honest about it. Some layers
describe a month or a quarter that cannot be identified: `bio13` is precipitation of the
*wettest* month, and which month that is varies from cell to cell. There is no single
calendar month to divide by, so a mean month is used and the name says so. Where the month
*is* known — the twelve slices of a monthly climate stack — each slice uses its own real
duration instead.

`January` through `December` are also available, as the month **numbers** `1` to `12`
(re-exported from `Dates`), for naming a slice of a monthly layer. They are ordinals, not
durations: `january_duration` is a length of time and will not index anything.

## Layers that change over time

A layer that varies in time is not a special type of layer. It is an ordinary layer holding
the values current now, plus a **change** that says how those values depend on elapsed time.
Declare one when the environment is built, by wrapping the layer's spec in [`Varying`](@ref):

```@example clock
warming = GridHabitat(regime = Varying(UniformSpec(285.0K, axis = Temperature),
                                             IncrementBy(0.02K / year)),
                            supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                                 axis = SolarRadiation),
                            area = area)
warming.regime.change
```

The same wrapper takes a data-backed spec — `Varying(SourceSpec(WorldClim{BioClim}, :bio1),
IncrementBy(0.02K / year))` warms a real temperature layer the same way.

A change is written as a **shape** inside a **recipe**. The shape says what the values are;
the recipe says how to read them:

| recipe | the shape's values are | accepts |
| --- | --- | --- |
| [`ReplaceWith`](@ref) | the layer's value itself | a series, or a pattern giving full values |
| [`OffsetBy`](@ref) | an offset from the values as they stood | a series or a pattern |
| [`IncrementBy`](@ref) | a rate, accumulated each step | a series, a pattern, or a constant |

The shapes are a bare constant (a steady drift, only with `IncrementBy`),
[`PatternedChange`](@ref) for an arbitrary function of elapsed time — one sinusoid per
`timescale` by default, but a ramp, step or sigmoid works equally — and
[`SeriesChange`](@ref) for a stack of stored slices such as a read climate series. Add
changes with `+` to drive one layer with several at once, such as a monthly cycle plus a
multi-year warming trend.

Only `IncrementBy` accepts a bare constant, and the reason is worth knowing: writing the same
absolute value or the same offset every step is idempotent after the first step. That is a
one-off adjustment to the ecosystem, not something changing over time, so it is rejected
here rather than silently doing nothing.

### How a series is indexed

A [`SeriesChange`](@ref) is indexed by **elapsed time**, never by a step counter. Each stored
slice has a time coordinate, and the rule is:

> A slice is current from its own coordinate until the next slice's coordinate.

The time looked up is `origin + elapsed`, where `origin` is the coordinate that elapsed time
zero corresponds to — by default the first stored slice, so a series starts at its own
beginning whatever its axis is anchored to.

Indexing by time rather than by step is what makes a series independent of the timestep: one
twelve-month step and twelve one-month steps land on the same slice, and a daily step
through a monthly series holds each slice for the whole of its own month. A cursor advanced
once per call could not do this.

`atend` decides what happens once elapsed time runs past the last slice:
[`ErrorAtEnd`](@ref) (the default) says so plainly, [`HoldAtEnd`](@ref) keeps the last slice
for the rest of the run, [`RepeatAtEnd`](@ref) cycles by a true modulus of the series' own
period — the right choice for driving a long run from a twelve-month climatology — and
[`RevertToLayer`](@ref) hands the layer back, so it returns to the values it had before the
series was attached.

There is no matching `atstart`, because before its first slice a series has only one
sensible reading: **it has not started, so it says nothing and the layer stands.** That
happens automatically whenever a run begins before its series does — give a series an
`origin` earlier than its first slice, or a run an epoch earlier than a dated series' start,
and the layer keeps its own values until the data begins. `RevertToLayer` is the same rule
made available at the far end, where the alternatives are genuine.

Both ends read the same way: **outside its own span a series contributes nothing**, and
the layer is whatever it would be without it. Each mode expresses that in its own terms —
`ReplaceWith` restores the layer's values, `OffsetBy` offsets by zero, and `IncrementBy`
accumulates nothing at all.

`simulate!` checks this before the first timestep rather than discovering it part-way
through, so a run that will outlast its data says so immediately. `ErrorAtEnd` refuses;
`HoldAtEnd` warns, and quantifies it — a series holding its last slice for 80% of a run is
rarely what was meant.

### What a series' coordinates mean

A series also records what its time coordinates *are*, as an
[`AbstractSeriesCalendar`](@ref). This cannot be read off the values: coordinates of one,
two and three months are equally "the first three months of the year" and "one, two and
three months into my experiment".

| calendar | the slices are | elapsed zero is |
|---|---|---|
| [`DatedSeries`](@ref) | real dates | the slice covering the epoch |
| [`MonthOfYearSeries`](@ref) | months of the year, year unknown | the slice for the epoch's own month |
| [`UndatedSeries`](@ref) | plain offsets | the first slice, or wherever `origin` says |

A source whose lookup holds real dates is a `DatedSeries` and anything else an
`UndatedSeries`, both worked out for you. A climatology is the case you have to declare —
`calendar = MonthOfYearSeries()` — because nothing about twelve monthly grids says whether
they are the twelve months of a year or twelve consecutive samples.

`origin` is accepted **only** for an `UndatedSeries`, the one case where nothing else can
set the zero point. For the other two the epoch fixes the phase, and a second knob for the
same thing could only contradict it. To say that an otherwise undated series really begins
at a particular date, give it real `times` rather than an `origin`: that states what every
slice is, not merely where zero sits.

## The epoch: giving the run a date

Elapsed time on its own says how *long* a run is, not *when* it happens. The epoch is the
date elapsed zero corresponds to, and it settles both questions that need one: what date to
report, and — more usefully — **which slice of a seasonal series the run begins on**.

Without an epoch a monthly climatology starts at slice one, so every run starts in January.
That is fine if January is what you meant:

```@setup epoch
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using DimensionalData
using Dates: DateTime      # not plain `using Dates`: it would shadow Unitful's `day`
using DimensionalData.Lookups: NoLookup

area = StudyArea(extent = (50km, 50km), cellsize = 10km, verbosity = :silent)
# A twelve-month temperature climatology whose slice k holds 279 + k K, so the value names
# the month it came from.
slices = cat((fill((279.0 + i) * K, 5, 5) for i in 1:12)...; dims = 3)
stack = DimArray(slices, (Y(NoLookup()), X(NoLookup()), Ti((1:12) .* month_mean_duration)))
seasonal() = Varying(UniformSpec(280.0K, axis = Temperature),
                     ReplaceWith(SeriesChange(stack, calendar = MonthOfYearSeries(),
                                              atend = RepeatAtEnd())))
makeenv() = GridHabitat(regime = seasonal(),
                              supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                                   axis = SolarRadiation),
                              area = area)
species = build_species(5, tolerance = (285.0K, 3.0K), toleranceaxis = Temperature, demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                        abundance = 1000, seed = 1)
```

```@example epoch
eco = build_ecosystem(species, makeenv())
eco.habitat.regime.matrix[1, 1]        # 280 K: the January slice
```

Give the run a July epoch and it starts in July instead:

```@example epoch
eco = build_ecosystem(species, makeenv(), epoch = DateTime(2015, 7, 1))
eco.habitat.regime.matrix[1, 1]        # 286 K: the July slice
```

Note that this is the value **at construction**, before any timestep has run:
`build_ecosystem` writes each series into its layer, so the very first step's births, deaths
and movement see the environment the series and epoch describe rather than reaching it a
step late. A *rate* is the exception, and deliberately so — `IncrementBy` accumulates, and
at the start of a run nothing has accumulated yet.

The run then has a calendar, so it can report a date as well as an elapsed time:

```@example epoch
simulate!(eco, 90day, 1day)
EcoSISTEM.simulationdate(eco)
```

### Where the epoch comes from

An explicit `epoch` always wins. Left out, it is resolved from the environment's own series
the same way [`StudyArea`](@ref) resolves a CRS — adopt it if it is unambiguous, ask if it
is not:

- **exactly one** dated series → its start date becomes the epoch, and every other series is
  phased to it;
- **several that disagree** → an error naming the candidates, because no default could be
  right;
- **none at all** → no epoch. `simulationdate` is `nothing`, and the run behaves exactly as
  one that never mentions dates.

An epoch *before* a dated series begins is an error rather than a clamp: `atend` says what
to do past a series' end, but there are no values before its beginning to hold or cycle. A
`MonthOfYearSeries` has no beginning to precede, so any date phases it.

## When the data accumulated

The catalogue answers two independent questions about a layer's time, in two separate
columns, and conflating them is the classic way to get a rate wrong:

- **Temporal resolution** — how often the layer is *sampled*. A monthly climate variable has
  twelve slices; an annual one has a single grid.
- **Accumulation period** — what interval each value was *measured over*. A month's rainfall
  total accumulated over that month; a degree-day sum accumulated over a year.

These come apart in both directions. Solar radiation is sampled monthly but is already
reported as a daily flux, so it accumulates over a **day** while being sampled by month.
`bio13` is emphatically monthly in character yet arrives as one grid, so it has no sampling
resolution at all — only an accumulation period.

An accumulation period takes one of three forms, and the difference matters to the modeller,
not just to the reader:

| form | example | what it means |
| --- | --- | --- |
| constant | `bio12`, over a year | a fixed scale factor; the layer and any tolerance convert by the same number |
| per slice | `prec`, over each slice's own calendar month | each of the twelve slices is divided by its own month's real length |
| per cell | `gsp`, over the growing season (`gsl`) | the divisor varies from cell to cell |

A **constant** period is a matter of presentation: dividing by it is a single multiplication,
so a stock and a rate carry the same information and nothing hinges on which you look at.

A **per-cell** period is not. "This species needs 500 mm over the growing season" and "this
species needs 5 mm a day during the growing season" are *different biological claims*: the
first is 5 mm/day in a hundred-day season and 2.5 mm/day in a two-hundred-day one. Which one
you mean has to be decided, because the data cannot decide it for you — which is why
[`GrowingSeasonPrecipitation`](@ref) reads as a total when used as a condition and as a daily
rate when used as a resource.

Two functions report the two readings, and it is worth knowing which you are looking at:
[`layerunit`](@ref EcoSISTEM.layerunit) gives the amount the table declares,
while [`layerrate`](@ref EcoSISTEM.layerrate) gives the unit a read actually
yields. Monthly precipitation is catalogued as `L m⁻²` and read as `L m⁻² d⁻¹`.

Note that having an accumulation period does not by itself make a layer a rate. The axis
decides that. Degree-days accumulate over a year, but the heat *sum* is the meaningful
reading — dividing it would give a mean daily temperature excess that nobody asked for — so
[`CumulativeHeat`](@ref) keeps its accumulated form. The period says what interval a value
covers; the axis says which reading means something.

## Layers store rates, not per-step amounts

One invariant ties this page together:

> No layer stores a value that depends on the timestep.

Layers store rates, and a per-step quantity is computed as `rate × timestep` where it is
used. This is why a monthly total is turned into a daily rate as it is read, rather than
into "the amount per simulated month".

The alternative fails immediately. Suppose a layer stored February's rainfall as an amount
per step: it would be right for one timestep and wrong for every other, and the timestep is
not known until [`simulate!`](@ref) is called — long after the environment was built. Storing
the rate makes the layer independent of how the run is stepped, and leaves the arithmetic to
the one place that knows.

## Worked examples

[`examples/VaryingClimate.jl`][varying] runs the three shapes side by side on one landscape and one
species pool, so the only thing differing between the runs is how the environment moved: a **static**
control, a **steady drift** (`Varying(spec, IncrementBy(rate))`) and a **seasonal cycle**
(`Varying(spec, OffsetBy(PatternedChange(...)))`).

It also asserts what separates them, and the assertion is not the obvious one. A cycle is **not**
distinguished by "ending where it started" — that holds only at whole periods, and
[`simulate!`](@ref) takes `length((0s):timestep:duration)` steps, an *inclusive* range from zero, so
a four-year run at monthly steps advances **49** months rather than 48. What actually separates them
is that a cycle stays bounded by its amplitude however long it runs, while `IncrementBy` grows
without limit.

[`notebooks/Introduction.jl`][notebook] builds a warming gradient and a warming peak the same way,
with plots, and is the gentler introduction.

A layer change is a pure function of elapsed time. For change that acts on the **ecosystem** —
abundances, the active mask — see [Interventions](@ref), which is a deliberately separate mechanism
and explains why.

[varying]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/VaryingClimate.jl
[notebook]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/notebooks/Introduction.jl
