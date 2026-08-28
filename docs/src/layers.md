```@meta
CurrentModule = EcoSISTEM
```

# Layers, conditions and resources

An EcoSISTEM environment is assembled from **layers** — gridded quantities that species
respond to. This page is about what a layer *means*: which quantities are conditions a
species is matched against, which are resources it competes for, and how to find out which
of those a given dataset can be.

For *when* a layer's values apply, and over what interval they were measured, see
[Time in EcoSISTEM](@ref).

## What a layer is

Every layer carries two type parameters that together say what it represents:

- a **role** — `Condition` or `Resource`. A condition is *matched* to a species' tolerance;
  a resource is *consumed* against a species' demand.
- a **niche axis** — an [`EcoSISTEM.NicheAxis`](@ref) naming the physical quantity measured:
  [`Temperature`](@ref), [`Precipitation`](@ref), [`SolarRadiation`](@ref),
  [`CarbonFlux`](@ref), [`LandCoverTypology`](@ref) and so on. Axes are grouped under
  abstract supertypes (`EcoSISTEM.TemperatureAxis`, `EcoSISTEM.WaterAxis`, …), so a method can
  apply to one leaf or to a whole family. The leaves are exported and the supertypes are
  `public` rather than exported, so a supertype needs qualifying (or naming in your `using`):
  that is the convention for every abstract type here, since it is the leaf you name day to day.

There are two concrete layer types, and a container for several of them:

| type | holds |
| --- | --- |
| [`ContinuousLayer`](@ref) | numeric values — a temperature, a rainfall rate, a supply |
| [`CategoricalLayer`](@ref) | class codes — land-cover or climate classes; always a `Condition` |
| [`LayerCollection`](@ref) | several layers used together, positionally or by name |

Every layer is a two-dimensional grid over `(Y, X)`, with real geographic or projected
coordinates when it came from data and a plain index grid when it was generated
synthetically. The role-specific aliases are the names you will normally write:

| alias | role | value |
| --- | --- | --- |
| [`ContinuousRegime`](@ref) | condition | any numeric quantity |
| [`CategoricalRegime`](@ref) | condition | class codes |
| [`Supply`](@ref)`{SolarRadiation}` | resource | `kJ/day` per cell |
| [`Supply`](@ref)`{Precipitation}` | resource | `L/day` per cell |
| [`Supply`](@ref)`{CarbonFlux}` | resource | `g/day` per cell |

A layer holds **one** grid of values: the ones current now. A layer that varies in time is
not a different type — it is one of these carrying a change rule that decides, from the
simulation clock, what "now" means. See
[Layers that change over time](@ref) for the vocabulary.

Layers reach an environment through [`GridHabitat`](@ref), which takes a condition
layer as `regime`, a resource layer as `supply`, and the [`StudyArea`](@ref) that decided
the grid:

```@example layers
using EcoSISTEM, EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols

area = StudyArea(extent = (500km, 500km), cellsize = 10km, verbosity = :silent)
env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                        supply = UniformSpec(1.0kJ / (m^2 * day),
                                             axis = SolarRadiation),
                        area = area)
```

Pass a tuple to either keyword for a multi-variable environment — two conditions, or two
resources, or both — and the pieces are held in a [`LayerCollection`](@ref). Naming them
makes them addressable:

```@example layers
env = GridHabitat(regime = (temperature = UniformSpec(285.0K, axis = Temperature),
                                  rain = UniformSpec(2.0mm / day,
                                                     axis = Precipitation)),
                        supply = (sunlight = UniformSpec(1.0kJ / (m^2 * day),
                                                         axis = SolarRadiation),
                                  water = UniformSpec(2.0mm / day,
                                                      axis = Precipitation)),
                        area = area)
(regimes = keys(env.regime), supplies = keys(env.supply))
```

## Conditions and resources

EcoSISTEM follows the Begon–Townsend–Harper distinction, and the whole vocabulary rests on
it:

|  | condition | resource |
| --- | --- | --- |
| the environment offers | a **regime** | a **supply** |
| the species brings | a **tolerance** | a **demand** |
| they are combined by | a suitability function | a ratio of supply to demand |
| example | temperature, land-cover class | light, water, carbon |
| is it used up? | no | see below — it is *rival*, which is not the same thing |

A condition is a state of the world that suits a species more or less well. A species
carries a [`NicheTolerance`](@ref) — a distribution over the axis — and a
[`NicheSuitability`](@ref) turns the distance between a cell's value and that tolerance into
a multiplier on the species' rates. A resource is something individuals draw on and
therefore compete for: a species carries a demand ([`Demand{SolarRadiation}`](@ref),
[`Demand{Precipitation}`](@ref), [`Demand{CarbonFlux}`](@ref)), and what matters is
how that demand compares to what the cell supplies.

### Deciding which one a quantity is

Modelling a condition as though it were a resource produces no error and no warning — just
wrong numbers — so it is worth being deliberate. A quantity is a resource only if **both**
of these hold:

1. **Each individual has a per-capita demand for it**, so that the demands of everyone
   present can be added up; and
2. **it is rival** — more individuals present means less of it to go round.

Light passes both: a plant intercepts photons another plant then cannot, and shading is
exactly that rivalry made visible. Water passes both. Temperature fails both — an individual
does not have a "demand" for 15 °C that can be summed across a population, and one organism
being warm does not make its neighbour colder. Temperature is the textbook condition, and
belongs in a regime with a tolerance.

Evaporative demand (potential evapotranspiration) is the instructive failure. It looks like
a resource — it has units, it varies over space, it clearly matters to plants — but it is
the atmosphere's demand *for* water, not a pool of anything. Nothing consumes it, and one
plant transpiring does not reduce its neighbour's PET. Modelled as a supply it would enter
as the resource available, so high evaporative demand would *increase* births, when in
reality it means drought stress. A condition modelled as a resource invents competition that
is not there; a *demand* modelled as a resource is backwards.

## Supply and demand

A supply is a **flow**, not a reservoir. Each cell's supply is a rate — energy, water or
carbon per day — and it is recomputed in full every timestep rather than being drawn down.
Nothing is depleted.

It is nonetheless rival, and the rivalry lives in the comparison. Each timestep, every cell
sums the demands of all the individuals in it, and the birth rate is scaled by the ratio of
what the cell supplies to what its occupants collectively ask for. That ratio falls as the
cell fills, which is the density dependence a resource is supposed to produce — it does not
require anything to be consumed.

Because both sides are in the **same** unit — a flow against a flow, or a stock against a
stock — the ratio is a dimensionless count, and independent of the timestep.
This is why supplies are stored as rates rather than as amounts; see
[When the data accumulated](@ref) for what that means for data that arrives as a total.

Where a species draws on several resources at once, they combine by **Liebig's law of the
minimum**: the scarcest resource sets the outcome, and one that is never scarce never binds.
A generous or irrelevant resource therefore costs nothing but is also doing nothing.

### The resource families

**Both** sides are decided by the **axis** you declare, never by units — which cannot tell
`m/s` from `mm/day`:

| family | areal input | per-cell supply | species demand |
| --- | --- | --- | --- |
| solar | `kJ/m²/day` | [`Supply`](@ref)`{SolarRadiation}`, `kJ/day` | `Demand{SolarRadiation}`, `kJ/day` |
| water | `mm/day` (that is, `L/m²/day`) | [`Supply`](@ref)`{Precipitation}`, `L/day` | `Demand{Precipitation}`, `L/day` |
| carbon | `g/m²/day` | [`Supply`](@ref)`{CarbonFlux}`, `g/day` | `Demand{CarbonFlux}`, `g/day` |
| space | a **fraction** of the cell, 0–1 | [`Supply`](@ref)`{SurfaceArea}`, `m²` | `Demand{SurfaceArea}`, `m²` |

An areal rate becomes an absolute per-cell one by multiplying by the cell's area, so a
coarser grid gives each cell more of everything, as it should. On the species side
[`build_species`](@ref) takes the axis explicitly, as `demandaxis`:

```julia
build_species(n, tolerance = (298.0K, 2.0K), toleranceaxis = Temperature,
                 demand = 10.0kJ / day,      demandaxis = SolarRadiation, …)
```

The unit is still **checked** against the axis — a `L/day` demand declared on
`SolarRadiation` is refused — it just no longer *decides*. A **bare number is refused** too: it
carries no unit to check, and there is no free/dimensionless resource left for it to mean.

**Space is the odd one out, deliberately.** It is the only resource that is a *stock* rather
than a flow: a fraction of ground, not a rate of anything. The ratio the model needs stays a
dimensionless count either way (`m² ÷ m²` as much as `kJ/day ÷ kJ/day`), and supplies are
recomputed in full each timestep rather than depleted, so a standing stock needs no change to
the loop. Ask for one with [`SurfaceSpec`](@ref) — `SurfaceSpec()` for the whole cell,
`SurfaceSpec(0.4)` for a partly-available one. The twelve EarthEnv land-cover bands are space
layers too: each is the proportion of a cell covered by one class.

Water uses `mm/day` for its areal form because a millimetre of rain over a square metre is a
litre — a depth per unit time *is* a volume flow per unit area. Solar radiation and carbon
have no such identity, so their areal forms name the area explicitly.

!!! note "Carbon is not independent of light and water"
    A carbon supply is normally built from net primary productivity, which is an *estimate
    of what the climate allows a plant to produce* — from light and water, among other
    things. Using all three as separate resources therefore counts one limitation twice.
    Nothing prevents it, and under Liebig's minimum a redundant resource simply never binds,
    but it should be a deliberate choice rather than an accident of loading three layers.

## Layers the model has no axis for

Not every layer you want on the grid is one the simulation consumes. A survey of somewhere,
a covariate you mean to compare output against — these have no niche axis, and inventing one
for them would be dishonest.

They are still first-class layers. Name [`EcoSISTEM.NicheAxis`](@ref) — the root of the
catalogue — and the layer is carried on the simulation's own grid, cell for cell, while
**accepting any unit at all**.

```julia
odd = UniformSpec(2.0u"kJ^2", axis = EcoSISTEM.NicheAxis)   # claiming nothing, deliberately
habitat = GridHabitat(regime = odd, supply = <a real supply>, area = area)
```

That builds. Nothing objects to `kJ²`, and that is deliberate rather than an oversight:
there is no axis to disagree with the unit. An axis is what gives a value meaning, so a
layer on the root constrains nothing.

`axis` is **required** — there is no default. Naming the root is how you say *"I am
claiming nothing about what this measures"*, and saying it out loud is the point: a layer
that silently defaulted could not be told from one whose author forgot.

Three rules follow, and they are the whole of it:

  - **It can be a regime, never a supply.** `EcoSISTEM.NicheAxis` declares no resource, so offering
    one as a supply is refused by name. "Regime or reference, never a supply" is the shape to
    remember.
  - **It is not a wildcard.** The root pairs with the root and with nothing else — matching is
    *identity*, so a species on `Temperature` is refused against a root-axis regime, and told so.
    That is deliberate: a layer that declines to say what it measures must not be silently read as
    saying whatever the species happens to need. If you want the pairing, declare a real axis on
    both sides.
  - **The species side must be in the same frame.** A tolerance built in bare numbers cannot
    match a regime carrying `kJ²`. The failure comes from the arity/unit alignment check at
    `build_ecosystem` — *"layer `tolerance` is Float64 in the species tolerance but
    Quantity{…} in the environment regime"* — which names neither `EcoSISTEM.NicheAxis` nor the unit
    as the cause, so it is worth recognising.

!!! note "Matching a united axis-less layer"
    `build_species` reads bare tolerance parameters in the axis's own frame, which for
    `EcoSISTEM.NicheAxis` is bare numbers — so a regime carrying a real unit needs a tolerance built
    in that unit. Construct it directly and pass it in:

    ```julia
    tol = NicheTolerance(EcoSISTEM.NicheAxis, Normal, params, support = u"kJ^2")
    species = build_species(n, tolerance = tol, demand = …, demandaxis = …)
    ```

    A **pre-built tolerance is used as given**, carrying its own axis and frame, so
    `toleranceaxis` is not needed alongside it. It must already cover every species.

**The useful case is composition.** A layer with no axis is an *ingredient*: a soil-type
map has no tolerance and no supply of its own, yet it legitimately changes how much of the
rain is available. Feed it to a [`ConstructedSpec`](@ref) together with a precipitation layer
and declare the **result** on `Precipitation`, and *that* is what carries a tolerance or
becomes a supply. The axis is a property of the dish, not of the ingredients.

### What a combine is handed, and what it must give back

A [`ConstructedSpec`](@ref)'s `combine` receives one raster per layer and **must return a
raster**. That is the whole contract, and it does not depend on how the spec is later used:
a *mask* is simply a raster whose element type is `Bool`.

Nothing has to be unwrapped to satisfy it. A raster behaves like the values it holds and stays
a raster, so the natural way to write a combine is also the correct one:

```julia
# a mask — Bool-valued, still a raster; it claims nothing, so it names the root axis
ConstructedSpec(EarthEnv{LandCover}, axis = EcoSISTEM.NicheAxis) do lc
    compress_landcover(lc) .!= landcoverclass(:open_water)
end

# a derived layer — several bands added together, its meaning declared by `axis`
ConstructedSpec(EarthEnv{LandCover}, [:shrubs, :herbaceous], axis = SurfaceArea) do bands...
    sum(bands)
end
```

**A combine never names an array type**, in or out. That is deliberate: the array a raster
holds is an implementation detail, and a combine is *your* code — it should not have to be
rewritten because ours changed.

A derived raster drops what it can no longer claim. It keeps no layer code and no dataset:
its source becomes `DerivedData`, recording what it was computed *from* without claiming to *be*
it, because a combine is free to change what the values are — summing eight land-cover bands gives
a quantity that is none of them, and multiplying them by an incident flux gives solar radiation.
Its meaning comes from the spec's `axis`, exactly as for any other layer.

## Data you already hold

Everything above names a *source* and lets EcoSISTEM read it. Sometimes you have the data already
— computed elsewhere, or read by hand — and there is a pathway for that, though it is not the one
to reach for first:

```julia

regime = in_memory_raster(my_raster, axis = Temperature)
```

**Prefer naming the source where you can.** A [`SourceSpec`](@ref) lets EcoSISTEM read only the
window your study area needs, cache it between layers, and take the unit, niche axis, accumulation
period and value type from the catalogue. An in-memory raster gives all of that up and describes
itself only by the `axis` you pass.

**And `axis` is not optional in spirit, even though it has a default.** A raster carries values
and possibly a layer code, but no niche axis — nothing about it says whether those numbers are a
temperature, a rainfall rate or a cover fraction. That is why a bare raster is refused as a regime
or a supply, and why this function exists: it is the place the declaration goes.

## Investigating the catalogue

The shipped catalogue describes every layer EcoSISTEM can read: its unit, its niche axis,
what kind of quantity it is, and how it is laid out in time.
[`layerinfo`](@ref EcoSISTEM.layerinfo) is the way in.

```@example catalogue
using EcoSISTEM, RasterDataSources

EcoSISTEM.layerinfo(WorldClim{Climate}, :prec)
```

(The catalogue query functions are public but not exported, so they are written
`EcoSISTEM.layerinfo` rather than bare — `layerunit` and `layeraxis` *are* exported and need
no prefix.)

Called with a code alone it searches every shipped table, since the same code can appear in
more than one dataset. To go the other way — from an axis to the layers on it — use
[`layersbyaxis`](@ref EcoSISTEM.layersbyaxis), which accepts a leaf axis, an
abstract group, or nothing at all:

```@example catalogue
(precipitation = length(EcoSISTEM.layersbyaxis(Precipitation)),  # just that axis
 water = length(EcoSISTEM.layersbyaxis(EcoSISTEM.WaterAxis)),    # all its leaves too
 everything = length(EcoSISTEM.layersbyaxis()),                  # the whole catalogue
 unclassified = length(EcoSISTEM.layersbyaxis(nothing)))         # no axis declared
```

None of this reads a raster: the catalogue is shipped with the package, so a layer can be
looked up before any data is downloaded.

[`layeraxes`](@ref EcoSISTEM.layeraxes) returns the axis hierarchy itself, which
is the quickest way to see what axes exist before drilling into one. And two functions
answer two different questions about units:
[`layerunit`](@ref EcoSISTEM.layerunit) reports what the table declares (an
amount, `L m⁻²`), while [`layerrate`](@ref EcoSISTEM.layerrate) reports what
reading the layer actually yields (a rate, `L m⁻² d⁻¹`). The difference between them is the
subject of [When the data accumulated](@ref).

The columns worth knowing:

| column | says |
| --- | --- |
| `unit` | the physical unit of the stored amount |
| `axis` | the [`EcoSISTEM.NicheAxis`](@ref) the layer is modelled on |
| `category` | what kind of quantity it is — `rate`, `stock`, `balance`, `count`, `range`, `instantaneous`, `categorical` |
| `valuetype` | whether values are continuous, discrete or categorical class codes — catalogue metadata only; what the *package* acts on is the layer's axis, which says so with `categorical = true` in its [`@nicheaxis`](@ref) declaration |
| `temporal` / `numslices` | how often it is sampled, and how many slices a read returns |
| `period` | what interval a value accumulated over — see [Time in EcoSISTEM](@ref) |
| `sources` | which datasets provide it |

### Can this layer be a supply?

The answer is **not** the `category` column. A layer can be catalogued `rate` and still be a
condition: degree-day sums and evaporative demand are both rates, and neither is a resource.

What decides it is whether the axis **declares a supply type** — a statement, in code, that
species compete for this. Four axes do:

| axis | as a resource |
| --- | --- |
| [`Precipitation`](@ref) | [`Supply`](@ref)`{Precipitation}`, `L/day` |
| [`SolarRadiation`](@ref) | [`Supply`](@ref)`{SolarRadiation}`, `kJ/day` |
| [`CarbonFlux`](@ref) | [`Supply`](@ref)`{CarbonFlux}`, `g/day` |
| [`SurfaceArea`](@ref) | [`Supply`](@ref)`{SurfaceArea}`, `m²` |

Everything else is a condition, and asking for it as a supply reports that clearly rather
than guessing a resource type. [`CumulativeHeat`](@ref) and [`Evapotranspiration`](@ref) are
the cases the two-part test above rules out — both are rates, and neither is consumed.

**`SurfaceArea` is the one that shows the `category` column really cannot decide this.** It is
catalogued `instantaneous` — the twelve land-cover bands are a snapshot of what covers a cell, not
a rate of anything — and it is nonetheless a resource, because *room* is exactly what species
compete for. An intensive physical **state** times an area is meaningless (`kPa × m²` is a force);
an areal **fraction** times an area is an area, and that is a resource.

[`GrowingSeasonPrecipitation`](@ref) is a fourth route to a supply and is **not** a fourth
resource axis: it declares no supply type. Because its accumulation period is another layer
(`gsl`, growing-season length) that varies by cell, a `gsp` spec used as a supply is rewritten
onto [`Precipitation`](@ref) and divided by that layer first — on the *source* grid, since
dividing is cell-wise but nonlinear and so does not commute with regridding. The result is an
ordinary water supply.

For the full picture of how these declarations relate to the catalogue's columns, see
[Axes, units and roles: how a layer is classified](@ref).

[`SnowWaterEquivalent`](@ref) is the one that looks like it should be on the list and is
not. A snowpack certainly contains water — but how much of it is *available*, and when,
depends on melt timing, which this model does not represent. A snow-water supply would
therefore be an amount with no defensible rate attached to it. Precipitation, growing-season
precipitation and productivity are the water and carbon inputs that do have one.

## Worked examples

Each of the distinctions on this page is demonstrated on its own, at a size that runs in seconds:

| | |
|---|---|
| **A categorical layer**, and why the *axis* rather than the values decides it — with `SimpleCategoricalTolerance` matching a set of classes and `CategoricalSuitability` inferred from it | [`examples/CategoricalLandCover.jl`][categorical] |
| **Composing layers into a new axis** — the same land-cover bands read once as available ground and once as usable sunlight, where only the declared `axis` says which is which | [`examples/AvailableGround.jl`][ground] |
| **A real raster and a real shapefile**, built entirely from lazy specs with no array passed by hand | [`examples/ScottishCultivatedLand.jl`][scot] |
| **The same layer in both roles** — CHELSA growing-season precipitation as a *condition* and as a *resource*, meaning something different each time | [`examples/other/gsp.jl`][gsp] |
| **A purely synthetic environment**, needing no data at all | [`examples/SimulatedEcosystem.jl`][sim] |

For layers that change as the run proceeds, see [Time in EcoSISTEM](@ref) and
[`examples/VaryingClimate.jl`][varying].

[categorical]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/CategoricalLandCover.jl
[ground]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/AvailableGround.jl
[scot]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/ScottishCultivatedLand.jl
[gsp]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/other/gsp.jl
[sim]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/SimulatedEcosystem.jl
[varying]: https://github.com/EcoJulia/EcoSISTEM.jl/blob/main/examples/VaryingClimate.jl
