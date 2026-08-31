```@meta
CurrentModule = EcoSISTEM
```

# Axes, units and roles: how a layer is classified

[Layers, conditions and resources](@ref) explains *what* the condition/resource distinction
means ecologically. This page is the mechanical companion: how EcoSISTEM actually decides,
for a given dataset, what unit its values end up in and whether it can be a regime, a supply,
both or neither - and which of those answers come from the shipped catalogue, which are
declared in code, and which are still open questions.

Read it if you are adding a layer, adding a niche axis, or trying to work out why a layer
will not build.

## The chain, end to end

A value travels from a CSV cell to a simulation through four stages. Each answers a
different question, and conflating them is the source of most of the surprises below.

| stage | question | answered by |
| --- | --- | --- |
| 1. stored | what unit does this *layer, from this source* hold? | the catalogue's `Units` column -> `EcoSISTEM.layerunit` |
| 2. read | what does reading it actually yield? | `Units` + `AccumulationPeriod` -> [`layerrate`](@ref EcoSISTEM.layerrate) |
| 3. canonical | what unit does this *quantity* mean, across every source? | `EcoSISTEM.canonicalunit(axis)` - declared in code |
| 4. role | is it matched against, or consumed? | `EcoSISTEM.supplytype` / `EcoSISTEM.demandtype` - declared in code |

Stages 1 and 2 are data. Stages 3 and 4 are **modelling decisions** that the data constrains
but cannot make.

## Stage 1-2: what the catalogue can tell you

The catalogue is shipped with the package, so all of this works before any raster is
downloaded.

```@example units
using EcoSISTEM, RasterDataSources

rec = EcoSISTEM.layerinfo(WorldClim{Climate}, :prec)
(unit = rec.unit, category = rec.category, period = rec.period, axis = rec.axis)
```

Two columns do the work here, and they answer different questions:

- **`AccumulationPeriod`** - the interval a value accumulated over. This is what turns a
  total into a rate.
- **`Temporal Resolution`** - how often the layer is *sampled*. Nothing to do with units.

They used to be one column, which was invisibly wrong for `prec` (both a month) and plainly
wrong for `srad` (sampled monthly, but already a per-day flux). See
[When the data accumulated](@ref).

So the read unit is not always the stored unit:

```@example units
(stored = rec.unit,
 read   = EcoSISTEM.layerrate(rec.unit, rec.period, rec.axis))
```

!!! warning "The division is decided by the axis, not by the period"
    Having an accumulation period does **not** mean a layer is read as a rate. The period says
    *what interval this accumulated over*; the axis says *which reading is canonical*.

    A consequence worth knowing: an axis that declares **no** canonical unit is never
    divided, so it is read as the accumulated amount. Two layers with identical catalogue
    entries can therefore be read in different units, decided by a declaration elsewhere.

### So when *is* an accumulated layer divided?

The test is not the `Category` column, and it is not whether the layer has a period. It is:

> **Is the accumulated quantity a meaningful mean when divided by its period?**

| layer | category | divided? | why |
| --- | --- | --- | --- |
| `prec` | rate | yes | a mean daily rainfall is a real quantity |
| `pet` | rate | yes | evaporative demand happens continuously |
| `cmi` | balance | yes | a mean daily surplus or deficit is as real as a monthly one |
| `gdd0` | - | no | the **sum** is the quantity - a species is matched against accumulated warmth, and dividing gives a mean daily excess nobody asked for |
| `swb` | balance | no | a **capped** cumulative: it stops accumulating when the soil is full, so the mean is not a rate of anything |
| `gsp` | stock | **depends on the role** | see below - it is the only row where the answer is not the same for a regime and a supply |

!!! note "`Category` and "has a time in the unit" are independent"
    `Category` classifies the *semantics* - `rate` accumulates and is non-negative, `balance`
    accumulates and is signed, `range` is a spread. Any of them may or may not carry a time
    denominator. The climate moisture index is `prec - pet`, and both of those are per month,
    so `cmi` is a **balance whose dimension is a rate**. Dividing it by its month does not
    stop it being a balance - it is still signed, and still never a supply.

    What the "no" rows have in common is not their category or their sign: it is that
    dividing would produce something meaningless. A capped total is not a flow, and a heat
    sum *is* the thing being measured.

!!! note "`gsp` is the one where the role decides"
    Growing-season precipitation accumulates over `gsl`, the growing-season length, which
    **varies by cell**. So there is no fixed denominator, and `layerrate` deliberately leaves
    it as the stored amount - but that is not the end of the story:

    - as a **regime**, it stays `L m^-2`. "How much water fell over the season" is what a
      species is matched against, and that is the quantity you want;
    - as a **supply**, it *is* divided - by the `gsl` layer, cell by cell - and the result is
      an ordinary water supply in `L/day`.

    So a stock reading and a rate reading of `gsp` are, in the words of
    [`EcoSISTEM.PerCellAccumulationPeriod`](@ref), *"genuinely different hypotheses rather than two
    spellings of one"* - which is exactly why the division is deferred until the role is known
    rather than being settled when the layer is read.

    And it must happen on the **source** grid: dividing is cell-wise but nonlinear, so it
    does not commute with regridding. Native cells of (100 mm, 50 d) and (100 mm, 100 d) give
    1.5 mm/d if divided early and 1.33 mm/d if divided late.

## Stage 3: canonical units

`EcoSISTEM.canonicalunit(axis)` is the unit every layer on an axis is normalised to when it
is materialised. It exists because the catalogue *cannot* answer the question:

- **One axis spans rows that legitimately disagree.** WorldClim's `srad` is `kJ m^-2`,
  CHELSA's `rsds` is `MJ m^-2`. Both cells are right; without a declared scale the two are
  never reconciled.
- **Synthetic layers have no row at all.** `UniformSpec(285.0K, axis = Temperature)` never
  touches a catalogue but must be comparable with a layer that did.
- **Tolerances have no row either.** A species' niche is not a layer, and
  [`NicheTolerance`](@ref) takes its support unit from the axis.
- **The canonical reading is sometimes deliberately not the stored one** - see the warning
  above.

It has exactly two consumers, and they are the two halves of one comparison: the layer's
values are converted into it, and a species' tolerance is *built* in it. That is what makes
`pdf(tolerance, value)` meaningful rather than an accident of scale.

### Condition and resource units differ, by exactly the cell area

An axis can have two canonical units, and they have different dimensions:

```@example units
using Unitful
ax = Precipitation
(condition = EcoSISTEM.canonicalunit(ax),
 resource  = EcoSISTEM.canonicalunit(EcoSISTEM.Resource, ax))
```

A **condition is intensive** - how wet a place is, per unit area, independent of how you drew
the grid. A **supply is extensive** - an actual pool in an actual cell, which is what species
compete over. So:

> **Resource unit = areal flux density × cell area.**

That conversion is `cancel`: `kJ/m^2/day × m^2 -> kJ/day`, `L/m^2/day × m^2 -> L/day`,
`g/m^2/day × m^2 -> g/day`. Draw bigger cells and the condition is unchanged while the supply
rises - which is why a coarser grid gives each cell more of everything, as it should.

### Which grids get a cell-area correction, and which do not

Three things look alike here and only one of them is corrected, so it is worth separating them
before the question arises.

**The rule is a single line, and nothing about the *layer* takes part in it:**

> A cell-area correction depends only on the **units of the study grid's coordinates**.
> Coordinates in **metres** (a projected grid) give one constant area; coordinates in
> **degrees** (a geographic grid) give one area per latitude.

So a `SurfaceSpec()` on a British National Grid study area is exactly `cellsize^2` in every cell,
and so is the cell area a solar or water supply is multiplied by. They are treated identically,
because it is the same grid.

1. **A lat/long grid's cells really do change size.** A 0.1° cell at 60°N covers less ground than
   one at 55°N, because meridians converge - 8-20% across a country-sized extent. That is
   geometry, not data, and it *is* corrected: each cell is multiplied by its own true area.
2. **A projected grid's cells do not.** That is what projecting is *for*: British National Grid
   defines a plane in metres, so a 1 km cell is 1 km × 1 km everywhere on it, by construction.
   Nothing is corrected because nothing varies.
3. **A projected plane is still not a perfect map of the ground.** Transverse Mercator - which
   BNG and UTM both are - is *conformal*: it preserves shape and pays for that by distorting
   area, by about **0.22% across Great Britain**, varying with **easting** (distance from the
   2°W central meridian) rather than with latitude. This is **not** corrected, and the remedy is
   not a factor but a different projection: an **equal-area** CRS (LAEA, Albers) has no such
   distortion by construction.

!!! note "The correction is for what you *inspect*, not for what you simulate"
    A geographic grid **cannot be simulated** - `build_ecosystem` refuses it, because dispersal
    assumes one uniform cell size. So the latitude correction only ever affects values you read
    from a built environment: its `supply` matrices, `materialise(spec, area, role = Resource)`,
    anything plotted. Every grid a simulation actually runs on is projected, and on those the
    correction is exactly 1.

    It is kept because those are supported operations that were wrong by up to 20%, and because
    if the dispersal restriction is ever lifted the area side is already correct.

The remaining approximation on a geographic grid is **dispersal**, not supply: the cell *side*
feeds one lookup table per species, and making it per-cell would mean one per cell. `_cellsize`
says so when it reports that a grid's east-west length varies with latitude.

**Stated as "areal flux density", not "the Condition unit", deliberately.** The two
*coincide* for [`Precipitation`](@ref) and [`SolarRadiation`](@ref), which is why the shorter
phrasing is tempting - but [`CarbonFlux`](@ref) has **no** Condition unit at all, so there is
nothing to multiply, and `EcoSISTEM.NicheAxis` declares a Condition unit (`NoUnits`, a bare index)
that describes a *different quantity* from the flow its Resource unit describes. The rule is
about the areal reading, which is not always the canonical Condition one.

This is also why a **supply's rate of change** is written in the layer's per-cell unit
(`kJ/day/s`) rather than the per-area unit its spec was written in (`kJ/km^2/day/s`): the spec
is intensive, the built layer is extensive, and `cancel` is the boundary between them.

## Stage 4: is it a condition, a resource, both, or neither?

These are **two independent declarations**, not one choice:

```@example units
axes = (Temperature, Precipitation, CarbonFlux, ClimateMoisture)
[(axis = nameof(typeof(a)),
  condition = EcoSISTEM.canonicalunit(a),
  resource = EcoSISTEM.supplytype(a)) for a in axes]
```

All four combinations are meaningful:

| | can be a regime | can be a supply | example |
| --- | --- | --- | --- |
| **condition only** | yes | no | [`Temperature`](@ref) - the common case |
| **both** | yes | yes | [`Precipitation`](@ref), [`SolarRadiation`](@ref) |
| **resource only** | no | yes | [`CarbonFlux`](@ref) |
| **neither** | no | no | a reference layer - see below |

[`CarbonFlux`](@ref) is resource-only **on purpose**: net primary productivity is something
plants consume, not something they are matched against. There is no meaningful "my optimum
NPP is 3 g/m^2/day and more is worse" - more is just more. Declaring no `Condition` unit is
what stops a tolerance being built against it.

### Resource-hood is declared, not derived

Being a resource means a concrete supply type exists for the axis:

```@example units
(EcoSISTEM.supplytype(SolarRadiation), EcoSISTEM.demandtype(SolarRadiation))
```

There are four supply families, chosen by the **dimension** of the per-cell value:

| family | areal input | per-cell supply | species demand |
| --- | --- | --- | --- |
| solar | `kJ/m^2/day` | [`Supply`](@ref)`{SolarRadiation}`, `kJ/day` | [`Demand{SolarRadiation}`](@ref) |
| water | `mm/day` (= `L/m^2/day`) | [`Supply`](@ref)`{Precipitation}`, `L/day` | [`Demand{Precipitation}`](@ref) |
| carbon | `g/m^2/day` | [`Supply`](@ref)`{CarbonFlux}`, `g/day` | [`Demand{CarbonFlux}`](@ref) |

!!! warning "The `category` column does not decide this"
    A layer can be catalogued `rate` and still be a condition - evaporative demand (`pet`) is
    a rate that nothing consumes. And a layer can be catalogued `stock` and still become a
    supply: `gsp` (growing-season precipitation) accumulates over `gsl` (growing-season
    length), which varies per cell, so as a supply it is divided by that other layer and the
    result is a *water* supply.

    [`GrowingSeasonPrecipitation`](@ref) therefore declares **no supply type of its own**.
    The spec is rewritten onto [`Precipitation`](@ref) before a supply is built, because the
    division has to happen on the source grid - it is cell-wise but nonlinear, so it does not
    commute with regridding.

    What the catalogue *can* rule out is `category = balance`: a sign-indefinite quantity
    such as the climate moisture index (precipitation minus evapotranspiration) can never be
    a supply, because a negative amount of a consumable has no meaning.

## Layers that are neither

A layer on an axis with no canonical unit and no supply type cannot take part in the
simulation - no tolerance can be built against it, and no species can consume it. That does
**not** make it useless. Such a layer is still worth having on the same grid, at the same
resolution, built by the same machinery: as ground truth to compare simulation output
against, as a backdrop for display, or as a covariate to reason about afterwards.

[`materialise`](@ref EcoSISTEM.materialise) is the route - it puts any spec on a decided grid
without requiring a tolerance or a role, handing back the layer it becomes:

```@example units
using Unitful, Unitful.DefaultSymbols
area = StudyArea(extent = (30.0km, 30.0km), cellsize = 10.0km, verbosity = :silent)
reference = EcoSISTEM.materialise(UniformSpec(5.0, axis = ClimateMoisture), area)
(size = size(reference.matrix), value = first(reference.matrix))
```

Such a layer cannot be a *member* of a built environment: `GridHabitat` requires one
species tolerance per regime layer, so a reference layer is held alongside the ecosystem
rather than inside it.

### And a reference layer can still reach the simulation - as an ingredient

It need not stay outside. A [`ConstructedSpec`](@ref) declares the axis of its **result**, not
of its inputs, so a layer that carries no axis can be *combined into* one that does:

> a soil-type layer has no tolerance and no supply of its own, yet it legitimately changes how
> much of the rainfall is available, through its drainage properties. Combine the two and
> declare the result on [`Precipitation`](@ref), and *that* layer is what carries a tolerance
> or becomes a supply.

This resolves an apparent contradiction. "An axis that declares nothing cannot be simulated"
and "a soil layer usefully modifies water availability" are both true, because they are
statements about **different layers**: the input has no axis, the result does. The reference
layer is an ingredient; the axis is a property of the dish.

## Declaring your own niche axis

[`@nicheaxis`](@ref) is the **only supported way** to add an axis - this package's own are declared
with it too. It emits the internal methods together, so they cannot disagree.

```@example units
@nicheaxis(SoilTemperature <: Temperature)
EcoSISTEM.canonicalunit(SoilTemperature)
```

That is the simplest case: a subtype that declares nothing **inherits its parent's answers**, so
`SoilTemperature` is a temperature in `K` and is matched to temperature layers. An axis is always an
`abstract type`, so one that is a leaf today can gain children tomorrow without becoming a breaking
change.

To declare an axis that measures something new, give it a unit - and say which role it plays:

| what you write | what it means |
| --- | --- |
| `condition = U` | a **condition**: layers are canonicalised to `U`, and a tolerance is built in it |
| `resource = U supply = S demand = D` | a **resource**: species consume it. All three, or none |
| both | a condition *and* a resource, as [`Precipitation`](@ref) is |
| `condition = nothing` | *not* a condition - with `resource`, a supply-only axis such as [`CarbonFlux`](@ref) |
| `reference` | **neither** - buildable and composable, never simulated |
| nothing at all | inherit from the parent |

`bounds = (lo, hi)` is optional on any of them and is stated in that axis's own canonical unit.

### Axes whose values are class labels

`categorical = true` says an axis holds **class codes** rather than measurements, so a layer on it is
resampled by nearest class instead of being interpolated - averaging two land-cover codes gives a
third code that means something else entirely.

```@example units
(EcoSISTEM.iscategorical(LandCoverTypology), EcoSISTEM.iscategorical(Temperature))
```

Declaring it on a group covers the whole branch, which is why this package declares it **once**, on
[`TypologyAxis`](@ref); `LandCoverTypology` and `ClimateTypology` declare nothing and inherit it. The
default is `false`, so an axis that says nothing anywhere up its chain of parents holds measurements.

`categorical` is independent of the roles above, and in particular **combines with `reference`**: a
land-cover layer carried purely as ground truth is still put on a grid, and still must not be
interpolated between its classes. This is where it differs from `bounds` - a bound is *stated in* a
unit, so an axis that declares none has nothing to state one in, whereas how values may be combined
matters whether or not a species responds to them.

### Two rules the macro enforces, and why

!!! warning "Declaring the abstract type by hand does not register it"
    ```julia
    abstract type MyAxis <: EcoSISTEM.NicheAxis end   # ← not enough
    ```
    This gives you a type, but none of the methods the rest of the package asks for. It used to
    answer `nothing` for its canonical unit - that is, *"no axis was named"* about an axis you had
    just named - and then accept any unit for it. It now **errors at first use**, naming the macro.

!!! warning "An axis may not contradict the unit it inherits"
    ```julia
    @nicheaxis(MyGroup <: EcoSISTEM.NicheAxis, condition = u"K")
    @nicheaxis(MyLeaf <: MyGroup, condition = u"K" * u"d")   # ← refused
    ```
    Axis matching is **subtype-aware**: a layer on a general axis can be promoted onto a more
    specific one when the two are paired. Promotion relabels the layer without touching its values -
    and there is no meaningful conversion it could do instead, because a temperature and an
    accumulation of temperature over time are different *quantities*, not one quantity in two units.

    So a group whose leaves genuinely measure different things should declare **no unit at all**
    and let each leaf state its own. That is what [`TemperatureAxis`](@ref) does: it groups
    *things to do with temperature*, and its leaves range over `K`, `K*day` and dimensionless
    ratios.

    *Shedding* a claim is different from contradicting one, and stays legal: `reference` and
    `condition = nothing` both say *"I assert no condition unit"*, which promotes safely.

## Known gaps

These are recorded here because they are the difference between what this page describes and
what a reader might reasonably expect.

**`nothing` means exactly two things.** An axis that declares no canonical unit is either
resource-only by design ([`CarbonFlux`](@ref)) or a pure grouping node whose leaves each carry their
own. It never means *"simply undecided"*: an axis that has not run [`@nicheaxis`](@ref) **errors**
rather than answering `nothing`, so "nobody declared this" cannot be spelled the same way as a
deliberate refusal.

**One leaf axis declares no condition unit**, and it is the deliberate case:
[`CarbonFlux`](@ref) is resource-only. Every other leaf declares one, so a rate-dimensioned axis such
as [`Evapotranspiration`](@ref) or [`ClimateMoisture`](@ref) is scaled by a stated unit rather than
escaping the question.

## How the rule was arrived at

The package once inferred a layer's *meaning* from its **unit**, and the clearest illustration of why
that cannot work is the defect it caused:

> A supply was built by looking at its value's dimension. `m/s` and `mm/day` are both `𝐋 𝐓^-1`, so a
> **wind-speed layer offered as a supply was silently built as a water supply** - 3 m/s over a 1 km^2
> cell becoming 2.6 * 10^14 L/day of rainfall.

Nothing about that is a missing check: a unit genuinely cannot say what a quantity *means*, because
two different quantities can share one. The fix was to make the axis the single declaration, on every
side:

  - a **supply**'s type comes from `supplytype(axis)`, and an axis declaring no resource is refused
    by name;
  - a **demand**'s type comes from `demandtype(axis)`, with the axis named at the call site
    (`demandaxis`), and the unit checked against it rather than consulted to choose it;
  - a **tolerance** is built in its axis's frame, as it always was;
  - the free/dimensionless family - the one entry that could not have an axis, because it meant
    "no axis" - was removed on both sides, which is what made deleting the tables possible at all.

Two families arrived after the rule and cost nothing to add, which is the point of it: **space**
(`SurfaceArea`, the first resource that is a stock rather than a flow) needed only a declaration, and
`pet`/`cmi` became flows by declaring their canonical units rather than by changing any code that
reads them.
