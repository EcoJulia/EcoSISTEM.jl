# Named regions

A study area has to be positioned somewhere. You can give a bounding box, or your own shapefile - or
you can name the ground: `"Scotland"`, `"Madagascar"`, `"GREAT BRITAIN"`. The names come from
[Natural Earth](https://www.naturalearthdata.com/), which is in the public domain, at 1:10m.

Two things use them. `EcoSISTEM.boundingbox` gives a name's extent, read from a table shipped with
the package at **no download**. `NaturalEarthSpec` gives its actual outline, fetching the polygons
the first time it is built. Both resolve a name by the same rule, so the box that one reports is the
box the other's shape has.

## A name needs its level

A name on its own often means more than one thing, so every region is identified by a **level** as
well as a name.

```@example regions
using EcoSISTEM, Unitful, Unitful.DefaultSymbols
length(EcoSISTEM.naturalearth_levels())
```

The levels divide into political ones (`ADMIN` a country, `SUBUNIT` its constituent countries),
statistical groupings (`REGION_UN`, `SUBREGION`), a code (`ISO_A3_EH`), and Natural Earth's physical
vocabulary, prefixed `Physical` - `Physical Island`, `Physical Continent`, `Physical Desert` and the
rest.

Where the levels agree you need not say which you mean:

```@example regions
EcoSISTEM.boundingbox("Scotland")
```

Where they disagree, naming one is required rather than guessed at, and the error tabulates the
choice:

```@example regions
try
    EcoSISTEM.boundingbox("Africa")
catch e
    println(sprint(showerror, e))
end
```

`EcoSISTEM.naturalearth_regions(level)` lists what exists at one, with each region's box and area:

```@example regions
EcoSISTEM.naturalearth_regions("Physical Tundra")
```

Passing a region where a level belongs is a common slip, and the error says which level that region
is actually defined at rather than only that no such level exists.

## How much of a name to take

A name almost never denotes one connected piece of ground. `coverage` says how much to take, and the
default is **everything the name covers**, because that is what the source means by it.

```@example regions
france = EcoSISTEM.boundingbox("France", level = "ADMIN")
```

That reaches from French Guiana in South America to Mayotte in the Indian Ocean. For the ground most people mean, ask for the principal
landmass:

```@example regions
EcoSISTEM.boundingbox("France", level = "ADMIN", coverage = LargestLandmass())
```

`LargestLandmass(count = n)` takes the largest `n` - New Zealand's two main islands, Japan's four -
and `LandmassesAbove(threshold)` takes everything clearing a threshold, which is how the United
Kingdom loses Rockall. The threshold is either an **area** or a **share of the region's own total**:

```julia
LandmassesAbove(1km^2)     # an absolute size
LandmassesAbove(5percent)  # ...or a fraction of the region, which travels better between regions
```

A share must be a percentage rather than a bare number. `0.05` and `5percent` are the same quantity,
but only one of them says which it means when read beside `1km^2`, so the bare number is refused.

!!! tip "Check the share before asking for the principal landmass"
    A region's *share* - how much of its area its largest component holds - is what says whether
    `LargestLandmass()` is a sensible answer for it. New Zealand's is 56%, so asking for its
    principal landmass silently returns South Island alone; the Solomon Islands' is 20%.
    `naturalearth_regions` prints the share whenever it is below 90%.

!!! warning "Some regions have no bounding box at all"
    54 of the shipped regions cross the antimeridian, and an interval of longitude cannot describe
    one. `boundingbox` refuses those and says so, naming `LargestLandmass()` as the remedy - a single
    landmass cannot cross the date line, Natural Earth having split its polygons there.

## Finding a region

`investigate_regions` asks which named regions relate to something you already have - a coordinate, a
raster, a layer, a study area:

```@example regions
EcoSISTEM.investigate_regions(LatLong(55.95°, -3.19°), kind = :political, limit = 5)
```

`Encloses` is the default. `Overlaps(x)` asks which regions your data reaches into, ordered by how
much ground they share; `Within(x)` asks which lie entirely inside it, which is the question *"what
can I simulate in full with the data I have?"*

A row of the report can be turned straight into a spec. Use `only` when the answer should be unique -
it refuses if it is not - and `first` to take the best by the report's own ordering:

```@example regions
match = only(EcoSISTEM.investigate_regions(LatLong(55.95°, -3.19°), level = "SUBUNIT"))
NaturalEarthSpec(match)
```

### Boxes, and then outlines

By default the query compares **bounding boxes**, which is what makes it free - and it is loose.
**Norway's box encloses Edinburgh**, running west to Jan Mayen and north to Svalbard, so a query at
that coordinate lists Norway beside the United Kingdom.

`exact = true` checks the survivors against the regions' real outlines instead. It downloads the
geometry, so it is not free - but it removes those false positives, and it reaches the 54 regions
that cross the antimeridian, which have no box a query can compare at all:

```julia
investigate_regions(LatLong(55.95°, -3.19°), kind = :political, exact = true)
# Norway is gone; CONTINENT Europe, which wraps, appears
```

Refinement is lazy and in box order, stopping as soon as the answer cannot change - refining only
ever removes a match or shrinks its overlap, so a confirmed `limit` cannot be displaced by anything
later. A continental query that matches 361 regions by box typically fetches a few dozen, and the
report says how many.

## Using a region as a study area

`NaturalEarthSpec` is a `within` mask, so it both restricts which cells are simulated and sets the
grid's extent:

```julia
area = StudyArea(regime = temperature, supply = rainfall,
                 within = NaturalEarthSpec("Scotland", coverage = LargestLandmass()),
                 crs = EPSG(27700), cellsize = 1.0km)
```

`outline = false` takes the region's bounding box instead of its coastline, which is the cheaper
thing to want when the region is only saying *where* to work.

## Combining regions

Regions compose as **geometry**, through `ConstructedShapeSpec` - the vector mirror of
`ConstructedRasterSpec`. The result is exact and carries no resolution of its own, so the study grid
is still free to be decided afterwards.

```@example regions
uk = NaturalEarthSpec("United Kingdom", level = "ADMIN")
ireland = NaturalEarthSpec("Ireland", level = "ADMIN")
man = NaturalEarthSpec("Isle of Man", level = "ADMIN")
ConstructedShapeSpec(ShapeUnion(), uk, ireland, man)
```

Members may be any shape spec, including a `ShapeSpec` of your own study area and another
`ConstructedShapeSpec`. Operations are `ShapeUnion`, `ShapeIntersection` and `ShapeDifference` for
combining two or more, and `ShapeBuffer`, `ShapeSimplify` and `ShapeConvexHull` for transforming
exactly one:

```@example regions
scotland = NaturalEarthSpec("Scotland", coverage = LargestLandmass())
ConstructedShapeSpec(ShapeBuffer(50km), scotland)
```

An **arbitrary function** is accepted too, exactly as `ConstructedRasterSpec` takes any `combine`, so
anything the geometry library offers is reachable without a new operation.

## If you want X, use...

| if you want | write |
|---|---|
| France, whole | `NaturalEarthSpec("France", level = "ADMIN")` |
| France, metropolitan | `NaturalEarthSpec("France", level = "GEOUNIT")` |
| France continentale | the same, `coverage = LargestLandmass()` |
| Scotland, all of it | `NaturalEarthSpec("Scotland", level = "SUBUNIT")` |
| the Scottish mainland | the same, `coverage = LargestLandmass()` |
| Great Britain | `NaturalEarthSpec("United Kingdom", level = "ADMIN", coverage = LargestLandmass())` |
| ...or as a named island | `NaturalEarthSpec("GREAT BRITAIN", level = "Physical Island")` |
| GB and Northern Ireland | `coverage = LargestLandmass(count = 2)` |
| New Zealand, both islands | `coverage = LargestLandmass(count = 2)` |
| Japan, the four main islands | `coverage = LargestLandmass(count = 4)` |
| Europe without Russia | `NaturalEarthSpec("EUROPE", level = "Physical Continent")` |
| a country by its code | `NaturalEarthSpec("FRA", level = "ISO_A3_EH")` |
| the British Isles **with Shetland** | `ConstructedShapeSpec(ShapeUnion(), uk, ireland, man)` |
| ...and without Rockall | the same, `coverage = LandmassesAbove(1km^2)` |
| your own study area plus a country | `ConstructedShapeSpec(ShapeUnion(), ShapeSpec("mysite.geojson"), ireland)` |
| a country minus one of its parts | `ConstructedShapeSpec(ShapeDifference(), uk, scotland)` |
| within 50 km of a coastline | `ConstructedShapeSpec(ShapeBuffer(50km), coast)` |
| a coarser outline for a coarse grid | `ConstructedShapeSpec(ShapeSimplify(0.1°), scotland)` |
| the area an archipelago occupies | `ConstructedShapeSpec(ShapeConvexHull(), islands)` |
| just the box, no coastline | any of the above with `outline = false` |

## These are cartographic outlines

Natural Earth draws its shapes to be printed at a scale, not to be authoritative boundaries, and this
bites in practice rather than in principle. Two measured examples:

- Its own **`BRITISH ISLES`** polygon stops at 59.80°N, so it **omits Shetland**. The union of the
  United Kingdom, Ireland and the Isle of Man reaches 60.85°N and does not.
- Its physical continents are **landmass** outlines, so **`EUROPE` does not contain the British
  Isles** at all - Paris is inside it, central England is not. An intersection with it is empty, and
  is refused rather than returned as a grid with nothing active.

A third to know about: the physical `Channel Islands` are the **Californian** ones.

Use a named region as a convenience. If your study needs a particular boundary, supply it as a
`ShapeSpec` of your own - and combine it with a named one if that helps.
