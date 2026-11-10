```@meta
CurrentModule = EcoSISTEM
```

# Species

The species of a run are declared before it is built, in a [`SpeciesPool`](@ref): every species it
may hold, each once, with the values known of each. A species list is then built from the pool one
spec per trait - a tolerance, a demand, a movement - by [`build_species`](@ref). Nothing takes
values per species by position: values reach the pool keyed by species, and a spec names their
column.

```@setup sp
using EcoSISTEM, EcoSISTEM.Units
using DataFrames, Distributions
using Unitful, Unitful.DefaultSymbols
```

## A pool

A pool is made from a table - a CSV file, its provenance record, or any table in memory - from a
list of names, or from a count. Every column but the names is a value of each species, `missing`
where it is unknown:

```@example sp
envelopes = DataFrame(species = ["Acacia tortilis", "Aloe vera", "Ficus sycomorus"],
                      tmin = [285.2, 288.1, missing],
                      ssr = [1.4e7, 1.6e7, 1.5e7])
pool = SpeciesPool(envelopes, verbosity = :silent)
```

Values computed later join it by species, from a table or from a function of each species' row, and
the pool says how each name of its data was matched:

```@example sp
widths = DataFrame(species = ["Ficus sycomorus", "Aloe vera"], tmin_sd = [2.5, 1.8])
add_trait!(pool, widths, verbosity = :silent)
add_trait!(pool, :warm) do row
    return coalesce(row.tmin, 0.0) > 286.0
end
DataFrame(pool)
```

A species may be deactivated - kept, but to hold no individuals until it is activated - or removed;
both take names, a group at a rank where the pool has a name authority (`rank = "genus"`), or a
function of the row. With a name authority, such as the [`WFOPlantList`](@ref), every name is
standardised to the authority's accepted species before anything else, and a name it does not
recognise is discarded with a warning:

```julia
pool = SpeciesPool("envelopes.csv", authority = WFOPlantList(release = "2026-06"),
                   merge = MergeWith(first))
```

A pool may be given a phylogeny with [`add_phylogeny!`](@ref): a tree of
[Phylo](https://github.com/EcoJulia/Phylo.jl)'s, read, built or drawn, whose tips are matched to the
pool's species by name as a table's rows are (a tip written `Quercus_robur` matching `Quercus robur`).
Tips naming no species of the pool are pruned. Every species needs a tip, inactive ones too, since a
species list is built from every species of its pool: `onmissing` says what becomes of one with
none, refused by default ([`RefuseMissingSpecies`](@ref)) or removed from the pool
([`RemoveMissingSpecies`](@ref)), while a copy is placed beside its original. A species list built from the pool is then phylogenetic, so a diversity
measure of its ecosystem is too. Traits evolved along the tree, or read with it, join the pool as any
other values do:

```julia
using Phylo
tree = open(parsenewick, "plants.tre")
add_phylogeny!(pool, tree)
```

A table of records - occurrences, say - often names species a run does not hold. Filtering it by
[`IsActiveSpecies`](@ref) keeps the rows naming one of the active species of a pool or a species list,
each name found as the pool finds names (through its authority, where it has one, so a synonym
counts), and reports what it removed, by kind and with how many rows each name took. `filter!`
changes a `DataFrame` in place; `filter` returns a new table:

```@example sp
records = DataFrame(species = ["Aloe vera", "Opuntia ficus-indica", "Ficus sycomorus",
                               "Aloe vera"])
filter!(IsActiveSpecies(pool, verbosity = :silent), records)
```

[`InActiveCell`](@ref) filters by position instead, keeping the rows whose latitude and longitude
(Darwin Core's `decimalLatitude` and `decimalLongitude` by default) lie in one of a study area's or a
habitat's active cells. A record falling just outside them - a coastal plant whose point lands in a sea
cell - is kept where an active cell lies within `snap`, a ground distance:

```julia
filter!(InActiveCell(habitat, snap = 10km), occurrences)
```

[`FromOccurrences`](@ref) places a run's starting individuals from such records, matching and placing
them as the filters do. Each record is one individual by default, so the species keep the
proportions their records give them; `per_record` places that many for each record instead, and
`total` spreads a number over all the records in proportion to them; and a `do` block makes any other
table from the records, each with its species, its `cell` and, where dated, its `time`:

```julia
eco = build_ecosystem(FromOccurrences(occurrences, snap = 10km, total = 100_000), species, habitat)
recent = FromOccurrences(occurrences) do rows, habitat, context
    return rows[rows.year .>= 2000, [:species, :cell]]    # one individual a record since 2000
end
```

To start each cell as full as its resources allow, wrap a placement in [`ToCapacity`](@ref): its
numbers become the shares of the species in each cell, and the cell is filled until what they demand
of its scarcest resource is what it supplies - `ToCapacity(recent)`, say, or
`ToCapacity(FromOccurrences(occurrences))`.

A `weight` weighs each record by a layer at its position, read at the layer's own resolution rather
than the habitat's. [`WeightPerRecord`](@ref) gives each record the layer's value in its cell, as
suits records that are equal samples, a systematic survey's; [`WeightPerCell`](@ref) gives each cell
its value once, shared among the records made there, as suits records whose numbers follow how hard
each place was searched. `scale` makes each unit of weight that many individuals; without one the
weights are shares, for `ToCapacity`, and a `do` block sees each record's `weight`:

```julia
# A vegetation index times each 0.05-degree cell's land area, from a file of your own.
region = StudyArea(within = NaturalEarthSpec("Kenya"), cellsize = 0.05°)
amount = EcoSISTEM.materialise(RasterFileSpec("evi_land_area.tif", axis = SurfaceArea,
                                              unit = km^2), region)
seeded = FromOccurrences(occurrences, weight = WeightPerCell(amount), scale = 2.0 / km^2)
```

## Species from a pool

Each trait is a spec: [`ToleranceSpec`](@ref) names a response family's parameters as the family
names them, [`DemandSpec`](@ref) a demand on one resource, [`MovementSpec`](@ref) a dispersal mode
and its kernel. Each parameter is a trait value - one value for every species, a
[`TraitColumn`](@ref) of the pool, a [`DrawFrom`](@ref) a distribution, or a function of the row. A
column of plain numbers says its unit, and is noted in it, so a wrong unit shows; a missing value is
filled by the column's `fill`:

```@example sp
niche = ToleranceSpec(Temperature,
                      μ = TraitColumn(:tmin, unit = K, fill = 286.0),
                      σ = TraitColumn(:tmin_sd, unit = K, fill = Uniform(1.0, 3.0)))
dispersal = MovementSpec(InsectDispersal(),
                         meandist = DrawFrom(Uniform(0.6, 2.4), unit = km))
nothing # hide
```

Before building anything, [`fill_traits`](@ref) shows what these specs will read from the pool. It
returns two tables, one row per species. `values` is the pool's own table with every missing value
the specs read filled in - here *Ficus sycomorus*' unknown `tmin` becomes 286.0, and *Acacia
tortilis*' unknown `tmin_sd` is drawn from `Uniform(1.0, 3.0)` - and a column for each value drawn
for every species, here each species' mean dispersal distance, `movement_meandist`:

```@example sp
filled = fill_traits(pool, niche, dispersal, seed = 1, verbosity = :silent)
filled.values
```

`filled` marks which of those values were filled in or drawn, rather than given in the pool:

```@example sp
filled.filled
```

[`build_species`](@ref) takes exactly this step first, with the same `seed` giving the same values,
and keeps both tables on the species list it builds:

```@example sp
species = build_species(pool, tolerance = niche,
                        demand = DemandSpec(SolarRadiation, TraitColumn(:ssr, unit = J / day)),
                        movement = dispersal, seed = 1, verbosity = :silent)
species.traits.movement_meandist == filled.values.movement_meandist
```

Every value drawn - a fill drawn from a distribution, a `DrawFrom` - comes from a stream of the
species' own, seeded from `seed` and the species' key in the pool, so a species draws the same
values however many others are drawn for.

Species of one pool may disperse by different modes, each with its own kernel, chosen by a key of
each species' row with [`ChosenBy`](@ref):

```@example sp
syndromes = SpeciesPool(DataFrame(species = ["a", "b", "c"],
                                  syndrome = ["wind", "insect", "wind"]),
                        verbosity = :silent)
dispersers = build_species(syndromes,
                           tolerance = ToleranceSpec(Temperature, μ = 10.0°C, σ = 2.0K),
                           demand = DemandSpec(SolarRadiation, 1.0e7J / day),
                           movement = ChosenBy(:syndrome,
                                               "wind" => MovementSpec(WindDispersal(),
                                                                      meandist = 5.0km,
                                                                      scale = 10.0km),
                                               "insect" => MovementSpec(InsectDispersal(),
                                                                        meandist = 1.0km)),
                           verbosity = :silent)
dispersers.movement
```
