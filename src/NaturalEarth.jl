# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The vocabulary of named geographic regions: how much of a name to take, and where each kind of
# name is defined.
#
# Natural Earth publishes the same world twice - a *cultural* set whose regions are countries and
# the groupings statisticians put them in, and a *physical* set whose regions are landmasses,
# ranges, deserts and basins. Neither is a superset of the other, and for the same word they can
# mean different things: the cultural `CONTINENT` "Europe" is a list of whole countries and so
# reaches the Pacific, while the physical `Continent` "EUROPE" is a coastline stopping at the Urals.
# A level therefore names both the file and the attribute a region is defined by, and a name is only
# meaningful together with its level.
#
# This file declares that vocabulary and nothing that reads it: the selection itself is in
# `rasters.jl`, beside the shapefile machinery it shares.

using Unitful

using Unitful.DefaultSymbols

# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------

"""
    AbstractCoverage

How much of what a region name covers to take - [`AllTerritories`](@ref) or
[`LargestLandmass`](@ref).

A name almost never denotes one connected piece of ground. "France" includes Guadeloupe and
Martinique, "Norway" includes Bouvet Island in the South Atlantic, and "Chile" includes Easter
Island, so the extent of everything a name covers can be many times the extent of the ground most people
mean by it. A coverage says which of the two is wanted, and is the same choice whether a bounding
box or an actual shape is being asked for.
"""
abstract type AbstractCoverage end

"""
    AllTerritories()

Take everything the name covers, however scattered.

This is faithful to the source and is what a political question usually wants: the territory of
France really does include French Guiana. It is also what makes a selection able to span the
antimeridian, since a country with territory either side of the date line has no single interval of
longitude containing it.
"""
struct AllTerritories <: AbstractCoverage end

"""
    LargestLandmass(count::Integer = 1)

Take only the `count` largest connected pieces of ground the name covers.

Dissolving a selection and keeping its largest connected component is what turns "France" into
France continentale, "United Kingdom" into Great Britain and "Scotland" into the Scottish mainland,
without an area threshold or a hand-written list of exceptions. It is the largest *component* rather
than the largest *part*: components are measured after neighbouring features have been merged, so a
continent comes out as its mainland rather than as its largest country.

A landmass selected this way cannot cross the antimeridian, because Natural Earth splits its
polygons at the date line and so no single connected component spans it.

The count is a keyword rather than a positional argument, unlike the sibling
[`RandomCells`](@ref) and [`SpreadingCells`](@ref). Those name a plural, so a bare number attaches to
it and `RandomCells(20)` reads as twenty cells; this one names a single landmass, so a bare number
would attach to nothing and could as easily be read as an ordinal. Naming it at the call site is
what says which was meant.

# Arguments

  - `count`: how many components to keep, largest first. A `count` beyond the number of components
    returns all of them, so `LargestLandmass(count = 2)` on a single-component region is the region.
"""
struct LargestLandmass <: AbstractCoverage
    count::Int

    function LargestLandmass(; count::Integer = 1)
        count >= 1 ||
            throw(ArgumentError("LargestLandmass needs at least one component, got count = $count"))
        return new(Int(count))
    end
end

# These nest inside a spec and a report, where the default `show` prints the whole type signature
# for what the caller wrote as `LargestLandmass()`.
Base.show(io::IO, ::AllTerritories) = print(io, "AllTerritories()")
function Base.show(io::IO, c::LargestLandmass)
    return print(io, "LargestLandmass(count = $(c.count))")
end

# ---------------------------------------------------------------------------
# Levels
# ---------------------------------------------------------------------------

# Where Natural Earth's files are published. The path is `<base>/<resolution>/<category>/<dataset>.zip`.
const _NE_BASE = "https://naciscdn.org/naturalearth"

# 1:10m is the only resolution shipped. The coarser scales are a different answer rather than a
# cheaper one - at 1:110m Malta is absent entirely and the United Kingdom loses Rockall, moving its
# western edge by over six degrees - so a box taken from one scale would not describe a shape taken
# from another. Simplify a geometry on demand instead of selecting a coarser source.
const _NE_RESOLUTION = "10m"

"""
    NaturalEarthLevel

One kind of named region: which Natural Earth file defines it, and which attribute of that file
carries the name.

Levels are what make a name unambiguous. "Africa" is a cultural `CONTINENT` of 55 countries and a
UN `REGION_UN` of 62, with different members and different extents, so a name is always looked up at
a stated level.

# Fields

  - `name`: the level's own name, as it is written in the shipped region table.
  - `dataset`: the Natural Earth file it comes from, without the `.zip` extension.
  - `category`: `"cultural"` or `"physical"`, the source's own division and the directory it
    publishes each file under.
  - `field`: the attribute of that file whose values name a region at this level.
  - `kind`: what sort of division it is - `:political`, `:statistical`, `:physical` or `:code`.
  - `description`: what the level means, and where it is likely to surprise.
  - `within`: an attribute and value restricting which features belong to this level, or `nothing`
    where the whole file does. The physical file holds every kind of landform in one layer, so its
    levels are distinguished this way.
"""
struct NaturalEarthLevel
    name::String
    dataset::String
    category::String
    field::String
    kind::Symbol
    description::String
    within::Union{Pair{String, String}, Nothing}
end

function Base.show(io::IO, l::NaturalEarthLevel)
    return print(io, "NaturalEarthLevel(\"", l.name, "\", :", l.kind, ")")
end

function Base.show(io::IO, ::MIME"text/plain", l::NaturalEarthLevel)
    println(io, "NaturalEarthLevel \"", l.name, "\" (", l.kind, ")")
    println(io, "  from: ", l.dataset, " (", l.category, "), field ", l.field)
    return print(io, "  ", l.description)
end

# The four admin-0 files, which differ in how finely they split a country rather than in their
# fields: `countries` maps France as one feature including its overseas departments, `map_units`
# separates metropolitan France, and `map_subunits` splits the United Kingdom into its four
# constituent countries. Several attributes name a region on the same file, so a file appears once
# per attribute worth selecting on.
const _NE_COUNTRIES = "ne_10m_admin_0_countries"
const _NE_MAP_UNITS = "ne_10m_admin_0_map_units"
const _NE_MAP_SUBUNITS = "ne_10m_admin_0_map_subunits"
const _NE_SOVEREIGNTY = "ne_10m_admin_0_sovereignty"
const _NE_PHYSICAL = "ne_10m_geography_regions_polys"

# The landform classes of the physical file, each of which becomes a level. Taken from the file
# itself rather than from its documentation.
#
# "Dragons-be-here" is deliberately absent: its single feature is Null Island at (0, 0), a
# cartographers' joke, and it would otherwise be offered as a selectable region.
const _NE_PHYSICAL_CLASSES = ("Continent", "Island", "Island group",
                              "Range/mtn",
                              "Plateau", "Desert", "Pen/cape", "Peninsula",
                              "Isthmus", "Geoarea", "Plain", "Lowland",
                              "Foothills", "Basin", "Valley", "Gorge",
                              "Depression", "Delta", "Coast", "Tundra",
                              "Wetlands", "Lake")

# A physical level: one landform class of the shared physical layer, selected on `FEATURECLA`.
#
# The level is named "Physical <class>" rather than by the class alone, because the physical
# `Continent` and the cultural `CONTINENT` would otherwise differ only in case - and they are
# exactly the pair a user must be able to tell apart, physical Europe stopping at the Urals where
# cultural Europe takes in the whole of Russia and reaches the Pacific.
function _physicallevel(class::AbstractString)
    return NaturalEarthLevel("Physical " * class, _NE_PHYSICAL, "physical",
                             "NAME",
                             :physical,
                             "Natural Earth's cartographic \"$class\" outlines. These are drawn " *
                             "to place a map label, not derived from data, so they carry no " *
                             "authority as a boundary.",
                             "FEATURECLA" => class)
end

"""
    NATURALEARTH_LEVELS

Every kind of named region the shipped region table covers, in the order they are reported.

The cultural levels come first, coarsest grouping to finest division, then the physical landform
classes. A level's `name` is what the table's own `Level` column holds.
"""
const NATURALEARTH_LEVELS = NaturalEarthLevel[NaturalEarthLevel("CONTINENT",
                                                                _NE_COUNTRIES,
                                                                "cultural",
                                                                "CONTINENT",
                                                                :political,
                                                                "Continental grouping of whole countries. It assigns each country entire, so " *
                                                                "all of Russia falls in Europe; the physical \"Continent\" level is the " *
                                                                "geographic alternative.",
                                                                nothing),
                                              NaturalEarthLevel("REGION_UN",
                                                                _NE_COUNTRIES,
                                                                "cultural",
                                                                "REGION_UN",
                                                                :statistical,
                                                                "Region of the UN geoscheme, the grouping used for UN statistics.",
                                                                nothing),
                                              NaturalEarthLevel("SUBREGION",
                                                                _NE_COUNTRIES,
                                                                "cultural",
                                                                "SUBREGION",
                                                                :statistical,
                                                                "Subregion of the UN geoscheme, a division of REGION_UN.",
                                                                nothing),
                                              NaturalEarthLevel("REGION_WB",
                                                                _NE_COUNTRIES,
                                                                "cultural",
                                                                "REGION_WB",
                                                                :statistical,
                                                                "Region as the World Bank groups countries.",
                                                                nothing),
                                              NaturalEarthLevel("SOVEREIGNT",
                                                                _NE_SOVEREIGNTY,
                                                                "cultural",
                                                                "SOVEREIGNT",
                                                                :political,
                                                                "Sovereign state, taking in every territory it is sovereign over.",
                                                                nothing),
                                              NaturalEarthLevel("ADMIN",
                                                                _NE_COUNTRIES,
                                                                "cultural",
                                                                "ADMIN",
                                                                :political,
                                                                "Country as Natural Earth maps it, including overseas departments and " *
                                                                "remote territories.",
                                                                nothing),
                                              NaturalEarthLevel("GEOUNIT",
                                                                _NE_MAP_UNITS,
                                                                "cultural",
                                                                "GEOUNIT",
                                                                :political,
                                                                "Map unit: a country split where its parts are mapped separately, which is " *
                                                                "what separates metropolitan France from its overseas departments.",
                                                                nothing),
                                              NaturalEarthLevel("SUBUNIT",
                                                                _NE_MAP_SUBUNITS,
                                                                "cultural",
                                                                "SUBUNIT",
                                                                :political,
                                                                "Map subunit, the finest split of a country Natural Earth offers, and what " *
                                                                "gives the constituent countries of the United Kingdom.",
                                                                nothing),
                                              NaturalEarthLevel("ISO_A3",
                                                                _NE_COUNTRIES,
                                                                "cultural",
                                                                "ISO_A3", :code,
                                                                "ISO 3166-1 alpha-3 country code, for looking a country up by code rather " *
                                                                "than by name.",
                                                                nothing),
                                              map(_physicallevel,
                                                  _NE_PHYSICAL_CLASSES)...]

# The level of this name, or `nothing`. Matched case-insensitively because the physical file's own
# names are inconsistently cased, and a caller cannot be expected to reproduce that.
function _findlevel(name::AbstractString)
    i = findfirst(l -> lowercase(l.name) == lowercase(name),
                  NATURALEARTH_LEVELS)
    return isnothing(i) ? nothing : NATURALEARTH_LEVELS[i]
end

# The cached download a level reads from. Several levels share one file, and `CachedAsset` keys the
# cache on the file name, so naming the same dataset twice costs one download.
function _nesource(level::NaturalEarthLevel)
    url = join((_NE_BASE, _NE_RESOLUTION, level.category,
                level.dataset * ".zip"),
               "/")
    return CachedAsset(NaturalEarthLevel, url)
end
