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

using CSV

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

"""
    LandmassesAbove(area)

Take every connected piece of ground the name covers that is at least `area` in size.

This is the coverage for "everything except the specks": the United Kingdom without Rockall, which
[`LargestLandmass`](@ref) can only express by counting components and so needs to know how many
there are. Rockall is 0.031 km2 and the next smallest British component is 2.536 km2, an eighty-fold
gap, so any threshold between them does the same thing.

Unlike the other coverages this one cannot be answered from the shipped region table, which records
only the largest few components' sizes - so [`boundingbox`](@ref) refuses it, while a spec built onto
a real grid has the geometry and can.

# Arguments

  - `area`: the smallest component to keep, as a `Unitful` area (`1km^2`, `500u"m^2"`).
"""
struct LandmassesAbove{A} <: AbstractCoverage
    area::A

    function LandmassesAbove(area)
        dimension(area) == dimension(1.0km^2) ||
            throw(ArgumentError("LandmassesAbove needs an area, got $(area) which is " *
                                "$(dimension(area))"))
        area > zero(area) ||
            throw(ArgumentError("LandmassesAbove needs a positive area, got $area"))
        return new{typeof(area)}(area)
    end
end

function Base.show(io::IO, c::LandmassesAbove)
    return print(io, "LandmassesAbove(", c.area, ")")
end

# Whether a coverage is the one a spec would have taken had none been named, so that `show` prints
# only what a caller would have to write. One method per coverage rather than a chain of `isa`, so a
# coverage added later fails here loudly instead of silently printing as the default.
#
# The default is `AllTerritories` because that is what Natural Earth means by a name: its "France"
# is the one including Guadeloupe. Taking only the principal landmass is a real choice about what
# ground is wanted, so it is the caller's to make and to see written down.
_isdefaultcoverage(::AllTerritories) = true
_isdefaultcoverage(::LargestLandmass) = false
_isdefaultcoverage(::LandmassesAbove) = false

# ---------------------------------------------------------------------------
# Combining regions
# ---------------------------------------------------------------------------

"""
    AbstractShapeOperation

How several named regions combine into one - [`ShapeUnion`](@ref), [`ShapeIntersection`](@ref) or
[`ShapeDifference`](@ref).

Regions compose as **geometries**, not as rasters. A union of outlines is exact and independent of
any grid, where combining rasterised masks would have to pick a resolution before the study grid was
decided and would carry two sets of edge effects into the answer.
"""
abstract type AbstractShapeOperation end

"""    ShapeUnion() <: AbstractShapeOperation - every member's ground, merged. """
struct ShapeUnion <: AbstractShapeOperation end

"""    ShapeIntersection() <: AbstractShapeOperation - only ground common to every member. """
struct ShapeIntersection <: AbstractShapeOperation end

"""
    ShapeDifference() <: AbstractShapeOperation

The first member's ground with every later member's removed.
"""
struct ShapeDifference <: AbstractShapeOperation end

"""
    ShapeBuffer(distance)

Grow the shape outwards by `distance`, or shrink it inwards where `distance` is negative.

This is how a study area becomes "within 50 km of this coastline", which no set operation can
express. It takes one member: a buffer transforms a shape rather than combining several.

⚠️ The distance is applied in the shape's own coordinates, which for Natural Earth are degrees of
latitude and longitude - so a buffer given as a length is converted at the equator and is narrower
in east-west terms the further from it the shape lies. Give an angle to say exactly what is meant.

# Arguments

  - `distance`: how far to grow, as a length (`50km`) or an angle (`0.5°`).
"""
struct ShapeBuffer{D} <: AbstractShapeOperation
    distance::D

    function ShapeBuffer(distance)
        dim = dimension(distance)
        dim in (dimension(1.0km), NoDims) ||
            throw(ArgumentError("ShapeBuffer needs a length or an angle, got $distance which is " *
                                "$dim"))
        return new{typeof(distance)}(distance)
    end
end

"""
    ShapeSimplify(tolerance)

Replace the outline with a coarser one no further than `tolerance` from it.

Only 1:10m polygons ship, that being the one scale at which a name's box and its shape can agree, so
this is how a coarser outline is had when a fine one is not wanted - a continent's coastline carries
tens of thousands of vertices that no continental-resolution grid can use. It takes one member.

# Arguments

  - `tolerance`: how far the simplified outline may depart from the original, as a length or an
    angle, in the same coordinates [`ShapeBuffer`](@ref) uses.
"""
struct ShapeSimplify{D} <: AbstractShapeOperation
    tolerance::D

    function ShapeSimplify(tolerance)
        dim = dimension(tolerance)
        dim in (dimension(1.0km), NoDims) ||
            throw(ArgumentError("ShapeSimplify needs a length or an angle, got $tolerance which " *
                                "is $dim"))
        tolerance > zero(tolerance) ||
            throw(ArgumentError("ShapeSimplify needs a positive tolerance, got $tolerance"))
        return new{typeof(tolerance)}(tolerance)
    end
end

"""
    ShapeConvexHull() <: AbstractShapeOperation

The smallest convex outline containing the shape, as a rubber band round it.

Takes one member. Useful for turning a scattered archipelago into the single area it occupies, and
for smoothing away a coastline's detail without choosing a tolerance.
"""
struct ShapeConvexHull <: AbstractShapeOperation end

Base.show(io::IO, o::ShapeBuffer) = print(io, "ShapeBuffer(", o.distance, ")")
function Base.show(io::IO, o::ShapeSimplify)
    return print(io, "ShapeSimplify(", o.tolerance, ")")
end

# How many members an operation takes. The set operations need at least two to mean anything; the
# transforms act on exactly one. Stated per operation so that the constructor can refuse a wrong
# count where it was written, rather than failing inside a reduce.
_minmembers(::AbstractShapeOperation) = 2
_minmembers(::ShapeBuffer) = 1
_minmembers(::ShapeSimplify) = 1
_minmembers(::ShapeConvexHull) = 1
# A bare function is the escape hatch, mirroring `ConstructedRasterSpec`'s `combine`: it is handed
# every member's geometry and may do anything with them, so no count can be required.
_minmembers(::Function) = 1

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
# "Dragons-be-here" is deliberately absent. Its single feature is Null Island at (0, 0), a
# cartographers' joke - and a real, valid polygon about a metre across, which is the practical
# objection rather than the joke: at the table's three decimal places its box rounds to zero width,
# so it would be the one region in 2 444 with no extent at all.
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

# The nine cultural levels and the physical ones, assembled in report order.
#
# Written as calls to a local `add` rather than as one array literal: the literal aligns every
# continuation line to its opening bracket, which put the descriptions past column 100.
function _buildlevels()
    levels = NaturalEarthLevel[]
    function add(name, dataset, field, kind, description)
        return push!(levels,
                     NaturalEarthLevel(name, dataset, "cultural", field, kind,
                                       description, nothing))
    end
    add("CONTINENT", _NE_COUNTRIES, "CONTINENT", :political,
        "Continental grouping of whole countries. It assigns each country entire, so all of " *
        "Russia falls in Europe; the physical \"Continent\" level is the geographic alternative.")
    add("REGION_UN", _NE_COUNTRIES, "REGION_UN", :statistical,
        "Region of the UN geoscheme, the grouping used for UN statistics.")
    add("SUBREGION", _NE_COUNTRIES, "SUBREGION", :statistical,
        "Subregion of the UN geoscheme, a division of REGION_UN.")
    add("REGION_WB", _NE_COUNTRIES, "REGION_WB", :statistical,
        "Region as the World Bank groups countries.")
    add("SOVEREIGNT", _NE_SOVEREIGNTY, "SOVEREIGNT", :political,
        "Sovereign state, taking in every territory it is sovereign over.")
    add("ADMIN", _NE_COUNTRIES, "ADMIN", :political,
        "Country as Natural Earth maps it, including overseas departments and remote territories.")
    add("GEOUNIT", _NE_MAP_UNITS, "GEOUNIT", :political,
        "Map unit: a country split where its parts are mapped separately, which is what separates " *
        "metropolitan France from its overseas departments.")
    add("SUBUNIT", _NE_MAP_SUBUNITS, "SUBUNIT", :political,
        "Map subunit, the finest split of a country Natural Earth offers, and what gives the " *
        "constituent countries of the United Kingdom.")
    add("ISO_A3_EH", _NE_COUNTRIES, "ISO_A3_EH", :code,
        "ISO 3166-1 alpha-3 country code, for looking a country up by code rather than by name. " *
        "Named for the source field it reads. Natural Earth ships both ISO_A3 and ISO_A3_EH; the " *
        "latter is unset for 14 countries against the former's 22, France and Norway among the " *
        "eight it recovers. The 14 that remain have no assigned code - Kosovo and Somaliland among " *
        "them - and so are absent from this level.")
    append!(levels, map(_physicallevel, _NE_PHYSICAL_CLASSES))
    return levels
end

"""
    NATURALEARTH_LEVELS

Every kind of named region the shipped region table covers, in the order they are reported.

The cultural levels come first, coarsest grouping to finest division, then the physical landform
classes. A level's `name` is what the table's own `Level` column holds.
"""
const NATURALEARTH_LEVELS = _buildlevels()

# The level of this name, or `nothing`. Matched case-insensitively because the physical file's own
# names are inconsistently cased, and a caller cannot be expected to reproduce that.
function _findlevel(name::AbstractString)
    i = findfirst(l -> lowercase(l.name) == lowercase(name),
                  NATURALEARTH_LEVELS)
    return isnothing(i) ? nothing : NATURALEARTH_LEVELS[i]
end

# The level of this name, or an error naming the closest matches.
#
# A wrong level otherwise surfaces as "no region named X at level Y", which blames the name when the
# level is what is wrong. Suggestions are by prefix and substring, which is what turns a reasonable
# guess at "ISO_A3" into the "ISO_A3_EH" this package actually ships.
function _checklevel(name::AbstractString)
    level = _findlevel(name)
    isnothing(level) || return level
    lower = lowercase(name)
    # Matched both ways round, so a guess longer than the real name is caught as well as one that is
    # shorter: "ISO_A3" finds "ISO_A3_EH", and "SUBUNITS" finds "SUBUNIT".
    near = [l.name
            for l in NATURALEARTH_LEVELS
            if occursin(lower, lowercase(l.name)) ||
        occursin(lowercase(l.name), lower)]
    return error("There is no region level called \"$name\"." *
                 (isempty(near) ? "" :
                  " Did you mean " * join(near, ", ") * "?") *
                 " `EcoSISTEM.NATURALEARTH_LEVELS` lists all $(length(NATURALEARTH_LEVELS)).")
end

# The cached download a level reads from. Several levels share one file, and `CachedAsset` keys the
# cache on the file name, so naming the same dataset twice costs one download.
function _nesource(level::NaturalEarthLevel)
    url = join((_NE_BASE, _NE_RESOLUTION, level.category,
                level.dataset * ".zip"),
               "/")
    return CachedAsset(NaturalEarthLevel, url)
end

# ---------------------------------------------------------------------------
# The shipped table
# ---------------------------------------------------------------------------

# How many of the largest components the table describes individually. Kept in step with the
# generator, `examples/scripts/naturalearth_regions.jl`, which writes that many columns.
const NPARTS = 5

_regionspath() = pkgdir(@__MODULE__, "data", "NaturalEarth", "regions.csv")

# The parsed table, indexed as level -> lowercased name -> row, and the levels each name appears at.
#
# Memoised because a bounding box is looked up per build and re-reading 2 444 rows each time is what
# would make the lookup cost milliseconds rather than microseconds - the same reason `_layertable`
# caches the layer catalogue.
const _REGION_INDEX = Ref{Any}(nothing)

function _regionindex()
    isnothing(_REGION_INDEX[]) || return _REGION_INDEX[]
    bylevel = Dict{String, Dict{String, Any}}()
    bynames = Dict{String, Vector{String}}()
    for row in CSV.File(_regionspath())
        get!(bylevel, row.Level, Dict{String, Any}())[lowercase(row.Name)] = row
        push!(get!(bynames, lowercase(row.Name)) do
                  return String[]
              end, row.Level)
    end
    _REGION_INDEX[] = (bylevel = bylevel, bynames = bynames)
    return _REGION_INDEX[]
end

# The row for a name at a stated level, or `nothing` if the name is not defined there.
function _regionrow(level::AbstractString, name::AbstractString)
    byname = get(_regionindex().bylevel, level, nothing)
    return isnothing(byname) ? nothing : get(byname, lowercase(name), nothing)
end

# The extent a coverage asks for, out of a table row.
#
# `AllTerritories` is the whole selection, which the row states directly. `LargestLandmass` unions
# the row's largest components, which it carries individually up to `NPARTS`; asking for more than
# that, but fewer than the region has, is the one question the table cannot answer, since the rest
# were never written down.
function _coveragebox(row, ::AllTerritories, region, lvl)
    return (west = row.West, south = row.South, east = row.East,
            north = row.North,
            wraps = row.Wraps)
end

# The shipped table records only the largest few components' sizes, so an area threshold cannot be
# applied to it - the components it would have to test were never written down. A spec built onto a
# real grid has the geometry and answers this fine, which is what the message points at.
function _coveragebox(row, c::LandmassesAbove, region, lvl)
    return error("`boundingbox` cannot answer `$c`: the shipped region table records the sizes of " *
                 "only the largest $NPARTS components, so it cannot say which of \"$region\" " *
                 "($lvl) clear $(c.area). Build the shape instead - " *
                 "`NaturalEarthSpec(\"$region\", level = \"$lvl\", coverage = $c)` - which reads " *
                 "the geometry.")
end

function _coveragebox(row, c::LargestLandmass, region, lvl)
    c.count >= row.Parts &&
        return (west = row.West, south = row.South, east = row.East,
                north = row.North, wraps = row.Wraps)
    c.count > NPARTS &&
        error("\"$region\" ($lvl) has $(row.Parts) separate components, and the shipped table " *
              "describes only the largest $NPARTS individually, so a box for the largest " *
              "$(c.count) cannot be given. Ask for $NPARTS or fewer, or for " *
              "`AllTerritories()`, or build the shape itself.")
    west = minimum(i -> row[Symbol("Part$(i)West")], 1:(c.count))
    south = minimum(i -> row[Symbol("Part$(i)South")], 1:(c.count))
    east = maximum(i -> row[Symbol("Part$(i)East")], 1:(c.count))
    north = maximum(i -> row[Symbol("Part$(i)North")], 1:(c.count))
    # A component cannot cross the date line, Natural Earth having split its polygons there, so a
    # union of components has an honest west-to-east interval however scattered they are.
    return (west = west, south = south, east = east, north = north,
            wraps = false)
end

# The coverage a *level comparison* is made with.
#
# Resolving a bare name asks whether the levels would give the same answer, which means reading the
# shipped table - and `LandmassesAbove` is the one coverage the table cannot answer, since it records
# the sizes of only the largest few components. Comparing the whole selections instead is the
# conservative substitute: it is what the threshold will be applied to, so levels that disagree
# there cannot agree afterwards, and a name whose full extents match at every level is selecting the
# same ground however it is later filtered.
_tablecoverage(c::AbstractCoverage) = c
_tablecoverage(::LandmassesAbove) = AllTerritories()

# Which level a bare name means, for the coverage being asked for.
#
# A name may be defined at several levels, and most that are give the *same* answer at all of them -
# "Scotland" is both a map unit and a map subunit of the same ground - so requiring a level there
# would be ceremony for nothing. Where they genuinely disagree there is no defensible choice to make
# on the caller's behalf, and it is an error showing what each would give.
#
# The comparison is of the boxes the *requested coverage* produces, not of the whole selections,
# because that is what the caller will receive. It matters: "Africa" is 55 countries as a continent
# and 62 as a UN region, whose full extents differ by 54 degrees of longitude - but both have the
# same principal landmass, so asking for that one is not ambiguous at all and is answered.
function _resolvelevel(name::AbstractString, coverage::AbstractCoverage)
    levels = get(_regionindex().bynames, lowercase(name), String[])
    isempty(levels) &&
        error("No region named \"$name\" at any level. Names are matched case-insensitively but " *
              "must otherwise be the source's own - \"United Kingdom\" rather than \"UK\", " *
              "\"South America\" rather than \"SouthAmerica\". " *
              "`EcoSISTEM.NATURALEARTH_LEVELS` lists the levels a name may be defined at.")
    length(levels) == 1 && return levels[1]
    rows = [_regionrow(l, name) for l in levels]
    compare = _tablecoverage(coverage)
    boxes = [_coveragebox(r, compare, name, l) for (r, l) in zip(rows, levels)]
    allequal((b.west, b.south, b.east, b.north) for b in boxes) &&
        return levels[1]
    return error("\"$name\" means different regions at different levels, and the extents they " *
                 "give for $coverage disagree:\n" *
                 _levelchoices(levels, rows, boxes) *
                 "Name one with `level = \"...\"`.")
end

# The candidate levels as a table, so that a caller can choose without looking anything up. Each
# line names the level, the source's own spelling of the name there, and the box it would give.
function _levelchoices(levels, rows, boxes)
    width = maximum(length, levels)
    return join(["    " * rpad(l, width) * "  " * rpad(r.Name, 18) *
                 "W " * lpad(b.west, 9) * "  S " * lpad(b.south, 8) *
                 "  E " * lpad(b.east, 9) * "  N " * lpad(b.north, 8) * "\n"
                 for (l, r, b) in zip(levels, rows, boxes)])
end
