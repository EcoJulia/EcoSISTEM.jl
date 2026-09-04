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

# Imported here rather than relied on from a later file: this file is included before
# `Coordinates.jl`, and a method signature resolves at definition time.
import Extents

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
    LandmassesAbove(threshold)

Take every connected piece of ground the name covers that clears `threshold`.

This is the coverage for "everything except the specks": the United Kingdom without Rockall, which
[`LargestLandmass`](@ref) can only express by counting components and so needs to know how many there
are. Rockall is 0.031 km2 and the next smallest British component is 2.536 km2, an eighty-fold gap,
so any threshold between them does the same thing.

Unlike the other coverages this one cannot be answered from the shipped region table, which records
only the largest few components' sizes - so [`boundingbox`](@ref) refuses it, while a spec built onto
a real grid has the geometry and can.

# Arguments

  - `threshold`: how big a component must be to be kept, either as an **area** (`1km^2`,
    `500u"m^2"`) or as a **share of the region's own total** (`5percent`, from `Unitful`). A share is
    the more portable of the two, since what counts as a speck depends on how big the region is: the
    same `1km^2` keeps every part of Great Britain and discards most of the Maldives.

    A share must be written as a percentage rather than as a bare number. `0.05` and `5percent` are
    the same quantity, but only one of them says which it means when read beside `1km^2`, and a bare
    number is refused for that reason.
"""
struct LandmassesAbove{T} <: AbstractCoverage
    threshold::T

    function LandmassesAbove(threshold)
        threshold isa Unitful.AbstractQuantity ||
            throw(ArgumentError("LandmassesAbove needs an area or a percentage, got the bare " *
                                "number $threshold. Write `$(threshold)percent` for a share of " *
                                "the region's own area, or give an area such as `1km^2`."))
        isarea = dimension(threshold) == dimension(1.0km^2)
        isshare = dimension(threshold) == NoDims
        (isarea || isshare) ||
            throw(ArgumentError("LandmassesAbove needs an area or a percentage, got $threshold " *
                                "which is $(dimension(threshold))."))
        threshold > zero(threshold) ||
            throw(ArgumentError("LandmassesAbove needs a positive threshold, got $threshold."))
        isshare && ustrip(NoUnits, threshold) > 1 &&
            throw(ArgumentError("LandmassesAbove was given a share of $threshold, which is more " *
                                "than the whole region; nothing can clear it."))
        return new{typeof(threshold)}(threshold)
    end
end

# A share is of the region's own total, so it only becomes an area once the components are known.
# One method per kind rather than a branch, so a threshold that is neither cannot reach here.
_thresholdarea(c::LandmassesAbove{<:Unitful.Area}, total) = c.threshold

function _thresholdarea(c::LandmassesAbove{<:Unitful.DimensionlessQuantity},
                        total)
    return total * ustrip(NoUnits, c.threshold)
end

function Base.show(io::IO, c::LandmassesAbove)
    return print(io, "LandmassesAbove(", c.threshold, ")")
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
    # The likeliest mistake is passing a region where a level belongs, so look it up as one: saying
    # "there is no such level" is true and useless when the answer is "that is a region, and here is
    # the level it lives at".
    at = get(_regionindex().bynames, lower, String[])
    isempty(at) ||
        return error("\"$name\" is a region, not a level - it is defined at " *
                     join(at, ", ") *
                     ". A level is the *kind* of region: try " *
                     "`EcoSISTEM.naturalearth_regions(\"$(at[1])\")` to list that kind, or " *
                     "`EcoSISTEM.boundingbox(\"$name\", level = \"$(at[1])\")` for this one.")
    # Matched both ways round, so a guess longer than the real name is caught as well as one that is
    # shorter: "ISO_A3" finds "ISO_A3_EH", and "SUBUNITS" finds "SUBUNIT".
    near = [l.name
            for l in NATURALEARTH_LEVELS
            if occursin(lower, lowercase(l.name)) ||
        occursin(lowercase(l.name), lower)]
    # A typo of a region name lands here too, having matched neither a level nor a region exactly,
    # so near region names are offered as well: "Yemn" is a misspelt country, not a misspelt level.
    nearnames = _nearnames(lower)
    return error("There is no region level called \"$name\"." *
                 (isempty(near) ? "" :
                  " Did you mean the level " * join(near, ", ") * "?") *
                 (isempty(nearnames) ? "" :
                  " Did you mean the region " * join(nearnames, ", ") * "?") *
                 " `EcoSISTEM.naturalearth_levels()` lists all " *
                 "$(length(NATURALEARTH_LEVELS)) levels.")
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
# generator, `data/src/naturalearth_regions.jl`, which writes that many columns.
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
              "`EcoSISTEM.naturalearth_levels()` lists the levels a name may be defined at.")
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

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

"""
    naturalearth_levels()

Return every kind of named region there is, as a vector of [`NaturalEarthLevel`](@ref)s.

A name means nothing without its level - "Africa" is a continent of 55 countries and a UN region of
62 - so this is where to look before naming one. Each entry says which file defines it, which
attribute carries the name, what sort of division it is, and where it is likely to surprise.
"""
naturalearth_levels() = copy(NATURALEARTH_LEVELS)

"""
    naturalearth_regions(level)

Return every region defined at `level`, sorted by name, with the extent and area of each.

Reads the shipped table, so it costs no download. What comes back is a [`RegionReport`](@ref), the
same thing [`investigate_regions`](@ref) returns: it displays as a table, it iterates and indexes,
and any row of it converts straight into a [`NaturalEarthSpec`](@ref). The names alone are
`[m.name for m in naturalearth_regions(level)]`.

Names are matched case-insensitively wherever they are used, so the spelling here is for display
rather than something to reproduce exactly.

# Arguments

  - `level`: the level to list, as a name (`"ADMIN"`) or a [`NaturalEarthLevel`](@ref). An unknown
    one is an error suggesting the closest matches.
"""
function naturalearth_regions(level::AbstractString)
    return naturalearth_regions(_checklevel(level))
end

function naturalearth_regions(level::NaturalEarthLevel)
    byname = get(_regionindex().bylevel, level.name, nothing)
    matches = RegionMatch[]
    isnothing(byname) || for row in values(byname)
        box = _rowextent(row)
        push!(matches,
              RegionMatch{typeof(box)}(level, row.Name, box,
                                       row.AreaKm2 * km^2,
                                       row.Parts, _rowshare(row), nothing))
    end
    # By name, because this is for browsing: a reader looking for one knows how it is spelled long
    # before they know how big it is.
    sort!(matches, by = m -> m.name)
    return RegionReport(nothing, matches)
end

# ---------------------------------------------------------------------------
# Spatial relations
# ---------------------------------------------------------------------------

"""
    AbstractSpatialRelation

How a named region may relate to something you have - [`Encloses`](@ref), [`Overlaps`](@ref) or
[`Within`](@ref).

A relation **carries the thing being asked about**, so `Encloses(mylayer)` reads as what it means and
there is no argument order to get the wrong way round. It is callable, taking a region's extent and
answering whether the relation holds, which also makes it usable directly as a filter.
"""
abstract type AbstractSpatialRelation end

"""
    Encloses(x)

Regions that completely contain `x`.

Answers *"which named regions is my data inside?"*, and is the relation
[`investigate_regions`](@ref) uses when none is named. It is the only one that accepts a **point**,
since nothing can lie within a point and nothing overlaps one.

# Arguments

  - `x`: what to ask about - a study area or its report, a raster, a layer, a habitat, an ecosystem,
    an `Extents.Extent`, a [`LatLong`](@ref), or a match from an earlier report.
"""
struct Encloses{E} <: AbstractSpatialRelation
    extent::E

    Encloses(x) = new{typeof(_wgsextent(x))}(_wgsextent(x))
end

"""
    Overlaps(x)

Regions that share real ground with `x`.

Answers *"which regions does my data reach into?"*. Sharing only a boundary does not count, which is
the whole difference from a bare intersection test: a region touching your data along an edge
contains none of it, and `boundingbox(..., round = ...)` snaps boxes onto a lattice, so that case is
reachable rather than theoretical.

# Arguments

  - `x`: as [`Encloses`](@ref), but not a point - a point has no area to share.
"""
struct Overlaps{E} <: AbstractSpatialRelation
    extent::E

    function Overlaps(x)
        return new{typeof(_areaextent(x, "Overlaps"))}(_areaextent(x,
                                                                   "Overlaps"))
    end
end

"""
    Within(x)

Regions that lie completely inside `x`.

Answers *"which regions can I simulate in full with the data I have?"* - the converse of
[`Encloses`](@ref), and what walks *down* a hierarchy where that walks up.

# Arguments

  - `x`: as [`Encloses`](@ref), but not a point - no region fits inside a point.
"""
struct Within{E} <: AbstractSpatialRelation
    extent::E

    Within(x) = new{typeof(_areaextent(x, "Within"))}(_areaextent(x, "Within"))
end

# The subject's extent, refused if it has no area.
#
# `Overlaps` and `Within` are both empty for a point - measured, a point does not even overlap
# itself - so an empty report would read as "nothing found" when the truth is that the question is
# malformed. A real grid always has area, so nothing legitimate is refused.
function _areaextent(x, what::AbstractString)
    extent = _wgsextent(x)
    (extent.X[1] == extent.X[2] || extent.Y[1] == extent.Y[2]) &&
        error("`$what` needs something with area, and $(typeof(x)) gives a zero-width extent. " *
              "Only `Encloses` is meaningful for a point: nothing lies within one, and nothing " *
              "overlaps one.")
    return extent
end

# Delegated to `Extents`, whose predicates are the documented API for exactly this.
#
# `Encloses` uses `covers`, NOT `contains`. The two agree on every pair of boxes, which is what makes
# the choice look free - but they differ on a **point**: `contains` is false for a point inside a
# box, `covers` is true. Since a point is the one subject only `Encloses` accepts, `contains` would
# make the commonest query - "which regions enclose this coordinate?" - silently answer nothing.
(r::Encloses)(region) = Extents.covers(region, r.extent)
(r::Overlaps)(region) = Extents.overlaps(region, r.extent)
(r::Within)(region) = Extents.within(region, r.extent)

# The same three questions asked of real geometry rather than of boxes, for `exact = true`.
#
# ArchGDAL has no `covers`, and its `contains` excludes the boundary - the same trap `Extents`
# sets - so containment is written as "no part of the one lies outside the other", which is exact
# and boundary-inclusive. ⚠️ And GEOS `overlaps` means *partial* overlap, excluding containment,
# which is not what `Overlaps` means here: sharing ground is a positive intersection area.
function _exactencloses(region, subject)
    return ArchGDAL.isempty(ArchGDAL.difference(subject,
                                                region))
end

function (r::Encloses)(region::ArchGDAL.AbstractGeometry)
    return _exactencloses(region,
                          _subjectgeometry(r))
end

function (r::Within)(region::ArchGDAL.AbstractGeometry)
    return _exactencloses(_subjectgeometry(r), region)
end

function (r::Overlaps)(region::ArchGDAL.AbstractGeometry)
    return ArchGDAL.geomarea(ArchGDAL.intersection(region, _subjectgeometry(r))) >
           0
end

# The subject as geometry: a rectangle, or a point where the extent has no width. Cached on first
# use, since every refinement in a query asks for the same one.
function _subjectgeometry(r::AbstractSpatialRelation)
    e = r.extent
    w, s, x, n = ustrip(°, e.X[1]), ustrip(°, e.Y[1]), ustrip(°, e.X[2]),
                 ustrip(°, e.Y[2])
    (w == x && s == n) && return ArchGDAL.createpoint(w, s)
    return ArchGDAL.fromWKT("POLYGON (($w $s, $x $s, $x $n, $w $n, $w $s))")
end

function Base.show(io::IO, r::AbstractSpatialRelation)
    return print(io, nameof(typeof(r)), "(",
                 r.extent, ")")
end

# ---------------------------------------------------------------------------
# The query
# ---------------------------------------------------------------------------

"""
    RegionMatch

One named region a query found, and enough of it to act on without looking anything else up.

# Fields

  - `level`: which kind of region it is, as [`NaturalEarthLevel`](@ref).
  - `name`: the source's own spelling of the name.
  - `extent`: its bounding box, in degrees.
  - `area`: the total area of the region, in `km^2`.
  - `parts`: how many separate components it has.
  - `share`: what fraction of the region's area its **largest** component holds. This is the number
    that says whether [`LargestLandmass`](@ref) is a sensible answer for the region at all: New
    Zealand's is 0.56, so asking for its principal landmass silently returns South Island alone, and
    the Solomon Islands' is 0.20.
  - `overlap`: the area of the box it shares with whatever was asked about, in `km^2`, which is what
    an [`Overlaps`](@ref) report orders by. `nothing` for a listing, which asked about nothing.
"""
struct RegionMatch{E}
    level::NaturalEarthLevel
    name::String
    extent::E
    area::typeof(1.0km^2)
    parts::Int
    share::Float64
    overlap::Union{typeof(1.0km^2), Nothing}
end

# A match is a subject in its own right, so a report's row feeds straight back into a new query:
# `Encloses(match)` walks up the hierarchy and `Within(match)` walks down.
_wgsextent(m::RegionMatch) = m.extent

"""
    RegionReport

What [`investigate_regions`](@ref) found: the matching regions, in order, and the question asked.

A report is a container - it iterates, indexes and has a length - so `only(report)` asserts that the
answer was unique, `first(report)` takes the best one and `report[i]` takes a chosen one. Any of
those gives a [`RegionMatch`](@ref), which is what [`NaturalEarthSpec`](@ref) accepts.

# Fields

  - `relation`: the question that was asked, including what it was asked about.
  - `matches`: the regions that answered it, ordered as described in [`investigate_regions`](@ref).
  - `exact`: whether the answer was checked against the regions' real outlines rather than their
    bounding boxes.
  - `refined`: how many candidates had their geometry fetched to answer it. Zero unless `exact`.
"""
struct RegionReport{R <: Union{AbstractSpatialRelation, Nothing}, M}
    relation::R
    matches::Vector{M}
    exact::Bool
    refined::Int
end

# A report that was answered from boxes alone, which is every report but an `exact = true` query.
RegionReport(relation, matches) = RegionReport(relation, matches, false, 0)

Base.length(r::RegionReport) = length(r.matches)
Base.getindex(r::RegionReport, i) = r.matches[i]
Base.iterate(r::RegionReport, s...) = iterate(r.matches, s...)
Base.isempty(r::RegionReport) = isempty(r.matches)
Base.eltype(::Type{RegionReport{R, M}}) where {R, M} = M
Base.lastindex(r::RegionReport) = lastindex(r.matches)

"""
    investigate_regions(x; level = nothing, kind = nothing, limit = 20)
    investigate_regions(relation::AbstractSpatialRelation; ...)

Find the named regions that relate to `x`, as [`investigate_study_area`](@ref) reports on a grid
before one is built.

Given anything with a position - a study area, a raster, a layer, a habitat, an ecosystem, an
`Extents.Extent`, a [`LatLong`](@ref) or an earlier match - this asks which of the 2 444 shipped
regions [`Encloses`](@ref) it. Name a relation instead to ask a different question:
[`Overlaps`](@ref) for regions your data reaches into, [`Within`](@ref) for regions your data covers
entirely.

!!! warning "The answer is about boxes, not outlines"
    Every region is compared by its bounding box, because that is what costs no download. A box can
    be far larger than the ground it names: Chile's spans 43 degrees of longitude because of Easter
    Island, so a query "overlapping" Chile may share no Chilean land at all. Pass `exact = true` to
    check against the real outlines instead, or build the shape with [`NaturalEarthSpec`](@ref).

# Arguments

  - `x`: what to ask about, or a relation carrying it.
  - `level`: restrict to one level, by name. `EcoSISTEM.naturalearth_levels()` lists them.
  - `kind`: restrict to levels of one sort - `:political`, `:statistical`, `:physical` or `:code`.
  - `limit`: how many matches to keep. A continental query can match hundreds, and the ordering puts
    the useful ones first.
  - `exact`: check the surviving candidates against the regions' **real outlines** instead of their
    boxes, which needs the geometry and so downloads. It removes the false positives a box tier
    cannot avoid - Norway's box encloses Edinburgh, its coastline does not - and reaches the 54
    regions that cross the antimeridian, which have no box a query can compare at all.

    Refinement is lazy and in box order, stopping as soon as the answer cannot change: refining only
    ever removes a match or shrinks its overlap, so a confirmed `limit` cannot be displaced by
    anything later. The report says how many regions it had to fetch.
"""
function investigate_regions(relation::AbstractSpatialRelation; level = nothing,
                             kind = nothing, limit::Integer = 20,
                             exact::Bool = false)
    wanted = isnothing(level) ? nothing : _checklevel(level).name
    candidates = RegionMatch[]
    for l in NATURALEARTH_LEVELS
        isnothing(wanted) || l.name == wanted || continue
        isnothing(kind) || l.kind === kind || continue
        byname = get(_regionindex().bylevel, l.name, nothing)
        isnothing(byname) && continue
        for row in values(byname)
            # A wrapping row has no box that is a single interval, so the box tier cannot compare it
            # at all - it is skipped there, as `boundingbox` skips the same rows. Geometry has no
            # such problem, Natural Earth having split its polygons at the date line, so `exact`
            # reaches 54 regions that are otherwise invisible.
            box = _rowextent(row)
            if row.Wraps
                exact || continue
            else
                relation(box) || continue
            end
            # A region matches itself under `Encloses` and `Within`, both being reflexive. Excluded
            # by identity rather than by extent: distinct regions legitimately share a box, and
            # dropping those would delete real answers.
            _issubject(relation, l, row) && continue
            push!(candidates,
                  RegionMatch{typeof(box)}(l, row.Name, box, row.AreaKm2 * km^2,
                                           row.Parts, _rowshare(row),
                                           _overlaparea(box, relation.extent,
                                                        row.Wraps)))
        end
    end
    _ordermatches!(candidates, relation)
    exact ||
        return RegionReport(relation,
                            candidates[1:min(limit, length(candidates))])
    return _refineexactly(relation, candidates, limit)
end

# Refine box matches against real geometry, lazily and in box order, stopping as soon as the answer
# cannot change.
#
# Refinement only ever **removes** a match or **shrinks** its overlap, never the reverse, which is
# what makes an early stop exact rather than an approximation:
#
#   - `Encloses` and `Within` order by the region's own area, which refinement does not touch. Once
#     `limit` candidates have survived, nothing later in the order can displace them.
#   - `Overlaps` orders by shared area, which refinement does change - but only downwards, since the
#     true shared ground lies inside both boxes. So a box overlap bounds the true one, and the scan
#     can stop once every confirmed match beats the best any remaining candidate could reach.
#
# Without this a continental `Overlaps` query would fetch 617 regions' geometry to report 20.
function _refineexactly(relation, candidates, limit)
    confirmed = RegionMatch[]
    refined = 0
    for m in candidates
        length(confirmed) >= limit && _canstop(relation, confirmed, m, limit) &&
            break
        refined += 1
        geom = _regiongeometryof(m)
        relation(geom) || continue
        push!(confirmed, _withexactoverlap(relation, m, geom))
        _ordermatches!(confirmed, relation)
    end
    return RegionReport(relation, confirmed[1:min(limit, length(confirmed))],
                        true, refined)
end

# Whether the scan can stop. For the two area-ordered relations the order is fixed, so having enough
# is enough; for `Overlaps` the weakest confirmed match must already beat the best the next
# candidate could possibly turn out to be, which is its box overlap.
_canstop(::Encloses, confirmed, next, limit) = true
_canstop(::Within, confirmed, next, limit) = true

function _canstop(::Overlaps, confirmed, next, limit)
    isnothing(next.overlap) && return false
    return confirmed[limit].overlap >= next.overlap
end

# The region's own outline, in WGS84, merged into one geometry.
#
# Deliberately not `_dissolve`: that merges, splits into components and measures each one's area on
# an equal-area projection, and a predicate reads none of that. On a continent it is a transform per
# component - 858 of them for Asia - for an answer thrown away.
function _regiongeometryof(m::RegionMatch)
    return _mergegeoms(_selectfeatures(_findlevel(m.level.name), m.name))
end

# The match with its overlap replaced by the real shared area, which is what `Overlaps` orders by.
# The other two order by the region's own area, so the intersection is computed only where it will
# be read - it costs a geometry operation per candidate.
_withexactoverlap(relation, m::RegionMatch, geom) = m

function _withexactoverlap(r::Overlaps, m::RegionMatch{E}, geom) where {E}
    shared = ArchGDAL.intersection(geom, _subjectgeometry(r))
    return RegionMatch{E}(m.level, m.name, m.extent, m.area, m.parts, m.share,
                          _equalarea(shared))
end

function investigate_regions(x; kw...)
    return investigate_regions(Encloses(x); kw...)
end

# What fraction of a region's area its largest component holds, straight from the shipped columns -
# nothing extra is stored for it. A low share is the warning that `LargestLandmass()` will throw most
# of the region away.
_rowshare(row) = row.Part1Area / row.AreaKm2

# The box a table row states, as an `Extents.Extent` of degrees.
function _rowextent(row)
    return Extents.Extent(Y = (row.South * °, row.North * °),
                          X = (row.West * °, row.East * °))
end

# Whether this row IS what was asked about, so that a region does not answer a question about itself.
function _issubject(relation, level, row)
    e = relation.extent
    return _rowextent(row) == e
end

# The area the two boxes share, as a true solid angle rather than a product of degree spans: a
# degree of longitude narrows towards the poles, and ordering by a flat product would rank high
# latitudes above equatorial ones for the same ground.
function _overlaparea(a::Extents.Extent, b::Extents.Extent)
    south = max(a.Y[1], b.Y[1])
    north = min(a.Y[2], b.Y[2])
    west = max(a.X[1], b.X[1])
    east = min(a.X[2], b.X[2])
    (north < south || east < west) && return 0.0km^2
    return _sphericalboxarea(south, north, west, east)
end

# A wrapping box is two intervals - west to the date line, and the date line to east - so it is
# comparable after all, as the sum of its halves.
#
# Worth the few lines: without it a wrapping candidate has no overlap to report, and the `exact`
# scan's stopping rule needs an upper bound on every remaining candidate. One unbounded candidate
# would force it to refine the rest.
function _overlaparea(box::Extents.Extent, subject::Extents.Extent, wraps::Bool)
    wraps || return _overlaparea(box, subject)
    west = Extents.Extent(Y = box.Y, X = (box.X[1], 180.0°))
    east = Extents.Extent(Y = box.Y, X = (-180.0°, box.X[2]))
    return _overlaparea(west, subject) + _overlaparea(east, subject)
end

# The area of a lat/long box on the sphere: proportional to the longitude span times the difference
# of the sines of the latitudes.
function _sphericalboxarea(south, north, west, east)
    r = 6371.0087714km
    dlong = ustrip(NoUnits, uconvert(°, east - west) / 1° * (pi / 180))
    return r^2 * dlong * (sin(uconvert(°, north)) - sin(uconvert(°, south)))
end

# Smallest enclosing first for `Encloses`, since the tightest region containing your data is the one
# you almost always want. The other two have no such natural "best", so they order by how much
# ground they share with the subject, largest first.
_ordermatches!(matches, ::Encloses) = sort!(matches, by = m -> m.area)
_ordermatches!(matches, ::Overlaps) = sort!(matches, by = m -> -m.overlap)
_ordermatches!(matches, ::Within) = sort!(matches, by = m -> -m.area)

function Base.show(io::IO, m::RegionMatch)
    return print(io, "RegionMatch(", m.level.name, " \"", m.name, "\")")
end

function Base.show(io::IO, ::MIME"text/plain", m::RegionMatch)
    println(io, "RegionMatch \"", m.name, "\" at level ", m.level.name, " (",
            m.level.kind, ")")
    println(io, "  extent: ", m.extent)
    return print(io, "  area: ", round(typeof(1.0km^2), m.area, digits = 0),
                 " in ",
                 m.parts, " part", m.parts == 1 ? "" : "s")
end

function Base.show(io::IO, r::RegionReport)
    return print(io, "RegionReport(", length(r), " matches)")
end

function Base.show(io::IO, ::MIME"text/plain", r::RegionReport)
    isempty(r) && return print(io, _emptydescription(r.relation))
    println(io, length(r), " named region", length(r) == 1 ? "" : "s", " ",
            _headline(r.relation, length(r) != 1), ", ",
            _orderdescription(r.relation), ":")
    # The level column is dropped where every row shares one, which is what a listing is: repeating
    # it down the page would say nothing and crowd out the extent.
    levels = unique(m -> m.level.name, r.matches)
    showlevel = length(levels) > 1
    lw = showlevel ? maximum(m -> length(m.level.name), r.matches) : 0
    nw = min(30, maximum(m -> length(m.name), r.matches))
    showlevel && print(io, "  ", rpad("level", lw))
    println(io, "  ", rpad("name", nw), "       W       S       E       N",
            lpad("area/km2", 12), lpad("parts", 6), lpad("largest", 8))
    for m in r.matches
        showlevel && print(io, "  ", rpad(m.level.name, lw))
        print(io, "  ", rpad(first(m.name, nw), nw))
        print(io, "  ", _boxcolumns(m.extent))
        # A whole number of square kilometres, never scientific: a column of areas is read by
        # comparing them, and `4.45995e6` beside `830486.0` cannot be.
        print(io, lpad(round(Int, ustrip(km^2, m.area)), 12))
        print(io, lpad(m.parts, 6))
        # The largest component's share, shown only where it is a warning: a region whose principal
        # landmass is most of it needs no annotation, and a column of "100%" would bury the cases
        # that matter.
        println(io,
                m.share < 0.9 ? lpad(string(round(Int, 100m.share), "%"), 7) :
                "")
    end
    isnothing(r.relation) && return
    # Said every time rather than left to the docstring, and inverted when it no longer applies: a
    # box can be far larger than the ground it names, and a reader who does not know that will trust
    # the list too far.
    r.exact &&
        return print(io, "\nChecked against the regions' real outlines, ",
                     r.refined,
                     " of which were fetched to answer this.")
    return print(io,
                 "\nCompared by bounding box, which costs no download but is coarse - Chile's box " *
                 "\nspans 43 degrees because of Easter Island. Pass `exact = true` to check " *
                 "against\nthe real outlines instead.")
end

# A region's box as four fixed-width columns of degrees, so a column of them lines up and can be
# read down.
function _boxcolumns(e)
    return string(lpad(round(ustrip(°, e.X[1]), digits = 2), 8),
                  lpad(round(ustrip(°, e.Y[1]), digits = 2), 8),
                  lpad(round(ustrip(°, e.X[2]), digits = 2), 8),
                  lpad(round(ustrip(°, e.Y[2]), digits = 2), 8))
end

_orderdescription(::Nothing) = "by name"
_orderdescription(::Encloses) = "smallest first"
_orderdescription(::Overlaps) = "most overlap first"
_orderdescription(::Within) = "largest first"

# The relation as a verb agreeing with its subject, so a report reads as a sentence rather than as a
# type name spelled in lower case.
# What the report is a report *of*. A listing asked no question, so it names the level instead of a
# relation.
_headline(rel, plural) = _relationverb(rel, plural) * " it"
_headline(::Nothing, plural) = "at this level"

_emptydescription(rel) = "No named region " * _relationverb(rel, false) * " it."
_emptydescription(::Nothing) = "No named regions at this level."

_relationverb(::Encloses, plural) = plural ? "enclose" : "encloses"
_relationverb(::Overlaps, plural) = plural ? "overlap" : "overlaps"
_relationverb(::Within, plural) = plural ? "lie within" : "lies within"

# Region names close enough to `lower` to be worth offering, by edit distance. Only ever reached on
# an error path, so walking every name is fine; capped so a vague query cannot print a page.
function _nearnames(lower::AbstractString, limit::Integer = 4)
    tol = length(lower) <= 4 ? 1 : 2
    hits = String[]
    for (key, levels) in _regionindex().bynames
        abs(length(key) - length(lower)) > tol && continue
        _editdistance(key, lower) <= tol || continue
        row = _regionrow(levels[1], key)
        push!(hits, "\"$(row.Name)\" (at $(levels[1]))")
        length(hits) >= limit && break
    end
    return sort!(hits)
end

# Levenshtein distance, iterative with one row of state.
function _editdistance(a::AbstractString, b::AbstractString)
    prev = collect(0:length(b))
    row = similar(prev)
    for (i, ca) in enumerate(a)
        row[1] = i
        for (j, cb) in enumerate(b)
            row[j + 1] = min(prev[j + 1] + 1, row[j] + 1, prev[j] + (ca != cb))
        end
        prev, row = row, prev
    end
    return prev[end]
end
