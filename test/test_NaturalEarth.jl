# SPDX-License-Identifier: LGPL-3.0-or-later

module TestNaturalEarth

using EcoSISTEM
using EcoSISTEM: NaturalEarthLevel, NATURALEARTH_LEVELS
using EcoSISTEM: _findlevel, _nesource, _dissolve, _coverageof, _equalarea
using EcoSISTEM: _regionbox, _regionindex, _regionspath, NPARTS
using CSV
using Unitful
using Unitful.DefaultSymbols
import ArchGDAL
using Test

# An axis-aligned lat/long rectangle, the fixture every geometry test here is built from.
function _box(w, s, e, n)
    return ArchGDAL.fromWKT("POLYGON (($w $s, $e $s, $e $n, $w $n, $w $s))")
end

@testset "coverages say how much of a name to take" begin
    @test AllTerritories() isa EcoSISTEM.AbstractCoverage
    @test LargestLandmass() isa EcoSISTEM.AbstractCoverage

    # The default is the single largest landmass, which is what turns a country into its mainland.
    @test LargestLandmass().count == 1
    @test LargestLandmass(count = 3).count == 3

    # The count must be named. A bare number would attach to a singular type name and so could be
    # read as an ordinal - "the third largest" rather than "the largest three".
    @test_throws MethodError LargestLandmass(3)

    # Zero components is not a coverage, and is refused at construction rather than returning an
    # empty selection that a caller would have to recognise.
    @test_throws ArgumentError LargestLandmass(count = 0)
    @test_throws ArgumentError LargestLandmass(count = -1)

    # These nest inside other values, so both need a compact `show` rather than the default, which
    # prints the whole type signature.
    @test sprint(show, AllTerritories()) == "AllTerritories()"
    @test sprint(show, LargestLandmass(count = 2)) ==
          "LargestLandmass(count = 2)"
end

@testset "the level table names where each kind of region is defined" begin
    @test !isempty(NATURALEARTH_LEVELS)
    @test all(l -> l isa NaturalEarthLevel, NATURALEARTH_LEVELS)

    # Nine cultural levels - the attributes worth selecting on across the four admin-0 files - and
    # one per landform class of the physical file. Named rather than counted, so that adding a level
    # has to be looked at and removing one shrinks the list instead of passing silently.
    cultural = filter(l -> l.category == "cultural", NATURALEARTH_LEVELS)
    @test [l.name for l in cultural] ==
          ["CONTINENT", "REGION_UN", "SUBREGION", "REGION_WB", "SOVEREIGNT",
        "ADMIN", "GEOUNIT", "SUBUNIT", "ISO_A3"]

    physical = filter(l -> l.category == "physical", NATURALEARTH_LEVELS)
    @test length(physical) == 22
    @test all(l -> startswith(l.name, "Physical "), physical)
    @test all(l -> l.kind === :physical, physical)

    # Every physical level reads the one shared layer, filtered to its own landform class.
    @test all(l -> l.field == "NAME", physical)
    @test all(l -> l.within isa Pair && l.within.first == "FEATURECLA",
              physical)
    @test all(l -> l.within.second == l.name[(length("Physical ") + 1):end],
              physical)

    # "Dragons-be-here" is Null Island at (0, 0), a cartographers' joke rather than a place.
    @test !any(l -> occursin("Dragons", l.name), NATURALEARTH_LEVELS)

    @test all(l -> l.category in ("cultural", "physical"), NATURALEARTH_LEVELS)
    @test all(l -> l.kind in (:political, :statistical, :physical, :code),
              NATURALEARTH_LEVELS)
    @test all(l -> !isempty(l.description), NATURALEARTH_LEVELS)
end

@testset "a level is found by name, case-insensitively but unambiguously" begin
    @test _findlevel("ADMIN").field == "ADMIN"
    @test _findlevel("admin").field == "ADMIN"
    @test isnothing(_findlevel("NOT_A_LEVEL"))

    # The lookup ignores case because the physical file's own region names are inconsistently cased,
    # so no two levels may differ only in case. The cultural `CONTINENT` and the physical landform
    # class `Continent` are exactly that pair, which is why the physical levels are prefixed: without
    # it, asking for the coastline of Europe silently returned the list of countries assigned to it,
    # a selection spanning the whole globe rather than stopping at the Urals.
    lowered = [lowercase(l.name) for l in NATURALEARTH_LEVELS]
    @test length(unique(lowered)) == length(lowered)

    @test _findlevel("CONTINENT").category == "cultural"
    @test _findlevel("Physical Continent").category == "physical"
    @test _findlevel("CONTINENT") !== _findlevel("Physical Continent")
end

@testset "each level resolves to its published file" begin
    # A level's download is derived from its own fields, so a mistyped dataset cannot point at a
    # file that happens to exist for another level.
    @test _nesource(_findlevel("ADMIN")).url ==
          "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"
    @test _nesource(_findlevel("Physical Island")).url ==
          "https://naciscdn.org/naturalearth/10m/physical/ne_10m_geography_regions_polys.zip"

    # The finer admin-0 splits are separate files, not separate attributes of one: a country's
    # overseas departments are dropped by reading `map_units`, not by selecting differently.
    @test occursin("map_units", _nesource(_findlevel("GEOUNIT")).url)
    @test occursin("map_subunits", _nesource(_findlevel("SUBUNIT")).url)
    @test occursin("sovereignty", _nesource(_findlevel("SOVEREIGNT")).url)

    # One resolution ships, so a box and the shape it describes come from the same geometry.
    @test all(NATURALEARTH_LEVELS) do l
        return occursin("/10m/", _nesource(l).url) &&
               occursin("ne_10m_", _nesource(l).url)
    end
end

@testset "areas are measured on an equal-area projection, in the right axis order" begin
    # A 10-degree by 1-degree box at 50 degrees north. Were the transform built so that GDAL read
    # each pair latitude-first, the same polygon would be measured at 5 degrees north and 50 east -
    # real ground, and no error, but 55% larger. That is the whole reason `_gdalcrs` exists.
    here = _equalarea(_box(0, 50, 10, 51))
    swapped = _equalarea(_box(50, 0, 51, 10))
    @test here ≈ 78922.19km^2 rtol=1e-4
    @test swapped ≈ 122483.23km^2 rtol=1e-4
    @test !isapprox(here, swapped, rtol = 0.1)

    # Equal-area means equal-area: the same angular box covers less ground the further it is from
    # the equator, and a degree square at the equator is about 12 300 square kilometres.
    @test _equalarea(_box(0, 0, 1, 1)) ≈ 12308.46km^2 rtol=1e-4
    @test _equalarea(_box(0, 60, 1, 61)) < _equalarea(_box(0, 0, 1, 1))
end

@testset "dissolving merges neighbours, then splits what is left apart" begin
    # Two polygons sharing an edge are one piece of ground, and must come back as one component -
    # this is what makes the largest component of a continent its mainland rather than its largest
    # country.
    touching = _dissolve([_box(0, 0, 1, 1), _box(1, 0, 2, 1)])
    @test length(touching) == 1
    @test touching[1].area ≈ 2 * _equalarea(_box(0, 0, 1, 1))
    @test touching[1].envelope.MinX ≈ 0 && touching[1].envelope.MaxX ≈ 2

    # Separated ones stay separate, ordered largest first so that a coverage can take a prefix.
    apart = _dissolve([_box(0, 0, 1, 1), _box(10, 10, 13, 13)])
    @test length(apart) == 2
    @test issorted([p.area for p in apart], rev = true)
    @test apart[1].envelope.MinX ≈ 10
    @test apart[2].envelope.MinX ≈ 0

    # An envelope belongs to its own component, not to the selection: this is what lets a bounding
    # box exclude a remote territory rather than merely reporting it separately.
    @test apart[1].envelope.MaxX ≈ 13
    @test apart[2].envelope.MaxX ≈ 1

    @test isempty(_dissolve(ArchGDAL.IGeometry[]))
end

@testset "a coverage takes the whole selection or its largest components" begin
    parts = _dissolve([
                          _box(0, 0, 1, 1),
                          _box(10, 10, 13, 13),
                          _box(20, 20, 22, 22)
                      ])
    @test length(parts) == 3

    @test length(_coverageof(parts, AllTerritories())) == 3
    @test _coverageof(parts, AllTerritories()) == parts

    @test length(_coverageof(parts, LargestLandmass())) == 1
    @test _coverageof(parts, LargestLandmass())[1] === parts[1]
    @test length(_coverageof(parts, LargestLandmass(count = 2))) == 2

    # Asking for more components than exist returns them all rather than raising: how many pieces a
    # name has is a property of the data, not something a caller can be expected to know.
    @test length(_coverageof(parts, LargestLandmass(count = 99))) == 3

    @test isempty(_coverageof(_dissolve(ArchGDAL.IGeometry[]),
                              LargestLandmass()))
end

@testset "the antimeridian is read off the widest gap, not assumed" begin
    # A component's envelope is all `_regionbox` sees, so synthetic ones pin the rule exactly.
    part(w, e) = (geometry = nothing,
                  envelope = ArchGDAL.GDAL.OGREnvelope(w, e, 0.0, 10.0),
                  area = 1.0km^2)

    # An ordinary region: the widest empty stretch is the one outside it, so the box is plain.
    plain = _regionbox([part(0.0, 10.0), part(20.0, 30.0)])
    @test (plain.west, plain.east) == (0.0, 30.0)
    @test !plain.wraps

    # One straddling the date line arrives as components either side. Smallest-to-largest would say
    # -180 to 180, which is true and useless; the tight box runs east from 170 through the line.
    wrapped = _regionbox([part(170.0, 180.0), part(-180.0, -170.0)])
    @test (wrapped.west, wrapped.east) == (170.0, -170.0)
    @test wrapped.wraps
    @test wrapped.west > wrapped.east      # how RFC 7946 writes a wrapping box

    # Latitude has no seam and is never reinterpreted.
    @test (plain.south, plain.north) == (0.0, 10.0)
    @test isempty(_dissolve(ArchGDAL.IGeometry[])) &&
          !_regionbox(EcoSISTEM._RegionPart[]).wraps
end

@testset "the shipped region table is well formed" begin
    rows = collect(CSV.File(_regionspath()))
    @test length(rows) > 2000

    # A region is identified by its level and its name together, so that pair must be unique: two
    # rows sharing it would make a lookup silently order-dependent.
    ids = [(r.Level, lowercase(r.Name)) for r in rows]
    @test length(unique(ids)) == length(ids)

    # The table and the level table must describe the same set of levels. A level added to the code
    # and not regenerated, or dropped from the code and left in the table, is a silent gap.
    inconst = Set(l.name for l in EcoSISTEM.NATURALEARTH_LEVELS)
    intable = Set(r.Level for r in rows)
    @test intable ⊆ inconst
    @test inconst ⊆ intable

    # Each rule collects the rows that break it rather than asserting per row: 2 444 rows would
    # otherwise contribute tens of thousands of assertions, and a failure that names the offender is
    # more use than one that only counts.
    label(r) = r.Level * "/" * r.Name
    @test isempty([label(r)
                   for r in rows
                   if !(-90 <= r.South <= 90 && -90 <= r.North <= 90 &&
                        -180 <= r.West <= 180 && -180 <= r.East <= 180)])
    # Latitude has no seam, so south is always below north.
    @test isempty([label(r) for r in rows if r.South > r.North])
    # `West > East` is exactly the wrapping case, and `Wraps` must say so rather than leaving a
    # consumer to infer it from the numbers.
    @test isempty([label(r) for r in rows if (r.West > r.East) != r.Wraps])
    @test isempty([label(r) for r in rows if r.Parts < 1 || r.AreaKm2 <= 0])
end

@testset "a table row describes its own largest components" begin
    rows = collect(CSV.File(_regionspath()))
    label(r) = r.Level * "/" * r.Name
    described(r) = min(NPARTS, r.Parts)
    areas(r) = [r[Symbol("Part$(i)Area")] for i in 1:described(r)]

    # Components are written largest first, which is what lets `LargestLandmass` take a prefix.
    @test isempty([label(r) for r in rows if !issorted(areas(r), rev = true)])
    @test isempty([label(r) for r in rows if !all(>(0), areas(r))])

    # They are components *of* the region, so each sits inside the region's own box - except where
    # that box wraps, when "inside" is not a comparison of numbers.
    @test isempty([label(r)
                   for r in rows
                   if !r.Wraps &&
                      any(i -> r[Symbol("Part$(i)West")] < r.West ||
                               r[Symbol("Part$(i)East")] > r.East,
                          1:described(r))])

    # Beyond the region's own component count the columns are blank, not zero: a zero would read as
    # a real component of no area.
    @test isempty([label(r)
                   for r in rows
                   if any(i -> !ismissing(r[Symbol("Part$(i)Area")]),
                          (described(r) + 1):NPARTS)])

    # The table would be pointless if every region had one component; it earns its width.
    @test count(r -> r.Parts > 1, rows) > 500
end

end
