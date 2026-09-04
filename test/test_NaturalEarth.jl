# SPDX-License-Identifier: LGPL-3.0-or-later

module TestNaturalEarth

using EcoSISTEM
using EcoSISTEM: NaturalEarthLevel, NATURALEARTH_LEVELS
using EcoSISTEM: _findlevel, _nesource, _dissolve, _coverageof, _equalarea
using EcoSISTEM: _regionbox, _regionindex, _regionspath, NPARTS
using EcoSISTEM: _checklevel, _indegrees, _asraster, _mergegeoms
using EcoSISTEM: _wgsextent, _overlaparea, investigate_regions,
                 naturalearth_regions
using EcoSISTEM: naturalearth_levels, RegionMatch, RegionReport
using EcoSISTEM: _shapecomponents, _applyshapeop, _naturalearthgeoms,
                 _preparemask
import Rasters
import Extents
using CSV
using EcoSISTEM.Units
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
        "ADMIN", "GEOUNIT", "SUBUNIT", "ISO_A3_EH"]

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

@testset "a region reports the share its largest component holds" begin
    admin = naturalearth_regions("ADMIN")
    share(n) = only(m for m in admin if m.name == n).share

    # The number that says whether `LargestLandmass()` is a sensible answer at all. New Zealand's
    # principal landmass is barely over half of it, so asking for that alone silently returns South
    # Island - which the share makes visible where the component count does not.
    @test share("New Zealand") ≈ 0.56 atol=0.01
    @test share("Solomon Islands") < 0.25
    @test share("United Kingdom") > 0.85
    @test all(m -> 0 < m.share <= 1, admin)
    # A single-component region is all of itself.
    @test share("Vatican") ≈ 1.0

    # Derived from the shipped columns rather than stored, so the table did not have to grow for it.
    @test !occursin("Share", first(readlines(_regionspath())))

    # Shown only where it is a warning: a column of "100%" would bury the cases that matter.
    shown = sprint(show, MIME"text/plain"(),
                   EcoSISTEM.RegionReport(nothing,
                                          [m
                                           for m in admin
                                           if m.name in
                                              ("New Zealand", "Vatican")]))
    @test occursin("56%", shown)
    @test !occursin("100%", shown)
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
          !_regionbox(EcoSISTEM._ShapeComponent[]).wraps
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

# Downloading Natural Earth's polygons is a network cost the shipped table exists to avoid, so the
# geometry tests are guarded the way the heavy raster reads are; `ECOSISTEM_HEAVY_DATA=true` forces
# them on. Everything above this point - the table, the coverages, the box arithmetic - runs anywhere.
function geometrytests()
    return haskey(ENV, "ECOSISTEM_HEAVY_DATA") ?
           ENV["ECOSISTEM_HEAVY_DATA"] == "true" : !haskey(ENV, "RUNNER_OS")
end

@testset "a region spec resolves its name when it is written, not when it is built" begin
    # Validation is against the shipped table, so it costs no download and lands at the call site.
    @test NaturalEarthSpec("Scotland") isa EcoSISTEM.AbstractShapeSpec
    @test_throws ErrorException NaturalEarthSpec("Atlantis")
    @test_throws ErrorException NaturalEarthSpec("Scotland",
                                                 level = "CONTINENT")
    # The same rule `boundingbox` uses, which is what makes its box describe this spec's shape.
    @test_throws ErrorException NaturalEarthSpec("Africa")
    @test NaturalEarthSpec("Africa", level = "CONTINENT").level == "CONTINENT"

    # The source's own spelling is stored, not the caller's.
    @test NaturalEarthSpec("great britain", level = "Physical Island").name ==
          "GREAT BRITAIN"

    # `show` has to round-trip: these nest inside a `StudyArea`, and a coverage silently dropped from
    # the display would make two different specs print identically.
    @test sprint(show, NaturalEarthSpec("Scotland")) ==
          "NaturalEarthSpec(\"Scotland\", level = \"GEOUNIT\")"
    # A coverage the caller had to write is shown; the default is not, so that what `show` prints is
    # what would reproduce the spec.
    for cov in (LargestLandmass(), LargestLandmass(count = 2), LandmassesAbove(1km^2))
        @test occursin(sprint(show, cov),
                       sprint(show,
                              NaturalEarthSpec("Scotland", coverage = cov)))
    end
    @test !occursin("coverage",
                    sprint(show,
                           NaturalEarthSpec("Scotland",
                                            coverage = AllTerritories())))
    @test occursin("outline = false",
                   sprint(show, NaturalEarthSpec("Scotland", outline = false)))
end

@testset "LandmassesAbove keeps the components that clear a threshold" begin
    @test LandmassesAbove(1km^2) isa EcoSISTEM.AbstractCoverage
    @test LandmassesAbove(5percent) isa EcoSISTEM.AbstractCoverage
    # An area or a percentage, and nothing else. A **bare number** is refused although `0.05` and
    # `5percent` are the same quantity: only one of them says which it means when read beside
    # `1km^2`.
    @test_throws ArgumentError LandmassesAbove(0.05)
    @test_throws ArgumentError LandmassesAbove(5km)
    @test_throws ArgumentError LandmassesAbove(0km^2)
    @test_throws ArgumentError LandmassesAbove(-1percent)
    # ...and a share larger than the whole region can never be cleared.
    @test_throws ArgumentError LandmassesAbove(150percent)

    # It filters a component list exactly as the other coverages take a prefix of one.
    parts = _dissolve([
                          _box(0, 0, 3, 3),
                          _box(10, 10, 11, 11),
                          _box(20, 20, 20.1, 20.1)
                      ])
    @test length(parts) == 3
    kept = _coverageof(parts, LandmassesAbove(parts[2].area))
    @test length(kept) == 2
    @test all(p -> p.area >= parts[2].area, kept)
    @test isempty(_coverageof(parts, LandmassesAbove(1e9km^2)))

    # A share is of the region's OWN total, so the same threshold means different areas for
    # different regions - which is the point of offering it, since what counts as a speck depends on
    # how big the region is.
    total = sum(p -> p.area, parts)
    @test EcoSISTEM._thresholdarea(LandmassesAbove(50percent), total) ≈
          total / 2
    @test EcoSISTEM._thresholdarea(LandmassesAbove(1km^2), total) == 1km^2
    # The largest part here is 9/(9+1+0.01) of the whole, so a 50% share keeps only it.
    @test length(_coverageof(parts, LandmassesAbove(50percent))) == 1
    @test length(_coverageof(parts, LandmassesAbove(1percent))) == 2
    @test isempty(_coverageof(EcoSISTEM._ShapeComponent[],
                              LandmassesAbove(5percent)))

    # The shipped table records only the largest few components' sizes, so it cannot answer this and
    # says so rather than guessing from what it has.
    @test_throws ErrorException EcoSISTEM.boundingbox("United Kingdom",
                                                      level = "ADMIN",
                                                      coverage = LandmassesAbove(1km^2))
end

@testset "a shape file is a shape spec, so it can join a combination" begin
    # A vector file and a named region are both geometry, so they combine on equal terms. Written to
    # a temporary GeoJSON rather than downloaded, so this needs no network.
    path = joinpath(mktempdir(), "square.geojson")
    write(path,
          """{"type":"FeatureCollection","features":[{"type":"Feature","properties":{},
             "geometry":{"type":"Polygon","coordinates":[[[-5,55],[-3,55],[-3,57],[-5,57],
             [-5,55]]]}}]}""")
    file = ShapeSpec(path)
    @test file isa EcoSISTEM.AbstractShapeSpec

    parts = _shapecomponents(file)
    @test length(parts) == 1
    @test parts[1].envelope.MinX ≈ -5.0 && parts[1].envelope.MaxX ≈ -3.0

    # ...and it is accepted as a member, which is the gap this closes.
    combined = ConstructedShapeSpec(ShapeUnion(), file,
                                    NaturalEarthSpec("Scotland",
                                                     coverage = LargestLandmass()))
    @test combined isa EcoSISTEM.AbstractShapeSpec
    # A `ConstructedShapeSpec` is itself a shape spec, so it can nest.
    @test ConstructedShapeSpec(ShapeConvexHull(), combined) isa
          EcoSISTEM.AbstractShapeSpec
end

@testset "an operation may be a function, mirroring the raster side" begin
    a, b = _box(0, 0, 2, 2), _box(1, 0, 3, 2)
    # The escape hatch: handed one geometry per member, returning one geometry. This is what
    # `ConstructedRasterSpec`'s `combine` is for rasters, and what the three named operations alone
    # could not express - ArchGDAL's `symdifference` among them.
    sym = _applyshapeop((x, y) -> ArchGDAL.symdifference(x, y), [a, b])
    @test ArchGDAL.geomarea(sym) ≈ 4.0

    ne = NaturalEarthSpec("Scotland")
    @test ConstructedShapeSpec((x, y) -> x, ne, ne) isa
          EcoSISTEM.AbstractShapeSpec
    # A function says nothing about how many shapes it wants, so one is enough.
    @test ConstructedShapeSpec(gs -> gs, ne) isa EcoSISTEM.AbstractShapeSpec
end

@testset "the unary operations transform one shape" begin
    sq = _box(0, 0, 2, 2)
    # A buffer grows the outline; a negative one shrinks it.
    @test ArchGDAL.geomarea(_applyshapeop(ShapeBuffer(0.5°), [sq])) > 4.0
    @test ArchGDAL.geomarea(_applyshapeop(ShapeBuffer(-0.5°), [sq])) < 4.0
    # A hull of a convex shape is itself; of a scattered pair, the area between them too.
    @test ArchGDAL.geomarea(_applyshapeop(ShapeConvexHull(), [sq])) ≈ 4.0
    @test ArchGDAL.geomarea(_applyshapeop(ShapeConvexHull(),
                                          [_mergegeoms([sq,
                                                           _box(10, 10, 12, 12)
                                                       ])])) > 8.0
    @test ArchGDAL.geomarea(_applyshapeop(ShapeSimplify(0.1°), [sq])) ≈ 4.0

    # A length is converted at the equator, where a degree of longitude is widest, so a buffer never
    # under-reaches what was asked for.
    @test _indegrees(0.5°) ≈ 0.5
    @test _indegrees(50km) ≈ 50 / 111.32 rtol=1e-3
    @test _indegrees(0.25) == 0.25

    # A distance is a length or an angle; anything else is refused where it was written.
    @test_throws ArgumentError ShapeBuffer(5km^2)
    @test_throws ArgumentError ShapeSimplify(5km^2)
    @test_throws ArgumentError ShapeSimplify(-1°)

    # The transforms take exactly one shape, the set operations at least two, and the constructor
    # refuses a wrong count rather than failing inside a reduce.
    ne = NaturalEarthSpec("Scotland")
    @test ConstructedShapeSpec(ShapeBuffer(0.1°), ne) isa
          EcoSISTEM.AbstractShapeSpec
    @test_throws ArgumentError ConstructedShapeSpec(ShapeUnion(), ne)
end

@testset "a level that does not exist blames the level, and suggests" begin
    # Otherwise a mistyped level surfaces as "no region named X", which points at the wrong argument.
    @test _checklevel("ADMIN").name == "ADMIN"
    @test_throws ErrorException _checklevel("nonsense")
    # Matched both ways round: the shipped level is `ISO_A3_EH`, and a guess either side finds it.
    @test occursin("ISO_A3_EH", sprint(showerror, try
                                           _checklevel("ISO_A3")
                                       catch e
                                           e
                                       end))
    @test occursin("SUBUNIT", sprint(showerror, try
                                         _checklevel("SUBUNITS")
                                     catch e
                                         e
                                     end))
end

@testset "geometry is not a raster, and says which question was meant" begin
    # A shape has no grid of its own, so it cannot be a `ConstructedRasterSpec` member: choosing a
    # resolution here would fix one before the study grid was decided. The refusal names the two
    # things that do want a grid rather than failing as a `MethodError` on a private function.
    for spec in (ShapeSpec("/tmp/nowhere.shp"), NaturalEarthSpec("Scotland"))
        err = try
            _asraster(spec)
        catch e
            sprint(showerror, e)
        end
        @test occursin("ConstructedShapeSpec", err)
        @test occursin("within", err)
    end

    # The released name still resolves, since renaming it was for the mirror rather than a change.
    @test ConstructedSpec === ConstructedRasterSpec
end

@testset "regions combine as geometry, before any grid exists" begin
    @test ShapeUnion() isa EcoSISTEM.AbstractShapeOperation
    uk = NaturalEarthSpec("United Kingdom", level = "ADMIN")
    @test_throws ArgumentError ConstructedShapeSpec(ShapeUnion(), uk)

    # The operations act on geometry, so synthetic squares pin them without any download.
    a, b = _box(0, 0, 2, 2), _box(1, 0, 3, 2)
    @test ArchGDAL.geomarea(_applyshapeop(ShapeUnion(), [a, b])) ≈ 6.0
    @test ArchGDAL.geomarea(_applyshapeop(ShapeIntersection(), [a, b])) ≈ 2.0
    @test ArchGDAL.geomarea(_applyshapeop(ShapeDifference(), [a, b])) ≈ 2.0
    # Difference takes every later member from the first, so it is not symmetric.
    @test ArchGDAL.geomarea(_applyshapeop(ShapeDifference(), [b, a])) ≈ 2.0
    @test ArchGDAL.isempty(_applyshapeop(ShapeIntersection(),
                                         [a, _box(10, 10, 11, 11)]))

    # An empty result is not a mask of nowhere: `_dissolve` refuses to make a component out of it,
    # which is what stops an extent being reported at the origin.
    @test isempty(_dissolve([_applyshapeop(ShapeIntersection(),
                                           [a, _box(10, 10, 11, 11)])]))
end

@testset "the exact predicates are boundary-inclusive where the box ones are" begin
    poly = ArchGDAL.fromWKT("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
    box = Extents.Extent(Y = (0.0°, 10.0°), X = (0.0°, 10.0°))

    # ArchGDAL has no `covers`, and its `contains` excludes the boundary - the same trap `Extents`
    # sets, and the same one that made `Encloses` silently answer nothing for a point. Containment is
    # written as "no part of the one lies outside the other", which has neither problem.
    @test Encloses(LatLong(5.0°, 5.0°))(poly)
    @test Encloses(LatLong(5.0°, 0.0°))(poly)      # exactly on the edge
    @test !Encloses(LatLong(20.0°, 20.0°))(poly)
    @test !ArchGDAL.contains(poly, ArchGDAL.createpoint(0.0, 5.0))

    inside = Extents.Extent(Y = (2.0°, 4.0°), X = (2.0°, 4.0°))
    straddling = Extents.Extent(Y = (2.0°, 4.0°), X = (8.0°, 14.0°))
    @test Encloses(inside)(poly) && !Encloses(straddling)(poly)
    @test Within(box)(ArchGDAL.fromWKT("POLYGON ((2 2, 4 2, 4 4, 2 4, 2 2))"))
    # `Overlaps` means sharing ground, which includes containment - unlike GEOS `overlaps`, which
    # means partial overlap and excludes it.
    @test Overlaps(inside)(poly) && Overlaps(straddling)(poly)
    @test !Overlaps(Extents.Extent(Y = (50.0°, 51.0°), X = (50.0°, 51.0°)))(poly)
end

@testset "a wrapping box is comparable once split at the date line" begin
    subject = Extents.Extent(Y = (0.0°, 10.0°), X = (170.0°, 179.0°))
    wrapping = Extents.Extent(Y = (0.0°, 10.0°), X = (175.0°, -175.0°))
    # Without splitting, a wrapping candidate has no overlap at all - and the exact scan's stopping
    # rule needs an upper bound on every remaining candidate, so one unbounded candidate would force
    # it to refine the rest.
    @test EcoSISTEM._overlaparea(wrapping, subject, true) > 0.0km^2
    # The two halves are what it sums: 175..180 overlaps the subject, -180..-175 does not.
    west = Extents.Extent(Y = (0.0°, 10.0°), X = (175.0°, 180.0°))
    @test EcoSISTEM._overlaparea(wrapping, subject, true) ≈
          EcoSISTEM._overlaparea(west, subject)
end

if geometrytests()
    @testset "a named region becomes a mask that agrees with its box" begin
        # The property the whole feature exists for: the table's box, which is rounded outwards,
        # must contain the extent of the shape the same name builds.
        for (name, level) in (("Scotland", nothing), ("Africa", "CONTINENT"),
            ("Madagascar", "ADMIN"))
            spec = isnothing(level) ? NaturalEarthSpec(name) :
                   NaturalEarthSpec(name, level = level)
            _, extent = _naturalearthgeoms(spec, Rasters.EPSG(4326))
            box = isnothing(level) ? EcoSISTEM.boundingbox(name) :
                  EcoSISTEM.boundingbox(name, level = level)
            # Both sides stripped to degrees explicitly. A bare number compared against a `°`
            # quantity is read by Unitful as radians, which compares silently and wrongly.
            deg(x) = ustrip(°, x)
            @test deg(minimum(box.X)) <= deg(minimum(extent.X))
            @test deg(maximum(box.X)) >= deg(maximum(extent.X))
            @test deg(minimum(box.Y)) <= deg(minimum(extent.Y))
            @test deg(maximum(box.Y)) >= deg(maximum(extent.Y))
        end
    end

    @testset "the British Isles recipe beats the polygon of that name" begin
        all_ = AllTerritories()
        uk = NaturalEarthSpec("United Kingdom", level = "ADMIN",
                              coverage = all_)
        ie = NaturalEarthSpec("Ireland", level = "ADMIN", coverage = all_)
        im = NaturalEarthSpec("Isle of Man", level = "ADMIN", coverage = all_)
        isles = ConstructedShapeSpec(ShapeUnion(), uk, ie, im)
        _, ext = _naturalearthgeoms(isles, Rasters.EPSG(4326))

        # Natural Earth's own BRITISH ISLES polygon stops at 59.80 and so drops Shetland; the union
        # of the three countries reaches 60.85 and does not. This is the cartographic-outline
        # caveat biting on the first composite anyone would reach for.
        @test ustrip(°, maximum(ext.Y)) > 60.8
        named = NaturalEarthSpec("BRITISH ISLES",
                                 level = "Physical Island group",
                                 coverage = all_)
        _, namedext = _naturalearthgeoms(named, Rasters.EPSG(4326))
        @test ustrip(°, maximum(namedext.Y)) < 59.9

        # ...and dropping components under a square kilometre drops exactly Rockall, which is what
        # moves the western edge from -13.69 to about -10.5.
        trimmed = ConstructedShapeSpec(ShapeUnion(), uk, ie, im,
                                       coverage = LandmassesAbove(1km^2))
        @test length(_shapecomponents(trimmed)) ==
              length(_shapecomponents(isles)) - 1
        _, trimext = _naturalearthgeoms(trimmed, Rasters.EPSG(4326))
        @test ustrip(°, minimum(ext.X)) ≈ -13.69 atol=0.01
        @test ustrip(°, minimum(trimext.X)) ≈ -10.48 atol=0.01
    end

    @testset "a spec selecting no ground says so" begin
        uk = NaturalEarthSpec("United Kingdom", level = "ADMIN")
        # Natural Earth's physical continents are landmass outlines, so EUROPE does not contain the
        # British Isles at all - measured, Paris is inside it and central England is not. An empty
        # intersection must be an error, not a grid with nothing active.
        eur = NaturalEarthSpec("EUROPE", level = "Physical Continent",
                               coverage = AllTerritories())
        @test_throws ErrorException _preparemask(ConstructedShapeSpec(ShapeIntersection(),
                                                                      uk, eur),
                                                 Rasters.EPSG(4326))
        @test_throws ErrorException _preparemask(NaturalEarthSpec("United Kingdom",
                                                                  level = "ADMIN",
                                                                  coverage = LandmassesAbove(1e9km^2)),
                                                 Rasters.EPSG(4326))
    end

    @testset "exact = true checks the real outlines, and drops the false positives" begin
        coarse = investigate_regions(LatLong(55.95°, -3.19°), kind = :political,
                                     limit = 5)
        fine = investigate_regions(LatLong(55.95°, -3.19°), kind = :political,
                                   limit = 5, exact = true)
        @test !coarse.exact && fine.exact
        @test coarse.refined == 0 && fine.refined > 0

        # Norway's *box* encloses Edinburgh - it runs west to Jan Mayen and north to Svalbard - and
        # its coastline does not. That is the whole point of the exact tier.
        @test "Norway" in [m.name for m in coarse]
        @test "Norway" ∉ [m.name for m in fine]

        # ...and every region the exact tier keeps really does contain the point.
        for m in fine
            @test Encloses(LatLong(55.95°, -3.19°))(EcoSISTEM._regiongeometryof(m))
        end

        # A wrapping region has no box that is a single interval, so the box tier cannot see it at
        # all; geometry can, Natural Earth having split its polygons at the date line.
        @test any(m -> m.name == "Europe", fine)
        @test !any(m -> m.name == "Europe", coarse)

        # The report says which tier answered it, rather than leaving the reader to assume.
        @test occursin("real outlines", sprint(show, MIME"text/plain"(), fine))
        @test occursin("bounding box", sprint(show, MIME"text/plain"(), coarse))
    end

    @testset "the exact scan stops as soon as the answer cannot change" begin
        # Refinement only ever removes a match or shrinks its overlap, so a scan in box order can
        # stop early and still be exact. Without that, a continental `Overlaps` query would fetch
        # every one of its 617 box candidates to report 20.
        africa = EcoSISTEM.boundingbox("Africa", level = "CONTINENT")
        all = investigate_regions(Overlaps(africa), kind = :political,
                                  limit = 10^6)
        @test length(all) > 300
        few = investigate_regions(Overlaps(africa), kind = :political,
                                  limit = 5,
                                  exact = true)
        @test length(few) == 5
        @test few.refined < length(all) / 4
        # Ordered by real shared ground, so each is at least as good as the next.
        @test issorted([m.overlap for m in few], rev = true)
    end

    @testset "outline = false gives the box, with every cell in it active" begin
        spec = NaturalEarthSpec("Scotland", outline = false)
        prep = _preparemask(spec, Rasters.EPSG(4326))
        # No payload is how `_preparemask` already says "restrict the extent, mask nothing".
        @test isnothing(prep.payload)
        @test !isnothing(prep.extent)
        outlined = _preparemask(NaturalEarthSpec("Scotland"),
                                Rasters.EPSG(4326))
        @test !isnothing(outlined.payload)
        # Both describe the same ground, so they agree on where it is.
        @test prep.extent == outlined.extent
    end
end

@testset "the listings answer from the table, with no download" begin
    @test length(naturalearth_levels()) == length(EcoSISTEM.NATURALEARTH_LEVELS)
    # A copy, so a caller cannot mutate the package's own table.
    @test naturalearth_levels() !== EcoSISTEM.NATURALEARTH_LEVELS

    admin = naturalearth_regions("ADMIN")
    @test length(admin) == 258
    # A report, so it iterates and indexes, and every row carries the box and area a reader browsing
    # a level actually wants - not just the name.
    @test admin isa RegionReport
    @test first(admin) isa RegionMatch
    @test issorted([m.name for m in admin])
    @test "United Kingdom" in [m.name for m in admin]
    @test all(m -> m.area > 0km^2 && m.parts >= 1, admin)
    # A listing asked no spatial question, so there is no overlap to report.
    @test all(m -> isnothing(m.overlap), admin)
    # ...and a row of it converts straight into a spec, as a query's row does.
    uk = only(m for m in admin if m.name == "United Kingdom")
    @test NaturalEarthSpec(uk).level == "ADMIN"

    # `_checklevel`'s error promises this function exists; it must accept what that suggests.
    @test "FRA" in [m.name for m in naturalearth_regions("ISO_A3_EH")]
    @test_throws ErrorException naturalearth_regions("ISO_A3")

    # Passing a REGION where a level belongs is the likeliest slip, and saying only "no such level"
    # is true and useless when the answer is "that is a region, and here is the level it lives at".
    err = try
        naturalearth_regions("Sabrina Coast")
    catch e
        sprint(showerror, e)
    end
    @test occursin("is a region, not a level", err)
    @test occursin("Physical Coast", err)

    # A typo of a region name matches neither a level nor a region exactly, and is the case that was
    # most useless before: "Yemn" is a misspelt country, not a misspelt level.
    typo = try
        naturalearth_regions("Yemn")
    catch e
        sprint(showerror, e)
    end
    @test occursin("Yemen", typo)
    # Every error here points at the *function*, which displays as a readable table, rather than at
    # the underlying `const`.
    @test occursin("naturalearth_levels()", typo)
    @test !occursin("NATURALEARTH_LEVELS", typo)

    # ...and a query resembling nothing says so without inventing a suggestion.
    nothinglike = try
        naturalearth_regions("qqqq")
    catch e
        sprint(showerror, e)
    end
    @test !occursin("Did you mean", nothinglike)

    @test EcoSISTEM._editdistance("Yemn", "Yemen") == 1
    @test EcoSISTEM._editdistance("abc", "abc") == 0
    @test EcoSISTEM._editdistance("", "abc") == 3
end

@testset "a subject's extent is asked polymorphically, and always in degrees" begin
    @test _wgsextent(Extents.Extent(Y = (1.0, 2.0), X = (3.0, 4.0))).Y ==
          (1.0, 2.0)
    # A point is a zero-area extent, which is what makes it answerable at all.
    pt = _wgsextent(LatLong(55.9°, -3.2°))
    @test pt.Y[1] == pt.Y[2] && pt.X[1] == pt.X[2]
    # Something with no grid cannot be asked, and says so rather than giving `nothing` onwards.
    @test_throws ErrorException _wgsextent(42)
end

@testset "the three relations, and which accept a point" begin
    box = Extents.Extent(Y = (54.0°, 58.0°), X = (-6.0°, -2.0°))
    inside = Extents.Extent(Y = (55.0°, 56.0°), X = (-5.0°, -4.0°))
    apart = Extents.Extent(Y = (10.0°, 11.0°), X = (10.0°, 11.0°))

    @test Encloses(inside)(box)          # the box encloses the smaller extent
    @test !Encloses(box)(inside)
    @test Within(box)(inside)            # the smaller extent lies within the box
    @test Overlaps(inside)(box)
    @test !Overlaps(apart)(box)
    @test Encloses(inside) isa EcoSISTEM.AbstractSpatialRelation

    # `Encloses` delegates to `covers`, not `contains`: the two agree on every pair of boxes but
    # differ on a point, and a point is the one subject only `Encloses` takes. Using `contains`
    # would make "which regions enclose this coordinate?" silently answer nothing.
    edinburgh = LatLong(55.95°, -3.19°)
    @test Encloses(edinburgh)(box)
    @test !Extents.contains(box, _wgsextent(edinburgh))

    # ...and the other two refuse a point at construction, where it was written, rather than
    # returning an empty report that reads as "nothing found".
    @test_throws ErrorException Within(edinburgh)
    @test_throws ErrorException Overlaps(edinburgh)
end

@testset "overlap area is a solid angle, not a product of degree spans" begin
    equator = Extents.Extent(Y = (0.0°, 1.0°), X = (0.0°, 1.0°))
    arctic = Extents.Extent(Y = (70.0°, 71.0°), X = (0.0°, 1.0°))
    # The same angular box covers far less ground near the pole, so ordering by a flat product would
    # rank high latitudes above equatorial ones for the same shared ground.
    @test _overlaparea(equator, equator) > _overlaparea(arctic, arctic)
    @test _overlaparea(equator, equator) ≈ 12308km^2 rtol=0.01
    # Disjoint boxes share nothing.
    @test _overlaparea(equator,
                       Extents.Extent(Y = (10.0°, 11.0°), X = (10.0°, 11.0°))) ==
          0.0km^2
end

@testset "investigate_regions finds what encloses a point, smallest first" begin
    report = investigate_regions(LatLong(55.95°, -3.19°), kind = :political)
    @test report isa RegionReport
    @test !isempty(report)
    # A container, so Base does the disambiguating: `only`, `first` and indexing all work.
    @test first(report) isa RegionMatch
    @test report[1] === first(report)
    @test length(collect(report)) == length(report)

    names = [m.name for m in report]
    @test "Scotland" in names && "United Kingdom" in names
    # Smallest enclosing first, which is the whole point of the ordering.
    @test issorted([m.area for m in report])
    @test first(report).name == "Scotland"

    # Filters narrow at source, which is the first tool for an ambiguous answer.
    @test all(m -> m.level.kind === :political, report)
    @test all(m -> m.level.name == "ADMIN",
              investigate_regions(LatLong(55.95°, -3.19°), level = "ADMIN"))
    @test length(investigate_regions(LatLong(55.95°, -3.19°), limit = 2)) == 2

    # A match is a subject in its own right, so a row feeds back in: `Encloses` walks up the
    # hierarchy, `Within` walks down.
    scotland = first(report)
    up = investigate_regions(Encloses(scotland), kind = :political)
    @test "United Kingdom" in [m.name for m in up]
    # ...and a region never answers a question about itself, `covers` and `within` both being
    # reflexive. Excluded by identity, never by extent: distinct regions share a box legitimately.
    @test !any(m -> m.extent == scotland.extent, up)
end

@testset "a match converts to a spec; a report does not" begin
    report = investigate_regions(LatLong(55.95°, -3.19°), level = "ADMIN")

    # TWO countries' boxes enclose Edinburgh: the United Kingdom, and **Norway**, whose box runs
    # west to Jan Mayen and north to Svalbard. That is the box tier being honest rather than wrong,
    # and it is why `only` is the safe way to take a single answer - it refuses here.
    @test length(report) == 2
    @test [m.name for m in report] == ["United Kingdom", "Norway"]
    @test_throws ArgumentError only(report)

    # `first` is meaningful under `Encloses`, whose order is smallest-enclosing-first.
    spec = NaturalEarthSpec(first(report))
    @test spec isa NaturalEarthSpec
    @test spec.name == "United Kingdom" && spec.level == "ADMIN"
    # The conversion re-derives nothing, so the spec's shape agrees with the box displayed.
    @test spec.name == first(report).name

    # `only` does work where the query really is unique.
    @test NaturalEarthSpec(only(investigate_regions(LatLong(55.95°, -3.19°),
                                                    level = "SUBUNIT"))).name ==
          "Scotland"

    # A report may hold several regions, so it names no single spec - and `first` is meaningful for
    # only one of the three orderings. The error names the three ways to choose.
    err = try
        NaturalEarthSpec(report)
    catch e
        sprint(showerror, e)
    end
    @test occursin("only(report)", err) && occursin("first(report)", err)
end

@testset "a report says its answer is coarse, every time" begin
    report = investigate_regions(LatLong(55.95°, -3.19°), limit = 3)
    shown = sprint(show, MIME"text/plain"(), report)
    # Stated in the display rather than left to the docstring: a box can be far larger than the
    # ground it names, and a reader who does not know that will trust the list too far.
    @test occursin("bounding box", shown)
    # The remedy it names is the one for a coarse *query* - the exact tier - rather than building a
    # shape, which answers a different question.
    @test occursin("exact = true", shown)
    # The relation reads as a verb agreeing with its subject.
    @test occursin("enclose it", shown)
    @test occursin("lies within it",
                   sprint(show, MIME"text/plain"(),
                          investigate_regions(Within(Extents.Extent(Y = (54.6°,
                                                                         58.7°),
                                                                    X = (-6.3°,
                                                                         -1.7°))),
                                              limit = 1)))
end

end
