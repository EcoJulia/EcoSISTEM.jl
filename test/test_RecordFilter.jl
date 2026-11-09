# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/RecordFilter.jl`: `IsActiveSpecies` keeping the rows of a table that name an active species
# of a pool or a species list, in place or into a new table, matching names through a name authority (a
# made-up WFO release), and its report - each kind of name dropped, its level and its detail; and
# `InActiveCell` keeping the rows whose position lies in an active cell of a study area or a habitat,
# or within `snap` of one, on geographic and projected grids. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_RecordFilter.jl"])'

module TestRecordFilter

using CSV
using DataFrames
using DimensionalData: DimensionalData, DimArray, X, Y
using EcoSISTEM
using EcoSISTEM: AbstractRecordFilter, ProblemWarning, ProblemError
using EcoSISTEM.Units
using Rasters
using RasterDataSources
using Test
using Unitful
using Unitful.DefaultSymbols

include("wfofixtures.jl")
include("rasterfixtures.jl")
include("buildfixtures.jl")

# The species of `pool`, every one alike.
function _species(pool::SpeciesPool)
    return build_species(pool,
                         tolerance = ToleranceSpec(Temperature, μ = 298.0K,
                                                   σ = 2.0K),
                         demand = DemandSpec(SolarRadiation, 1.0e-4kJ / day),
                         seed = 1, verbosity = :silent)
end

# A 5 x 7 grid of 1-degree cells over 50-55 N and 0-7 E, with `inactive` cells off: across the middle
# row a degree of longitude is about 68 km, and a degree of latitude is 111 km everywhere.
function _geohabitat(inactive = [CartesianIndex(3, 4), CartesianIndex(3, 3),
                         CartesianIndex(2, 4)])
    temp = _testraster(WorldClim{BioClim}, fill(298.0K, 5, 7),
                       lat = (50.0:1.0:54.0) .* °, long = (0.0:1.0:6.0) .* °)
    mask = trues(5, 7)
    mask[inactive] .= false
    regime = _reg(temp, axis = Temperature)
    return GridHabitat(regime = regime, supply = SUP,
                       area = _area(regime = regime, within = mask))
end

# A pool of three species, the third not active, and a copy of the first.
function _pool()
    pool = SpeciesPool(["a", "b", "c"], verbosity = :silent)
    deactivate_species!(pool, "c", verbosity = :silent)
    return copy_species!(pool, "a", as = "a copy")
end

@testset "a filter keeps the rows naming an active species, and reports the rest" begin
    pool = _pool()
    records = DataFrame(species = ["a", "b", "c", "zz", missing, "a copy", "zz",
                            "a"],
                        n = 1:8)
    kept = @test_logs((:info,
                       r"^2 records of 1 name of no species in the pool were removed: `zz` \(2\)\.$"),
                      (:info,
                       r"^1 record of 1 species not active was removed: `c` \(1\)\. `activate_species!` activates them\.$"),
                      (:info, r"^1 record naming no species was removed\.$"),
                      filter!(IsActiveSpecies(pool), records))
    @test kept === records
    @test records.n == [1, 2, 6, 8]
    @test IsActiveSpecies(pool) isa AbstractRecordFilter
    # A species list filters as the pool it was built from does, naming the way to activate a species.
    species = _species(pool)
    records = DataFrame(species = ["a", "c", "a copy", "zz"])
    @test_logs((:info, r"of no species in the species list was removed"),
               (:info,
                r"Activate them in their pool with `activate_species!` and build the species again"),
               filter!(IsActiveSpecies(species), records))
    @test records.species == ["a", "a copy"]
end

@testset "filter returns a new table, of the kind it was given" begin
    pool = _pool()
    quiet = IsActiveSpecies(pool, verbosity = :silent)
    records = DataFrame(species = ["a", "zz"], n = [1, 2])
    new = @test_logs filter(quiet, records)
    @test new.n == [1] && records.n == [1, 2]
    @test filter(quiet, (species = ["zz", "b"], n = [1, 2])) ==
          (species = ["b"], n = [2])
    mktempdir() do dir
        path = joinpath(dir, "records.csv")
        CSV.write(path, records)
        file = CSV.File(path)
        filtered = filter(quiet, file)
        @test filtered isa DataFrame
        @test filtered.species == ["a"] && filtered.n == [1]
        @test_throws "cannot be changed in place" filter!(quiet, file)
    end
end

@testset "a report is raised at its level, in more detail at `:verbose`" begin
    pool = _pool()
    records = DataFrame(species = vcat(["x$i" for i in 1:12], "a"))
    # The first ten names at `:normal`, the heaviest first and then by name.
    @test_logs((:warn,
                r"^12 records of 12 names .*: `x1` \(1\), `x10` \(1\), .*`x7` \(1\) and 2 more\.$"),
               filter(IsActiveSpecies(pool, level = ProblemWarning()), records))
    @test_logs((:info, r"`x8` \(1\), `x9` \(1\)\.$"),
               filter(IsActiveSpecies(pool, verbosity = :verbose), records))
    heavy = DataFrame(species = ["y", "x", "y"])
    @test_logs((:info, r": `y` \(2\), `x` \(1\)\.$"),
               filter(IsActiveSpecies(pool), heavy))
    # An error stops the filter before the table changes.
    @test_throws "12 records of 12 names" filter!(IsActiveSpecies(pool,
                                                                  level = ProblemError()),
                                                  records)
    @test nrow(records) == 13
    # The column naming the species, and what is refused.
    @test filter(IsActiveSpecies(pool, column = :taxon, verbosity = :silent),
                 DataFrame(taxon = ["a", "zz"])).taxon == ["a"]
    @test_throws "names no species in a `species` column; it has `taxon`" filter(IsActiveSpecies(pool),
                                                                                 DataFrame(taxon = ["a"]))
    @test_throws "row 1 of a table names its species as `3`" filter(IsActiveSpecies(pool),
                                                                    DataFrame(species = [3]))
    @test_throws "`verbosity` must be" IsActiveSpecies(pool, verbosity = :loud)
end

@testset "names are found through the pool's name authority" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        pool = SpeciesPool(["Acacia alpha", "Acacia zeta"], authority = wfo,
                           verbosity = :silent)
        # A synonym reaches its accepted species, an identifier its species; a name the authority
        # recognises but the pool does not hold is reported apart from one it does not recognise.
        given = ["Mimosa alpha", "Acacia alpha", "wfo-0000000010",
            "Acacia theta", "Acacia theta", "Algae x"]
        for holder in (pool, _species(pool))
            records = DataFrame(species = given)
            @test_logs((:info,
                        r"^1 record of 1 name not recognised by the WFO Plant List 2026-06 was removed: `Algae x` \(1\)\.$"),
                       (:info,
                        r"^2 records of 1 name recognised by the WFO Plant List 2026-06 but of no species in the (pool|species list) were removed: `Acacia theta` \(2\)\.$"),
                       filter!(IsActiveSpecies(holder), records))
            @test records.species == given[1:3]
        end
        # Names all found directly - an own name, an identifier - never read the authority's names.
        empty!(EcoSISTEM._WFOLOADED)
        @test filter(IsActiveSpecies(pool),
                     DataFrame(species = ["Acacia alpha", "wfo-0000000010"])).species ==
              ["Acacia alpha", "wfo-0000000010"]
        @test isempty(EcoSISTEM._WFOLOADED)
        filter(IsActiveSpecies(pool), DataFrame(species = ["Mimosa alpha"]))
        @test !isempty(EcoSISTEM._WFOLOADED)
    end
end

@testset "a filter by position keeps the records in active cells, and reports the rest" begin
    habitat = _geohabitat()
    records = DataFrame(decimalLatitude = [52.5, 52.5, 49.0, missing, 50.5],
                        decimalLongitude = [3.5, 0.5, 3.5, 1.0, 3.5], n = 1:5)
    kept = @test_logs((:info, r"^1 record off the grid was removed\.$"),
                      (:info, r"^1 record in an inactive cell was removed\.$"),
                      (:info, r"^1 record with no position was removed\.$"),
                      filter(InActiveCell(habitat), records))
    @test kept.n == [2, 5] && nrow(records) == 5
    @test InActiveCell(habitat) isa AbstractRecordFilter
    # A habitat's cells as they stand when it is called, an area's as it was decided.
    habitat.active[1, 1] = false
    corner = DataFrame(decimalLatitude = [50.5], decimalLongitude = [0.5])
    @test nrow(filter(InActiveCell(habitat, verbosity = :silent), corner)) == 0
    @test nrow(filter(InActiveCell(habitat.area, verbosity = :silent), corner)) ==
          1
    @test nrow(filter(InActiveCell(habitat.area.report, verbosity = :silent),
                      corner)) == 1
    # Other columns, in place, and the level the removals are reported at.
    renamed = DataFrame(lat = [52.5, 50.5], long = [3.5, 3.5])
    @test_logs((:warn, r"^1 record in an inactive cell was removed\.$"),
               filter!(InActiveCell(habitat, latitude = :lat, longitude = :long,
                                    level = ProblemWarning()), renamed))
    @test renamed.lat == [50.5]
end

@testset "`snap` keeps a record within a ground distance of an active cell" begin
    habitat = _geohabitat()
    records = DataFrame(decimalLatitude = [52.5], decimalLongitude = [3.5])
    # Its cell's western, northern and southern neighbours are inactive or 111 km away; its eastern
    # neighbour's centre is 68 km away.
    @test_logs((:info,
                r"^1 record in an inactive cell, with no active cell within 60 km, was removed\.$"),
               filter(InActiveCell(habitat, snap = 60km), records))
    @test_logs((:info,
                r"^1 record in an inactive cell within 70 km of an active cell was kept\.$"),
               filter(InActiveCell(habitat, snap = 70km), records))
    # The nearest active cell, by the ground distance between centres, ties to the lower `(y, x)`.
    yx = EcoSISTEM._gridyx(habitat)
    cell = CartesianIndex(3, 4)
    @test EcoSISTEM._celldistance(yx, cell, CartesianIndex(3, 5)) ≈ 67.77km rtol = 1e-3
    @test EcoSISTEM._celldistance(yx, cell, CartesianIndex(4, 4)) ≈ 111.32km
    nearest(h, snap) = EcoSISTEM._nearestactive(yx, parent(h.active), cell,
                                                snap)
    @test nearest(habitat, 70km) == CartesianIndex(3, 5)
    @test isnothing(nearest(habitat, 60km)) && isnothing(nearest(habitat, 0km))
    @test nearest(_geohabitat([cell]), 70km) == CartesianIndex(3, 3)
    @test nearest(_geohabitat([
                                  cell,
                                  CartesianIndex(3, 3),
                                  CartesianIndex(3, 5)
                              ]),
                  120km) == CartesianIndex(2, 4)
    # A projected grid's distances are straight, in its own units.
    temp = _bngraster(WorldClim{BioClim}, fill(298.0K, 9, 9))
    bng = _env(_reg(temp, axis = Temperature), SUP)
    byx = EcoSISTEM._gridyx(bng)
    @test EcoSISTEM._celldistance(byx, CartesianIndex(4, 4),
                                  CartesianIndex(5, 5)) ≈ hypot(2.5, 2.5) * km
end

@testset "a filter by position refuses what it cannot place" begin
    habitat = _geohabitat()
    @test_throws "gives no position in a `decimalLatitude` column; it has `lat`" filter(InActiveCell(habitat),
                                                                                        DataFrame(lat = [1.0],
                                                                                                  decimalLongitude = [1.0]))
    @test_throws "row 1 of a table gives its latitude as `\"north\"`" filter(InActiveCell(habitat),
                                                                             DataFrame(decimalLatitude = ["north"],
                                                                                       decimalLongitude = [1.0]))
    @test_throws "a finite, non-negative length; got -1 km" InActiveCell(habitat,
                                                                         snap = -1km)
    @test_throws TypeError InActiveCell(habitat, snap = 3)
    # A grid with no real-world position cannot place a latitude and a longitude.
    synthetic = StudyArea(extent = (7.0km, 5.0km), cellsize = 1.0km,
                          verbosity = :silent)
    @test_throws "the place given is in °" filter(InActiveCell(synthetic),
                                                  DataFrame(decimalLatitude = [1.0],
                                                            decimalLongitude = [1.0]))
end

@testset "an ecosystem filters as its species list and its habitat do" begin
    # A projected grid, since an ecosystem cannot be simulated on a geographic one, and two places in
    # different cells of it.
    habitat = _env(_reg(_bngraster(WorldClim{BioClim}, fill(298.0K, 9, 9)),
                        axis = Temperature), SUP)
    lattice = [(lat, long) for lat in 55.6:0.02:56.0 for long in -4.6:0.02:-3.8]
    cells = EcoSISTEM.getcellat(habitat,
                                [LatLong(lat * °, long * °)
                                 for (lat, long) in lattice])
    here = findfirst(!isnothing, cells)
    there = findfirst(c -> !isnothing(c) && c != cells[here], cells)
    (lat1, long1), (lat2, long2) = lattice[here], lattice[there]
    species = _species(_pool())
    eco = build_ecosystem(NoAbundance(), species, habitat, verbosity = :silent)
    records = DataFrame(species = ["a", "c", "a copy", "zz", "b"],
                        decimalLatitude = [lat1, lat1, lat2, lat2, 0.0],
                        decimalLongitude = [long1, long1, long2, long2, 0.0])
    bylist = filter(IsActiveSpecies(species, verbosity = :silent), records)
    @test filter(IsActiveSpecies(eco, verbosity = :silent), records) == bylist
    @test bylist.species == ["a", "a copy", "b"]
    # Its habitat's cells as they stand: one deactivated is filtered out, as is a place off the grid.
    eco.habitat.active[cells[there]] = false
    bycells = filter(InActiveCell(eco.habitat, verbosity = :silent), records)
    @test filter(InActiveCell(eco, verbosity = :silent), records) == bycells
    @test bycells.species == ["a", "c"]
end

@testset "a filter shows as it is written" begin
    pool = _pool()
    shown = sprint(show, pool)
    @test sprint(show, IsActiveSpecies(pool)) == "IsActiveSpecies($shown)"
    @test sprint(show,
                 IsActiveSpecies(pool, column = :taxon,
                                 level = ProblemWarning(),
                                 verbosity = :verbose)) ==
          "IsActiveSpecies($shown, column = :taxon, level = ProblemWarning(), verbosity = :verbose)"
    plain = sprint(show, MIME"text/plain"(), IsActiveSpecies(pool))
    @test startswith(plain, "IsActiveSpecies\n") &&
          occursin("level     ProblemNotice()", plain)
    habitat = _geohabitat()
    @test sprint(show, InActiveCell(habitat)) ==
          "InActiveCell($(sprint(show, habitat)))"
    @test sprint(show,
                 InActiveCell(habitat.area, latitude = :lat, snap = 10km,
                              level = ProblemWarning())) ==
          "InActiveCell($(sprint(show, habitat.area)), latitude = :lat, snap = 10 km, level = ProblemWarning())"
    plain = sprint(show, MIME"text/plain"(), InActiveCell(habitat, snap = 5km))
    @test startswith(plain, "InActiveCell\n") &&
          occursin("snap      5 km", plain) &&
          occursin("cells     GridHabitat", plain)
end

end
