# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/Placement.jl` and the placing in `build_ecosystem` and `repopulate!`: the default placement
# and the one it equals, what each shipped placement places, a table's rows, the rows it drops and the
# names it finds through a name authority, a function as a placement, and what the ecosystem keeps of
# how it was placed. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_Placement.jl"])'

module TestPlacement

using ArchGDAL
using CSV
using DataFrames
using Dates: Date
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using Distributions
using EcoSISTEM
using EcoSISTEM: getcellat, materialise
using EcoSISTEM.Units
using Extents: Extent
using Rasters
using RasterDataSources
using Test
using Unitful
using Unitful.DefaultSymbols

include("rasterfixtures.jl")
include("buildfixtures.jl")
include("wfofixtures.jl")

# A 7 x 5 grid of 1 km cells with its first cell inactive; the supply rises across it, so supply
# weights differ between cells.
const MASK = let mask = trues(7, 5)
    mask[1, 1] = false
    mask
end

const SOLAR = GradientSpec(1.0e5kJ / (km^2 * day), 9.0e5kJ / (km^2 * day),
                           axis = SolarRadiation)

function _habitat(; regime = UniformSpec(298.0K, axis = Temperature),
                  supply = SOLAR)
    area = StudyArea(extent = (7.0km, 5.0km), cellsize = 1.0km, within = MASK,
                     verbosity = :silent)
    return build_habitat(regime = regime, supply = supply, area = area,
                         verbosity = :silent)
end

# The species of `pool`, every one alike.
function _species(pool::SpeciesPool;
                  demand = DemandSpec(SolarRadiation, 1.0e-4kJ / day))
    return build_species(pool,
                         tolerance = ToleranceSpec(Temperature, μ = 298.0K,
                                                   σ = 2.0K),
                         demand = demand, seed = 1, verbosity = :silent)
end

# `n` species from a pool, every one alike.
function _species(n = 4; demand = DemandSpec(SolarRadiation, 1.0e-4kJ / day))
    return _species(SpeciesPool(n, verbosity = :silent), demand = demand)
end

# An ecosystem built quietly by a placement.
function _placed(p, species, habitat; seed = 1)
    return build_ecosystem(p, species, habitat,
                           seed = seed,
                           verbosity = :silent)
end

# Each species' individuals, in species order.
_perspecies(eco) = vec(sum(eco.abundances.matrix, dims = 2))

# The storage column of the grid's inactive cell.
const INACTIVE = LinearIndices(size(MASK))[1, 1]

# A 9 x 9 grid of 2.5 km cells over central Scotland in British National Grid, its centre cell
# inactive, where occurrence records place, since an ecosystem refuses a geographic grid.
function _bnghabitat()
    mask = trues(9, 9)
    mask[5, 5] = false
    return _env(_reg(_bngraster(WorldClim{BioClim}, fill(298.0K, 9, 9)),
                     axis = Temperature), SUP, within = mask)
end

# A latitude and longitude in each cell of `habitat`'s grid, from a lattice over central Scotland.
function _positions(habitat)
    lattice = [(lat, long) for lat in 55.6:0.005:56.0
               for long in -4.7:0.005:-3.8]
    cells = getcellat(habitat,
                      [LatLong(lat * 1.0°, long * 1.0°)
                       for (lat, long) in lattice])
    found = Dict{CartesianIndex{2}, Tuple{Float64, Float64}}()
    for (place, cell) in zip(lattice, cells)
        isnothing(cell) || get!(found, cell, place)
    end
    return found
end

# Occurrence records of the four species of `_species()`, the fourth to be made inactive: species 1
# three times in cell `A` and once in `B`, species 2 twice in `B`, species 3 once in the inactive cell
# `C`; then records of the inactive species, of a name the list does not hold, off the grid and with
# no position; `n` numbers them.
const A, B, C = CartesianIndex(2, 3), CartesianIndex(6, 7), CartesianIndex(5,
                                                                           5)

function _records(habitat)
    at = _positions(habitat)
    where = [A, A, A, B, B, B, C, A, A]
    return DataFrame(species = ["1", "1", "1", "1", "2", "2", "3", "4", "zz",
                         "1",
                         "1"],
                     decimalLatitude = vcat([at[c][1] for c in where], 10.0,
                                            missing),
                     decimalLongitude = vcat([at[c][2] for c in where], 0.0,
                                             1.0),
                     n = 1:11)
end

@testset "the default placement places 1000 of each active species by supply, and says so" begin
    for (habitat, species) in ((_habitat(), _species()),
                               (_habitat(supply = (SOLAR,
                                                   GradientSpec(1.0mm / day,
                                                                5.0mm / day,
                                                                axis = Precipitation))),
                                _species(demand = (DemandSpec(SolarRadiation,
                                                              1.0e-4kJ / day),
                                                   DemandSpec(Precipitation,
                                                              1.0e-4Unitful.L /
                                                              day)))))
        placed = @test_logs (:info,
                             r"`placement` defaulted to `DefaultEcosystem\(\)`, 1000 individuals per active species") build_ecosystem(species,
                                                                                                                                      habitat,
                                                                                                                                      seed = 7)
        @test placed.placement === DefaultEcosystem()
        @test sum(placed.abundances.matrix) == 4000
        # The community total `BySupply` is given, draw for draw.
        bysupply = _placed(BySupply(total = 4000), species, habitat, seed = 7)
        @test placed.abundances.matrix == bysupply.abundances.matrix
        @test placed.rngs == bysupply.rngs
    end
    # A species that is not active is placed none, and the rest 1000 each.
    species = _species()
    species.traits.active[2] = false
    eco = _placed(DefaultEcosystem(), species, _habitat())
    @test _perspecies(eco)[2] == 0 && sum(eco.abundances.matrix) == 3000
end

@testset "each placement places what it says" begin
    habitat = _habitat()
    species = _species()
    active = count(MASK)
    # A community total, or each species' own; none in the inactive cell.
    total = _placed(BySupply(total = 500), species, habitat)
    @test sum(total.abundances.matrix) == 500
    @test all(iszero, total.abundances.matrix[:, INACTIVE])
    @test _perspecies(_placed(BySupply(each = 30), species, habitat)) ==
          fill(30, 4)
    # A density is of the cells it may place in, rounded half to even.
    @test sum(_placed(BySupply(total = 2.0 / km^2), species, habitat).abundances.matrix) ==
          2 * active
    @test _perspecies(_placed(BySupply(each = 0.5 / km^2), species, habitat)) ==
          fill(round(Int, 0.5 * active), 4)
    # Each species' own number from its traits, by a function or a draw from its own stream.
    @test _perspecies(_placed(BySupply(each = row -> 10 * row.key), species,
                              habitat)) == [10, 20, 30, 40]
    drawn = BySupply(each = DrawFrom(Poisson(50)))
    @test _placed(drawn, species, habitat).abundances.matrix ==
          _placed(drawn, species, habitat).abundances.matrix
    @test _placed(drawn, species, habitat).abundances.matrix !=
          _placed(drawn, species, habitat, seed = 2).abundances.matrix
    # The supply weights the cells: the richest half holds most.
    weighted = _placed(BySupply(total = 20_000), species, habitat)
    bycell = reshape(vec(sum(weighted.abundances.matrix, dims = 1)), size(MASK))
    @test sum(bycell[:, 4:5]) > sum(bycell[:, 1:2]) ||
          sum(bycell[5:7, :]) > sum(bycell[1:3, :])
    # Given no numbers, 1000 per active species, said.
    @test_logs (:info,
                r"places 1000 individuals per active species, 4000 in all") build_ecosystem(BySupply(),
                                                                                            species,
                                                                                            habitat,
                                                                                            verbosity = :normal)
    @test sum(_placed(BySuitability(), species, habitat).abundances.matrix) ==
          4000
    @test sum(_placed(NoAbundance(), species, habitat).abundances.matrix) == 0
    # Placements add, and weights multiply: a layer weight of zero on half the grid keeps everything
    # off that half.
    added = _placed(BySupply(total = 300) + BySuitability(each = 5), species,
                    habitat)
    @test sum(added.abundances.matrix) == 320
    halved = materialise(UniformSpec(1.0, axis = EcoSISTEM.NicheAxis),
                         habitat.area)
    halved.matrix[:, 1:2] .= 0.0
    onhalf = _placed(Weighted(BySupply(), ByLayer(halved), total = 1_000),
                     species, habitat)
    @test sum(onhalf.abundances.matrix) == 1_000
    @test all(iszero,
              reshape(vec(sum(onhalf.abundances.matrix, dims = 1)),
                      size(MASK))[:, 1:2])
end

@testset "suitability places each species where it suits" begin
    habitat = _habitat(regime = GradientSpec(290.0K, 306.0K,
                                             axis = Temperature))
    pool = SpeciesPool(DataFrame(species = ["1", "2"],
                                 optimum = [291.0, 305.0]),
                       verbosity = :silent)
    species = build_species(pool,
                            tolerance = ToleranceSpec(Temperature,
                                                      μ = TraitColumn(:optimum,
                                                                      unit = K),
                                                      σ = 1.0K),
                            demand = DemandSpec(SolarRadiation, 1.0e-4kJ / day),
                            seed = 1, verbosity = :silent)
    eco = _placed(BySuitability(each = 2_000), species, habitat)
    grids = [reshape(eco.abundances.matrix[sp, :], size(MASK)) for sp in 1:2]
    # The two species' warm and cool ends, whichever way the gradient runs.
    @test grids[1] != grids[2]
    rows = [sum(g, dims = 2) for g in grids]
    cols = [sum(g, dims = 1) for g in grids]
    @test argmax(vec(rows[1])) != argmax(vec(rows[2])) ||
          argmax(vec(cols[1])) != argmax(vec(cols[2]))
    # By its axis, or refused naming the members there are.
    @test sum(_placed(BySuitability(Temperature, each = 10), species, habitat).abundances.matrix) ==
          20
    @test_throws r"no member Precipitation; it has `Temperature` \(Temperature\)" _placed(BySuitability(Precipitation),
                                                                                          species,
                                                                                          habitat)
    @test_throws r"no member :summer" _placed(BySupply(:summer), species,
                                              habitat)
end

@testset "what a placement cannot place is said, or stops it" begin
    habitat = _habitat()
    species = _species()
    @test_throws "not both" BySupply(total = 10, each = 1)
    @test_throws "whole number" BySupply(total = 2.5)
    @test_throws "density per area" BySupply(each = 3.0km)
    @test_throws "give its `total` or `each` to the `Weighted`" Weighted(BySupply(total = 1))
    @test_throws "multiplies weights" Weighted(NoAbundance())
    @test_throws "a placement says how many" build_ecosystem(5, species,
                                                             habitat)
    # A weight that is zero everywhere places nobody, and says so.
    nowhere = ByLayer(UniformSpec(0.0, axis = EcoSISTEM.NicheAxis), total = 100)
    eco = @test_logs (:warn,
                      r"100 individuals were not placed, their species having no weight") build_ecosystem(nowhere,
                                                                                                          species,
                                                                                                          habitat,
                                                                                                          verbosity = :silent)
    @test sum(eco.abundances.matrix) == 0
    @test_throws "finite and non-negative" _placed(ByLayer(UniformSpec(-1.0,
                                                                       axis = EcoSISTEM.NicheAxis)),
                                                   species, habitat)
end

@testset "a table places its rows, and drops what it cannot place" begin
    habitat = _habitat()
    species = _species()
    table = (species = ["1", "1", "2", "zz", "3"],
             cell = [CartesianIndex(2, 3), CartesianIndex(2, 3),
                 CartesianIndex(7, 5), CartesianIndex(3, 3),
                 CartesianIndex(1, 1)],
             count = [4, 6, 7, 9, 11])
    eco = @test_logs (:warn,
                      r"9 individuals of 1 name of no species in the species list were not placed: `zz` \(9\)\.") (:warn,
                                                                                                                   r"11 individuals in 1 row naming a cell the placement may not place in") build_ecosystem(FromTable(table),
                                                                                                                                                                                                            species,
                                                                                                                                                                                                            habitat,
                                                                                                                                                                                                            verbosity = :silent)
    linear = LinearIndices(size(MASK))
    @test eco.abundances.matrix[1, linear[2, 3]] == 10
    @test eco.abundances.matrix[2, linear[7, 5]] == 7
    @test sum(eco.abundances.matrix) == 17
    # A place given as a location is the cell containing it; one off the grid is dropped.
    located = DataFrame(species = ["2", "4"],
                        location = [EcoSISTEM.SpatialLocation(3.5km, 1.5km),
                            EcoSISTEM.SpatialLocation(30.0km, 1.5km)])
    eco = @test_logs (:warn,
                      r"1 individual in 1 row naming a place off the grid") build_ecosystem(FromTable(located),
                                                                                            species,
                                                                                            habitat,
                                                                                            verbosity = :silent)
    @test eco.abundances.matrix[2,
                                linear[getcellat(habitat,
                                                 EcoSISTEM.SpatialLocation(3.5km,
                                                                           1.5km))]] ==
          1
    @test sum(eco.abundances.matrix) == 1
    # A species not active is dropped, naming the way to activate it.
    species.traits.active[3] = false
    @test_logs (:warn,
                r"5 individuals of 1 species not active were not placed: `3` \(5\)\. Activate them in their pool with `activate_species!`") build_ecosystem(FromTable((species = ["3"],
                                                                                                                                                                       cell = [CartesianIndex(4,
                                                                                                                                                                                              4)],
                                                                                                                                                                       count = [5])),
                                                                                                                                                            species,
                                                                                                                                                            habitat,
                                                                                                                                                            verbosity = :silent)
    species.traits.active[3] = true
    # A species may be named by its name authority's identifier.
    species.traits.id[4] = "id-4"
    eco = _placed(FromTable((species = ["id-4"], cell = [CartesianIndex(4, 4)])),
                  species, habitat)
    @test eco.abundances.matrix[4, linear[4, 4]] == 1
    # What stops it: a count that is not whole, a row after the start, a date with no epoch, no
    # species column, no place.
    onecell = [CartesianIndex(4, 4)]
    @test_throws "a whole, non-negative number" _placed(FromTable((species = ["1"],
                                                                   cell = onecell,
                                                                   count = [-1])),
                                                        species, habitat)
    @test_throws "places the starting state" _placed(FromTable((species = ["1"],
                                                                cell = onecell,
                                                                time = [1.0year])),
                                                     species, habitat)
    @test _placed(FromTable((species = ["1"], cell = onecell,
                             time = [-1.0year])), species, habitat).abundances.matrix[1,
                                                                                      linear[4,
                                                                                             4]] ==
          1
    @test_throws "no epoch" _placed(FromTable((species = ["1"], cell = onecell,
                                               time = [Date(2020)])),
                                    species, habitat)
    @test_throws "in a `species` column" _placed(FromTable((name = ["1"],
                                                            cell = onecell)),
                                                 species, habitat)
    @test_throws "by `latitude` and `longitude`" _placed(FromTable((species = ["1"],)),
                                                         species, habitat)
    @test_throws "the grid's own key" _placed(FromTable((species = ["1"],
                                                         cell = [7])),
                                              species, habitat)
    @test_throws "no file there" FromTable(joinpath(tempdir(),
                                                    "no-such-table.csv"))
end

@testset "a table's species are found through the pool's name authority" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        species = _species(SpeciesPool(["Acacia alpha", "Acacia zeta"],
                                       authority = wfo, verbosity = :silent))
        # A synonym places as its accepted species; a name the authority recognises but the list
        # does not hold is reported apart from one it does not recognise.
        table = (species = ["Mimosa alpha", "Acacia theta", "Algae x"],
                 cell = fill(CartesianIndex(4, 4), 3), count = [3, 4, 5])
        eco = @test_logs((:warn,
                          r"^5 individuals of 1 name not recognised by the WFO Plant List 2026-06 were not placed: `Algae x` \(5\)\.$"),
                         (:warn,
                          r"^4 individuals of 1 name recognised by the WFO Plant List 2026-06 but of no species in the species list were not placed: `Acacia theta` \(4\)\.$"),
                         build_ecosystem(FromTable(table), species, _habitat(),
                                         verbosity = :silent))
        @test _perspecies(eco) == [3, 0]
    end
end

@testset "a file's latitudes and longitudes place in a projected grid, and its record is kept" begin
    temp = _bngraster(WorldClim{BioClim}, fill(298.0K, 7, 9),
                      north = (640000.0:2500.0:655000.0) .* m)
    habitat = _env(_reg(temp, axis = Temperature),
                   UniformSpec(1.0e5kJ / (km^2 * day), axis = SolarRadiation))
    species = _species(2)
    # Places on a lattice over the grid's part of Scotland, the ones in it found as the grid finds
    # them.
    lattice = [LatLong(lat * 1.0°, long * 1.0°)
               for lat in 55.6:0.02:56.0 for long in -4.6:0.02:-3.8]
    cells = getcellat(habitat, lattice)
    inside = findall(!isnothing, cells)
    @test length(inside) > 10
    rows = DataFrame(species = ["1", "2"],
                     latitude = [ustrip(lattice[i].y) for i in inside[1:2]],
                     longitude = [ustrip(lattice[i].x) for i in inside[1:2]],
                     count = [3, 5])
    mktempdir() do dir
        path = joinpath(dir, "founders.csv")
        CSV.write(path, rows)
        eco = _placed(FromTable(path), species, habitat)
        bykey = _placed(FromTable((species = ["1", "2"],
                                   cell = cells[inside[1:2]],
                                   count = [3, 5])), species, habitat)
        @test eco.abundances.matrix == bykey.abundances.matrix
        @test sum(eco.abundances.matrix) == 8
        @test any(r -> r.role == :abundance && r.path == "founders.csv",
                  eco.inputs)
        @test provenance(eco).context.run.placement ==
              "FromTable($(repr(path)))"
        # A table that can be read only once places once.
        once = FromTable(CSV.Rows(path,
                                  types = [String, Float64, Float64, Int]))
        eco = _placed(once, species, habitat)
        @test sum(eco.abundances.matrix) == 8
        @test_throws "can be read only once" repopulate!(eco,
                                                         verbosity = :silent)
    end
end

@testset "occurrence records place an individual each, by their names and places" begin
    habitat = _bnghabitat()
    species = _species()
    species.traits.active[4] = false
    records = _records(habitat)
    linear = LinearIndices(EcoSISTEM.getgridshape(habitat))
    eco = @test_logs((:warn,
                      r"^1 record of 1 name of no species in the species list was not placed: `zz` \(1\)\.$"),
                     (:warn,
                      r"^1 record of 1 species not active was not placed: `4` \(1\)\..*activate_species!"),
                     (:warn, r"^1 record off the grid was not placed\.$"),
                     (:warn,
                      r"^1 record in an inactive cell was not placed\.$"),
                     (:warn, r"^1 record with no position was not placed\.$"),
                     _placed(FromOccurrences(records), species, habitat))
    @test eco.abundances.matrix[1, linear[A]] == 3
    @test eco.abundances.matrix[1, linear[B]] == 1
    @test eco.abundances.matrix[2, linear[B]] == 2
    @test sum(eco.abundances.matrix) == 6
    # Filtered first, only what is left is reported; `snap` moves a record into an active cell, and
    # a placement at a lower level says so as it does.
    kept = filter(IsActiveSpecies(species, verbosity = :silent), records)
    kept = filter(InActiveCell(habitat, snap = 3km, verbosity = :silent), kept)
    snapped = @test_logs((:info,
                          r"^1 record in an inactive cell within 3 km of an active cell was moved into it\.$"),
                         build_ecosystem(FromOccurrences(kept, snap = 3km),
                                         species,
                                         habitat, seed = 1))
    @test sum(snapped.abundances.matrix[3, :]) == 1
    @test snapped.abundances.matrix[3, linear[C]] == 0
    @test_logs((:info, r"^1 record off the grid was not placed\.$"),
               build_ecosystem(FromOccurrences(records[10:10, :],
                                               level = EcoSISTEM.ProblemNotice()),
                               species, habitat, seed = 1))
    # From a CSV file, its record kept; and placed again by `repopulate!`.
    mktempdir() do dir
        path = joinpath(dir, "occurrences.csv")
        CSV.write(path, records[1:6, :])
        eco = _placed(FromOccurrences(path), species, habitat)
        @test sum(eco.abundances.matrix) == 6
        @test any(r -> r.role == :abundance && r.path == "occurrences.csv",
                  eco.inputs)
        repopulate!(eco, verbosity = :silent)
        @test sum(eco.abundances.matrix) == 6
    end
end

@testset "occurrence records give `per_record` each, or share a `total` in proportion" begin
    habitat = _bnghabitat()
    species = _species()
    records = _records(habitat)[1:6, :]
    linear = LinearIndices(EcoSISTEM.getgridshape(habitat))
    # A number for each record: three of species 1's in `A`, one in `B`, two of species 2's in `B`.
    tripled = _placed(FromOccurrences(records, per_record = 3), species,
                      habitat)
    @test tripled.abundances.matrix[1, linear[A]] == 9
    @test tripled.abundances.matrix[1, linear[B]] == 3
    @test tripled.abundances.matrix[2, linear[B]] == 6
    @test _placed(FromOccurrences(records, per_record = 1), species,
                  habitat).abundances.matrix ==
          _placed(FromOccurrences(records), species, habitat).abundances.matrix
    # The community's number over every record, so species 1, with twice species 2's records, gets
    # about twice as many, and each cell about its share of the records - not equally over the species,
    # nor over the cells they occupy.
    big = _placed(FromOccurrences(records, total = 60_000), species, habitat)
    @test sum(big.abundances.matrix) == 60_000
    @test _perspecies(big)[1] ≈ 40_000 rtol = 0.02
    @test big.abundances.matrix[1, linear[A]] ≈ 30_000 rtol = 0.02
    @test big.abundances.matrix[1, linear[B]] ≈ 10_000 rtol = 0.05
    @test big.abundances.matrix[2, linear[B]] ≈ 20_000 rtol = 0.03
    @test big.abundances.matrix ==
          _placed(FromOccurrences(records, total = 60_000), species,
                  habitat).abundances.matrix
    @test big.abundances.matrix !=
          _placed(FromOccurrences(records, total = 60_000), species, habitat,
                  seed = 2).abundances.matrix
    # What is refused: a density, which records have no area for; both numbers; a number that is not
    # whole; a number with a rule, which makes its own.
    @test_throws "no area for a density to multiply" FromOccurrences(records,
                                                                     total = 2.0 /
                                                                             km^2)
    @test_throws "a `total` or a number `per_record`, not both" FromOccurrences(records,
                                                                                total = 10,
                                                                                per_record = 2)
    @test_throws "`per_record` is a whole number" FromOccurrences(records,
                                                                  per_record = 1.5)
    @test_throws "not both" FromOccurrences(records, per_record = 2
                                            ) do rows, _, _
        return rows
    end
    @test_throws "`FromOccurrences` was given the path" FromOccurrences("nowhere.csv")
end

@testset "a rule places what it makes of the occurrence records" begin
    habitat = _bnghabitat()
    species = _species()
    records = _records(habitat)[1:6, :]
    linear = LinearIndices(EcoSISTEM.getgridshape(habitat))
    seen = Ref{Any}()
    byn = FromOccurrences(records) do rows, habitat, context
        seen[] = rows
        return DataFrames.combine(DataFrames.groupby(rows, [:species, :cell]),
                                  :n => sum => :count)
    end
    eco = _placed(byn, species, habitat)
    @test eco.abundances.matrix[1, linear[A]] == 1 + 2 + 3
    @test eco.abundances.matrix[1, linear[B]] == 4
    @test eco.abundances.matrix[2, linear[B]] == 5 + 6
    @test Set(propertynames(seen[])) ==
          Set([:species, :decimalLatitude, :decimalLongitude, :n, :cell])
    @test seen[].cell[1] == A
    @test_throws "returned a value of type Int64" _placed(FromOccurrences((_, _,
                                                                           _) -> 3,
                                                                          records),
                                                          species, habitat)
end

@testset "a weight reads a layer at each record, for each record or shared in its cell" begin
    habitat = _bnghabitat()
    species = _species()
    records = _records(habitat)[1:6, :]
    linear = LinearIndices(EcoSISTEM.getgridshape(habitat))
    # An amount on the habitat's own grid, 30 in `A` and 60 in `B`, where species 1 has three
    # records in `A` and one in `B`, and species 2 two in `B`.
    function amounts()
        layer = materialise(UniformSpec(1.0km^2, axis = SurfaceArea),
                            habitat.area)
        layer.matrix[A] = 30.0km^2
        layer.matrix[B] = 60.0km^2
        return layer
    end
    amount = amounts()
    placed(weight; kw...) = _placed(FromOccurrences(records, weight = weight,
                                                    scale = 1.0 / km^2; kw...),
                                    species, habitat)
    # Every record its cell's value; or each cell's value once, shared among its records.
    perrecord = placed(WeightPerRecord(amount)).abundances.matrix
    @test perrecord[1, linear[A]] == 90
    @test perrecord[1, linear[B]] == 60
    @test perrecord[2, linear[B]] == 120
    percell = placed(WeightPerCell(amount)).abundances.matrix
    @test percell[1, linear[A]] == 30
    @test percell[1, linear[B]] == 20
    @test percell[2, linear[B]] == 40
    @test sum(percell) == 90
    # A record of the data shares its cell whether or not it is placed: a fourth in `A`, of a name
    # the list does not hold, leaves species 1 three quarters of 30, rounded to even.
    unknown = vcat(records, _records(habitat)[9:9, :])
    shared = @test_logs((:warn, r"`zz`"),
                        _placed(FromOccurrences(unknown,
                                                weight = WeightPerCell(amount),
                                                scale = 1.0 / km^2), species,
                                habitat)).abundances.matrix
    @test shared[1, linear[A]] == 22
    @test shared[1, linear[B]] == 20

    # Shared among the records in each cell of the layer's own grid, not the habitat's: two places
    # in `A` that lie in different cells of a grid twice as fine.
    finespec = _reg(_bngraster(WorldClim{BioClim}, fill(1.0, 18, 18),
                               east = (245000.0:1250.0:266250.0) .* m,
                               north = (640000.0:1250.0:661250.0) .* m))
    fine = materialise(finespec, _area(regime = finespec))
    lattice = [LatLong(lat * 1.0°, long * 1.0°)
               for lat in 55.6:0.002:56.0 for long in -4.7:0.002:-3.8]
    coarsecells = getcellat(habitat, lattice)
    finecells = getcellat(fine, lattice)
    ina = [i
           for i in eachindex(lattice)
           if coarsecells[i] == A && !isnothing(finecells[i])]
    here = first(ina)
    there = first(i for i in ina if finecells[i] != finecells[here])
    split = DataFrame(species = ["1", "1", "1", "2"],
                      decimalLatitude = [ustrip(lattice[i].y)
                                         for i in (here, here, here, there)],
                      decimalLongitude = [ustrip(lattice[i].x)
                                          for i in (here, here, here, there)])
    onfine = _placed(FromOccurrences(split, weight = WeightPerCell(fine),
                                     scale = 12), species,
                     habitat).abundances.matrix
    @test onfine[1, linear[A]] == 12
    @test onfine[2, linear[A]] == 12
    coarse = materialise(UniformSpec(1.0, axis = EcoSISTEM.NicheAxis),
                         habitat.area)
    oncoarse = _placed(FromOccurrences(split, weight = WeightPerCell(coarse),
                                       scale = 12), species,
                       habitat).abundances.matrix
    @test oncoarse[1, linear[A]] == 9
    @test oncoarse[2, linear[A]] == 3

    # Without a scale the weights are shares, which `ToCapacity` fills each cell by and nothing else
    # takes.
    filled = _placed(ToCapacity(FromOccurrences(records,
                                                weight = WeightPerCell(amount))),
                     species, habitat).abundances.matrix
    @test filled[2, linear[B]] / filled[1, linear[B]] ≈ 2 rtol = 0.01
    @test filled[2, linear[A]] == 0
    @test_throws "only `ToCapacity` takes" _placed(FromOccurrences(records,
                                                                   weight = WeightPerCell(amount)),
                                                   species, habitat)
    # A rule sees each record's weight, to make any other number from it.
    seen = Ref{Any}()
    ruled = FromOccurrences(records, weight = WeightPerRecord(amount)
                            ) do rows, _, _
        seen[] = rows.weight
        return (species = rows.species, cell = rows.cell,
                count = round.(Int, max.(ustrip.(km^2, rows.weight) ./ 30, 1)))
    end
    @test sum(_placed(ruled, species, habitat).abundances.matrix) == 9
    @test seen[] == [30.0, 30.0, 30.0, 60.0, 60.0, 60.0] .* km^2
    @test occursin("weight = WeightPerCell(",
                   sprint(show,
                          FromOccurrences(records,
                                          weight = WeightPerCell(amount),
                                          scale = 1.0 / km^2)))

    # Where the layer has no value a record takes the weight's `novalue`, and is refused without
    # one; a negative value is refused.
    gap = amounts()
    gap.matrix[B] = NaN * km^2
    @test_throws "has no value" placed(WeightPerCell(gap))
    nob = placed(WeightPerCell(gap, novalue = 0.0km^2)).abundances.matrix
    @test nob[1, linear[A]] == 30
    @test iszero(nob[:, linear[B]])
    below = amounts()
    below.matrix[A] = -1.0km^2
    @test_throws "negative" placed(WeightPerRecord(below))
    # Refused: a scale whose units do not cancel the weight's, a scale with no weight, numbers of
    # records beside a weight, a scale beside a rule, and a `novalue` in other units.
    @test_throws "give `scale` in km⁻²" FromOccurrences(records,
                                                        weight = WeightPerCell(amount),
                                                        scale = 2.0)
    @test_throws "give a `weight` too" FromOccurrences(records, scale = 2.0)
    @test_throws "not from `total` or `per_record`" FromOccurrences(records,
                                                                    weight = WeightPerCell(amount),
                                                                    total = 10)
    @test_throws "not both" FromOccurrences(records,
                                            weight = WeightPerCell(amount),
                                            scale = 1.0 / km^2) do rows, _, _
        return rows
    end
    @test_throws "in the layer's units, km²" WeightPerCell(amount, novalue = 3)
end

@testset "occurrence records are dated through the run's epoch and calendar" begin
    habitat = _bnghabitat()
    species = _species()
    at = _positions(habitat)[A]
    function records(column, values)
        n = length(values)
        table = DataFrame(species = fill("1", n),
                          decimalLatitude = fill(at[1], n),
                          decimalLongitude = fill(at[2], n))
        table[!, column] = values
        return table
    end
    seen = Ref{Any}()
    keep = (rows, _, _) -> (seen[] = rows; rows)
    epoch = Date(2020, 1, 1)
    # A date, a month, a year and a range each at its start; what is not a date dropped in a notice.
    dated = records(:eventDate,
                    [
                        "2019-05-03",
                        "2019-05",
                        "2019",
                        "2019-03/2019-05",
                        "",
                        "nope"
                    ])
    @test_logs((:info,
                r"^2 records with no date that could be read were not placed\.$"),
               build_ecosystem(FromOccurrences(keep, dated, date = :eventDate),
                               species, habitat, seed = 1, epoch = epoch))
    elapsed(d) = EcoSISTEM._elapsedat(epoch, d, EcoSISTEM.ExactDates())
    @test seen[].time ≈
          elapsed.([Date(2019, 5, 3), Date(2019, 5, 1), Date(2019, 1, 1),
                       Date(2019, 3, 1)])
    # A year's, a month's and a day's columns, a month left out taking the year's start.
    parts = DataFrame(species = ["1", "1"], decimalLatitude = [at[1], at[1]],
                      decimalLongitude = [at[2], at[2]], year = [2018, 2019],
                      month = [6, missing])
    @test_throws "no epoch" _placed(FromOccurrences(parts,
                                                    date = (:year, :month)),
                                    species, habitat)
    eco = build_ecosystem(FromOccurrences(keep, parts, date = (:year, :month)),
                          species, habitat, seed = 1, epoch = epoch,
                          verbosity = :silent)
    @test seen[].time ≈ elapsed.([Date(2018, 6, 1), Date(2019, 1, 1)])
    # Elapsed times as given; and at build, nothing after the start.
    timed = records(:elapsed, [-1.0year, -2.0year])
    eco = _placed(FromOccurrences(timed, time = :elapsed), species, habitat)
    @test sum(eco.abundances.matrix) == 2
    @test_throws "give `date = nothing`" build_ecosystem(FromOccurrences(records(:eventDate,
                                                                                 ["2021-01-01"]),
                                                                         date = :eventDate),
                                                         species, habitat,
                                                         epoch = epoch,
                                                         verbosity = :silent)
    @test_throws "give `time = nothing`" _placed(FromOccurrences(records(:elapsed,
                                                                         [1.0year]),
                                                                 time = :elapsed),
                                                 species, habitat)
    @test_throws "by `date` or by `time`, not both" FromOccurrences(timed,
                                                                    date = :eventDate,
                                                                    time = :elapsed)
    @test_throws "no `eventDate` column" _placed(FromOccurrences(timed,
                                                                 date = :eventDate),
                                                 species, habitat)
end

# How far each occupied cell of `eco` is filled - the largest, over the resources, of what its
# individuals demand over what it supplies - and how far rounding each species to whole individuals can
# move that: half of every species' demand, over the supply.
function _fullness(eco)
    supplies = map(b -> b.matrix, values(EcoSISTEM.getsupply(eco.habitat)))
    demands = map(d -> d.resource, values(eco.spplist.demand))
    linear = LinearIndices(EcoSISTEM.getgridshape(eco.habitat))
    abundances = eco.abundances.matrix
    filled = Dict{CartesianIndex{2}, NamedTuple}()
    for cell in CartesianIndices(size(MASK))
        (MASK[cell] && sum(abundances[:, linear[cell]]) > 0) || continue
        ratios = [ustrip(NoUnits,
                         sum(abundances[sp, linear[cell]] * d[sp]
                             for sp in axes(abundances, 1)) / b[cell])
                  for (b, d) in zip(supplies, demands)]
        slack = maximum(ustrip(NoUnits, sum(d) / 2 / b[cell])
                        for (b, d) in zip(supplies, demands))
        filled[cell] = (fullness = maximum(ratios), slack = slack)
    end
    return filled
end

# Whether every occupied cell of `eco` is full to within rounding.
_isfull(eco) = all(f -> abs(f.fullness - 1) <= f.slack, values(_fullness(eco)))

@testset "ToCapacity fills each cell to what its scarcest resource supports" begin
    habitat = _habitat()
    demand = 1.0e3kJ / day
    species = _species(demand = DemandSpec(SolarRadiation, demand))
    supply = EcoSISTEM.getsupply(habitat).matrix
    linear = LinearIndices(size(MASK))
    cell = CartesianIndex(3, 4)
    capacity = ustrip(NoUnits, supply[cell] / demand)
    shares = (species = ["1", "2"], cell = [cell, cell], count = [1, 3])
    eco = _placed(ToCapacity(FromTable(shares)), species, habitat)
    @test eco.abundances.matrix[1, linear[cell]] == round(Int, capacity / 4)
    @test eco.abundances.matrix[2, linear[cell]] == round(Int, 3 * capacity / 4)
    @test sum(eco.abundances.matrix) ==
          round(Int, capacity / 4) + round(Int, 3 * capacity / 4)
    # A fraction of capacity; weights for shares; a time ignored, with no epoch to place it.
    half = _placed(ToCapacity(FromTable(shares), fraction = 0.5), species,
                   habitat)
    @test half.abundances.matrix[2, linear[cell]] ==
          round(Int, 3 * capacity / 8)
    weighted = _placed(ToCapacity((_, _, _) -> (species = ["1", "2"],
                                                cell = [cell, cell],
                                                count = [0.25, 0.75],
                                                time = [Date(2030), Date(2031)])),
                       species, habitat)
    @test weighted.abundances.matrix == eco.abundances.matrix
    # Species demanding more take more of the supply each: the cell holds fewer of them.
    species.demand.resource[2] = 3 * demand
    mixed = _placed(ToCapacity(FromTable(shares)), species, habitat)
    total = capacity / ((1 + 3 * 3) / 4)
    @test mixed.abundances.matrix[1, linear[cell]] == round(Int, total / 4)
    @test mixed.abundances.matrix[2, linear[cell]] == round(Int, 3 * total / 4)
    @test _isfull(mixed)
    # A species' name fills every active cell with that species alone.
    wheat = _placed(ToCapacity("2"), species, habitat)
    @test _perspecies(wheat)[[1, 3, 4]] == [0, 0, 0]
    @test count(>(0), wheat.abundances.matrix[2, :]) == count(MASK)
    @test _isfull(wheat)
    # Every active cell filled by a placement putting some of each species everywhere.
    full = _placed(ToCapacity(BySupply()), species, habitat)
    @test length(_fullness(full)) == count(MASK)
    @test _isfull(full)
    @test sum(full.abundances.matrix[:, INACTIVE]) == 0
end

@testset "ToCapacity is bound by the scarcest resource" begin
    water = GradientSpec(1.0mm / day, 5.0mm / day, axis = Precipitation)
    habitat = _habitat(supply = (SOLAR, water))
    species = _species(demand = (DemandSpec(SolarRadiation, 1.0e3kJ / day),
                                 DemandSpec(Precipitation,
                                            7.0e3Unitful.L / day)))
    eco = _placed(ToCapacity(BySupply()), species, habitat)
    @test _isfull(eco)
    # Each resource alone would hold more somewhere than the two together: each binds somewhere.
    supplies = map(b -> b.matrix, values(EcoSISTEM.getsupply(habitat)))
    ratio = [ustrip(NoUnits, supplies[1][c] / 1.0e3kJ * day) /
             ustrip(NoUnits, supplies[2][c] / (7.0e3Unitful.L) * day)
             for c in CartesianIndices(size(MASK)) if MASK[c]]
    @test minimum(ratio) < 1 < maximum(ratio)
end

@testset "ToCapacity reports what it cannot fill, and refuses what is not a placement" begin
    habitat = _habitat()
    species = _species(demand = DemandSpec(SolarRadiation, 1.0e3kJ / day))
    shares = (species = ["1", "zz", "2"],
              cell = [
                  CartesianIndex(3, 4),
                  CartesianIndex(3, 4),
                  CartesianIndex(1, 1)
              ],
              count = [1, 1, 1])
    eco = @test_logs((:warn,
                      r"^1 row of 1 name of no species in the species list was not placed: `zz` \(1\)\.$"),
                     (:warn,
                      r"^1 row naming a cell the placement may not place in was not placed\.$"),
                     _placed(ToCapacity(FromTable(shares)), species, habitat))
    @test _perspecies(eco)[[2, 3, 4]] == [0, 0, 0] && _perspecies(eco)[1] > 0
    # Species demanding nothing are bounded by no supply.
    free = _species(demand = DemandSpec(SolarRadiation, 0.0kJ / day))
    empty = @test_logs((:warn,
                        r"^1 cell whose species demand nothing of any resource"),
                       _placed(ToCapacity(FromTable((species = ["1"],
                                                     cell = [CartesianIndex(3,
                                                                            4)],
                                                     count = [1]))), free,
                               habitat))
    @test sum(empty.abundances.matrix) == 0
    @test_throws "finite, non-negative number; got -1" ToCapacity(BySupply(),
                                                                  fraction = -1)
    @test_throws "a placement says how many individuals start where" ToCapacity(3)
    @test_throws "a share is a finite, non-negative number" _placed(ToCapacity((_,
                                                                                _,
                                                                                _) -> (species = ["1"],
                                                                                       cell = [CartesianIndex(3,
                                                                                                              4)],
                                                                                       count = [-1.0])),
                                                                    species,
                                                                    habitat)
    # Placed again by `repopulate!`, as it was.
    eco = _placed(ToCapacity(BySupply(), fraction = 0.5), species, habitat)
    start = copy(eco.abundances.matrix)
    repopulate!(eco, verbosity = :silent)
    @test eco.abundances.matrix == start
end

@testset "a function places what it returns, given the context" begin
    habitat = _habitat()
    species = _species()
    seen = Ref{Any}(nothing)
    eco = build_ecosystem(species, habitat, seed = 1, verbosity = :silent
                          ) do species,
                               habitat,
                               context
        seen[] = context
        return (species = ["2"], cell = [CartesianIndex(3, 4)], count = [12])
    end
    @test sum(eco.abundances.matrix) == 12
    context = seen[]
    @test context isa PlacementContext
    @test sort(context.cells) == sort(findall(MASK))
    @test sort(collect(keys(context.streams))) == species.names
    @test context.time == 0.0s && isnothing(context.epoch)
    @test provenance(eco).context.run.placement == "an anonymous function"
    # What it found to say is said; anything but a table or a result is refused.
    told = EcoSISTEM.Problem(EcoSISTEM.ProblemWarning(), :noted,
                             "a note from the placement")
    @test_logs (:warn, "a note from the placement") build_ecosystem(species,
                                                                    habitat,
                                                                    verbosity = :silent
                                                                    ) do _,
                                                                         _,
                                                                         _
        return PlacedAbundances((species = String[],
                                 cell = CartesianIndex{2}[]),
                                problems = [told])
    end
    @test_throws "returned a value of type Int64" build_ecosystem((_, _,
                                                                   _) -> 3,
                                                                  species,
                                                                  habitat,
                                                                  verbosity = :silent)
end

@testset "each species is placed from a stream of its own, keyed on its key in the pool" begin
    habitat = _habitat()
    placement = BySupply(each = 50) + BySuitability(each = 10)
    eco = _placed(placement, _species(), habitat)
    # Placing leaves the streams the dynamics draw from at their seeded start.
    @test eco.rngs == EcoSISTEM.makerngs(1, 4)
    # A species more, or one made inactive, leaves every other species' placement as it is.
    more = _placed(placement, _species(5), habitat)
    @test more.abundances.matrix[1:4, :] == eco.abundances.matrix
    fewer = _species()
    fewer.traits.active[2] = false
    lessone = _placed(placement, fewer, habitat)
    @test _perspecies(lessone)[2] == 0
    @test lessone.abundances.matrix[[1, 3, 4], :] ==
          eco.abundances.matrix[[1, 3, 4], :]
    # A species removed from the pool moves the rest up a position but leaves their keys, and so their
    # placements.
    pool = SpeciesPool(4, verbosity = :silent)
    remove_species!(pool, "2", verbosity = :silent)
    removed = _placed(placement, _species(pool), habitat)
    @test removed.abundances.matrix == eco.abundances.matrix[[1, 3, 4], :]
    # A copy is placed from its own key's stream, after its original's.
    pool = SpeciesPool(4, verbosity = :silent)
    copy_species!(pool, "2", as = "2 copy")
    copied = _placed(placement, _species(pool), habitat)
    @test copied.abundances.matrix[5, :] == more.abundances.matrix[5, :]
    @test copied.abundances.matrix[5, :] != copied.abundances.matrix[2, :]
end

@testset "repopulate! places again, by the ecosystem's placement or another" begin
    habitat = _habitat()
    species = _species()
    placement = BySupply(total = 400) + BySuitability(each = 10)
    eco = _placed(placement, species, habitat)
    start = copy(eco.abundances.matrix)
    simulate!(eco, 2month_mean_duration, 1month_mean_duration)
    # The habitat is unchanged, so its own placement and seed place the start again, leaving the
    # clock and the dynamics' streams where the run left them.
    streams = copy.(eco.rngs)
    repopulate!(eco, verbosity = :silent)
    @test eco.abundances.matrix == start
    @test eco.rngs == streams && eco.elapsed > 0.0s
    other = copy(repopulate!(eco, seed = 99, verbosity = :silent).abundances.matrix)
    @test other != start && sum(other) == sum(start)
    repopulate!(eco, NoAbundance(), verbosity = :silent)
    @test sum(eco.abundances.matrix) == 0 && eco.placement === NoAbundance()
    @test provenance(eco).context.run.placement == "NoAbundance()"
    eco.placement = nothing
    @test_throws "records no placement" repopulate!(eco)
end

# The bytes placing `p` into `eco` allocates, inside a function so that nothing global is counted.
function _placingbytes(eco, p, seed)
    context = EcoSISTEM._placementcontext(eco, seed)
    return @allocated EcoSISTEM._placeinto!(eco, p, context, :silent,
                                            clear = true)
end

# A weighted placement draws straight into the abundances, with no table. Its draws are the ones its
# table holds, which a sum with `NoAbundance()` places, at build, by `repopulate!`, at an
# `AddAbundance` firing over part of the grid, and for an `AddSpecies` arrival. So a build holds
# little more than its abundances, where a table holds a row for each species in each cell it
# reaches.
@testset "a weighted placement draws straight into the abundances, as its table would place them" begin
    habitat = _habitat()
    # Four species, the second not active; built again for each run an arrival changes.
    function fresh()
        species = _species()
        species.traits.active[2] = false
        return species
    end
    species = fresh()
    layer = UniformSpec(2.0, axis = EcoSISTEM.NicheAxis)
    for p in (BySupply(total = 10^6), BySupply(each = 200),
        BySupply(each = DrawFrom(Poisson(50))), BySuitability(each = 50),
        ByLayer(layer, total = 3000),
        Weighted(BySupply(), BySuitability(), total = 2000))
        direct = _placed(p, species, habitat, seed = 4)
        table = _placed(p + NoAbundance(), species, habitat, seed = 4)
        @test direct.abundances.matrix == table.abundances.matrix
        @test _perspecies(direct)[2] == 0 && sum(direct.abundances.matrix) > 0
        repopulate!(direct, seed = 9, verbosity = :silent)
        repopulate!(table, seed = 9, verbosity = :silent)
        @test direct.abundances.matrix == table.abundances.matrix
    end
    @test _placed(DefaultEcosystem(), species, habitat, seed = 7).abundances.matrix ==
          _placed(BySupply(total = 3000) + NoAbundance(), species, habitat,
                  seed = 7).abundances.matrix
    # A species with no weight on any cell is said to be placed nowhere, as the table path says it.
    nowhere = ByLayer(UniformSpec(0.0, axis = EcoSISTEM.NicheAxis), total = 100)
    @test_logs (:warn, r"100 individuals were not placed") _placed(nowhere,
                                                                   species,
                                                                   habitat)
    @test_logs (:warn, r"100 individuals were not placed") _placed(nowhere +
                                                                   NoAbundance(),
                                                                   species,
                                                                   habitat)
    # At a firing over the grid's western cells, and for the inactive species arriving there.
    west = CellMask([x <= 2 for y in 1:7, x in 1:5])
    for (op, viatable) in ((AddAbundance(BySupply(total = 500)),
                            AddAbundance(BySupply(total = 500) + NoAbundance())),
                           (AddSpecies("2", BySupply(each = 40)),
                            AddSpecies("2",
                                       BySupply(each = 40) + NoAbundance())))
        runs = map((op, viatable)) do operation
            eco = _placed(BySupply(total = 2000), fresh(), habitat, seed = 4)
            simulate!(eco, 2month_mean_duration, 1month_mean_duration,
                      intervention = Intervention(AtTime(1month_mean_duration),
                                                  west, operation))
            return eco.abundances.matrix
        end
        @test runs[1] == runs[2]
    end
    # Many species, each reaching every cell: the table path allocates a row for each of them, and
    # the direct path little more than one draw at a time.
    many = _species(200)
    eco = _placed(NoAbundance(), many, habitat)
    p = BySupply(each = 10^6)
    _placingbytes(eco, p, 4)
    _placingbytes(eco, p + NoAbundance(), 4)
    direct = _placingbytes(eco, p, 4)
    table = _placingbytes(eco, p + NoAbundance(), 4)
    @test direct < table / 10
    @test direct < 2 * sizeof(eco.abundances.matrix)
end

@testset "a placement shows as it is written" begin
    @test sprint(show, BySupply()) == "BySupply()"
    @test sprint(show, BySupply(SolarRadiation, :summer, total = 5)) ==
          "BySupply(SolarRadiation, :summer, total = 5)"
    @test sprint(show, BySuitability(each = TraitColumn(:n))) ==
          "BySuitability(each = TraitColumn(:n))"
    @test sprint(show, Weighted(BySupply(), BySuitability(), total = 9)) ==
          "Weighted(BySupply(), BySuitability(), total = 9)"
    @test sprint(show, ToCapacity(BySupply())) == "ToCapacity(BySupply())"
    @test sprint(show, ToCapacity(BySupply(), fraction = 0.5)) ==
          "ToCapacity(BySupply(), fraction = 0.5)"
    @test sprint(show, BySupply() + NoAbundance() + BySuitability()) ==
          "BySupply() + NoAbundance() + BySuitability()"
    @test BySupply() + (NoAbundance() + BySuitability()) isa
          SumOfPlacements{<:NTuple{3, Any}}
    mktempdir() do dir
        path = joinpath(dir, "records.csv")
        write(path, "species\n")
        @test sprint(show, FromOccurrences(path)) ==
              "FromOccurrences($(repr(path)))"
        @test sprint(show,
                     FromOccurrences(path, species = :scientificName,
                                     date = (:year, :month), snap = 3km,
                                     per_record = 10,
                                     level = EcoSISTEM.ProblemNotice())) ==
              "FromOccurrences($(repr(path)), species = :scientificName, date = (:year, :month), snap = 3 km, per_record = 10, level = ProblemNotice())"
    end
    @test startswith(sprint(show,
                            FromOccurrences((_, _, _) -> nothing,
                                            (species = ["a"],))),
                     "FromOccurrences(an anonymous function, a NamedTuple")
    @test sprint(show,
                 FromTable((species = ["a", "b"],
                            cell = [CartesianIndex(1, 1),
                                CartesianIndex(1, 2)]))) ==
          "FromTable(a NamedTuple, 2 rows)"
    @test sprint(show,
                 PlacedAbundances((species = ["a"],
                                   cell = [CartesianIndex(1, 1)]))) ==
          "PlacedAbundances(1 row, 0 records, 0 problems)"
end

end
