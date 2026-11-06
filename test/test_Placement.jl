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

@testset "a placement shows as it is written" begin
    @test sprint(show, BySupply()) == "BySupply()"
    @test sprint(show, BySupply(SolarRadiation, :summer, total = 5)) ==
          "BySupply(SolarRadiation, :summer, total = 5)"
    @test sprint(show, BySuitability(each = TraitColumn(:n))) ==
          "BySuitability(each = TraitColumn(:n))"
    @test sprint(show, Weighted(BySupply(), BySuitability(), total = 9)) ==
          "Weighted(BySupply(), BySuitability(), total = 9)"
    @test sprint(show, BySupply() + NoAbundance() + BySuitability()) ==
          "BySupply() + NoAbundance() + BySuitability()"
    @test BySupply() + (NoAbundance() + BySuitability()) isa
          SumOfPlacements{<:NTuple{3, Any}}
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
