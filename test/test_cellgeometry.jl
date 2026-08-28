# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for the grid-inspection accessors in `src/cellgeometry.jl` — `getcellsizes` /
# `getcellareas` / `getgridshape` / `getcellat` / `getspeciesstorage`.
#
# Must go through `Pkg.test`: the fixtures name `WorldClim{BioClim}`, a weak dependency.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_cellgeometry.jl"])'

module TestCellGeometry

using EcoSISTEM
using EcoSISTEM: getspeciesstorage, SpatialLocation,
                 in_memory_raster, SyntheticData
using EcoSISTEM.Units    # `day` — and `Dates` also exports one, so it must be named
using Unitful, Unitful.DefaultSymbols
using DimensionalData
using DimensionalData: Y, X, dims
using Rasters
using RasterDataSources
using Statistics: mean
using Test

const L = DimensionalData.Lookups

# A geographic fixture, `Intervals(Start)` like everything this package builds.
#
# **It spans 0°–85° deliberately.** Over a few degrees `sin` is near-linear, so the mean solid angle
# and the solid angle *at* the mean latitude agree to several digits — and the assertion that they
# differ would pass against an implementation that computed the wrong one. The fixture has to cover
# enough latitude for the nonlinearity to bite.
function _geo(lats, longs)
    mk(D, v) = D(Rasters.Projected(collect(v),
                                   sampling = L.Intervals(L.Start()),
                                   crs = Rasters.EPSG(4326),
                                   order = L.ForwardOrdered(),
                                   span = L.Regular(v[2] - v[1])))
    r = ClimateRaster(SyntheticData,
                      DimArray(fill(290.0K, length(lats), length(longs)),
                               (mk(Y, lats), mk(X, longs))))
    return StudyArea(regime = in_memory_raster(r, axis = Temperature),
                     verbosity = :silent)
end

# The closed form, written out independently of the implementation so the test is a check rather than
# a restatement: a cell spanning `φ1…φ2` by `dlong` subtends `Δλ · (sin φ₂ − sin φ₁)`.
_Ω(φ1, φ2, dlong) = ustrip(u"rad", dlong) * (sin(φ2) - sin(φ1)) * u"sr"

# A **projected** fixture carrying a real CRS (British National Grid), so the inverse transform has
# something to transform through. `syn` deliberately has none — a synthetic grid is plain cell
# indices with no real-world position — which is a different case and tested as one.
function _proj(ys, xs)
    mk(D, v) = D(Rasters.Projected(collect(v),
                                   sampling = L.Intervals(L.Start()),
                                   crs = Rasters.EPSG(27700),
                                   order = L.ForwardOrdered(),
                                   span = L.Regular(v[2] - v[1])))
    r = ClimateRaster(EcoSISTEM.SyntheticData,
                      DimArray(fill(290.0K, length(ys), length(xs)),
                               (mk(Y, ys), mk(X, xs))))
    return StudyArea(regime = in_memory_raster(r, axis = Temperature),
                     verbosity = :silent)
end

@testset "CellGeometry" begin
    lats = (0.0:5.0:80.0) .* °          # 17 rows, cells [φ, φ+5°)
    longs = (0.0:5.0:20.0) .* °         # 5 columns
    geo = _geo(lats, longs)
    syn = StudyArea(extent = (60km, 60km), cellsize = 10km,
                    verbosity = :silent)

    @testset "sizes: the grid's own units, one value per cell" begin
        at(place) = EcoSISTEM.getcellat(geo, place)

        @test EcoSISTEM.getcellsizes(geo).y[1, 1] == 5.0°
        @test EcoSISTEM.getcellsizes(°, geo).x[1, 1] == 5.0°
        @test EcoSISTEM.getcellsizes(syn).y[1, 1] == 10.0km
        @test EcoSISTEM.getcellsizes(m, syn).x[1, 1] == 10000.0m
        # Uniform extents compress to a `Fill`, so the answer costs the same however large the grid.
        @test EcoSISTEM.getcellsizes(syn).y isa EcoSISTEM.Fill
        @test size(EcoSISTEM.getcellsizes(syn).y) == EcoSISTEM.getgridshape(syn)

        # A metric size of a geographic grid: north-south is the same everywhere, east-west shrinks
        # as the cosine of latitude. The scalar API could express neither and refused outright; a
        # field says both, which is why it is materialised rather than compressed.
        metric = EcoSISTEM.getcellsizes(km, geo)
        @test size(metric.y) == EcoSISTEM.getgridshape(geo)
        @test allequal(metric.y)
        @test metric.x[1, 1] > metric.x[end, 1]
        @test !(metric.x isa EcoSISTEM.Fill)

        # Angular from a projected grid is now COMPUTED, by transforming each cell back through
        # its CRS — see the inverse-transform testset below. `syn` has no CRS at all (a synthetic
        # grid is plain cell indices), so there is genuinely nothing to transform through and
        # `nothing` is the answer rather than a gap. It says so silently, as an answer should.
        @test_logs isnothing(EcoSISTEM.getcellsizes(°, syn))
        @test isnothing(EcoSISTEM.getcellsizes(UniformSpec(290.0K,
                                                           axis = Temperature)))
    end

    @testset "areas: a solid angle, NOT the product of the extents" begin
        s = EcoSISTEM.getcellsizes(°, geo)
        angular = EcoSISTEM.getcellareas(°^2, geo)
        # The defect this testset exists for: `5° × 5°` is `25 °²` as a coordinate product, and
        # that is neither an area nor a solid angle.
        @test mean(angular) != s.y[1, 1] * s.x[1, 1]
        # The point value, against the closed form — the cell labelled 0° spans 0°–5° under
        # `Intervals(Start)`, so it is `sin(5°) − sin(0°)`, *not* the centred form.
        @test angular[EcoSISTEM.getcellat(geo, SpatialLocation(2.5°, 2.5°))] ≈
              uconvert(°^2, _Ω(0.0°, 5.0°, 5.0°))
        @test angular[EcoSISTEM.getcellat(geo, SpatialLocation(62.5°, 2.5°))] ≈
              uconvert(°^2, _Ω(60.0°, 65.0°, 5.0°))
        # It shrinks towards the poles — the property that makes it a solid angle at all — and the
        # field says so directly, where the scalar API needed a location to be asked twice.
        @test issorted(angular[:, 1], rev = true)
        @test !allequal(angular)
        # Every column of a row is the same: a solid angle depends on latitude, not longitude.
        @test allequal(angular[1, :])
    end

    @testset "areas: the mean is total ÷ count, not the value at the mean latitude" begin
        # The total telescopes — `Σ(sin φᵢ₊₁ − sin φᵢ) = sin φ_top − sin φ_bottom` — so the average
        # needs only the outer edges and the row count. Computed here from the fixture's own extent,
        # independently of how the implementation gets there.
        total = _Ω(0.0°, 85.0°, 5.0°)
        angular = EcoSISTEM.getcellareas(°^2, geo)
        @test mean(angular) ≈ uconvert(°^2, total / length(lats))
        # `sin` is nonlinear, so the mean is NOT the solid angle at the mean latitude. This is the
        # assertion the 0°–85° span exists for: over a few degrees the two agree.
        @test !isapprox(mean(angular),
                        angular[EcoSISTEM.getcellat(geo,
                                                    SpatialLocation(42.5°,
                                                                    2.5°))])
        # `mean` is Base's, not a keyword. That is the whole reason the accessor returns a field.
    end

    @testset "units: an area needs an AREA unit, in both spellings" begin
        @test first(EcoSISTEM.getcellareas(km^2, syn)) == 100.0km^2
        @test first(EcoSISTEM.getcellareas(m^2, syn)) == 1.0e8 * m^2
        @test first(EcoSISTEM.getcellareas(syn)) == 100.0km^2          # native
        @test unit(mean(EcoSISTEM.getcellareas(u"sr", geo))) == u"sr"
        @test mean(EcoSISTEM.getcellareas(u"sr", geo)) ≈
              uconvert(u"sr", mean(EcoSISTEM.getcellareas(°^2, geo)))
        # A length or an angle is a caller mistake, not an unanswerable question — so it errors
        # rather than returning `nothing`, and the message names the fix. This is what keeps the
        # first argument meaning the *same thing* in both functions: the unit of the answer.
        @test_throws ErrorException EcoSISTEM.getcellareas(km, syn)
        @test_throws ErrorException EcoSISTEM.getcellareas(°, geo)
        # …while `getcellsizes` is the one that takes a length or an angle.
        @test EcoSISTEM.getcellsizes(km, syn).y[1, 1] == 10.0km
    end

    @testset "everything that knows the grid gives the same answer" begin
        env = GridHabitat(regime = UniformSpec(290.0K,
                                               axis = Temperature),
                          supply = UniformSpec(1.0kJ / (km^2 * day),
                                               axis = SolarRadiation),
                          area = syn)
        for x in (env, env.regime, env.supply, syn.report, syn)
            # The supply reads the **grid**, not its own `size` field — which reports `1.0 m` on
            # every grid, the bug this accessor is here to make unreachable.
            @test EcoSISTEM.getcellsizes(x).y[1, 1] == 10.0km
            @test EcoSISTEM.getgridshape(x) == (6, 6)
        end
        # A spec has no grid of its own and says so, which is what lets callers filter on it.
        @test isnothing(EcoSISTEM.getcellsizes(UniformSpec(290.0K,
                                                           axis = Temperature)))
    end

    @testset "a place must be in the grid's own frame" begin
        # There are no `lat`/`long`/`mean` keyword combinations to refuse: a location is
        # `getcellat`'s single argument, and an average is Base's `mean` over the field. What can go
        # wrong is asking in the wrong frame.
        @test_throws ErrorException EcoSISTEM.getcellat(syn,
                                                        SpatialLocation(1.0°,
                                                                        1.0°))
        @test_throws ErrorException EcoSISTEM.getcellat(geo,
                                                        SpatialLocation(1.0km,
                                                                        1.0km))
        # And a place off the grid throws rather than answering: not on this grid is a caller
        # mistake, not a missing value.
        @test_throws Exception EcoSISTEM.getcellat(syn,
                                                   SpatialLocation(500.0km,
                                                                   0.0km))
    end

    # The point of this accessor is that a run can be *sized before it is built* — so the report
    # from `investigate_study_area` must give the same figure as the habitat eventually does.
    @testset "per-species storage, asked before and after building" begin
        env = GridHabitat(regime = UniformSpec(290.0K,
                                               axis = Temperature),
                          supply = UniformSpec(1.0kJ / (km^2 * day),
                                               axis = SolarRadiation),
                          area = syn)
        # Assert the arithmetic, not merely that the answers agree with each other: 60 km of 10 km
        # cells each way is 6 × 6 = 36 cells of `Int64`. A y/x mix-up cannot show here and does
        # not need to — this is a product over the whole grid, so it is symmetric in the two.
        @test getspeciesstorage(syn) == 36 * sizeof(Int64)
        @test getspeciesstorage(syn.report) == getspeciesstorage(syn)
        @test getspeciesstorage(env) == getspeciesstorage(syn)
        @test getspeciesstorage(env.regime) == getspeciesstorage(syn)

        # It must not drift from the report's own footprint — they are the same question, and the
        # accessor exists so callers need not reach into that field.
        @test getspeciesstorage(syn.report) == syn.report.footprint.perspecies

        # A rule rather than data has no grid to answer for, matching `getcellsizes`.
        @test isnothing(getspeciesstorage(UniformSpec(290.0K,
                                                      axis = Temperature)))
    end

    @testset "nothing, missing and a grid are three different answers" begin
        rule = UniformSpec(290.0K, axis = Temperature)
        lazy = ConstructedSpec(() -> nothing, axis = EcoSISTEM.NicheAxis)

        # The distinction `force` exists to let a caller act on. A rule adopts whatever grid it is
        # placed on and so has none — `nothing`. A data-backed spec HAS one, but only reading it,
        # potentially a download, would say what — `missing`. Collapsing the two would make "we did
        # not look" indistinguishable from "there is nothing to look at".
        for f in (EcoSISTEM.getgridshape, EcoSISTEM.getcellareas,
            EcoSISTEM.getcellsizes, EcoSISTEM.getcellcount,
            EcoSISTEM.getgridarea)
            @test isnothing(f(rule))
            @test ismissing(f(lazy))
        end
        @test EcoSISTEM.getgridshape(syn) == (6, 6)

        # `missing` as a *return* does not collide with `missing` as the "no unit given" argument:
        # they never occupy the same position, and each is read by the predicate at its own site.
        @test ismissing(EcoSISTEM.getcellareas(km^2, lazy))

        # A spec describing a computation cannot be read on its own — its grid is decided by the
        # area it is built with — so `force` says so rather than guessing.
        @test_throws ErrorException EcoSISTEM.getgridshape(lazy, force = true)
    end

    @testset "angular extents of a projected grid, by inverse CRS transform" begin
        # 100 km cells on the British National Grid, spanning enough of the country for the
        # transform's variation to be visible.
        bng = _proj((0.0:1.0e5:4.0e5)m, (0.0:1.0e5:3.0e5)m)

        ang = EcoSISTEM.getcellsizes(°, bng)
        @test size(ang.y) == EcoSISTEM.getgridshape(bng)
        @test unit(ang.y[1, 1]) === °

        # The property that makes this worth computing rather than approximating: a projected
        # cell's ANGULAR extent varies across the grid, even where its projected extent does not.
        # A single answer would be wrong everywhere it was not taken — which is why the old API
        # refused instead of averaging.
        @test !allequal(ang.x)
        # A degree of longitude spans more ground the further south you are, so a fixed-width cell
        # spans FEWER degrees there: the east-west angular extent grows northwards.
        @test ang.x[end, 1] > ang.x[1, 1]

        # Sanity against the geographic reading: 100 km is on the order of a degree, not a
        # thousandth or a hundred of one. This is a check that the units came back right, not a
        # restatement of the arithmetic.
        @test 0.5° < ang.y[1, 1] < 2.0°

        # The metric answer is untouched by any of this.
        @test EcoSISTEM.getcellsizes(bng).y[1, 1] == 1.0e5m
    end

    @testset "one route to a cell's area, and one step from a place to a cell" begin
        # A metric area of a geographic grid goes through the same `_cellareas` the **supply
        # constructor** uses. `Supply{A}(::ClimateRaster{WorldClim{BioClim}})` multiplies an areal
        # rate by each cell's true area, so a second computation here could silently disagree with
        # the one the model actually runs on.
        @test EcoSISTEM.getcellareas(km^2, geo)[:, 1] ≈
              vec(EcoSISTEM._cellareas(geo.report.active))

        # `getcellat` is the only step that turns a position into a grid position; asking one cell
        # is then indexing the field, which is why no other accessor takes a location.
        here = SpatialLocation(2.5°, 2.5°)
        @test EcoSISTEM.getcellarea(geo, here) ==
              EcoSISTEM.getcellareas(geo)[EcoSISTEM.getcellat(geo, here)]
        @test EcoSISTEM.getcellsize(geo, here).y ==
              EcoSISTEM.getcellsizes(geo).y[EcoSISTEM.getcellat(geo, here)]

        # The whole grid's area is the sum of its cells', which is correct on every projection —
        # multiplying one cell by a count would not be.
        @test EcoSISTEM.getgridarea(geo) ≈ sum(EcoSISTEM.getcellareas(geo))
        @test EcoSISTEM.getcellcount(syn) == prod(EcoSISTEM.getgridshape(syn))
    end
end

end
