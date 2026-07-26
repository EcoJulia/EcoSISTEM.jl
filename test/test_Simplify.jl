# SPDX-License-Identifier: LGPL-3.0-or-later

module TestSimplify

using EcoSISTEM
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using Unitful, Unitful.DefaultSymbols
using AxisArrays
using RasterDataSources
using ArchGDAL
using Test

# Small synthetic climate rasters (no downloads) for the data-driven builders.
function _testraster(T, values; lat = (0.0:1.0:4.0) .* °,
                     long = (0.0:1.0:4.0) .* °)
    return ClimateRaster(T,
                         AxisArray(values, Axis{:latitude}(lat),
                                   Axis{:longitude}(long)))
end

# A synthetic single-polygon shapefile (a rectangle over `(x1,y1)-(x2,y2)` in `sr`, or with no CRS
# metadata at all if `sr = nothing`), written to a temp dir for `shapemask` tests — no checked-in
# binary fixture needed.
function _testshapefile(x1, y1, x2, y2; sr = nothing)
    path = joinpath(mktempdir(), "test.shp")
    poly = ArchGDAL.createpolygon([[
                                      (x1, y1),
                                      (x2, y1),
                                      (x2, y2),
                                      (x1, y2),
                                      (x1, y1)
                                  ]])
    ArchGDAL.create(path; driver = ArchGDAL.getdriver("ESRI Shapefile")
                    ) do dataset
        makelayer = layer -> begin
            ArchGDAL.addfeature(layer) do feature
                return ArchGDAL.setgeom!(feature, poly)
            end
            ArchGDAL.copy(layer, dataset = dataset)
        end
        if sr === nothing
            ArchGDAL.createlayer(makelayer; name = "test",
                                 geom = ArchGDAL.wkbPolygon)
        else
            ArchGDAL.createlayer(makelayer; name = "test",
                                 geom = ArchGDAL.wkbPolygon,
                                 spatialref = sr)
        end
    end
    return path
end

@testset "build_environment" begin
    extent = (5km, 5km)     # with cellsize 1km -> a 5×5 grid
    cellsize = 1km

    @testset "cell-count derivation" begin
        u = build_environment(extent, cellsize, 298.0K)
        @test size(u.regime.matrix) == (5, 5)
        # extent not a whole number of cells warns
        @test_logs (:warn,) match_mode=:any build_environment((5km, 5km), 2km,
                                                              298.0K)

        # A square extent can't distinguish x from y (both dimensions have the
        # same count), which is exactly why a historical X/Y mixup between this
        # synthetic path and the data-driven build_environment path went
        # undetected. `extent` is documented as (x, y); with a non-square
        # extent the resulting matrix must be (ny, nx) — dimension 1 = y —
        # matching the convention used throughout Generate.jl/Habitats.jl, not
        # (nx, ny).
        nonsquare = build_environment((12km, 4km), cellsize, 298.0K)
        @test size(nonsquare.regime.matrix) == (4, 12)
        @test size(nonsquare.active) == (4, 12)
    end

    @testset "dispatch on the regime argument type" begin
        # scalar temperature -> flat temperature regime
        u = build_environment(extent, cellsize, 298.0K)
        @test u isa GridHabitat
        @test EcoSISTEM.iscontinuous(u.regime)
        @test u.supply isa SolarSupply

        # temperature pair -> gradient (linear or peaked)
        g = build_environment(extent, cellsize, (274.0K, 303.0K))
        @test EcoSISTEM.iscontinuous(g.regime)
        gp = build_environment(extent, cellsize, (274.0K, 303.0K);
                               peaked = true)
        @test EcoSISTEM.iscontinuous(gp.regime)

        # integer -> discrete niches
        niche = build_environment(extent, cellsize, 3)
        @test !EcoSISTEM.iscontinuous(niche.regime)

        # length pair with no supply -> rainfall used as a WaterSupply
        rain = build_environment(extent, cellsize, (0.0mm, 100.0mm))
        @test rain.supply isa WaterSupply
    end

    @testset "supply unit picks type; active honoured" begin
        water = build_environment(extent, cellsize, 298.0K;
                                  supply = 1.0mm / day)
        @test water.supply isa WaterSupply

        active = fill(true, (5, 5))
        active[1, 1] = false
        env = build_environment(extent, cellsize, 298.0K; active = active)
        @test env.active == active
    end

    @testset "unsupported regime argument has no method" begin
        @test_throws MethodError build_environment(extent, cellsize, :banana)
    end
end

@testset "build_environment from data" begin
    temp = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))

    @testset "continuous regime + supply from quantity/raster" begin
        env = build_environment(temp)
        @test env isa GridHabitat
        @test EcoSISTEM.iscontinuous(env.regime)
        @test size(env.regime.matrix) == (5, 5)
        @test env.supply isa SolarSupply                    # default 1kJ/km^2/day

        water = _testraster(WorldClim{BioClim}, fill(50.0mm / day, 5, 5))
        @test build_environment(temp; supply = water).supply isa WaterSupply
    end

    @testset "discrete (land cover) regime" begin
        landcover = _testraster(EarthEnv{LandCover},
                                Float64.(repeat(1:5, 1, 5)))
        env = build_environment(landcover)
        @test !EcoSISTEM.iscontinuous(env.regime)
    end

    @testset "grid: default no resample, size= resamples with a warning" begin
        @test size(build_environment(temp).regime.matrix) == (5, 5)  # native
        @test_logs (:warn,) match_mode=:any build_environment(temp;
                                                              size = 300.0km)
    end

    @testset "masks" begin
        withnan = fill(290.0K, 5, 5)
        withnan[1, 1] = NaN * K
        dm = datamask(_testraster(WorldClim{BioClim}, withnan))
        @test dm isa AxisArray
        @test !dm[1, 1] && dm[2, 2]

        landcover = _testraster(EarthEnv{LandCover},
                                Float64.(repeat(1:5, 1, 5)))
        lm = landmask(landcover; water = [1])
        @test count(lm) == 20                                # class 1 excluded (one row)
        @test build_environment(landcover; active = lm) isa GridHabitat

        env = build_environment(temp; active = circlemask(radius = 100.0km))
        @test env.active isa Matrix{Bool}
    end

    @testset "shapemask" begin
        # No CRS metadata (no .prj) — treated as already WGS84 lat/long. Polygon covers lat/long
        # [-0.5,2.5] x [-0.5,2.5], so grid points 0,1,2 (of `temp`'s native 0:1:4 grid) are inside.
        path = _testshapefile(-0.5, -0.5, 2.5, 2.5)
        env = build_environment(temp; active = shapemask(path))
        @test env.active isa Matrix{Bool}
        @test env.active == [i <= 3 && j <= 3 for i in 1:5, j in 1:5]

        # British National Grid (EPSG:27700) — reprojects to ~lat 55.72-55.81, long -4.39 to -4.23.
        bng = ArchGDAL.importEPSG(27700)
        bngpath = _testshapefile(250000.0, 650000.0, 260000.0, 660000.0;
                                 sr = bng)
        scotland = _testraster(WorldClim{BioClim}, fill(290.0K, 2, 2);
                               lat = [55.75, 56.5] .* °,
                               long = [-4.6, -4.3] .* °)
        scotenv = build_environment(scotland; active = shapemask(bngpath))
        @test scotenv.active == [false true; false false]
    end

    @testset "geodesy: area-preserving cell size + top/bottom report" begin
        spanning = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5);
                               lat = (0.0:10.0:40.0) .* °,
                               long = (0.0:10.0:40.0) .* °)
        @test_logs (:info,) match_mode=:any build_environment(spanning)
    end

    @testset "cell size dispatches on coordinate units" begin
        deg = (0.0:1.0:2.0) .* °
        km_axis = (0.0:1.0:2.0) .* km
        # projected (length) coordinates -> planar, exactly sqrt(Δlat × Δlong)
        @test EcoSISTEM._cellsize(km_axis, km_axis) ≈ 1.0km
        # geographic (degree) coordinates -> spherical, shrinking toward the poles
        @test EcoSISTEM._cellsize((58.0:1.0:60.0) .* °, deg) <
              EcoSISTEM._cellsize(deg, deg)
    end

    @testset "attach units to a raster layer" begin
        # a bare (unitless) layer gets its unit, preserving axes and source type
        bare = _testraster(WorldClim{BioClim}, fill(290.0, 5, 5))
        united = EcoSISTEM._attach_unit(bare, u"K")
        @test eltype(united.array) <: Unitful.Temperature
        @test united isa ClimateRaster{WorldClim{BioClim}}
        @test AxisArrays.axisnames(united.array) == (:latitude, :longitude)
        # NoUnits is a no-op -> stays bare numbers (what discrete land cover wants)
        landcover = _testraster(EarthEnv{LandCover},
                                Float64.(repeat(1:5, 1, 5)))
        @test eltype(EcoSISTEM._attach_unit(landcover, NoUnits).array) ==
              Float64
    end
end

@testset "multi-raster environments" begin
    temp = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))
    rain = _testraster(WorldClim{BioClim}, fill(50.0mm, 5, 5))
    landcover = _testraster(EarthEnv{LandCover}, Float64.(repeat(1:5, 1, 5)))

    @testset "RegimeCollection2 + SupplyCollection2" begin
        env = build_environment((temp, rain);
                                supply = (1.0kJ / (km^2 * day), 1.0mm / day))
        @test env.regime isa EcoSISTEM.RegimeCollection2
        @test env.supply isa EcoSISTEM.SupplyCollection2
        @test Base.size(env.regime.one.matrix) == (5, 5)
    end

    @testset "mixed continuous + discrete, and a 3-tuple" begin
        env = build_environment((temp, landcover))
        @test EcoSISTEM.iscontinuous(env.regime) == [true, false]
        @test build_environment((temp, rain, landcover)).regime isa
              EcoSISTEM.RegimeCollection3
    end

    @testset "arity errors" begin
        @test_throws ErrorException build_environment((temp, rain, landcover,
                                                       temp))
        @test_throws ErrorException build_environment((temp,))
    end

    @testset "default active = AND of non-missing cells" begin
        withnan = fill(290.0K, 5, 5)
        withnan[1, 1] = NaN * K
        env = build_environment((_testraster(WorldClim{BioClim}, withnan),
                                 rain))
        @test env.active isa Matrix{Bool}
        @test !env.active[1, 1] && env.active[2, 2]
    end
end

@testset "_default_suitability derives NF from the tolerance for every tolerance kind" begin
    # Generalises the `NicheTolerance` fix: `DiscreteTolerance`/`LandCoverTolerance` must also take their
    # nichefit's `NF` from the tolerance, not `eltype(regime)` — otherwise nichefit trivially mirrors
    # whatever the regime happens to be, and a genuine tolerance/regime disagreement goes uncaught.
    disc = DiscreteTolerance(fill(1, 5))
    regime = simplenichehabitat(3, (5, 5), 10000.0kJ / km^2 / day, 25.0km^2).regime
    @test EcoSISTEM._default_suitability(disc, regime) isa
          MatchSuitability{eltype(disc)}

    lc = LandCoverTolerance([[0.2, 0.8] for _ in 1:5])
    @test EcoSISTEM._default_suitability(lc, regime) isa
          LandCoverSuitability{eltype(lc)}
end

@testset "niche-axis threading through the data path" begin
    # the shipped layer table resolves a source+code to its niche axis (no download)
    @test EcoSISTEM._specaxis((WorldClim{BioClim}, 1)) === Temperature
    @test EcoSISTEM._specaxis((WorldClim{BioClim}, 12)) === Precipitation
    # a bare raster has dropped its code, so no axis
    temp = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))
    @test EcoSISTEM._specaxis(temp) === EcoSISTEM.Unclassified
    # an explicit axis tags the data regime; the default leaves it unclassified
    @test EcoSISTEM.axisof(build_environment(temp; axis = Temperature).regime) ===
          Temperature
    @test EcoSISTEM.axisof(build_environment(temp).regime) ===
          EcoSISTEM.Unclassified
    # a matched-axis species assembles; a mismatched one is a clear build-time error
    env = build_environment(temp; axis = Temperature)
end

end
