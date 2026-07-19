# SPDX-License-Identifier: LGPL-3.0-or-later

module TestLayerUnits

using EcoSISTEM
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using RasterDataSources
using Unitful
using Unitful.DefaultSymbols
using Test

@testset "layerunit" begin
    @testset "temperature and precipitation (BioClim)" begin
        @test layerunit(WorldClim{BioClim}, 1) == K            # annual mean temperature
        @test layerunit(WorldClim{BioClim}, 12) == mm          # annual precipitation
        # dimensionless indices come back as NoUnits
        @test layerunit(WorldClim{BioClim}, 3) == NoUnits      # isothermality
    end

    @testset "compound units and Symbol/String keys (Climate)" begin
        @test layerunit(WorldClim{Climate}, :srad) == kJ * m^-2 * day^-1
        @test layerunit(WorldClim{Climate}, "wind") == m * s^-1
        @test dimension(layerunit(WorldClim{Climate}, :vapr)) ==
              dimension(u"kPa")
    end

    @testset "land cover is dimensionless" begin
        @test layerunit(EarthEnv{LandCover}, 1) == NoUnits
        @test layerunit(EarthEnv{LandCover}, 12) == NoUnits
    end

    @testset "other tables resolve (BOM, month, elevation)" begin
        @test layerunit(WorldClim{Elevation}, :elev) == m
        @test layerunit(WorldClim{BioClimPlus}, :cmi_max) ==
              kg * m^-2 * month^-1
    end

    @testset "_datasettype extracts the wrapped dataset" begin
        @test ClimatePref._datasettype(WorldClim{BioClim}) === BioClim
        @test ClimatePref._datasettype(EarthEnv{LandCover}) === LandCover
        @test ClimatePref._datasettype(CHELSA{BioClim}) === BioClim
    end

    @testset "dataset type resolves through any wrapper" begin
        # CHELSA{BioClim} shares the BioClim layer table
        @test layerunit(CHELSA{BioClim}, 1) == K
    end

    @testset "unknown layer / source errors" begin
        @test_throws ErrorException layerunit(WorldClim{BioClim}, 99)
        @test_throws ErrorException layerunit(RasterDataSources.AWAP, 1)
    end
end

@testset "layeraxis / _resolve_axis" begin
    # populated Axis columns resolve to the right axis (incl. the hierarchy leaves)
    @test layeraxis(WorldClim{BioClim}, 1) === MeanTemperature
    @test layeraxis(WorldClim{BioClim}, 12) === Precipitation
    @test layeraxis(WorldClim{Climate}, :srad) === SolarRadiation
    @test layeraxis(EarthEnv{LandCover}, 1) === LandType
    @test layeraxis(WorldClim{Elevation}, :elev) === Altitude
    @test layeraxis(EarthEnv{HabitatHeterogeneity}, :shannon) === Heterogeneity
    # BioClimPlus has no Axis column populated yet ⇒ unclassified
    @test layeraxis(WorldClim{BioClimPlus}, :cmi_max) === nothing
    # a name resolves to its concrete NicheAxis by autodiscovery (through abstract groups)
    @test ClimatePref._resolve_axis("MeanTemperature") === MeanTemperature
    @test ClimatePref._resolve_axis("Heterogeneity") === Heterogeneity
    @test ClimatePref._resolve_axis("Altitude") === Altitude
    # a name that isn't a loaded NicheAxis errors clearly
    @test_throws ErrorException ClimatePref._resolve_axis("NotAnAxis")
end

@testset "guard: shipped data/ tables are well-formed" begin
    datadir = pkgdir(EcoSISTEM, "data", "RasterDataSources")
    csvs = filter(endswith(".csv"), readdir(datadir))
    @test !isempty(csvs)
    for f in csvs
        base = first(splitext(f))
        # every table's basename names a real RasterDataSources dataset type
        @test isdefined(RasterDataSources, Symbol(base))
        rows = ClimatePref._layertable(joinpath(datadir, f))
        for (_, cell) in rows
            # every non-blank Units parses …
            isempty(cell.units) ||
                @test uparse(cell.units, unit_context = [Unitful, Units]) isa
                      Unitful.Units
            # … and every non-blank Axis resolves to a loaded NicheAxis
            isempty(cell.axis) ||
                @test ClimatePref._resolve_axis(cell.axis) <: NicheAxis
        end
    end
end

end
