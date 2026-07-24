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
        @test layerunit(WorldClim{BioClim}, 1) == °C        # annual mean temperature
        @test layerunit(WorldClim{BioClim}, 2) == K            # mean diurnal range
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
        @test layerunit(CHELSA{BioClim}, 1) == °C
    end

    @testset "code aliases + reconciled codes resolve" begin
        # BioClim layers take both an int and a `bioN` symbol (`;`-separated `Code` aliases)
        @test layerunit(WorldClim{BioClim}, :bio1) ==
              layerunit(WorldClim{BioClim}, 1) == °C
        @test layerunit(WorldClim{BioClim}, :bio12) == mm
        # LandCover takes an int or a descriptive symbol
        @test layerunit(EarthEnv{LandCover}, :open_water) ==
              layerunit(EarthEnv{LandCover}, 12)
        # BioClimPlus documents bio1-19
        @test layerunit(CHELSA{BioClimPlus}, :bio1) == °C
        # CHELSA{Climate} monthly layers live in Climate.csv; `pr` is now `mm` (not kg m⁻² month⁻¹)
        @test layerunit(CHELSA{Climate}, :cmi) == kg * m^-2 * month^-1
        @test layerunit(CHELSA{Climate}, :pr) == mm
        @test layerunit(CHELSA{Climate}, :tas) == °C
    end

    @testset "unknown layer / source errors" begin
        @test_throws ErrorException layerunit(WorldClim{BioClim}, 99)
        @test_throws ErrorException layerunit(RasterDataSources.AWAP, 1)
    end
end

@testset "layeraxis / _resolve_axis" begin
    # populated Axis columns resolve to the right axis (incl. the hierarchy leaves)
    @test layeraxis(WorldClim{BioClim}, 1) === Temperature
    @test layeraxis(WorldClim{BioClim}, 2) === TemperatureRange
    @test layeraxis(WorldClim{BioClim}, 12) === Precipitation
    @test layeraxis(WorldClim{Climate}, :srad) === SolarRadiation
    @test layeraxis(EarthEnv{LandCover}, 1) === LandCoverTypology
    @test layeraxis(WorldClim{Elevation}, :elev) === Altitude
    @test layeraxis(EarthEnv{HabitatHeterogeneity}, :shannon) === Heterogeneity
    # the max/min/mean collapse: quarter temps all fold into Temperature
    @test layeraxis(CHELSA{BioClimPlus}, :bio8) === Temperature
    # the aggregate families: level vs range; and the new axes
    @test layeraxis(CHELSA{BioClimPlus}, :cmi_max) === ClimateMoisture
    @test layeraxis(CHELSA{BioClimPlus}, :cmi_range) === ClimateMoistureRange
    @test layeraxis(CHELSA{Climate}, :vpd) === VaporPressureDeficit
    @test layeraxis(CHELSA{Climate}, :hurs) === RelativeHumidity
    @test layeraxis(CHELSA{Climate}, :pet) === Evapotranspiration
    @test layeraxis(CHELSA{Climate}, :clt) === CloudCover
    @test layeraxis(CHELSA{BioClimPlus}, :gdd0) === CumulativeHeat
    @test layeraxis(CHELSA{BioClimPlus}, :gst) === Temperature
    @test layeraxis(CHELSA{BioClimPlus}, :swe) === SnowWaterEquivalent
    @test layeraxis(CHELSA{BioClimPlus}, :gsp) === GrowingSeasonPrecipitation
    @test layeraxis(CHELSA{BioClimPlus}, :fgd) === DayOfYear
    @test layeraxis(CHELSA{BioClimPlus}, :gsl) === DayRange
    @test layeraxis(CHELSA{BioClimPlus}, :npp) === CarbonFlux
    @test layeraxis(CHELSA{BioClimPlus}, :fcf) === FrostChangeFrequency
    @test layeraxis(CHELSA{BioClimPlus}, :kg0) === ClimateTypology
    # a name resolves to its concrete NicheAxis by autodiscovery (through abstract groups)
    @test ClimatePref._resolve_axis("Temperature") === Temperature
    @test ClimatePref._resolve_axis("Heterogeneity") === Heterogeneity
    @test ClimatePref._resolve_axis("Altitude") === Altitude
    # an abstract grouping supertype is NOT a leaf, so it does not resolve as a layer's axis
    @test_throws ErrorException ClimatePref._resolve_axis("TemperatureAxis")
    @test_throws ErrorException ClimatePref._resolve_axis("WaterAxis")
    @test_throws ErrorException ClimatePref._resolve_axis("DayAxis")
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

@testset "guard: tables reconcile to RasterDataSources layers(T)" begin
    RDS = RasterDataSources
    # each shipped table ↔ the RDS source type(s) that map to it (via `_datasettype` → CSV name)
    tablesources = Dict("BioClim" => Any[WorldClim{BioClim}, CHELSA{BioClim}],
                        "BioClimPlus" => Any[CHELSA{BioClimPlus}],
                        "Climate" => Any[WorldClim{Climate}, CHELSA{Climate}],
                        "Elevation" => Any[WorldClim{Elevation}],
                        "HabitatHeterogeneity" =>
                            Any[EarthEnv{HabitatHeterogeneity}],
                        "LandCover" => Any[EarthEnv{LandCover}])
    # `ncdf` is listed by `layers(CHELSA{Climate})` but is a spurious RDS entry (no such CHELSA
    # variable / file), so it is deliberately absent from the table — carve it out here.
    known_spurious = Set(["ncdf"])
    datadir = pkgdir(EcoSISTEM, "data", "RasterDataSources")
    # every code any source accepts: `layers(T)` plus its `layerkeys(T)` aliases
    accepted(T) = union(Set(string(c) for c in RDS.layers(T)),
                        Set(string(c)
                            for c in (try
                                          RDS.layerkeys(T)
                                      catch
                                          RDS.layers(T)
                                      end)))
    for (base, Ts) in tablesources
        acc = setdiff(union((accepted(T) for T in Ts)...), known_spurious)
        csvcodes = Set(keys(ClimatePref._layertable(joinpath(datadir,
                                                             "$base.csv"))))
        # exactly reconciled: every fetchable code is documented, and no phantom rows exist
        @test acc == csvcodes
    end
end

@testset "catalogue query helpers (layerinfo / layersbyaxis / layeraxes)" begin
    CP = ClimatePref
    # layerinfo(T, code): one detailed record
    r = CP.layerinfo(WorldClim{BioClim}, 1)
    @test r isa CP.LayerRecord
    @test r.name == "Annual Mean Temperature"
    @test r.unit == °C
    @test r.axis === Temperature
    @test "bio1" in r.aliases && "1" in r.aliases
    @test r.sources == ["WorldClim", "CHELSA"]     # full list of supporting sources
    @test_throws ErrorException CP.layerinfo(WorldClim{BioClim}, 99)

    # layerinfo(code): every dataset that has the code
    @test Set(x.dataset for x in CP.layerinfo(:bio1)) ==
          Set([:BioClim, :BioClimPlus])
    @test only(CP.layerinfo(:tas)).sources == ["CHELSA"]   # a CHELSA-only Climate code

    # temporal metadata surfaced from the Layers / Temporal Resolution columns
    @test CP.layerinfo(WorldClim{Climate}, :tmin).temporal == "Monthly"
    @test CP.layerinfo(WorldClim{Climate}, :tmin).numlayers == 12
    @test CP.layerinfo(WorldClim{BioClim}, 1).numlayers == 1       # static single layer
    @test CP.layerinfo(WorldClim{BioClim}, 1).temporal == ""

    # layersbyaxis: concrete leaf, and abstract group spanning several axes
    @test Set(x.dataset for x in CP.layersbyaxis(Temperature)) ⊇
          Set([:BioClim, :Climate])
    temp = CP.layersbyaxis(TemperatureAxis)
    @test all(x.axis <: TemperatureAxis for x in temp)
    @test length(temp) > length(CP.layersbyaxis(Temperature))
    # the water umbrella spans precipitation + humidity + fluxes/stocks
    water = CP.layersbyaxis(WaterAxis)
    @test all(x.axis <: WaterAxis for x in water)
    @test Set(x.axis for x in water) ⊇
          Set([Precipitation, RelativeHumidity, Evapotranspiration,
                  ClimateMoisture])

    # layeraxes: nested tree of NicheAxis subtypes down to the concrete leaves, each leaf
    # carrying the names of the shipped layers that use it
    tree = CP.layeraxes()
    @test tree.axis === NicheAxis
    @test isempty(tree.names)                       # NicheAxis itself is never a layer's axis
    water = only(c for c in tree.children if c.axis === WaterAxis)
    precip = only(c for c in water.children if c.axis === PrecipitationAxis)
    precipleaf = only(c for c in precip.children if c.axis === Precipitation)
    @test isempty(precip.names)                     # abstract group: no layers of its own
    @test "Annual Mean Temperature" ∉ precipleaf.names # sanity: wrong branch has no bio1
    @test "Annual Precipitation" in precipleaf.names
    @test isempty(precipleaf.children)               # a concrete leaf has no subtypes

    # scoping to a subtree with the optional argument
    tempaxis = CP.layeraxes(TemperatureAxis)
    @test tempaxis.axis === TemperatureAxis
    @test any(c -> c.axis === Temperature, tempaxis.children)
    @test "Annual Mean Temperature" in
          only(c for c in tempaxis.children if c.axis === Temperature).names

    # public but NOT exported
    @test :layerinfo in names(ClimatePref)
    @test !Base.isexported(ClimatePref, :layerinfo)
end

end
