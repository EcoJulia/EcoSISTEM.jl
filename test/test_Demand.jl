# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDemand

using EcoSISTEM
# `[C7-VIS]` C: these are `public` rather than exported — a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: SeriesLayerChange
using Diversity
using InteractiveUtils
using Test
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using RasterDataSources
using DimensionalData: DimensionalData, DimArray, X, Y, Ti
using Rasters   # `rasterfixtures.jl` builds genuine `Projected` lookups

@testset "Demand and supply types" begin
    numspecies = 10
    abun = fill(10, numspecies)
    # **The free/dimensionless demands are gone**, and for one reason: with the free *supply*
    # removed, nothing was left for them to pair with, so neither could reach a simulation.
    # `SimpleDemand` was deleted outright; `SizeDemand` is **commented out** in `Demand.jl`, parked
    # against a metabolic reading that is not implemented (see the note there). Neither is defined.
    @test !isdefined(EcoSISTEM, :SimpleDemand)
    @test !isdefined(EcoSISTEM, :SizeDemand)
    # The generic accessors are checked on a real resource demand instead — they are
    # demand-agnostic, so they need no demand type of their own to exercise them.
    resource_vec = Demand{SolarRadiation}(fill(2.0kJ / day, numspecies))
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test counttypes(resource_vec) == length(resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1

    numSpecies = 10
    # Test Demand{SolarRadiation}
    resource_vec = Demand{SolarRadiation}(fill(0.2 * kJ / day, numSpecies))
    @test_nowarn resource_vec = Demand{SolarRadiation}(fill(0.2 * kJ / day,
                                                            numSpecies))
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test counttypes(resource_vec) == length(resource_vec.resource)

    # Test Demand{Precipitation}
    resource_vec = Demand{Precipitation}(fill(0.2 * Unitful.L / day,
                                              numSpecies))
    @test_nowarn resource_vec = Demand{Precipitation}(fill(0.2 * Unitful.L /
                                                           day,
                                                           numSpecies))
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test counttypes(resource_vec) == length(resource_vec.resource)

    # Test Demand{CarbonFlux} — the carbon family's demand side, matched against an `npp`-derived
    # supply. `mean(resource)` is the default exchange rate, exactly as Solar/Water.
    resource_vec = Demand{CarbonFlux}(fill(0.2 * g / day, numSpecies))
    @test_nowarn resource_vec = Demand{CarbonFlux}(fill(0.2 * g / day,
                                                        numSpecies))
    @test EcoSISTEM._getdemand(abun, resource_vec) ==
          sum(abun .* resource_vec.resource)
    @test EcoSISTEM.numdemands(typeof(resource_vec)) == 1
    @test eltype(resource_vec) == typeof(resource_vec.resource[1])
    @test counttypes(resource_vec) == length(resource_vec.resource)
    # a demand given in another mass-per-time unit converts to the canonical `g/day`
    @test Demand{CarbonFlux}(fill(1.0e-3 * kg / day, numSpecies)).resource[1] ==
          1.0g / day

    resource_vec1 = Demand{SolarRadiation}(fill(0.2 * kJ / day, numSpecies))
    resource_vec2 = Demand{Precipitation}(fill(0.2 * Unitful.L / day,
                                               numSpecies))
    demand = SpeciesRequirementCollection((resource_vec1, resource_vec2))
    @test EcoSISTEM.numdemands(typeof(demand)) == 2
    @test map(eltype, values(demand)) ==
          (typeof(resource_vec1.resource[1]), typeof(resource_vec2.resource[1]))
    @test counttypes(demand) == length(resource_vec1.resource)

    # The free/dimensionless supply is gone (its meaning could only come from its unit), so the
    # per-cell accessors are exercised on a solar supply instead — they are axis-agnostic.
    supply = fill(100.0 * EcoSISTEM.canonicalunit(EcoSISTEM.Resource,
                                          SolarRadiation), 2, 2)
    @test_nowarn Supply{SolarRadiation}(supply)
    supply = Supply{SolarRadiation}(supply)
    @test countsubcommunities(supply) == 4
    @test EcoSISTEM._getsupply(supply) == supply.matrix
    @test eltype(supply) == typeof(supply.matrix[1])

    # Test Supply{SolarRadiation}
    sol = fill(200.0 * kJ / day, 100, 100)
    @test_nowarn Supply{SolarRadiation}(sol)
    supply1 = Supply{SolarRadiation}(sol)
    @test countsubcommunities(supply1) == 100 * 100
    @test EcoSISTEM._getsupply(supply1) == supply1.matrix
    @test eltype(supply1) == typeof(supply1.matrix[1])

    # Test Supply{Precipitation}
    water = fill(2000.0 * Unitful.L / day, 100, 100)
    @test_nowarn Supply{Precipitation}(water)
    supply2 = Supply{Precipitation}(water)
    @test countsubcommunities(supply2) == 100 * 100
    @test EcoSISTEM._getsupply(supply2) == supply2.matrix
    @test eltype(supply2) == typeof(supply2.matrix[1])

    # Test Supply{CarbonFlux}
    carbon = fill(5.0 * g / day, 100, 100)
    @test_nowarn Supply{CarbonFlux}(carbon)
    supply3 = Supply{CarbonFlux}(carbon)
    @test countsubcommunities(supply3) == 100 * 100
    @test EcoSISTEM._getsupply(supply3) == supply3.matrix
    @test eltype(supply3) == typeof(supply3.matrix[1])
    # The axis is recorded on the type — `npp` is a resource in its own right, not water or light.
    @test supply3 isa ContinuousLayer{EcoSISTEM.Resource, CarbonFlux}

    # **`Supply{A}` leaves its value type free, so the unit guarantee lives in the constructor.**
    # The old per-resource aliases pinned the element type in the signature and got this from
    # dispatch; these four are what replaced that, and each is a different mistake.
    #
    # 1. A dimensionally-correct value at another scale is **converted**, not refused — which is
    #    what lets CHELSA's `MJ·m⁻²·d⁻¹` meet WorldClim's `kJ·m⁻²·d⁻¹` on one axis.
    megajoules = Supply{SolarRadiation}(fill(0.2u"MJ/d", 3, 3))
    @test all(==(200.0kJ / day), megajoules.matrix)
    # …and it is the *same concrete type* as one built canonically, not a parallel MJ-flavoured one.
    @test typeof(megajoules) ===
          typeof(Supply{SolarRadiation}(fill(200.0kJ / day,
                                             3, 3)))

    # 2. A wrong dimension is refused **here**, naming the axis and the unit it declares. Left to
    #    build, it would surface as a `DimensionError` inside the threaded hot loop instead.
    @test_throws "measured in `kJ d⁻¹`" Supply{SolarRadiation}(fill(2000.0Unitful.L /
                                                                    day, 3, 3))
    # 3. A bare number cannot be checked against anything, so it is refused too.
    @test_throws "carry no unit at all" Supply{SolarRadiation}(fill(3.0, 3, 3))
    # 4. The wind-speed case, at the constructor rather than through `GridHabitat`: an axis
    #    that declares no resource is refused by name, not silently given a supply type.
    @test_throws "not a consumable resource" Supply{WindSpeed}(fill(3.0m / s, 3,
                                                                    3))

    # Test a time-varying solar supply: the same type, carrying a `SeriesLayerChange` over the stack
    sol = fill(200.0 * kJ / day, 100, 100, 10)
    supply1 = EcoSISTEM._setseries!(Supply{SolarRadiation}(sol[:, :, 1]), sol)
    @test supply1 isa Supply{SolarRadiation}
    @test supply1.change isa SeriesLayerChange
    @test countsubcommunities(supply1) == 100 * 100
    @test EcoSISTEM._getsupply(supply1) == supply1.matrix
    @test eltype(supply1) == typeof(supply1.matrix[1])

    # Test a time-varying water supply
    water = fill(2000.0 * Unitful.L / day, 100, 100, 10)
    supply2 = EcoSISTEM._setseries!(Supply{Precipitation}(water[:, :, 1]),
                                    water)
    @test supply2 isa Supply{Precipitation}
    @test countsubcommunities(supply2) == 100 * 100
    @test EcoSISTEM._getsupply(supply2) == supply2.matrix
    @test eltype(supply2) == typeof(supply2.matrix[1])

    # Test SupplyCollection
    supply = LayerCollection((supply1, supply2))
    @test_nowarn LayerCollection((supply1, supply2))
    @test countsubcommunities(supply) == 100 * 100
    # Named by axis, the two supplies being on distinguishable ones.
    @test keys(supply) == (:SolarRadiation, :Precipitation)
    @test EcoSISTEM._getsupply(supply, :SolarRadiation) == supply1.matrix
    @test map(eltype, values(supply)) ==
          (typeof(supply1.matrix[1]), typeof(supply2.matrix[1]))
    [sum(supply1.matrix), sum(supply2.matrix)]
end

# **Placed before `Worldclim/Bioclim supplies`, though that testset no longer depends on the order.**
# **Kept as written because the failure mode it names is permanent**: a testset that errors on its
# first line silently skips every assertion after it — which is how the *only* tests of
# `Supply{A}(::ClimateRaster)` can go unrun — and a file that errors mid-way still emits a
# `Test Summary`, so neither a failure-marker grep nor a summary-presence check can see the loss.
# Only reading the pass count against what the file contains can.
include("rasterfixtures.jl")

@testset "supplies scale by each cell's true area on a geographic raster" begin
    # `Supply{A}(::ClimateRaster)` derives its cell area from the raster's *own* grid — there is no
    # `StudyArea` in this constructor — so it is the one supply path where the latitude correction
    # has to come from the data itself. The pre-existing coverage could not see this: its
    # latitudes were in **metres** (a projected grid), where the correction is exactly 1.0.
    rain = fill(2.0mm / day, 5, 4)                       # non-square, so a Y/X swap shows up

    # Degrees: a real geographic grid, where cells shrink towards the pole.
    geo = _testraster(WorldClim{BioClim}, rain,
                      lat = (50.0:2.0:58.0) .* °, long = (0.0:1.0:3.0) .* °)
    gsup = Supply{Precipitation}(geo)
    gvals = Array(gsup.matrix)
    @test size(gvals) == (5, 4)
    @test all(≈(gvals[1, 1]), gvals[1, :])          # constant across X…
    @test issorted(gvals[:, 1], rev = true)         # …and falling northwards down Y
    @test !isapprox(gvals[1, 1], gvals[end, 1], rtol = 0.05)

    # Metres: a projected grid is exact already, so every cell gets the same supply.
    proj = _testraster(WorldClim{BioClim}, rain,
                       lat = (0.0:2000.0:8000.0) .* m,
                       long = (0.0:1000.0:3000.0) .* m,
                       crs = Rasters.EPSG(27700))
    pvals = Array(Supply{Precipitation}(proj).matrix)
    @test all(≈(pvals[1, 1]), pvals)
end

# guard: every axis that declares a demand is measured in something the model can divide.
#
# **Restated at step 8, because the space resource is the first demand that is NOT a rate.**
# Asserting `dimension(elt) / _basedimension(elt) == 𝐓⁻¹` of *every* demand holds only while every
# resource is a flow — energy, water, carbon. A space demand is an **area** (m² per individual) with
# no time in it at all.
# **The rule is narrowed, not dropped, and the exception is named rather than tolerated.** What the
# model actually requires is that supply ÷ demand be a dimensionless, timestep-independent count —
# which holds for a **flow ÷ flow** (`kJ/day ÷ kJ/day`) and equally for a **stock ÷ stock**
# (`m² ÷ m²`). What it would *not* survive is mixing the two on one axis, and that is what this now
# checks: each axis is one or the other, consistently, on both sides.
# Supplies are never depleted — they are recomputed in full each step — so a standing stock needs
# no change to the loop. That is why space fits without one.
@testset "guard: every declared resource is a flow or a stock, consistently" begin
    allaxes(T = EcoSISTEM.NicheAxis, acc = Type[]) = begin
        push!(acc, T)
        for S in subtypes(T)
            allaxes(S, acc)
        end
        acc
    end
    resourceaxes = filter(A -> !isnothing(EcoSISTEM.demandtype(A)), allaxes())
    @test !isempty(resourceaxes)

    # The axes whose resource is a **stock** rather than a flow, named explicitly so that adding
    # another is a deliberate act rather than a silently weakened guard.
    stocks = (SurfaceArea,)

    for A in resourceaxes
        u = EcoSISTEM.canonicalunit(EcoSISTEM.Resource, A)
        d = Demand{A}(fill(1.0 * u, 3))
        @test EcoSISTEM.axisof(d) === A
        elt = eltype(d)
        # Both sides of an axis share one unit, so the ratio is dimensionless by construction —
        # that is what `Supply{A}`/`Demand{A}` buy. What is worth asserting is which *kind* of
        # quantity it is, and that the axis is honest about it.
        if any(==(A), stocks)
            @test dimension(elt) == dimension(m^2)          # a stock: an area, no time
            @test !hasmethod(EcoSISTEM._basedimension, Tuple{Type{elt}})
        else
            # a genuine rate: `substance × 𝐓⁻¹` exactly, not merely "some negative time exponent"
            # (energy already embeds 𝐓⁻² in its own definition).
            @test dimension(elt) / EcoSISTEM._basedimension(elt) == Unitful.𝐓^-1
        end
    end
end

@testset "Worldclim/Bioclim supplies" begin
    # The hand constructors take a raw ClimateRaster whose element type is the areal rate
    # (e.g. `mm/day`) at its own native resolution — `cancel` (against the axis's canonical
    # resource unit) converts it to the absolute `L/day`/`kJ/day` supply value against the cell area,
    # so the resulting `supply.matrix` values are scaled by that area, not equal to the raw
    # input; the assertions below compare the built supply against itself, not a hardcoded
    # number, so they hold regardless of the exact scaling.
    water = DimArray(fill(1.0mm / day, 10, 10, 12),
                     (Y(collect(1:10) .* m), X(collect(1:10) .* m),
                      Ti(collect(1:12) .* month_mean_duration)))
    worldclim = ClimateRaster(WorldClim{Climate}, water)
    # The monthly `ClimateRaster` constructors are deprecated (the modern path reads a stack
    # through `GridHabitat`), so they warn; the shim itself is covered in
    # `test_deprecations.jl`. What is being checked here is the area scaling they still do.
    supply = @test_deprecated WaterTimeBudget(worldclim, 1)
    @test countsubcommunities(supply) == 100
    @test EcoSISTEM._getsupply(supply) == supply.matrix
    @test eltype(supply) == typeof(supply.matrix[1])

    solar = DimArray(fill(1.0kJ / m^2 / day, 10, 10, 12),
                     (Y(collect(1:10) .* m), X(collect(1:10) .* m),
                      Ti(collect(1:12) .* month_mean_duration)))
    worldclim = ClimateRaster(WorldClim{Climate}, solar)
    supply = @test_deprecated SolarTimeBudget(worldclim, 1)
    @test countsubcommunities(supply) == 100
    @test EcoSISTEM._getsupply(supply) == supply.matrix
    @test eltype(supply) == eltype(supply.matrix)

    water = DimArray(fill(1.0mm / day, 10, 10),
                     (Y(collect(1:10) .* m), X(collect(1:10) .* m)))
    worldclim = ClimateRaster(WorldClim{BioClim}, water)
    supply = Supply{Precipitation}(worldclim)
    @test_nowarn Supply{Precipitation}(worldclim)
    @test countsubcommunities(supply) == 100
    @test EcoSISTEM._getsupply(supply) == supply.matrix
    @test eltype(supply) == eltype(supply.matrix)
end

end
