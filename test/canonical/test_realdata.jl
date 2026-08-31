# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Canonical results from **real** input: WorldClim rasters read through the ordinary spec/read path.
# The companion to `test_simulated.jl`, and the pair is the point - see README.md.
#
# **This is the file that would have caught the February bug.** A monthly total divided by a fixed
# 30.4375-day month is 7.7% wrong for February, and until this file existed *nothing in the suite read
# a monthly layer at all*: the error could be introduced, or corrected, without a single test moving.
# The blessed rates below pin the real per-month divisors.
#
# Downloads on first run (cached into `EcoSISTEM.assetdir()` afterwards), so it is deliberately not
# part of the unit suite.

module CanonicalRealData

using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using RasterDataSources
using DimensionalData

include("canonical.jl")
using .Canonical

# A land cell chosen by **coordinate**, not by array index: an index is meaningless if the source's
# grid ever changes resolution, and index (400, 300) turned out to be ocean - every value `NaN`, every
# assertion vacuous. Edinburgh is unambiguously land and unambiguously wet.
const SITE = (lat = 55.95°, long = -3.19°)
cellat(a) = a[Y(Near(SITE.lat)), X(Near(SITE.long))]

@testset "canonical: real monthly climate" begin
    prec = EcoSISTEM._read(SourceSpec(WorldClim{Climate}, :prec))
    rateunit = Unitful.L / Unitful.m^2 / Unitful.d

    # The blessed rates, month by month. Twelve numbers rather than a summary precisely because the
    # bug this guards against was *per month*: a single annual figure would have hidden it entirely,
    # since dividing every month by 30.4375 conserves the annual total exactly.
    monthly = ustrip.(rateunit, collect(cellat(prec.array)))
    canonical("realdata/prec_monthly_rate", monthly)

    # The February/January ratio, which is the bug's own signature: it depends *only* on the two
    # divisors, so it is identical for every cell on Earth and cannot drift with the data. If someone
    # reinstates a fixed month, this is the number that moves - from 31/28.25 back to 1.
    raw = read(WorldClim{Climate}, :prec)
    rawm = collect(cellat(raw.array))
    febjan = (monthly[2] / monthly[1]) / (rawm[2] / rawm[1])
    canonical("realdata/feb_jan_divisor_ratio", febjan)
    @test febjan ≈ 31 / 28.25 rtol=1e-6

    # The unit a read yields - the observable half of option C, and a plain regression guard on the
    # catalogue's `Units`/`AccumulationPeriod` split.
    @test unit(eltype(prec.array)) == rateunit
    @test layerunit(WorldClim{Climate}, :prec) == Unitful.L * m^-2

    # Properties that hold whatever the blessed numbers are, so a re-blessing cannot quietly accept
    # nonsense: rates are positive, finite, and of a plausible magnitude for terrestrial rainfall.
    @test all(>(0), monthly)
    @test all(isfinite, monthly)
    @test maximum(monthly) < 1000        # L m^-2 d^-1; a metre of rain a day would be wrong

    # A constant-period layer, which must convert by exactly the declared period - the contrast that
    # shows the per-slice divisor applies only where the catalogue says it should.
    bio13 = EcoSISTEM._read(SourceSpec(WorldClim{BioClim}, :bio13))
    rawb = read(WorldClim{BioClim}, :bio13)
    b, rb = ustrip(rateunit, cellat(bio13.array)), cellat(rawb.array)
    canonical("realdata/bio13_rate", b)
    @test b ≈ rb / 30.4375 rtol=1e-8
end

# **The discriminator, read from real data - the two halves of "is an accumulated layer divided?"
# proved on actual rasters rather than on the catalogue alone.**
#
# **Local only.** These are CHELSA `BioClimPlus` layers, far larger than the WorldClim rasters
# above, and the GitHub runners already download close to their time budget. `heavydata()` is `false`
# on a runner; set `ECOSISTEM_HEAVY_DATA=true` to force it.
#
# **Skipping cannot lose the blessed values**: `writereference` merges rather than replaces, so a
# re-blessing on CI leaves these keys exactly as they were blessed locally.
#
# And skipping does **not** leave the underlying fact unchecked - which of these layers divides is
# a catalogue property, asserted for all five in `test_LayerCatalogue.jl` on every build. What is gated
# is only whether the *read values* match, which is the part that needs the raster.
if !heavydata()
    @info "Skipping the BioClimPlus canonical reads: CHELSA layers are large and this is a CI " *
          "runner. Set ECOSISTEM_HEAVY_DATA=true to run them. The divide/do-not-divide rule they " *
          "illustrate is checked without downloading in test_LayerCatalogue.jl."
else
    @testset "canonical: the accumulation-period discriminator, on real rasters" begin
        rateunit = Unitful.L / Unitful.m^2 / Unitful.d

        # **Divides**: potential evapotranspiration is a flow - evaporative demand happens
        # continuously, so the accumulated total over its period is a meaningful mean daily rate.
        pet = EcoSISTEM._read(SourceSpec(CHELSA{BioClimPlus}, :pet_penman_mean))
        canonical("realdata/pet_rate", ustrip(rateunit, cellat(pet.array)))
        @test unit(eltype(pet.array)) == rateunit

        # **Does not divide**: a heat sum, where the accumulated total *is* the quantity a species
        # is matched against. Dividing would give a mean daily excess nobody asked for.
        gdd0 = EcoSISTEM._read(SourceSpec(CHELSA{BioClimPlus}, :gdd0))
        canonical("realdata/gdd0_sum",
                  ustrip(Unitful.K * Unitful.d, cellat(gdd0.array)))
        @test unit(eltype(gdd0.array)) == Unitful.K * Unitful.d

        # **The pair that proves `Category` does not decide it**: `cmi` is a `balance` and divides,
        # `swb` is a `balance` and does not - a capped cumulative is not a flow.
        cmi = EcoSISTEM._read(SourceSpec(CHELSA{BioClimPlus}, :cmi_mean))
        swb = EcoSISTEM._read(SourceSpec(CHELSA{BioClimPlus}, :swb))
        canonical("realdata/cmi_rate", ustrip(rateunit, cellat(cmi.array)))
        canonical("realdata/swb_amount",
                  ustrip(Unitful.L / Unitful.m^2, cellat(swb.array)))
        @test unit(eltype(cmi.array)) == rateunit
        @test unit(eltype(swb.array)) == Unitful.L / Unitful.m^2

        # Properties that hold whatever the blessed numbers are: a heat sum and an evaporative
        # demand are non-negative, while both balances may legitimately be either sign.
        @test ustrip(gdd0.array[Y(Near(SITE.lat)), X(Near(SITE.long))]) >= 0
        @test ustrip(pet.array[Y(Near(SITE.lat)), X(Near(SITE.long))]) >= 0
    end
end

end
