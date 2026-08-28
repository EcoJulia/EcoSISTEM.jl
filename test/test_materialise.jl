# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `src/materialise.jl` — what inspection shows against what the builder builds.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_materialise.jl"])'

module TestMaterialise

using Test
using EcoSISTEM
using EcoSISTEM: hasdata, landcoverclass
using EcoSISTEM: materialise
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using Rasters
using RasterDataSources
using ArchGDAL
using Distributions: Normal
using Extents: Extent
include("rasterfixtures.jl")
include("buildfixtures.jl")

# **`[ONE-PATH]`: `GridHabitat` puts in the habitat exactly what `materialise` shows.**
# This is a *structural* guard, not a numerical coincidence. A builder running its own near-copy of
# the inspection chain (`_materialiseon` → `_assemble` → `_resolve_regime`/`_resolve_supply`) drifts
# from it: three times over, in un-united dims, a hardcoded `Intervals(Start())`, and a `NicheSpec`
# that `materialise` built and the builder could not. Each
# was found by accident. The builder now calls `materialise`, so the assertions below can only fail
# if the two are pulled apart again.
#
# Both roles, both kinds of spec (data-backed and synthetic), on both kinds of positioned area —
# because the two paths differed *per kind*, so a single case would prove almost nothing.
# **`NicheSpec` is deliberately absent**: it is stochastic and unseeded (A19), so two
# materialisations of one spec disagree with each other, never mind with the builder.
@testset "what `materialise` shows is what `GridHabitat` builds" begin
    data = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)),
                axis = Temperature)
    watersrc = _reg(_bngraster(WorldClim{BioClim}, fill(4.0mm / day, 9, 9)),
                    axis = Precipitation)
    sun = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
    warm = UniformSpec(291.0K, axis = Temperature)

    for (label, area) in ("projected" => _area(regime = data),
                          "synthetic" => _area(extent = (40.0km, 40.0km),
                                cellsize = 10.0km))
        # A synthetic area cannot take a data-backed layer at all, so it tests the synthetic pair.
        regime = label == "projected" ? data : warm
        supply = label == "projected" ? watersrc : sun
        env = GridHabitat(regime = regime, supply = supply, area = area)
        seenreg = materialise(regime, area, role = EcoSISTEM.Condition)
        seensup = materialise(supply, area, role = EcoSISTEM.Resource)

        @test env.regime.matrix == seenreg.matrix
        @test env.supply.matrix == seensup.matrix
        # The **dims**, not just the values: A17 was two identical value arrays whose coordinates
        # were described differently, which is exactly what a value-only check misses.
        @test dims(env.regime.matrix, (Y, X)) == dims(seenreg.matrix, (Y, X))
        @test dims(env.supply.matrix, (Y, X)) == dims(seensup.matrix, (Y, X))
        @test env.regime.size == seenreg.size
        @test typeof(env.supply) === typeof(seensup)
    end

    # …including the *mixed* multi-layer regime, where one member is generated at the grid's shape
    # and the other sampled onto it — the arity and the names both survive the round trip.
    area = _area(regime = data)
    env = GridHabitat(regime = (temperature = data, extra = warm),
                      supply = sun, area = area)
    seen = materialise((temperature = data, extra = warm), area,
                       role = EcoSISTEM.Condition)
    @test keys(EcoSISTEM.NamedTuple(env.regime)) ==
          keys(EcoSISTEM.NamedTuple(seen))
    @test env.regime.temperature.matrix == seen.temperature.matrix
    @test env.regime.extra.matrix == seen.extra.matrix
end

end
