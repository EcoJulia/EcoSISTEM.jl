# SPDX-License-Identifier: LGPL-3.0-or-later

module TestSimplify

using EcoSISTEM
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using Unitful, Unitful.DefaultSymbols
using AxisArrays
using RasterDataSources
using Test

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

end
