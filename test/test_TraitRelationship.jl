# SPDX-License-Identifier: LGPL-3.0-or-later

module TestTraitRelationship

using EcoSISTEM
using Test
using Distributions
using Unitful
using Unitful.DefaultSymbols

@testset "Trait relationships" begin
    @test_nowarn NicheSuitability{Unitful.Temperature}()
    @test_nowarn MatchSuitability{Int64}()
    @test_nowarn NoFitContinuous{Int64}()
    @test_nowarn NoFitDiscrete{Int64}()

    # `NicheSuitability` evaluates a distribution's density at the (unit-stripped) regime value.
    @test NicheSuitability{Unitful.Temperature}()(Normal(1.0, 0.01), 1.0K) > 0.0
    @test NicheSuitability{typeof(1.0mm)}()(Uniform(1, 2), 1.0mm) == 1.0
    @test NicheSuitability{Int64}()(Trapezoid(1, 2, 3, 4), 1) == 0.0

    # A concrete `NF` forces a real `uconvert`, not a bare strip: a dimensionally-compatible but
    # differently-scaled `current` is corrected rather than silently misread (1000.0μm ≡ 1.0mm).
    @test NicheSuitability{typeof(1.0mm)}()(Uniform(1, 2), 1000.0Unitful.μm) ==
          1.0
    @test_throws Unitful.DimensionError NicheSuitability{typeof(1.0mm)}()(Uniform(1,
                                                                                  2),
                                                                          1.0K)
    @test EcoSISTEM.iscontinuous(NicheSuitability{Unitful.Temperature}()) ==
          true
    @test eltype(NicheSuitability{Unitful.Temperature}()) == Unitful.Temperature

    @test MatchSuitability{Int64}()(1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(MatchSuitability{Int64}()) == false
    @test eltype(MatchSuitability{Int64}()) == Int64

    @test NoFitContinuous{Int64}()(1, 1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(NoFitContinuous{Int64}()) == true
    @test eltype(NoFitContinuous{Int64}()) == Int64

    @test NoFitDiscrete{Int64}()(1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(NoFitDiscrete{Int64}()) == false
    @test eltype(NoFitDiscrete{Int64}()) == Int64

    tr2 = multiplicativeFit2(NoFitContinuous{Int64}(), NoFitDiscrete{Int64}())
    @test EcoSISTEM.iscontinuous(tr2) == [true, false]
    @test eltype(tr2) == [Int64, Int64]
    tr3 = multiplicativeFit3(NoFitContinuous{Int64}(),
                             NoFitDiscrete{Int64}(),
                             NicheSuitability{Unitful.Temperature}())
    @test EcoSISTEM.iscontinuous(tr3) == [true, false, true]
    @test eltype(tr3) == [Int64, Int64, Unitful.Temperature]

    @test EcoSISTEM.combinefit(tr2) == *
    @test EcoSISTEM.combinefit(tr3) == *

    tr2 = additiveFit2(NoFitContinuous{Int64}(), NoFitDiscrete{Int64}())
    @test EcoSISTEM.iscontinuous(tr2) == [true, false]
    @test eltype(tr2) == [Int64, Int64]
    tr3 = additiveFit3(NoFitContinuous{Int64}(),
                       NoFitDiscrete{Int64}(),
                       NicheSuitability{Unitful.Temperature}())
    @test EcoSISTEM.iscontinuous(tr3) == [true, false, true]
    @test eltype(tr3) == [Int64, Int64, Unitful.Temperature]

    @test EcoSISTEM.combinefit(tr2) == +
    @test EcoSISTEM.combinefit(tr3) == +
end

end
