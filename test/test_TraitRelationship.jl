# SPDX-License-Identifier: LGPL-3.0-or-later

module TestTraitRelationship

using EcoSISTEM
using Test
using Distributions
using Unitful
using Unitful.DefaultSymbols

@testset "Trait relationships" begin
    @test_nowarn Gauss{Unitful.Temperature}()
    @test_nowarn Match{Int64}()
    @test_nowarn Trapeze{Int64}()
    @test_nowarn Unif{Int64}()
    @test_nowarn NoRelContinuous{Int64}()
    @test_nowarn NoRelDiscrete{Int64}()

    @test Gauss{Unitful.Temperature}()(1.0K, 1.0K, 0.01K) > 0.0
    @test EcoSISTEM.iscontinuous(Gauss{Unitful.Temperature}()) == true
    @test eltype(Gauss{Unitful.Temperature}()) == Unitful.Temperature

    @test Match{Int64}()(1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(Match{Int64}()) == false
    @test eltype(Match{Int64}()) == Int64

    @test Trapeze{Int64}()(Trapezoid(1, 2, 3, 4), 1) == 0.0
    @test EcoSISTEM.iscontinuous(Trapeze{Int64}()) == true
    @test eltype(Trapeze{Int64}()) == Int64

    @test Unif{typeof(1.0mm)}()(Uniform(1, 2), 1.0mm) == 1.0
    @test EcoSISTEM.iscontinuous(Unif{typeof(1.0mm)}()) == true
    @test eltype(Unif{typeof(1.0mm)}()) == typeof(1.0mm)

    @test NoRelContinuous{Int64}()(1, 1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(NoRelContinuous{Int64}()) == true
    @test eltype(NoRelContinuous{Int64}()) == Int64

    @test NoRelDiscrete{Int64}()(1, 1) == 1.0
    @test EcoSISTEM.iscontinuous(NoRelDiscrete{Int64}()) == false
    @test eltype(NoRelDiscrete{Int64}()) == Int64

    tr2 = multiplicativeTR2(NoRelContinuous{Int64}(), NoRelDiscrete{Int64}())
    @test EcoSISTEM.iscontinuous(tr2) == [true, false]
    @test eltype(tr2) == [Int64, Int64]
    tr3 = multiplicativeTR3(NoRelContinuous{Int64}(), NoRelDiscrete{Int64}(),
                            Gauss{Unitful.Temperature}())
    @test EcoSISTEM.iscontinuous(tr3) == [true, false, true]
    @test eltype(tr3) == [Int64, Int64, Unitful.Temperature]

    @test EcoSISTEM.combineTR(tr2) == *
    @test EcoSISTEM.combineTR(tr3) == *

    tr2 = additiveTR2(NoRelContinuous{Int64}(), NoRelDiscrete{Int64}())
    @test EcoSISTEM.iscontinuous(tr2) == [true, false]
    @test eltype(tr2) == [Int64, Int64]
    tr3 = additiveTR3(NoRelContinuous{Int64}(), NoRelDiscrete{Int64}(),
                      Gauss{Unitful.Temperature}())
    @test EcoSISTEM.iscontinuous(tr3) == [true, false, true]
    @test eltype(tr3) == [Int64, Int64, Unitful.Temperature]

    @test EcoSISTEM.combineTR(tr2) == +
    @test EcoSISTEM.combineTR(tr3) == +
end

end
