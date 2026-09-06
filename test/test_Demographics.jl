# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDemographics

using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using Unitful
import EcoSISTEM.equalpop

@testset "Params" begin
    birth = 0.6 / month_mean_duration
    death = 0.6 / month_mean_duration
    longevity = 1.0
    survival = 0.0
    numSpecies = 10

    param = EqualPop(birth, death, longevity, survival)
    @test_nowarn EqualPop(birth, death, longevity, survival)
    equalparams = equalpop(param, numSpecies)
    @test length(equalparams.birth) == numSpecies
    @test all(equalparams.birth .== birth)
    @test all(equalparams.death .== death)
    @test_nowarn param = PopGrowth{typeof(unit(0.0 / month_mean_duration))}(fill(birth,
                                                                                 5),
                                                                            fill(death,
                                                                                 5),
                                                                            longevity,
                                                                            survival)
    @test_nowarn param = NoGrowth{typeof(unit(0.0 / month_mean_duration))}(fill(birth,
                                                                                5),
                                                                           fill(death,
                                                                                5),
                                                                           longevity,
                                                                           survival)

    param = PopGrowth{typeof(unit(0.0 / month_mean_duration))}(fill(birth,
                                                                    numSpecies),
                                                               fill(death,
                                                                    numSpecies),
                                                               longevity,
                                                               survival)
    equalparams = equalpop(param, numSpecies)
    @test length(equalparams.birth) == numSpecies
    @test all(equalparams.birth .== birth)
    @test all(equalparams.death .== death)

    param = NoGrowth{typeof(unit(0.0 / month_mean_duration))}(fill(birth,
                                                                   numSpecies),
                                                              fill(death,
                                                                   numSpecies),
                                                              longevity,
                                                              survival)
    equalparams = equalpop(param, numSpecies)
    @test length(equalparams.birth) == numSpecies
    @test all(equalparams.birth .== birth)
    @test all(equalparams.death .== death)
end

# A struct dump would print both rate vectors in full. Bounded size and the identifying facts are
# what is pinned.
@testset "per-species params print bounded, whatever the species count" begin
    u = typeof(unit(0.6 / year))
    same = PopGrowth{u}(fill(0.6 / year, 1000), fill(0.6 / year, 1000), 1.0,
                        0.2)
    spread = NoGrowth{u}(collect(range(0.1 / year, 0.9 / year, length = 1000)),
                         fill(0.5 / year, 1000), 1.0, 0.0)
    for p in (same, spread)
        @test length(repr(p)) < 160
        @test length(repr("text/plain", p)) < 200
        @test occursin(string(nameof(typeof(p))), repr(p))
        @test occursin("1000 species", repr(p))
        @test count(==('\n'), repr("text/plain", p)) == 5
    end
    @test occursin("birth 0.6 yr⁻¹, death 0.6 yr⁻¹", repr(same))
    @test occursin("birth 0.1 yr⁻¹ to 0.9 yr⁻¹", repr(spread))
    @test occursin("survival 0.0", repr(spread))
    # `EqualPop` holds scalars and keeps the default, which is short.
    @test length(repr(EqualPop(0.6 / year, 0.6 / year, 1.0, 0.2))) < 120
end

end
