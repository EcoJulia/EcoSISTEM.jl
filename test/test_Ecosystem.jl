# SPDX-License-Identifier: LGPL-3.0-or-later

module TestEcosystem

using EcoSISTEM
using Test
using Unitful, Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Diversity

include("TestCases.jl")

@testset "Ecosystem" begin
    @test_nowarn eco = Test1Ecosystem()
    eco = Test1Ecosystem()
    @test sum(eco.abundances.matrix, dims = 2)[:, 1] == eco.spplist.abun
    @test EcoSISTEM.tematch(eco.spplist, eco.abenv) == true
    @test EcoSISTEM.trmatch(eco.spplist, eco.relationship) == true

    sppl = SpeciesList{typeof(eco.spplist.traits),
                       typeof(eco.spplist.demand),
                       typeof(eco.spplist.movement), UniqueTypes,
                       typeof(eco.spplist.params)}(eco.spplist.names,
                                                   eco.spplist.traits,
                                                   eco.spplist.abun,
                                                   eco.spplist.demand,
                                                   UniqueTypes(length(eco.spplist.names)),
                                                   eco.spplist.movement,
                                                   eco.spplist.params,
                                                   eco.spplist.native)
    eco = Ecosystem(sppl, eco.abenv, eco.relationship)
    @test_nowarn addspecies!(eco, 10)

    @testset "get functions" begin
        # Test Simulation get functions
        @test_nowarn gettraitrel(eco)
        @test gettraitrel(eco) == eco.relationship
        @test_nowarn gethabitat(eco)
        @test gethabitat(eco) == eco.abenv.habitat
        @test_nowarn getsize(eco)
        @test getsize(eco) ==
              size(eco.abundances.matrix, 2) .* eco.abenv.habitat.size^2
        @test_nowarn getgridsize(eco)
        @test getgridsize(eco) == eco.abenv.habitat.size
        @test_nowarn getdispersaldist(eco, 1)
        @test getdispersaldist(eco, 1) == eco.spplist.movement.kernels[1].dist
        @test_nowarn getdispersaldist(eco, "1")
        @test getdispersaldist(eco, "1") == eco.spplist.movement.kernels[1].dist
        @test_nowarn getdispersalvar(eco, 1)
        @test_nowarn getdispersalvar(eco, "1")
        @test EcoSISTEM.getlookup(eco, "1") == eco.lookup[1]
        @test EcoSISTEM.getlookup(eco, 1) == eco.lookup[1]
        @test_nowarn resetrate!(eco, 0.1 / s)
        @test eco.abenv.habitat.dynamics.rate == 0.1 / s
        @test_throws MethodError resettime!(eco)
    end
    @testset "movement types" begin
        # Test other movement types
        eco = Test1Ecosystem()
        mov = AlwaysMovement(fill(LongTailKernel(10.0km, 10.0, 1e-10),
                                  length(eco.spplist.names)),
                             eco.spplist.movement.boundary)
        sppl = SpeciesList{typeof(eco.spplist.traits),
                           typeof(eco.spplist.demand), typeof(mov),
                           typeof(eco.spplist.types),
                           typeof(eco.spplist.params)}(eco.spplist.names,
                                                       eco.spplist.traits,
                                                       eco.spplist.abun,
                                                       eco.spplist.demand,
                                                       eco.spplist.types, mov,
                                                       eco.spplist.params,
                                                       eco.spplist.native)
        @test_nowarn Ecosystem(sppl, eco.abenv, eco.relationship)
    end
    @testset "diversity" begin
        # Test Diversity get functions
        @test_nowarn getmetaabundance(eco)
        @test_nowarn getpartition(eco)
        @test getpartition(eco) == eco.abenv
        @test_nowarn gettypes(eco)
        @test gettypes(eco) == eco.spplist
        @test_nowarn getordinariness!(eco)
        @test getordinariness!(eco) == eco.ordinariness
    end
    @testset "collection-trait incompatibility message" begin
        # `iscontinuous` of a collection is a `Vector{Bool}`; the constructor's incompatibility errors must
        # format that rather than crash on a `Vector{Bool}` in a `?:` (the old bug threw a `TypeError`).
        @test EcoSISTEM._kindlabel(true) == "continuous"
        @test EcoSISTEM._kindlabel(false) == "discrete"
        @test EcoSISTEM._kindlabel([true, false]) == "[continuous, discrete]"

        # Build a valid two-variable (temperature + rainfall) collection environment + species, then pass a
        # mismatched single relationship so `trmatch` fails on the collection path.
        numSpecies = 4
        grid = (5, 5)
        area = 10000.0km^2
        env1 = simplehabitatAE(10.0K, grid, 1000.0kJ / km^2, area;
                               axis = MeanTemperature)
        env2 = simplehabitatAE(10.0mm, grid, 100.0mm / km^2, area;
                               axis = Precipitation)
        habitat = HabitatCollection2(env1.habitat, env2.habitat)   # eltype [K, mm]
        supply = SupplyCollection2(env1.supply, env2.supply)
        abenv = GridAbioticEnv{typeof(habitat), typeof(supply)}(habitat,
                                                                env1.active,
                                                                supply,
                                                                env1.names)
        traits = TraitCollection2(Bin(MeanTemperature, Normal,
                                      fill(10.0K, numSpecies),
                                      fill(0.1K, numSpecies)),
                                  Bin(Precipitation, Uniform,
                                      fill(1.0mm, numSpecies),
                                      fill(5.0mm, numSpecies)))
        abun = rand(Multinomial(1000, numSpecies))
        movement = BirthOnlyMovement(GaussianKernel.(fill(1.0km, numSpecies),
                                                     10e-4))
        native = fill(true, numSpecies)
        resource = DemandCollection2(SolarDemand(fill(2.0kJ, numSpecies)),
                                     WaterDemand(fill(2.0mm, numSpecies)))
        param = EqualPop(0.6 / month, 0.6 / month, 1.0, 0.0, 1000.0)
        sppl = SpeciesList(numSpecies, traits, abun, resource, movement, param,
                           native)
        # the matching relationship is `multiplicativeTR2(...)`; a single `DistRel` mismatches
        badrel = DistRel{typeof(1.0K)}()
        err = try
            Ecosystem(sppl, abenv, badrel)
            nothing
        catch e
            e
        end
        @test err isa ErrorException                       # not a TypeError from the message itself
        @test occursin("incompatible", err.msg)
        @test !occursin("non-boolean", err.msg)
    end
end

end
