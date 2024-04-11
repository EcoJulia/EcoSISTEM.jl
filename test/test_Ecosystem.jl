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

    sppl = SpeciesList{typeof(eco.spplist.traits), typeof(eco.spplist.requirement), typeof(eco.spplist.movement), UniqueTypes, typeof(eco.spplist.params)}(eco.spplist.names, eco.spplist.traits, eco.spplist.abun, eco.spplist.requirement, UniqueTypes(length(eco.spplist.names)), eco.spplist.movement, eco.spplist.params, eco.spplist.native)
    eco = Ecosystem(sppl, eco.abenv, eco.relationship)
    @test_nowarn addspecies!(eco, 10)

    @testset "get functions" begin
        # Test Simulation get functions
        @test_nowarn gettraitrel(eco)
        @test gettraitrel(eco) == eco.relationship
        @test_nowarn gethabitat(eco)
        @test gethabitat(eco) == eco.abenv.habitat
        @test_nowarn getsize(eco)
        @test getsize(eco) == size(eco.abundances.matrix, 2) .* eco.abenv.habitat.size^2
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
        @test_nowarn resetrate!(eco, 0.1/s)
        @test eco.abenv.habitat.change.rate == 0.1/s
        @test_throws MethodError resettime!(eco)
    end
    @testset "movement types" begin
        # Test other movement types
        eco = Test1Ecosystem()
        mov = AlwaysMovement(fill(LongTailKernel(10.0km, 10.0, 1e-10), length(eco.spplist.names)), eco.spplist.movement.boundary)
        sppl = SpeciesList{typeof(eco.spplist.traits), typeof(eco.spplist.requirement), typeof(mov), typeof(eco.spplist.types), typeof(eco.spplist.params)}(eco.spplist.names, eco.spplist.traits, eco.spplist.abun, eco.spplist.requirement, eco.spplist.types, mov, eco.spplist.params, eco.spplist.native)
        @test_nowarn Ecosystem(sppl, eco.abenv, eco.relationship)
    end
    @testset "diversity"  begin
        # Test Diversity get functions
        @test_nowarn getmetaabundance(eco)
        @test_nowarn getpartition(eco)
        @test getpartition(eco) == eco.abenv
        @test_nowarn gettypes(eco)
        @test gettypes(eco) == eco.spplist
        @test_nowarn getordinariness!(eco)
        @test getordinariness!(eco) == eco.ordinariness
    end
end

end
