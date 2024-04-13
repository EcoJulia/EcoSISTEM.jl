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
    @test sum(eco.abundances.matrix, dims = 2)[:, 1] == eco.spplist.species.abun
    @test EcoSISTEM.tematch(eco.spplist, eco.abenv) == true
    @test EcoSISTEM.trmatch(eco.spplist, eco.relationship) == true

    eco = makeunique(eco)
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
        @test getdispersaldist(eco, 1) == eco.spplist.species.movement.kernels[1].dist
        @test_nowarn getdispersaldist(eco, "1")
        @test getdispersaldist(eco, "1") == eco.spplist.species.movement.kernels[1].dist
        @test_nowarn getdispersalvar(eco, 1)
        @test_nowarn getdispersalvar(eco, "1")
        @test EcoSISTEM.getlookup(eco, 1) == eco.lookup.species[1]
        @test_nowarn resetrate!(eco, 0.1/s)
        @test eco.abenv.habitat.change.rate == 0.1/s
        @test_throws MethodError resettime!(eco)
    end
    @testset "movement types" begin
        # Test other movement types
        eco = Test1Ecosystem()
        mov = AlwaysMovement(fill(LongTailKernel(10.0km, 10.0, 1e-10), length(eco.spplist.species.names)), eco.spplist.species.movement.boundary)
        species = SpeciesTypes{typeof(eco.spplist.species.traits), typeof(eco.spplist.species.requirement), typeof(mov), typeof(eco.spplist.species.types)}(eco.spplist.species.names, eco.spplist.species.traits, eco.spplist.species.abun, eco.spplist.species.requirement, eco.spplist.species.types, mov, eco.spplist.species.native)
        sppl = SpeciesList{typeof(species), NoPathogen, typeof(eco.spplist.params)}(species, NoPathogen(), eco.spplist.params)
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
