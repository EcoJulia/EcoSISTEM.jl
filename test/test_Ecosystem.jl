using Simulation
using Test
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units
using Diversity

include("TestCases.jl")

@test_nowarn eco = TestEcosystem()
eco = TestEcosystem()

# Test Simulation get functions
@test_nowarn gettraitrel(eco)
@test gettraitrel(eco) == eco.relationship
@test_nowarn gethabitat(eco)
@test gethabitat(eco) == eco.abenv.habitat
@test_nowarn getsize(eco)
@test getsize(eco) == size(eco.abundances.matrix, 2) .* eco.abenv.habitat.size^2
@test_nowarn getgridsize(eco)
@test getgridsize(eco) == eco.abenv.habitat.size
@test_nowarn getdispersaldist(eco)
@test getdispersaldist(eco) == eco.spplist.movement.kernel.dist
@test_nowarn getdispersaldist(eco, 1)
@test getdispersaldist(eco, 1) == eco.spplist.movement.kernel.dist[1]
@test_nowarn getdispersaldist(eco, "1")
@test getdispersaldist(eco, "1") == eco.spplist.movement.kernel.dist[1]
@test_nowarn getdispersalvar(eco)
@test getdispersalvar(eco) == (eco.spplist.movement.kernel.dist).^2 .* pi ./ 4
@test_nowarn getdispersalvar(eco, 1)
@test_nowarn getdispersalvar(eco, "1")


# Test Diversity get functions
@test_nowarn getmetaabundance(eco)
@test_nowarn getpartition(eco)
@test getpartition(eco) == eco.abenv
@test_nowarn gettypes(eco)
@test gettypes(eco) == eco.spplist
@test_nowarn getordinariness!(eco)
@test getordinariness!(eco) == eco.ordinariness
