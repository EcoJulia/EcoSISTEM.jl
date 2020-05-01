using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

epi = TestEpiSystem()

@test sum(epi.abundances.matrix, dims = 2)[:, 1] == epi.epilist.abun
@test_nowarn gettraitrel(epi)
@test gettraitrel(epi) == epi.relationship
@test_nowarn gethabitat(epi)
@test gethabitat(epi) == epi.epienv.habitat
@test_nowarn getsize(epi)
@test getsize(epi) == size(epi.abundances.matrix, 2) .* epi.epienv.habitat.size^2
@test_nowarn getgridsize(epi)
@test getgridsize(epi) == epi.epienv.habitat.size
@test_nowarn getdispersaldist(epi, 1)
@test getdispersaldist(epi, 1) == epi.epilist.movement.kernels[1].dist
@test_nowarn getdispersaldist(epi, "Virus")
@test getdispersaldist(epi, "Virus") == epi.epilist.movement.kernels[1].dist
@test_nowarn getdispersalvar(epi, 1)
@test_nowarn getdispersalvar(epi, "Virus")

# Test saving Function
Simulation.save("testepi.jlso", epi)
loaded_epi = Simulation.load("testepi.jlso", EpiSystem)
# Ideally we'd compare epi and loaded_epi, but `EpiSystem` still does not support comparisons
@test loaded_epi isa EpiSystem
rm("testepi.jlso") # Delete temporary file
