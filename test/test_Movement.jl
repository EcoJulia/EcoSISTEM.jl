using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

numSpecies = 10

@test_nowarn GaussianKernel.(fill(0.2km, numSpecies), 10e-4)
@test_nowarn LongTailKernel.(fill(0.2km, numSpecies), 1.0, 10e-4)
kernel = GaussianKernel.(fill(0.2m, numSpecies), 10e-4)
@test_nowarn AlwaysMovement(kernel, NoBoundary())
@test_nowarn AlwaysMovement(kernel, Cylinder())
@test_nowarn AlwaysMovement(kernel, Torus())
@test_nowarn BirthOnlyMovement(kernel)
@test_nowarn NoMovement(kernel)


kernel = GaussianKernel.(fill(0.2m, numSpecies), 10e-4)
mov = AlwaysMovement(kernel, NoBoundary())
@test Simulation.getkernels(mov) == mov.kernels
@test Simulation.getboundary(mov) == mov.boundary
mov = BirthOnlyMovement(kernel)
@test Simulation.getkernels(mov) == mov.kernels
@test Simulation.getboundary(mov) == mov.boundary
mov =  NoMovement(kernel)
@test Simulation.getkernels(mov) == mov.kernels
kernel = LongTailKernel.(fill(0.2m, numSpecies), 1.0, 10e-4)
mov = AlwaysMovement(kernel, NoBoundary())
@test Simulation.getkernels(mov) == mov.kernels
@test Simulation.getboundary(mov) == mov.boundary
