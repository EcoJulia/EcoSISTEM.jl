using Simulation
using Compat.Test
using Distributions
using Unitful.DefaultSymbols
using MyUnitful

numSpecies = 10

@test_nowarn GaussianKernel(0.2km, numSpecies, 10e-4)
kernel = GaussianKernel(0.2m, numSpecies, 10e-4)
@test_nowarn AlwaysMovement(kernel, NoBoundary())
@test_nowarn AlwaysMovement(kernel, Cylinder())
@test_nowarn AlwaysMovement(kernel, Torus())
@test_nowarn BirthOnlyMovement(kernel)
@test_nowarn NoMovement(kernel)
