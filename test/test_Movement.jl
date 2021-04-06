using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units

numSpecies = 10

@test_nowarn GaussianKernel.(fill(0.2km, numSpecies), 10e-4)
kernel = GaussianKernel.(fill(0.2m, numSpecies), 10e-4)
@test_nowarn AlwaysMovement(kernel, NoBoundary())
@test_nowarn AlwaysMovement(kernel, Cylinder())
@test_nowarn AlwaysMovement(kernel, Torus())
@test_nowarn BirthOnlyMovement(kernel)
@test_nowarn NoMovement(kernel)
