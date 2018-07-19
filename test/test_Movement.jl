using Simulation
using Base.Test
using Distributions
using Unitful.DefaultSymbols
using MyUnitful

@test_nowarn GaussianKernel(0.2km, numSpecies, 10e-4)
kernel = GaussianKernel(0.2m, numSpecies, 10e-4)
@test_nowarn AlwaysMovement(kernel)
@test_nowarn BirthOnlyMovement(kernel)
@test_nowarn NoMovement(kernel)
