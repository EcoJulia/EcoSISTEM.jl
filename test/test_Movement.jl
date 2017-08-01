using Simulation
using Base.Test
using Distributions

@test_nowarn GaussianKernel(0.2, numSpecies, 10e-4)
kernel = GaussianKernel(0.2, numSpecies, 10e-4)
@test_nowarn AlwaysMovement(kernel)
@test_nowarn BirthOnlyMovement(kernel)
@test_nowarn NoMovement(kernel)
