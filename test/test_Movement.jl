using EcoSISTEM
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units

@testset "Movement" begin
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
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    @test EcoSISTEM.getboundary(mov) == mov.boundary
    @test EcoSISTEM.getdispersaldist(mov, 1) == mov.kernels[1].dist
    @test EcoSISTEM.getdispersalvar(mov, 1) == (mov.kernels[1].dist)^2 * pi / 4
    mov = BirthOnlyMovement(kernel)
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    @test EcoSISTEM.getboundary(mov) == mov.boundary
    @test EcoSISTEM.getdispersaldist(mov, 1) == mov.kernels[1].dist
    @test EcoSISTEM.getdispersalvar(mov, 1) == (mov.kernels[1].dist)^2 * pi / 4
    mov =  NoMovement(kernel)
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    @test EcoSISTEM.getdispersaldist(mov, 1) == "No movement takes place"
    @test EcoSISTEM.getdispersalvar(mov, 1) == "No movement takes place"
    kernel = LongTailKernel.(fill(0.2m, numSpecies), 1.0, 10e-4)
    mov = AlwaysMovement(kernel, NoBoundary())
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    @test EcoSISTEM.getboundary(mov) == mov.boundary
end
