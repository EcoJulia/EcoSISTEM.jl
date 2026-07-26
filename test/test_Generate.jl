# SPDX-License-Identifier: LGPL-3.0-or-later

module TestGenerate

using EcoSISTEM
using Test
using Unitful
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using JLD2

include("TestCases.jl")

@testset "Update functions" begin
    eco = Test1Ecosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)) ==
          Vector{Int64}
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.habitat,
                                                        eco.spplist),
                           eco.spplist, eco.habitat, eco.nichefit, eco.rngs)
    @test_nowarn repopulate!(eco)

    # Test Cylinder
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Cylinder(), 1, 1, 1, eco, 10)
    # Test Torus
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Torus(), 1, 1, 1, eco, 10)

    # Test ecosystem with multiple supplies
    eco = TestMultiEcosystem()
    @test_nowarn update!(eco, 1month)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                              1, 1, 1, eco, 10)) ==
          Vector{Int64}
    @test_nowarn populate!(EcoSISTEM.emptygridlandscape(eco.habitat,
                                                        eco.spplist),
                           eco.spplist, eco.habitat, eco.nichefit, eco.rngs)
    @test_nowarn repopulate!(eco)
end

@testset "NoGrowth freezes with one or two supplies" begin
    # NoGrowth must zero the birth/death adjustment regardless of how many resource
    # supplies the environment has (the two-supply path previously skipped it).
    N = 8
    grid = (4, 4)
    area = 100.0km^2
    birth = 0.1 / month
    nogrowth = NoGrowth{typeof(unit(birth))}(fill(birth, N), fill(birth, N),
                                             1.0, 0.0, 1.0)
    tolerance = NicheTolerance(Temperature, Normal, fill(274.0K, N),
                               fill(0.5K, N))
    movement = BirthOnlyMovement(fill(GaussianKernel(1.0km, 1.0e-3), N))
    native = fill(true, N)
    abun = fill(10, N)
    nichefit = NicheSuitability{typeof(1.0K)}()

    # single supply
    sppl1 = SpeciesList(N, tolerance, abun, SolarDemand(fill(10.0kJ / day, N)),
                        movement, nogrowth, native)
    habitat1 = simplehabitat(274.0K, grid, 10000.0kJ / km^2 / day, area)
    eco1 = Ecosystem(sppl1, habitat1, nichefit)
    @test EcoSISTEM.resource_adjustment(eco1, eco1.habitat.supply, 1, 1) ==
          (0.0, 0.0)

    # two supplies (the previously-buggy path)
    resource2 = DemandCollection2(SolarDemand(fill(10.0kJ / day, N)),
                                  WaterDemand(fill(2.0Unitful.L / day, N)))
    sppl2 = SpeciesList(N, tolerance, abun, resource2, movement, nogrowth,
                        native)
    habitat_solar = simplehabitat(274.0K, grid, 10000.0kJ / km^2 / day, area)
    habitat_water = simplehabitat(274.0K, grid, 10.0mm / day, area)
    supply = SupplyCollection2(habitat_solar.supply, habitat_water.supply)
    habitat2 = GridHabitat{typeof(habitat_solar.regime), typeof(supply)}(habitat_solar.regime,
                                                                         habitat_solar.active,
                                                                         supply,
                                                                         habitat_solar.names)
    eco2 = Ecosystem(sppl2, habitat2, nichefit)
    @test EcoSISTEM.resource_adjustment(eco2, eco2.habitat.supply, 1, 1) ==
          (0.0, 0.0)
end

@testset "non-square grid: y/x dimension order regression" begin
    # A square grid cannot detect an X/Y mixup (sizes and indices coincide by
    # accident). This grid is deliberately asymmetric, so confused code either
    # throws a BoundsError (indexing a 2-row array with a column-range value up
    # to 6) or silently reads the wrong cells.
    N = 4
    ny, nx = 2, 6
    grid = (ny, nx)
    area = 100.0km^2

    # A partial active mask — only column x=1 is active, every row — to check
    # the active-mask indexing itself corresponds to (y, x), not (x, y):
    # deliberately asymmetric (not a full row/column/diagonal) so a swap would
    # show up as reading the wrong cells, not coincidentally the same ones.
    partial_active = fill(false, ny, nx)
    partial_active[:, 1] .= true
    partial_habitat = simplehabitat(274.0K, grid, 10000.0kJ / km^2 / day, area,
                                    partial_active)
    @test size(partial_habitat.regime.matrix) == grid
    @test size(partial_habitat.active) == grid
    @test all(partial_habitat.active[y, 1] for y in 1:ny)
    @test all(!partial_habitat.active[y, x] for y in 1:ny, x in 2:nx)

    tolerance = NicheTolerance(Temperature, Normal, fill(274.0K, N),
                               fill(0.5K, N))
    param = EqualPop(0.6 / month, 0.6 / month, 1.0, 0.0, 1000.0)
    kernel = GaussianKernel.(fill(1.0km, N), 1.0e-3)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, N)
    abun = fill(10, N)
    resource = SolarDemand(fill(2.0kJ / day, N))
    sppl = SpeciesList(N, tolerance, abun, resource, movement, param, native)

    # A fully-active habitat for exercising calc_lookup_moves!/update! at
    # extreme y/x coordinates — every cell reachable, so no cell's move
    # probabilities can go to zero regardless of dimension order.
    habitat = simplehabitat(274.0K, grid, 10000.0kJ / km^2 / day, area)
    @test size(habitat.regime.matrix) == grid
    @test size(habitat.active) == grid

    nichefit = NicheSuitability{eltype(habitat.regime)}()
    eco = Ecosystem(sppl, habitat, nichefit)

    # getdimension must report (ny, nx) — the actual array shape — not (nx, ny).
    @test EcoSISTEM.getdimension(eco) == grid

    # Calling calc_lookup_moves! at every combination of extreme y/x coordinates
    # must not throw — a BoundsError here is exactly what an X/Y mixup produces
    # on a non-square grid (e.g. checking a y=2-valid offset against nx=6's
    # bound instead of ny=2's, then indexing a nonexistent row of `active`).
    for (y, x) in ((1, 1), (ny, 1), (1, nx), (ny, nx))
        @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.spplist.movement.boundary,
                                                  y, x, 1, eco, 5)
        @test_nowarn EcoSISTEM.calc_lookup_moves!(Cylinder(), y, x, 1, eco, 5)
        @test_nowarn EcoSISTEM.calc_lookup_moves!(Torus(), y, x, 1, eco, 5)
    end

    # convert_coords must round-trip correctly through the (y,x) convention —
    # height (dimension 1) is ny, not nx.
    for i in 1:(ny * nx)
        (y, x) = EcoSISTEM.convert_coords(eco, i, ny)
        @test 1 <= y <= ny && 1 <= x <= nx
        @test EcoSISTEM.convert_coords(y, x, ny) == i
    end

    # A full timestep must run without error on this asymmetric grid.
    @test_nowarn update!(eco, 1month)
end

@testset "Multithreaded reproducibility" begin
    # A seeded run must give identical results regardless of the number of
    # threads. Run the same seeded simulation in child processes forced to use
    # 1, 2 and 4 threads (so the parallel `update!` path is genuinely exercised)
    # and check that all three final abundance matrices are identical. Running
    # the 2-thread case twice also confirms same-thread-count repeatability, so a
    # failure distinguishes "not reproducible at all" from "not reproducible
    # across thread counts".
    script = pkgdir(EcoSISTEM, "test", "threading_reproducibility.jl")
    dir = mktempdir()
    function run_with(nthreads, tag)
        out = joinpath(dir, "repro_$tag.jld2")
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t $nthreads $script $out`
        @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
        return load(out, "matrix")
    end
    m1 = run_with(1, "t1")
    m2 = run_with(2, "t2")
    m4 = run_with(4, "t4")
    m2b = run_with(2, "t2b")
    @test m1 == m2 == m4
    @test m2 == m2b
end

end
