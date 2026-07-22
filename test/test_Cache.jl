# SPDX-License-Identifier: LGPL-3.0-or-later

module TestCache

using EcoSISTEM
using Test
using Distributions
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units

# Build a small, fully deterministic ecosystem. With a fixed `seed` the initial
# abundances (a fixed vector) and the per-species RNG streams are identical every
# time, so runs are byte-for-byte reproducible.
function cache_test_eco(seed)
    numSpecies = 8
    grid = (4, 4)
    area = 100.0km^2
    individuals = 1_000
    totalK = 10000.0kJ / km^2

    resource = SolarDemand(fill(10.0kJ, numSpecies))
    param = EqualPop(0.2 / year, 0.2 / year, 1.0, 0.0, 1.0)
    kernel = fill(GaussianKernel(2.0km, 1.0e-3), numSpecies)
    movement = BirthOnlyMovement(kernel, NoBoundary())
    tolerance = NicheTolerance(MeanTemperature, Normal,
                               fill(274.0K, numSpecies),
                               fill(0.5K, numSpecies))
    native = fill(true, numSpecies)
    abun = fill(div(individuals, numSpecies), numSpecies)

    sppl = SpeciesList(numSpecies, tolerance, abun, resource, movement, param,
                       native)
    habitat = simplehabitatAE(274.0K, grid, totalK, area)
    rel = DistRel{typeof(1.0K)}()
    return Ecosystem(sppl, habitat, rel; seed = seed)
end

@testset "CachedEcosystem save/load/resume" begin
    seed = 20250708
    times = (0.0year):(1.0year):(5.0year)

    # Reference run: never cached, five one-year steps
    ref = cache_test_eco(seed)
    for _ in 1:5
        update!(ref, 1.0year)
    end
    expected = copy(ref.abundances.matrix)

    # Cached run, only simulated (and saved to disk) as far as year 3
    dir = mktempdir()
    eco1 = cache_test_eco(seed)
    cache1 = CachedEcosystem(eco1, dir, times)
    abundances(cache1, 3.0year)
    @test !isempty(filter(f -> endswith(f, ".jld2"), readdir(dir)))

    # Fresh cache over the same folder: its only in-memory state is year 0, so
    # asking for year 5 must load the year-3 snapshot from disk (restoring the
    # per-species RNG streams) and simulate the remaining two years.
    eco2 = cache_test_eco(seed)
    cache2 = CachedEcosystem(eco2, dir, times)
    resumed = abundances(cache2, 5.0year)

    # Resuming from the cache reproduces the uncached run exactly
    @test resumed.matrix == expected

    # clearcache removes the saved files
    @test_nowarn clearcache(cache2)
    @test isempty(filter(f -> endswith(f, ".jld2"), readdir(dir)))
end

@testset "CachedEcosystem matches Ecosystem" begin
    seed = 424242
    times = (0.0year):(1.0year):(5.0year)

    # Plain ecosystem, five one-year steps
    eco = cache_test_eco(seed)
    for _ in 1:5
        update!(eco, 1.0year)
    end

    # Same configuration driven as a CachedEcosystem, run to year 5. Its results
    # come from its own in-memory cache; the disk-reload path is covered above.
    cached = CachedEcosystem(cache_test_eco(seed), mktempdir(), times)
    result = abundances(cached, 5.0year)

    @test result.matrix == eco.abundances.matrix
end

@testset "CachedEcosystem caching interval independence" begin
    seed = 987654
    times = (0.0year):(1.0year):(5.0year)

    # Two caches, same seed and timestep but different disk-save cadence
    dirA = mktempdir()
    dirB = mktempdir()
    cacheA = CachedEcosystem(cache_test_eco(seed), dirA, times;
                             saveinterval = 1.0year)
    cacheB = CachedEcosystem(cache_test_eco(seed), dirB, times;
                             saveinterval = 2.0year)

    rA = abundances(cacheA, 4.0year)
    rB = abundances(cacheB, 4.0year)

    # The caching interval no longer affects the simulation result
    @test rA.matrix == rB.matrix

    # ...but the two caches genuinely saved at different frequencies
    countfiles(d) = length(filter(f -> endswith(f, ".jld2"), readdir(d)))
    @test countfiles(dirA) != countfiles(dirB)
end

end
