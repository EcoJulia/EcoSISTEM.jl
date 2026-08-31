# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Canonical results from **simulated** input: no data is read, so nothing here can change because a
# download changed. Every number is a pure function of the code plus the seed, which makes this the
# file that isolates *model* change from *data* change - see README.md.

module CanonicalSimulated

using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Distributions
using Random

include("canonical.jl")
using .Canonical

# One fixed seed for the whole file. `Ecosystem(...; seed)` builds one deterministic RNG stream per
# species, so a result is reproducible regardless of thread or MPI-rank count - that independence is
# what makes a blessed number meaningful rather than an artefact of the machine it was recorded on.
const SEED = 20260804

# A small, fully specified ecosystem. Deliberately modest (10×10, 20 species, 2 years) so the whole
# canonical suite stays runnable on demand; a canonical test earns its keep by being *run*, and one
# that takes an hour will not be.
function testecosystem(; numspecies = 20, grd = (10, 10), area = 100.0km^2,
                       temperature = 298.0K, seed = SEED)
    totalK = (4.5e11kJ / km^2 / day, 192.0mm / day)
    # One `GridHabitat` in place of two deprecated `simplehabitat` calls whose regimes were
    # thrown away. Verified bit-identical - regime, both supplies and the active mask - which is
    # why `reference.toml` does not move.
    side = sqrt(area)
    studyarea = StudyArea(extent = (side, side), cellsize = side / grd[1],
                          verbosity = :silent)
    habitat = GridHabitat(regime = UniformSpec(temperature,
                                               axis = Temperature),
                          supply = (UniformSpec(totalK[1],
                                                axis = SolarRadiation),
                                    UniformSpec(totalK[2],
                                                axis = Precipitation)),
                          area = studyarea,
                          topology = Torus())
    vars = fill(2.0K, numspecies)
    opts = temperature .+ vars .* range(-3, stop = 3, length = numspecies)
    demand = SpeciesRequirementCollection((Demand{SolarRadiation}(fill(450000.0kJ /
                                                                       m^2 /
                                                                       day *
                                                                       1.0m^2,
                                                                       numspecies)),
                                           Demand{Precipitation}(fill(192.0Unitful.L /
                                                                      m^2 /
                                                                      day *
                                                                      1.0m^2,
                                                                      numspecies))))
    movement = BirthOnlyMovement(GaussianKernel.(fill(2.4km, numspecies),
                                                 10e-10))
    param = EqualPop(0.15 / year, 0.15 / year, 1.0, 0.1, 1.0)
    # The initial abundances are drawn from the *seeded* stream too, not from the global RNG, or
    # the starting state would differ between runs and nothing downstream could be blessed.
    rng = Random.Xoshiro(seed)
    abun = rand(rng, Multinomial(1_000_000, numspecies))
    tolerance = NicheTolerance(Temperature, Normal, opts, vars)
    sppl = SpeciesList(numspecies, tolerance, abun, demand, movement, param,
                       fill(true, numspecies))
    return Ecosystem(sppl, habitat, NicheSuitability(tolerance), seed = seed)
end

@testset "canonical: simulated ecosystem" begin
    eco = testecosystem()
    simulate!(eco, 2year, 1month_mean_duration)
    abun = eco.abundances.matrix

    # Blessed totals. These are the numbers that move when the *model* changes - birth/death
    # arithmetic, dispersal, the resource ratio - and cannot move because a raster was re-released.
    canonical("simulated/total_abundance", sum(abun))
    canonical("simulated/occupied_cells", count(>(0), sum(abun, dims = 1)))
    canonical("simulated/species_surviving", count(>(0), sum(abun, dims = 2)))
    # Per-species totals, so a change that preserves the grand total but redistributes it is still
    # caught - the failure mode a single summary number hides.
    canonical("simulated/abundance_by_species", vec(sum(abun, dims = 2)))

    # The properties that must hold whatever the blessed numbers are. Kept alongside them
    # deliberately: a blessed number tells you *something changed*, a property tells you *the model is
    # still wrong or right*. Re-blessing silences the first and must never silence the second.
    @test sum(abun) > 0
    @test all(>=(0), abun)
    # Species nearest the habitat's own temperature should end up commonest.
    bysp = vec(sum(abun, dims = 2))
    @test argmax(bysp) in (9, 10, 11, 12)      # the middle of a 20-species optimum gradient

    # Reproducibility is the premise of the whole file, so it is asserted rather than assumed.
    second = testecosystem()
    simulate!(second, 2year, 1month_mean_duration)
    @test sum(second.abundances.matrix) == sum(abun)
    @test second.abundances.matrix == abun
end

end
