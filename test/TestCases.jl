# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoSISTEM
# **`Phylo` is a weak dependency**, so the `SpeciesList` constructor below — the one that *builds*
# a phylogeny rather than taking tolerances directly — lives in `EcoSISTEMPhyloExt` and only exists
# once `Phylo` is loaded. Loading it here is what activates that extension for every test file that
# includes this one.
# **`import`, not `using`, and that is load-bearing**: extension loading is triggered by the
# package being *loaded*, either way — but `using Phylo` would also drag its exports into every
# includer's namespace, where `Phylo.DiscreteTrait` shadows EcoSISTEM's deprecated binding of the
# same name. Measured: it broke `test_deprecations`' `DiscreteTrait === SimpleCategoricalTolerance`.
import Phylo
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Unitful
using Random

"""
    Test1Ecosystem(seed = nothing)

Build a small test ecosystem. If `seed` is supplied, the initial per-species
abundance totals and the per-species simulation RNGs are both made deterministic,
so the whole run is reproducible regardless of the number of threads.
"""
function Test1Ecosystem(; seed = nothing)
    # **Shrunk 2026-08-13**: 150 species on a 10 × 10 grid became 15 on 5 × 7. Work per timestep is
    # roughly species × cells, so this is ~40× less of it, and this fixture is shared by seven test
    # files.
    #
    # **Non-square, and that is not incidental.** The old grid was 10 × 10; a square fixture
    # cannot see a y/x transposition, which is exactly how `plot(eco)`'s `BoundsError` went unnoticed
    # (found 2026-08-13 the moment a notebook grid stopped being square) and how `get_neighbours`
    # once kept x-first parameter names after the index order switched. 5 × 7 costs 35 cells against
    # a square 5 × 5's 25 — ten cells to keep the property.
    numSpecies = 15
    numNiches = 2

    birth = 0.6 / month_mean_duration
    death = 0.6 / month_mean_duration
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month_mean_duration
    param = EqualPop(birth, death, long, surv, boost)

    # `(ny, nx)` — rows then columns, as everywhere in the package.
    grid = (5, 7)
    cellsize = 2.0km
    individuals = 2000 * numSpecies
    totalK = 1000000.0 * kJ / km^2 / day * numSpecies

    # The grid is decided first and only then built on. `extent` and `grid` are both `(y, x)`, so
    # they are written in the same order -- stated explicitly rather than derived from a total area,
    # which is what forced the old builder's `sqrt(area / prod(dimension))`.
    studyarea = StudyArea(extent = (grid[1] * cellsize, grid[2] * cellsize),
                          cellsize = cellsize, verbosity = :silent)
    # `NicheSpec` draws its pattern from the **global** RNG and takes no seed of its own (see the
    # deferred `NicheSpec` seed question). That is unchanged here: the old builder did the same, and
    # `Random.seed!` below still runs *after* the environment is built, so the niche layout was
    # never part of what `seed` pins. It does not need to be — a run is reproducible because the
    # habitat is built once and then shared.
    habitat = GridHabitat(regime = NicheSpec(numNiches,
                                             axis = EcoSISTEM.TypologyAxis),
                          supply = UniformSpec(totalK,
                                               axis = SolarRadiation),
                          area = studyarea)

    # Seed the global RNG so the initial abundance totals are also deterministic
    isnothing(seed) || Random.seed!(seed)
    abun = rand(Multinomial(individuals, numSpecies))

    kernel = GaussianKernel.(fill(1.0km, numSpecies), 10e-04)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, numSpecies)
    resource = Demand{SolarRadiation}(fill(2.0kJ / day, numSpecies))
    sppl = SpeciesList(numSpecies, numNiches, abun, resource, movement, param,
                       native)

    nichefit = CategoricalSuitability{EcoSISTEM.axisof(habitat.regime),
                                      eltype(habitat.regime)}()
    eco = isnothing(seed) ? Ecosystem(sppl, habitat, nichefit) :
          Ecosystem(sppl, habitat, nichefit, seed = seed)
    return eco
end

function TestMultiEcosystem()
    numSpecies = 150

    birth = 0.6 / month_mean_duration
    death = 0.6 / month_mean_duration
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month_mean_duration
    param = EqualPop(birth, death, long, surv, boost)

    grid = (10, 10)
    area = 10000.0km^2
    individuals = 20000 * numSpecies
    totalK1 = 1000000.0 * kJ / km^2 / day * numSpecies
    totalK2 = 100.0 * mm / day * numSpecies
    side = sqrt(area)
    studyarea = StudyArea(extent = (side, side), cellsize = side / grid[1],
                          verbosity = :silent)
    habitat = GridHabitat(regime = UniformSpec(10.0K,
                                               axis = Temperature),
                          supply = (UniformSpec(totalK1,
                                                axis = SolarRadiation),
                                    UniformSpec(totalK2,
                                                axis = Precipitation)),
                          area = studyarea)

    abun = rand(Multinomial(individuals, numSpecies))

    kernel = GaussianKernel.(fill(1.0km, numSpecies), 10e-04)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, numSpecies)
    resource1 = Demand{SolarRadiation}(fill(2.0kJ / day, numSpecies))
    resource2 = Demand{Precipitation}(fill(2.0Unitful.L / day, numSpecies))
    resource = SpeciesRequirementCollection((resource1, resource2))
    tolerance = NicheTolerance(Temperature, Normal, fill(10.0K, numSpecies),
                               fill(0.1K, numSpecies))
    sppl = SpeciesList(numSpecies, tolerance, abun, resource, movement, param,
                       native)

    nichefit = NicheSuitability{EcoSISTEM.axisof(habitat.regime),
                                eltype(habitat.regime)}()
    eco = Ecosystem(sppl, habitat, nichefit)
    return eco
end
