# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoSISTEM
# **`Phylo` is a weak dependency**: the species below are given a tree drawn with it, which
# `add_phylogeny!` takes through `EcoSISTEMPhyloExt`, the extension loading `Phylo` activates for
# every test file that includes this one.
# **`import`, not `using`, and that is load-bearing**: extension loading is triggered by the
# package being *loaded*, either way - but `using Phylo` would also drag its exports into every
# includer's namespace, where `Phylo.DiscreteTrait` shadows EcoSISTEM's deprecated binding of the
# same name. Measured: it broke `test_deprecations`' `DiscreteTrait === SimpleCategoricalTolerance`.
import Phylo
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Unitful
using Random
using DataFrames: DataFrame

"""
    Test1Ecosystem(seed = nothing, grid = (5, 7), numspecies = 15, mixed = false)

Build a small test ecosystem. If `seed` is supplied, the initial per-species
abundance totals and the per-species simulation RNGs are both made deterministic,
so the whole run is reproducible regardless of the number of threads. With
`mixed`, the species take four kinds of movement in turn - dispersing at birth in
two modes, always, and not at all - the most the step loop splits.
"""
function Test1Ecosystem(; seed = nothing, grid = (5, 7), numspecies = 15,
                        mixed = false)
    # **Non-square, and that is not incidental.** The old grid was 10 × 10; a square fixture
    # cannot see a y/x transposition.
    #
    # `numspecies` is overridable for the same reason `grid` is: the hot loop's allocation has to be
    # shown flat in *both* directions, and a fixture that can only vary one of them can only prove
    # half of it.
    numSpecies = numspecies
    numNiches = 2

    birth = 0.6 / month_mean_duration
    death = 0.6 / month_mean_duration
    long = 1.0
    surv = 0.0
    timestep = 1.0month_mean_duration

    # `(ny, nx)` - rows then columns, as everywhere in the package. Overridable so a test can
    # compare two grid sizes, which is how the hot loop's per-cell allocation is measured.
    cellsize = 2.0km
    individuals = 2000 * numSpecies
    totalK = 1000000.0 * kJ / km^2 / day * numSpecies

    # The grid is decided first and only then built on. `extent` and `grid` are both `(y, x)`, so
    # they are written in the same order -- stated explicitly rather than derived from a total area.
    studyarea = StudyArea(extent = (grid[1] * cellsize, grid[2] * cellsize),
                          cellsize = cellsize, verbosity = :silent)
    # This `NicheSpec` is unseeded, so it draws its pattern from the **global** RNG, and
    # `Random.seed!` below runs *after* the environment is built, so the niche layout is not part of
    # what `seed` pins. It does not need to be - a run is reproducible because the habitat is built
    # once and then shared.
    habitat = GridHabitat(regime = NicheSpec(numNiches,
                                             axis = EcoSISTEM.TypologyAxis),
                          supply = UniformSpec(totalK,
                                               axis = SolarRadiation),
                          area = studyarea)

    # Seed the global RNG so the initial abundance totals are also deterministic
    isnothing(seed) || Random.seed!(seed)
    abun = rand(Multinomial(individuals, numSpecies))

    # Each species' movement by its `movement` column: with `mixed`, every fourth stays put and the
    # rest take three modes in turn.
    kernel = (meandist = 1.0km, thresh = 10e-04)
    movement = ChosenBy(:movement,
                        :insect => MovementSpec(InsectDispersal(); kernel...),
                        :always => MovementSpec(AlwaysDispersal(); kernel...),
                        :wind => MovementSpec(WindDispersal(GaussianKernel);
                                     kernel...),
                        :none => MovementSpec(NoDispersal()))
    resource = Demand{SolarRadiation}(fill(2.0kJ / day, numSpecies))
    # A random tree, and niche classes evolved along it, drawn from the global generator after the
    # counts above; the pool takes both, and its species are phylogenetic.
    names = string.(1:numSpecies)
    tree = rand(Random.default_rng(),
                Phylo.Ultrametric{Phylo.BinaryTree{Phylo.OneRoot, DataFrame,
                                                   DataFrame}}(names))
    EcoSISTEM.assigntraits!(Random.default_rng(), tree, [0.5],
                            DataFrame(trait1 = collect(1:numNiches)))
    classes = EcoSISTEM.gettraits(tree, true)
    niche = Dict(zip(classes.species, classes.trait1))
    modes = mixed ? (:insect, :always, :wind) : (:insect,)
    pool = SpeciesPool(names, verbosity = :silent)
    add_phylogeny!(pool, tree, verbosity = :silent)
    add_trait!(pool,
               DataFrame(species = names, niche = [niche[n] for n in names],
                         movement = [mixed && sp % 4 == 0 ? :none :
                                     modes[mod1(sp, length(modes))]
                                     for sp in eachindex(names)]),
               verbosity = :silent)
    # `penalty = 0.5` is soft exclusion: a species does worse outside its niche class rather than
    # being unable to live there.
    tolerance = SimpleCategoricalTolerance(Int.(pool.species.niche),
                                           axis = EcoSISTEM.TypologyAxis,
                                           penalty = 0.5)
    sppl = build_species(pool, tolerance = tolerance, demand = resource,
                         movement = movement, birth = birth, death = death,
                         longevity = long, survival = surv,
                         verbosity = :silent)

    nichefit = CategoricalSuitability{EcoSISTEM.axisof(habitat.regime),
                                      eltype(habitat.regime)}()
    # Each species' drawn count placed by supply; serial even under MPI.
    return build_ecosystem(BySupply(each = row -> abun[row.key]), sppl, habitat,
                           nichefit = nichefit, seed = seed,
                           distributed = false)
end

function TestMultiEcosystem()
    numSpecies = 150

    birth = 0.6 / month_mean_duration
    death = 0.6 / month_mean_duration
    long = 1.0
    surv = 0.0
    timestep = 1.0month_mean_duration
    param = EqualPop(birth, death, long, surv)

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

    kernel = fill(GaussianKernel(meandist = 1.0km, thresh = 10e-04), numSpecies)
    movement = InsectDispersal().(kernel)
    resource1 = Demand{SolarRadiation}(fill(2.0kJ / day, numSpecies))
    resource2 = Demand{Precipitation}(fill(2.0Unitful.L / day, numSpecies))
    resource = SpeciesRequirementCollection((resource1, resource2))
    tolerance = NicheTolerance(Temperature, Normal, fill(10.0K, numSpecies),
                               fill(0.1K, numSpecies))
    sppl = SpeciesList(numSpecies, tolerance, resource, movement, param)

    nichefit = NicheSuitability{EcoSISTEM.axisof(habitat.regime),
                                eltype(habitat.regime)}()
    return build_ecosystem(BySupply(each = row -> abun[row.key]), sppl, habitat,
                           nichefit = nichefit, distributed = false)
end
