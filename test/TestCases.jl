using EcoSISTEM
using EcoSISTEM.ClimatePref
using JLD
using AxisArrays
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Unitful

function TestEcosystem()
    numSpecies = 150
    numNiches = 2

    birth = 0.6/month
    death = 0.6/month
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month
    param = EqualPop(birth, death, long, surv, boost)

    grid = (10, 10)
    area = 10000.0km^2
    individuals=20000 * numSpecies
    totalK = 1000000.0 * kJ/km^2 * numSpecies
    abenv = simplenicheAE(numNiches, grid, totalK, area)

    abun = rand(Multinomial(individuals, numSpecies))

    kernel = GaussianKernel.(fill(1.0km, numSpecies), 10e-04)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, numSpecies)
    energy = SolarRequirement(fill(2.0kJ, numSpecies))
    sppl = SpeciesList(numSpecies, numNiches, abun, energy, movement, param, native)

    rel = Match{eltype(abenv.habitat)}()
    eco = Ecosystem(sppl, abenv, rel)
    return eco
end

function TestMultiEcosystem()
    numSpecies = 150

    birth = 0.6/month
    death = 0.6/month
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month
    param = EqualPop(birth, death, long, surv, boost)

    grid = (10, 10)
    area = 10000.0km^2
    individuals=20000 * numSpecies
    totalK1 = 1000000.0 * kJ/km^2 * numSpecies
    totalK2 = 100.0 * mm/km^2 * numSpecies
    abenv1 = simplehabitatAE(10.0K, grid, totalK1, area)
    abenv2 = simplehabitatAE(10.0K, grid, totalK2, area)
    budget = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(budget)}(abenv1.habitat, abenv1.active, budget, abenv1.names)

    abun = rand(Multinomial(individuals, numSpecies))

    kernel = GaussianKernel.(fill(1.0km, numSpecies), 10e-04)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, numSpecies)
    energy1 = SolarRequirement(fill(2.0kJ, numSpecies))
    energy2 = WaterRequirement(fill(2.0mm, numSpecies))
    energy = ReqCollection2(energy1, energy2)
    traits = GaussTrait(fill(10.0K, numSpecies), fill(0.1K, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, energy, movement, param, native)

    rel = Gauss{eltype(abenv.habitat)}()
    eco = Ecosystem(sppl, abenv, rel)
    return eco
end
