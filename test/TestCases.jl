using Simulation
using Simulation.ClimatePref
using JLD
using AxisArrays
using Unitful.DefaultSymbols
using Distributions
using Simulation.Units
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

    grid = (50, 50)
    area = 10000.0km^2
    individuals=20000 * numSpecies
    totalK = 1000000.0 * numSpecies
    abenv = simplenicheAE(numNiches, grid, totalK, area)

    abun = rand(Multinomial(individuals, numSpecies))

    kernel = GaussianKernel(1.0km, numSpecies, 10e-04)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, numSpecies)
    energy = SimpleRequirement(fill(2.0, numSpecies))
    sppl = SpeciesList(numSpecies, numNiches, abun, energy, movement, param, native)

    rel = Match{eltype(abenv.habitat)}()
    eco = Ecosystem(trait_populate!, sppl, abenv, rel)
    return eco
end
