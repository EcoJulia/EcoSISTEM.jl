using Simulation
using Simulation.ClimatePref
using JLD
using AxisArrays
using Unitful.DefaultSymbols
using Distributions
using SimulationData.Units
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

    grid = (50, 50)
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

function TestEpiSystem()
    numvirus = 2
    numclasses = 4
    birth = [fill(1e-5/day, numclasses - 1); 0.0/day]
    death = [fill(1e-5/day, numclasses - 1); 0.0/day]
    beta_force = 5.0/day
    beta_env = 0.5/day
    sigma = 0.05/day
    virus_growth = 0.0001/day
    virus_decay = 0.07/day
    param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
    param = transition(param)

    grid = (2, 2)
    area = 10.0km^2
    epienv = simplehabitatAE(298.0K, grid, area, NoControl())

    abun_h = (
        Susceptible = 1000,
        Infected = 1,
        Recovered = 0,
        Dead = 0
    )
    disease_classes = (
        susceptible = ["Susceptible"],
        infectious = ["Infected"]
    )
    abun_v = (Environment = 0, Force = 0)

    dispersal_dists = fill(2.0km, grid[1] * grid[2])
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)

    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    return epi
end
function TestEpiSystemFromPopulation(
    initial_pop::AbstractMatrix;
    epienv_active=fill(true, size(initial_pop))
)
    numclasses = 4
    numvirus = 2
    birth = [fill(1e-5/day, numclasses - 1); 0.0/day]
    death = [fill(1e-5/day, numclasses - 1); 0.0/day]
    beta_force = 5.0/day
    beta_env = 0.5/day
    sigma = 0.05/day
    virus_growth = 0.0001/day
    virus_decay = 0.07/day
    param = SIRGrowth{typeof(unit(beta_force))}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
    param = transition(param)

    area = 10.0km^2
    epienv = simplehabitatAE(298.0K, size(initial_pop), area, epienv_active, NoControl())

    # Zero susceptible so we can test the specified initial_pop
    abun_h = (
        Susceptible = 0,
        Infected = 1,
        Recovered = 0,
        Dead = 0
    )
    disease_classes = (
        susceptible = ["Susceptible"],
        infectious = ["Infected"]
    )
    abun_v = (Environment = 0, Force = 0)

    dispersal_dists = fill(2.0km, size(initial_pop, 1) * size(initial_pop, 2))
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = EpiMovement(kernel)

    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = EpiList(traits, abun_v, abun_h, disease_classes, movement, param)

    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel, initial_pop)

    return epi
end


function TestCache()
    numSpecies = 3
    Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/test/Testdata/testTempBin.jld",
     "Temperature")
    energy_vec = SolarRequirement(fill(0.2*day^-1*kJ*m^-2, numSpecies))


    birth = 0.6/month
    death = 0.6/month
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month
    param = EqualPop(birth, death, long, surv, boost)

    file = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/World.tif"
    world = extractfile(file)
    europe = world[-10° .. 60°, 35° .. 80°]
    eu = ustrip.(europe)
    active = Array{Bool, 2}(.!isnan.(eu))

    grid = (94, 60)
    individuals=1000000
    era = ClimatePref.TestERA()
    wc = ClimatePref.TestWorldclim()
    wc = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3}, wc.array[-10° .. 60°, 35° .. 80°,:])
    bud = SolarTimeBudget(wc, 1)
    #active = Array(bud.matrix[:,:,1] .> 0*day^-1*kJ*m^-2)
    kernel = GaussianKernel.(fill(1.0km, numSpecies), 10e-04)
    movement = BirthOnlyMovement(kernel)
    common_species = ["Trifolium repens", "Urtica dioica", "Achillea millefolium"]
    native = fill(true, numSpecies)
    traits = Array(transpose(Temp[common_species,:]))
    traits = TempBin(traits)
    abun = rand(Multinomial(individuals, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
    abenv = eraAE(era, bud, active)
    rel = Trapeze{eltype(abenv.habitat)}()
    eco = Ecosystem(populate!, sppl, abenv, rel)
    times = 1month:1month:10years
    return eco
end
