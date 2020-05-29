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
function TestEpiSystem()
    numvirus = 1
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

    abun_h = [1000, 1, 0, 0]
    abun_v = [10]

    dispersal_dists = [fill(2.0km, numclasses - 1); 1e-2km]
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = SIR(traits, abun_v, abun_h, movement, param)

    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

    return epi
end
function TestEpiSystemFromPopulation(initial_pop::Matrix{<:Real})
    numclasses = 4
    numvirus = 1
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
    epienv = simplehabitatAE(298.0K, area, NoControl(), initial_pop)

    # Zero susceptible so we can test the specified initial_pop
    abun_h = [0, 1, 0, 0]
    abun_v = [10]

    dispersal_dists = [fill(2.0km, numclasses - 1); 1e-2km]
    kernel = GaussianKernel.(dispersal_dists, 1e-10)
    movement = AlwaysMovement(kernel)

    traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
    epilist = SIR(traits, abun_v, abun_h, movement, param)

    rel = Gauss{eltype(epienv.habitat)}()
    epi = EpiSystem(epilist, epienv, rel)

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
