include("CreateEco.jl")
include("RunSim.jl")
include("Scenarios.jl")

## Gradient habitat, 1000 species

numInvasive = 1; numSpecies = 1000; grd = (40,25); req=(450000.0kJ/m^2, 192.0nm/m^2); individuals = 10^10; area = 10_000.0*km^2; totalK = (4.5e11kJ/km^2, 192.0mm/km^2)

habloss = 1.0 /10year
scenario = [SimpleScenario(TempIncrease!, 0.0K/10year),
    SimpleScenario(TempIncrease!, 0.0K/10year),
    SimpleScenario(GeneralistInvasive, 200.0/year),
    SimpleScenario(SpecialistInvasive, 200.0/year),
    SimpleScenario(TempIncrease!, 1.0K/10year),
    SimpleScenario(TempIncrease!, 2.0K/10year),
    SimpleScenario(TempIncrease!, 3.0K/10year),
    SimpleScenario(RandHabitatLoss!, habloss),
    SimpleScenario(ClustHabitatLoss!, habloss)]

divfuns = [sorenson, meta_speciesrichness, meta_shannon, meta_simpson, mean_abun, geom_mean_abun, norm_meta_alpha, raw_meta_alpha, norm_meta_beta, raw_meta_beta, norm_meta_rho, raw_meta_rho, meta_gamma]
q = 1.0

simDict = Dict("times" => 50years, "burnin" => 10year, "interval" => 1month, "timestep" => 1month, "scenarios" => scenario, "divfuns" => divfuns, "q" => q, "reps" => 100)
lensim = length(0month:simDict["interval"]:simDict["times"])

for i in 1:simDict["reps"]
    abenv1 = tempgradAE(293.0K, 303.0K, grd, totalK[1], area, 0.0K/month)
    abenv2 = tempgradAE(293.0K, 303.0K, grd, totalK[2], area, 0.0K/month)
    bud = BudgetCollection2(abenv1.budget, abenv2.budget)
    abenv = GridAbioticEnv{typeof(abenv1.habitat), typeof(bud)}(abenv1.habitat, abenv1.active, bud, abenv1.names)

    vars = rand(Uniform(0.75, 2.5), numSpecies + numInvasive) .* K
    opts = rand(Uniform(293.0, 303.0), numSpecies + numInvasive) .* K

    av_dist = rand(Uniform(0.6, 2.4), numSpecies + numInvasive) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies + numInvasive)) ./year
    birth_rates = death_rates

    paramDict = Dict("numSpecies" => numSpecies, "numInvasive" => numInvasive, "numIndiv" => individuals, "reqs" => req, "opts" => opts, "vars" => vars, "birth" => birth_rates, "death" => death_rates, "s" => 1e-3, "boost" => 1.0, "kernel" => kernel, "totalK" => totalK, "bound" => Cylinder())

    diver = zeros(length(divfuns), lensim, length(scenario))
    runsim!(diver, abenv, paramDict, simDict, i, "/home/claireh/Documents/Chapter3/Region/Grad/")
end
