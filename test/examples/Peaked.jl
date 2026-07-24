# SPDX-License-Identifier: LGPL-3.0-or-later

include("CreateEco.jl")
include("RunSim.jl")
include("Scenarios.jl")

## Peaked regime, 10_000 species

numInvasive = 1;
numSpecies = 10_000;
grd = (50, 50);
demand = (450000.0kJ / m^2 / day, 192.0Unitful.L / m^2 / day);
individuals = 10^10;
area = 200_000.0 * km^2;
totalK = (4.5e11kJ / km^2 / day, 192.0mm / day);

habloss = 1.0 / 10year
scenario = [
    SimpleScenario(TempIncrease!, 0.0K / 10year),
    SimpleScenario(TempIncrease!, 0.0K / 10year),
    SimpleScenario(GeneralistInvasive, 200.0 / year),
    SimpleScenario(SpecialistInvasive, 200.0 / year),
    SimpleScenario(TempIncrease!, 1.0K / 10year),
    SimpleScenario(TempIncrease!, 2.0K / 10year),
    SimpleScenario(TempIncrease!, 3.0K / 10year),
    SimpleScenario(RandHabitatLoss!, habloss),
    SimpleScenario(ClustHabitatLoss!, habloss)
]
scenario_names = [
    "Neutral1",
    "Neutral2",
    "GeneralistInvasive",
    "SpecialistInvasive",
    "TempIncrease1",
    "TempIncrease2",
    "TempIncrease3",
    "RandHabitatLoss",
    "ClustHabitatLoss"
]

divfuns = [
    sorenson,
    meta_speciesrichness,
    meta_shannon,
    meta_simpson,
    mean_abun,
    geom_mean_abun,
    norm_meta_alpha,
    raw_meta_alpha,
    norm_meta_beta,
    raw_meta_beta,
    norm_meta_rho,
    raw_meta_rho,
    meta_gamma
]
q = 1.0

simDict = Dict("times" => 50years,
               "burnin" => 10year,
               "interval" => 1month,
               "timestep" => 1month,
               "scenarios" => scenario,
               "divfuns" => divfuns,
               "q" => q,
               "reps" => 55,
               "cacheInterval" => 1year,
               "scenario_names" => scenario_names)
lensim = length((0month):simDict["interval"]:simDict["times"])

for i in 1:simDict["reps"]
    habitat1 = peakedgradhabitat(293.0K, 303.0K, grd, totalK[1], area,
                                 0.0K / month)
    habitat2 = peakedgradhabitat(293.0K, 303.0K, grd, totalK[2], area,
                                 0.0K / month)
    supply = SupplyCollection2(habitat1.supply, habitat2.supply)
    habitat = GridHabitat{typeof(habitat1.regime),
                          typeof(supply)}(habitat1.regime,
                                          habitat1.active,
                                          supply,
                                          habitat1.names)

    vars = rand(Uniform(0.75, 2.5), numSpecies + numInvasive) .* K
    opts = rand(Uniform(293.0, 303.0), numSpecies + numInvasive) .* K

    av_dist = rand(Uniform(0.6, 2.4), numSpecies + numInvasive) .* km
    kernel = GaussianKernel.(av_dist, 10e-10)

    death_rates = abs.(rand(Normal(0.15, 0.135), numSpecies + numInvasive)) ./
                  year
    birth_rates = death_rates

    paramDict = Dict("numSpecies" => numSpecies,
                     "numInvasive" => numInvasive,
                     "numIndiv" => individuals,
                     "demand" => demand,
                     "opts" => opts,
                     "vars" => vars,
                     "birth" => birth_rates,
                     "death" => death_rates,
                     "s" => 1e-3,
                     "boost" => 1.0,
                     "kernel" => kernel,
                     "totalK" => totalK,
                     "bound" => Torus())

    diver = zeros(length(divfuns), lensim, length(scenario))
    runsim!(diver, habitat, paramDict, simDict, i,
            "/media/storage/Chapter3/Continent/Peaked/")
end
