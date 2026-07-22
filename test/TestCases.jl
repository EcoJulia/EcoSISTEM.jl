# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoSISTEM
using EcoSISTEM.ClimatePref
using AxisArrays
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Unitful
using Random

"""
    Test1Ecosystem(; seed = nothing)

Build a small test ecosystem. If `seed` is supplied, the initial per-species
abundance totals and the per-species simulation RNGs are both made deterministic,
so the whole run is reproducible regardless of the number of threads.
"""
function Test1Ecosystem(; seed = nothing)
    numSpecies = 150
    numNiches = 2

    birth = 0.6 / month
    death = 0.6 / month
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month
    param = EqualPop(birth, death, long, surv, boost)

    grid = (10, 10)
    area = 10000.0km^2
    individuals = 20000 * numSpecies
    totalK = 1000000.0 * kJ / km^2 * numSpecies
    habitat = simplenicheAE(numNiches, grid, totalK, area)

    # Seed the global RNG so the initial abundance totals are also deterministic
    isnothing(seed) || Random.seed!(seed)
    abun = rand(Multinomial(individuals, numSpecies))

    kernel = GaussianKernel.(fill(1.0km, numSpecies), 10e-04)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, numSpecies)
    resource = SolarDemand(fill(2.0kJ, numSpecies))
    sppl = SpeciesList(numSpecies, numNiches, abun, resource, movement, param,
                       native)

    rel = Match{eltype(habitat.regime)}()
    eco = isnothing(seed) ? Ecosystem(sppl, habitat, rel) :
          Ecosystem(sppl, habitat, rel; seed = seed)
    return eco
end

function TestMultiEcosystem()
    numSpecies = 150

    birth = 0.6 / month
    death = 0.6 / month
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month
    param = EqualPop(birth, death, long, surv, boost)

    grid = (10, 10)
    area = 10000.0km^2
    individuals = 20000 * numSpecies
    totalK1 = 1000000.0 * kJ / km^2 * numSpecies
    totalK2 = 100.0 * mm / km^2 * numSpecies
    # `axis = MeanTemperature` so the regime carries `TempChange` dynamics (the `TempIncrease!`/`TempFluct!`
    # scenarios reset its rate); dynamics now comes from the declared axis, not the value's K unit.
    habitat1 = simplehabitatAE(10.0K, grid, totalK1, area;
                               axis = MeanTemperature)
    habitat2 = simplehabitatAE(10.0K, grid, totalK2, area)
    supply = SupplyCollection2(habitat1.supply, habitat2.supply)
    habitat = GridHabitat{typeof(habitat1.regime), typeof(supply)}(habitat1.regime,
                                                                   habitat1.active,
                                                                   supply,
                                                                   habitat1.names)

    abun = rand(Multinomial(individuals, numSpecies))

    kernel = GaussianKernel.(fill(1.0km, numSpecies), 10e-04)
    movement = BirthOnlyMovement(kernel)
    native = fill(true, numSpecies)
    resource1 = SolarDemand(fill(2.0kJ, numSpecies))
    resource2 = WaterDemand(fill(2.0mm, numSpecies))
    resource = DemandCollection2(resource1, resource2)
    traits = Bin(MeanTemperature, Normal, fill(10.0K, numSpecies),
                 fill(0.1K, numSpecies))
    sppl = SpeciesList(numSpecies, traits, abun, resource, movement, param,
                       native)

    rel = DistRel{eltype(habitat.regime)}()
    eco = Ecosystem(sppl, habitat, rel)
    return eco
end
