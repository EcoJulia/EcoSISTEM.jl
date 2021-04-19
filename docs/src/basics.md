# The basics of EcoSISTEM.jl

## Install

EcoSISTEM has yet to be officially released, but you can download the package through github:

```julia
]add https://github.com/boydorr/EcoSISTEM.jl.git
```

## Setting up an ecosystem

The package runs on an `Ecosystem` containing information about species, the `SpeciesList`, their environment, `AbioticEnvironment` and the relationship between the two, `TraitRelationship`.

### A simple example

Load required packages
```julia
using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Distributions
using Diversity
```
Set up initial parameters for ecosystem
```julia
numSpecies = 10; grid = (5, 5); req= 10.0kJ; individuals=1000; area = 1000.0*km^2; totalK = 1.0kJ/km^2
```

Set up how much energy each species consumes
```julias
energy_vec = SolarRequirement(fill(req, numSpecies))
```

Set rates for birth and death
```julia
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.2
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)
```

Create kernel for movement
```julia
kernel = fill(GaussianKernel(10.0km, 10e-10), numSpecies)
movement = BirthOnlyMovement(kernel, Torus())
```

Create species list, including their temperature preferences, seed abundance and native status
```julia
opts = fill(274.0K, numSpecies)
vars = fill(0.5K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
```

Create abiotic environment - even grid of one temperature
```julia
abenv = simplehabitatAE(274.0K, grid, totalK, area)
```

Set relationship between species and environment (gaussian)
```julia
rel = Gauss{typeof(1.0K)}()
```

Create ecosystem
```julia
eco = Ecosystem(sppl, abenv, rel)
```

Run simulation
```julia
# EcoSISTEM Parameters
burnin = 5years; times = 50years; timestep = 1month; record_interval = 3months; repeats = 1
lensim = length(0years:record_interval:times)
# Burnin
@time simulate!(eco, burnin, timestep)
```

Plot using SpatialEcology
```julia
using SpatialEcology
plot(eco)
```
