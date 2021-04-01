
using Simulation
using Simulation.ClimatePref
using Simulation.Units
using Unitful
using Unitful.DefaultSymbols
using JLD
using Printf

file = "Documents/Chapter5/data/Africa.tif"
africa = readfile(file, -25°, 50°, -35°, 40°)
active =  Array{Bool, 2}(.!isnan.(africa'))
# Set up initial parameters for ecosystem
numSpecies = 50_000; grid = size(africa); req= 10.0kJ; individuals=3*10^8; area = 64e6km^2; totalK = 1000.0kJ/km^2

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(req, numSpecies))


# Set rates for birth and death
birth = 0.6/year
death = 0.6/year
longevity = 1.0
survival = 0.1
boost = 1.0
# Collect model parameters together
param = EqualPop(birth, death, longevity, survival, boost)

# Create kernel for movement
kernel = fill(GaussianKernel(15.0km, 10e-10), numSpecies)
movement = AlwaysMovement(kernel, Torus())


# Create species list, including their temperature preferences, seed abundance and native status
opts = fill(274.0K, numSpecies)
vars = fill(50.0K, numSpecies)
traits = GaussTrait(opts, vars)
native = fill(true, numSpecies)
# abun = rand(Multinomial(individuals, numSpecies))
abun = fill(div(individuals, numSpecies), numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
    movement, param, native)
sppl.params.birth


dir1 = "Documents/CERA"
temp = readCERA(dir1, "cera_20c_temp2m", "t2m")
africa_temp = ERA(temp[-25°.. 50°, -35°.. 40°, :])
# Create abiotic environment - even grid of one temperature
abenv = simplehabitatAE(274.0K, grid, totalK, area, active)
abenv = eraAE(africa_temp, africa_solar, active)

# Set relationship between species and environment (gaussian)
rel = Gauss{typeof(1.0K)}()

#Create ecosystem
eco = Ecosystem(sppl, abenv, rel)

import Simulation.simulate!
function simulate!(eco::Ecosystem, times::Unitful.Time, timestep::Unitful.Time, cacheInterval::Unitful.Time, cacheFolder::String, scenario_name::String)
  time_seq = 0s:timestep:times
  counting = 0
  for i in 1:length(time_seq)
      update!(eco, timestep);
      # Save cache of abundances
      if mod(time_seq[i], cacheInterval) == 0year
          JLD.save(joinpath(cacheFolder, scenario_name * (@sprintf "%02d.jld" uconvert(NoUnits,time_seq[i]/cacheInterval))), "abun", eco.abundances.matrix)
      end
  end
end

# Simulation Parameters
burnin = 10years; times = 100years; timestep = 1month; record_interval = 12months;
lensim = length(0years:record_interval:times)
@time simulate!(eco, burnin, timestep)
@time simulate!(eco, times, timestep, record_interval, "sdc/Africa", "Africa_run_coexist");
