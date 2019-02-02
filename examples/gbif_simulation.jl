using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Temperature.jld",
 "Temperature")
## Run simulation over a grid and plot
numSpecies=1

# Set up how much energy each species consumes
energy_vec = SolarRequirement(fill(2000.0*kJ, numSpecies))

# Set probabilities
birth = 0.6/month
death = 0.6/month
long = 1.0
surv = 0.5
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, long, surv, boost)

grid = (94, 60)
area = 10000.0km^2
totalK = 1000000.0
individuals=1000

# Load data for land cover
file = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/World.tif"
world = extractfile(file)
europe = world[-10° .. 60°, 35° .. 80°]
eu = ustrip.(europe)
@rput eu
R"image(t(eu))"
dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1990"
tempax2 = extractERA(dir2, "t2m", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2000"
tempax3 = extractERA(dir3, "t2m", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2010"
tempax4 = extractERA(dir4, "t2m", collect(361month:1month:38year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3}, srad.array[-10° .. 60°, 35° .. 80°,:])
srad = uconvert.(kJ, srad * month * (area/(grid[1]*grid[2])))
srad = SolarTimeBudget(srad, 1)
active = Array{Bool, 2}(.!isnan.(eu))

testtemp = tempax1
testtemp.array = tempax1.array[-10° .. 60°, 35° .. 80°, :]
# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

native = repmat([true], numSpecies)
traits = Array(transpose(Temp[1,:]))
traits = TempBin(traits)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, traits, rand(abun), energy_vec,
    movement, param, native)
abenv = eraAE(testtemp, srad, active)
#abenv = tempgradAE(-10.0°C, 30.0°C,  grid, totalK, area, 0.0°C/month, active)
#abenv.habitat.matrix = transpose(tempax1.array[-10° .. 60°, 35° .. 80°, 1])
rel = Trapeze{eltype(abenv.habitat)}()
eco = Ecosystem(sppl, abenv, rel)

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, timestep)
    eco.abenv.habitat.time = 1
    #resetrate!(eco, 0.1°C/month)
    simulate_record!(abun, eco, times, interval, timestep)
end
times = 10year
abun = runsim(eco, times)
plot_abun(abun, numSpecies, grid, 1)
plot_mean(abun, numSpecies, grid)

x = ustrip.(europe.axes[1].val)
y = ustrip.(europe.axes[2].val)
@rput x; @rput y
abunmat = reshape(abun[1,:,:,1], (94, 60, 121))
@rput abunmat
hab = ustrip.(eco.abenv.habitat.matrix)
@rput hab
bud = ustrip.(eco.abenv.budget.matrix)
@rput bud
R"par(mfrow=c(1,2));library(fields); library(viridis)
image.plot(t(hab[,,120]), main = 'Habitat', col=viridis(50))
image.plot(t(abunmat[,,120]), main = 'Abundances', col=magma(50))
   "
  R"pdf(file='plots/gbif_sim_era_sol.pdf', paper = 'a4r', height= 8.27, width=11.69 )
  library(rgdal)
  library(fields)
  months = c('January', 'February', 'March', 'April', 'May', 'June', 'July',
  'August', 'September', 'October', 'November', 'December')
  world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
  for (i in c(1:12)){
     par(mfrow=c(1,2))
     image.plot(x,y, hab[,,i],  col=viridis(50), main = months[i],
     zlim = c(-30, 40))
     plot(world,  add = T)
     image.plot(x, y, abunmat[,,(i+1)], main = 'Abundances', col=magma(50),
     zlim = c(0, 25))
     plot(world, add = T, border = 'white')}
     dev.off()"
 R"library(viridis)
 library(rgdal)
 library(fields)
 months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
 'August', 'September', 'October', 'November', 'December'), 10)
 years = c(1980:1989)
 world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
 for (i in c(1:120)){
    jpeg(paste('plots/gbif_sim/gbif_simulation_sol', i,'.jpg', sep=''), quality=100, width =1000)
    par(mfrow=c(1,3))
    image.plot(x,y, hab[,,i],  col=viridis(50), main = months[i],
    zlim = c(-30, 40), border = 'white')
    plot(world,  add = T)
    image.plot(x,y, bud[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = months[i],
    zlim = c(0, 30000))
    plot(world,  add = T, border = 'white')
    image.plot(x, y, abunmat[,,(i+1)], main = years[i%/%13 +1], col=magma(50),
    zlim = c(0, 2000))
    plot(world, add = T, border = 'white')
    dev.off()}"

 R"pdf(file='plots/gbif_sim.pdf', paper = 'a4r', height= 8.27, width=11.69 )
library(rgdal)
library(fields)
par(mfrow=c(1,2))
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
image.plot(x,y, t(hab), main = 'Habitat', col=viridis(50))
plot(world,  add = T)
image.plot(x, y, t(abunmat), main = 'Abundances', col=magma(50))
plot(world, add = T, border = 'white')
dev.off()"

############################################################
### COMBINE TEMPERATURE AND RAINFALL
############################################################
using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Temperature.jld",
 "Temperature")
 Rain = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Rainfall.jld",
  "Rainfall")

## Run simulation over a grid and plot
numSpecies=1

# Set up how much energy each species consumes
energy_vec = SolarRequirement(repmat([2*day^-1*kJ*m^-2], numSpecies))

# Set probabilities
birth = 0.6/month
death = 0.6/month
long = 1.0
surv = 0.5
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, long, surv, boost)

grid = (94, 60)
area = 10000.0km^2
totalK = 1000000.0
individuals=1000

# Load data for land cover
file = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/World.tif"
world = extractfile(file)
europe = world[-10° .. 60°, 35° .. 80°]
eu = ustrip.(europe)
@rput eu
R"image(t(eu))"
dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1990"
tempax2 = extractERA(dir2, "t2m", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2000"
tempax3 = extractERA(dir3, "t2m", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2010"
tempax4 = extractERA(dir4, "t2m", collect(361month:1month:38year))

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_1980"
precax1 = extractERA(dir1, "tp", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_1990"
precax2 = extractERA(dir2, "tp", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_2000"
precax3 = extractERA(dir3, "tp", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_2010"
precax4 = extractERA(dir4, "tp", collect(361month:1month:38year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
prec = extractworldclim(joinpath(dir, "wc2.0_5m_prec"))
prec.array = prec.array[-10° .. 60°, 35° .. 80°,:]
x = prec.array.axes[1]
y = prec.array.axes[2]
t = prec.array.axes[3]
prec.array = AxisArray(1.0.*(prec.array),
    Axis{:longitude}(x), Axis{:latitude}(y), Axis{:time}(t))


dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3},
            srad.array[-10° .. 60°, 35° .. 80°,:])
srad = SolarTimeBudget(srad, 1)
active = Array{Bool, 2}(.!isnan.(eu))

testtemp = tempax1
testtemp.array = tempax1.array[-10° .. 60°, 35° .. 80°, :]

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

native = repmat([true], numSpecies)
traits1 = Array(transpose(Temp["Trifolium repens",:]))
traits1 = TempBin(traits1)
traits2 = Array(transpose(Rain["Trifolium repens",:]))
traits2 = RainBin(traits2)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, TraitCollection2(traits1, traits2), abun, energy_vec,
movement, param, native)
abenv1 = eraAE(testtemp, srad, active)
abenv2 = worldclimAE(prec, srad, active)
hab = HabitatCollection2(abenv1.habitat, abenv2.habitat)
abenv = GridAbioticEnv{typeof(hab), typeof(abenv1.budget)}(hab, abenv1.active,
    abenv1.budget, abenv1.names)
rel1 = Trapeze{eltype(abenv.habitat.h1)}()
rel2 = Unif{eltype(abenv.habitat.h2)}()
rel = additiveTR2(rel1, rel2)
eco = Ecosystem(sppl, abenv, rel)

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, timestep)
    eco.abenv.habitat.h1.time = 1
    eco.abenv.habitat.h2.time = 1
    #resetrate!(eco, 0.1°C/month)
    simulate_record!(abun, eco, times, interval, timestep)
end
times = 10year
abun = runsim(eco, times)
plot_abun(abun, numSpecies, grid, 1)
plot_mean(abun, numSpecies, grid)

x = ustrip.(europe.axes[1].val)
y = ustrip.(europe.axes[2].val)
@rput x; @rput y
abunmat = reshape(abun[1,:,:,1], (94, 60, 121))
@rput abunmat
hab1 = ustrip.(eco.abenv.habitat.h1.matrix)
@rput hab1
hab2 = ustrip.(eco.abenv.habitat.h2.matrix)
@rput hab2
bud = ustrip.(eco.abenv.budget.matrix)
@rput bud
R"par(mfrow=c(1,2));library(fields); library(viridis)
image.plot(t(hab[,,120]), main = 'Habitat', col=viridis(50))
image.plot(t(abunmat), main = 'Abundances', col=magma(50))
   "
  R"pdf(file='plots/gbif_sim_era_rain.pdf', paper = 'a4r', height= 8.27, width=11.69 )
  library(rgdal)
  library(fields)
  months = c('January', 'February', 'March', 'April', 'May', 'June', 'July',
  'August', 'September', 'October', 'November', 'December')
  world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
  for (i in c(1:120)){
     par(mfrow=c(1,2))
     image.plot(x,y, hab[,,i],  col=viridis(50), main = months[i],
     zlim = c(-30, 40))
     plot(world,  add = T)
     image.plot(x, y, abunmat[,,(i+1)], main = 'Abundances', col=magma(50),
     zlim = c(0, 25))
     plot(world, add = T, border = 'white')}
     dev.off()"
 R"library(viridis)
 library(rgdal)
 library(fields)
 months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
 'August', 'September', 'October', 'November', 'December'), 10)
 years = c(1980:1989)
 world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
 for (i in c(1:120)){
    jpeg(paste('plots/gbif_sim/gbif_simulation_rain', i,'.jpg', sep=''),
    quality=100, width =1000, height =1000)
    par(mfrow=c(2,2))
    image.plot(x,y, hab1[,,i],  col=viridis(50), main = months[i],
    zlim = c(-30, 40), border = 'white')
    plot(world,  add = T)
    image.plot(x,y, hab2[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = months[i],
    zlim = c(0, 400), border = 'white')
    plot(world,  add = T)
    image.plot(x,y, bud[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = months[i],
    zlim = c(0, 30000))
    plot(world,  add = T, border = 'white')
    image.plot(x, y, abunmat[,,(i+1)], main = years[i%/%13 +1], col=magma(50),
    zlim = c(0, 500))
    plot(world, add = T, border = 'white')
    dev.off()}"


############################################################
### SEVERAL COMMON SPECIES
############################################################

using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Temperature.jld",
 "Temperature")
 Rain = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Rainfall.jld",
  "Rainfall")

## Run simulation over a grid and plot
numSpecies=3

# Set up how much energy each species consumes
energy_vec = SolarRequirement(repmat([2*day^-1*kJ*m^-2], numSpecies))

# Set probabilities
birth = 0.6/month
death = 0.6/month
long = 1.0
surv = 0.5
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, long, surv, boost)

grid = (94, 60)
area = 10000.0km^2
totalK = 1000000.0
individuals=1000

# Load data for land cover
file = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/World.tif"
world = extractfile(file)
europe = world[-10° .. 60°, 35° .. 80°]
eu = ustrip.(europe)

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1990"
tempax2 = extractERA(dir2, "t2m", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2000"
tempax3 = extractERA(dir3, "t2m", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_2010"
tempax4 = extractERA(dir4, "t2m", collect(361month:1month:38year))

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_1980"
precax1 = extractERA(dir1, "tp", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_1990"
precax2 = extractERA(dir2, "tp", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_2000"
precax3 = extractERA(dir3, "tp", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_mdfa_prec_2010"
precax4 = extractERA(dir4, "tp", collect(361month:1month:38year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
prec = extractworldclim(joinpath(dir, "wc2.0_5m_prec"))
prec.array = prec.array[-10° .. 60°, 35° .. 80°,:]
x = prec.array.axes[1]
y = prec.array.axes[2]
t = prec.array.axes[3]
prec.array = AxisArray(1.0.*(prec.array),
    Axis{:longitude}(x), Axis{:latitude}(y), Axis{:time}(t))


dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3},
            srad.array[-10° .. 60°, 35° .. 80°,:])
srad = SolarTimeBudget(srad, 1)
active = Array{Bool, 2}(.!isnan.(eu))

testtemp = tempax1
testtemp.array = tempax1.array[-10° .. 60°, 35° .. 80°, :]

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

common_species = ["Trifolium repens", "Urtica dioica", "Achillea millefolium"]
native = repmat([true], numSpecies)
traits1 = Array(transpose(Temp[common_species,:]))
traits1 = TempBin(traits1)
traits2 = Array(transpose(Rain[common_species,:]))
traits2 = RainBin(traits2)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, TraitCollection2(traits1, traits2), abun, energy_vec,
movement, param, native)
abenv1 = eraAE(testtemp, srad, active)
abenv2 = worldclimAE(prec, srad, active)
hab = HabitatCollection2(abenv1.habitat, abenv2.habitat)
abenv = GridAbioticEnv{typeof(hab), typeof(abenv1.budget)}(hab, abenv1.active,
    abenv1.budget, abenv1.names)
rel1 = Trapeze{eltype(abenv.habitat.h1)}()
rel2 = Unif{eltype(abenv.habitat.h2)}()
rel = additiveTR2(rel1, rel2)
eco = Ecosystem(sppl, abenv, rel)

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, timestep)
    eco.abenv.habitat.h1.time = 1
    eco.abenv.habitat.h2.time = 1
    #resetrate!(eco, 0.1°C/month)
    simulate_record!(abun, eco, times, interval, timestep)
end
times = 10year
abun = runsim(eco, times)

divtimes = collect(1:1:120)
alphas = zeros(94*60, 11, length(divtimes))
for i in 1:length(divtimes)
  met = Metacommunity(abun[:,:,divtimes[i],1], UniqueTypes(numSpecies))
  alphas[:,:,i] = norm_sub_alpha(met, 0:10)[:diversity]
end
alphas = reshape(alphas, 94, 60, 11, length(divtimes))

hab1 = ustrip.(eco.abenv.habitat.h1.matrix)
@rput hab1
hab2 = ustrip.(eco.abenv.habitat.h2.matrix)
@rput hab2
bud = ustrip.(eco.abenv.budget.matrix)
@rput bud
x = ustrip.(europe.axes[1].val)
y = ustrip.(europe.axes[2].val)
@rput x; @rput y
@rput alphas; @rput hab1; @rput hab2; @rput bud
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
for (i in c(1:120)){
jpeg(paste('plots/gbif_sim/gbif_simulation_multiple', i, '.jpg'), quality=1000,
    width =1000, height =1000)
im = alphas[ , , 2, i]
im[is.na(im)] = 0
par(mfrow=c(2,2))
image.plot(x,y, hab1[,,i],  col=viridis(50), main = months[i],
zlim = c(-30, 40))
plot(world,  add = T, border = 'white')
image.plot(x,y, hab2[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = years[i%/%13 +1],
zlim = c(0, 400))
plot(world,  add = T, border = 'white')
image.plot(x,y, bud[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = '',
zlim = c(0, 30000))
plot(world,  add = T, border = 'white')
image.plot(x, y, im, col=magma(30),zlim=c(0,3),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
}
"

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, timestep)
    eco.abenv.habitat.h1.time = 1
    eco.abenv.habitat.h2.time = 1
    #resetrate!(eco, 0.1°C/month)
    simulate_record_diversity!(abun, eco, times, interval, timestep,
        norm_sub_alpha, collect(1.0:10))
end
times = 10year
abun = runsim(eco, times)
############################################################
### WATER AS SECOND BUDGET
############################################################

using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Temperature.jld",
 "Temperature")
 Rain = load("/Users/claireh/Documents/PhD/GIT/ClimatePref/data/Rainfall.jld",
  "Rainfall")

## Run simulation over a grid and plot
numSpecies=3

# Set up how much energy each species consumes
energy_vec1 = SolarRequirement(repmat([2*day^-1*kJ*m^-2], numSpecies))
energy_vec2 = WaterRequirement(repmat([0.5*mm], numSpecies))
energy_vec = ReqCollection2(energy_vec1, energy_vec2)

# Set probabilities
birth = 0.6/month
death = 0.6/month
long = 1.0
surv = 0.5
boost = 1000.0
timestep = 1.0month

# Collect model parameters together (in this order!!)
param = EqualPop(birth, death, long, surv, boost)

grid = (94, 60)
area = 10000.0km^2
totalK = 1000000.0
individuals=1000

# Load data for land cover
file = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/World.tif"
world = extractfile(file)
europe = world[-10° .. 60°, 35° .. 80°]
eu = ustrip.(europe)

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
prec = extractworldclim(joinpath(dir, "wc2.0_5m_prec"))
prec.array = prec.array[-10° .. 60°, 35° .. 80°,:]
x = prec.array.axes[1]
y = prec.array.axes[2]
t = prec.array.axes[3]
prec.array = AxisArray(1.0.*(prec.array),
    Axis{:longitude}(x), Axis{:latitude}(y), Axis{:time}(t))
water = WaterTimeBudget(Array{typeof(1.0*mm), 3}(prec.array), 1)


dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3},
            srad.array[-10° .. 60°, 35° .. 80°,:])
srad = SolarTimeBudget(srad, 1)
active = Array{Bool, 2}(.!isnan.(eu))
bud = BudgetCollection2(srad, water)

testtemp = tempax1
testtemp.array = tempax1.array[-10° .. 60°, 35° .. 80°, :]

# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

common_species = ["Trifolium repens", "Urtica dioica", "Achillea millefolium"]
native = repmat([true], numSpecies)
traits1 = Array(transpose(Temp[common_species,:]))
traits1 = TempBin(traits1)
traits2 = Array(transpose(Rain[common_species,:]))
traits2 = RainBin(traits2)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, TraitCollection2(traits1, traits2), abun, energy_vec,
    movement, param, native)
abenv1 = eraAE(testtemp, srad, active)
abenv2 = worldclimAE(prec, srad, active)
hab = HabitatCollection2(abenv1.habitat, abenv2.habitat)
abenv = GridAbioticEnv{typeof(hab), typeof(bud)}(hab, abenv1.active,
    bud, abenv1.names)
rel1 = Trapeze{eltype(abenv.habitat.h1)}()
rel2 = Unif{eltype(abenv.habitat.h2)}()
rel = multiplicativeTR2(rel1, rel2)
eco = Ecosystem(sppl, abenv, rel)

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, lensim, 1)
    simulate!(eco, burnin, timestep)
    eco.abenv.habitat.h1.time = 1
    eco.abenv.habitat.h2.time = 1
    #resetrate!(eco, 0.1°C/month)
    simulate_record!(abun, eco, times, interval, timestep)
end

function runsim(eco::Ecosystem, times::Unitful.Time)
    burnin = 1year; interval = 1month
    lensim = length(0month:interval:times)
    abun = generate_storage(eco, 1, lensim, 1)
    simulate!(eco, burnin, timestep)
    eco.abenv.habitat.h1.time = 1
    eco.abenv.habitat.h2.time = 1
    #resetrate!(eco, 0.1°C/month)
    simulate_record_diversity!(abun, eco, times, interval, timestep,
        norm_sub_alpha, 1.0)
end

times = 1year
abun = runsim(eco, times)

hab1 = ustrip.(eco.abenv.habitat.h1.matrix)
@rput hab1
hab2 = ustrip.(eco.abenv.habitat.h2.matrix)
@rput hab2
bud1 = ustrip.(eco.abenv.budget.b1.matrix)
@rput bud1
bud2 = ustrip.(eco.abenv.budget.b2.matrix)
@rput bud2
x = ustrip.(europe.axes[1].val)
y = ustrip.(europe.axes[2].val)
@rput x; @rput y
@rput abun; @rput hab1; @rput hab2; @rput bud1; @rput bud2
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
for (i in c(1:120)){
jpeg(paste('plots/gbif_sim/gbif_simulation_water', i, '.jpg'), quality=1000,
    width =1000, height =1000)
im = alphas[ , , 2, i]
im[is.na(im)] = 0
par(mfrow=c(2,2))
image.plot(x,y, hab1[,,i],  col=viridis(50), main = months[i],
zlim = c(-30, 40))
plot(world,  add = T, border = 'white')
image.plot(x,y, hab2[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = years[i%/%13 +1],
zlim = c(0, 400))
plot(world,  add = T, border = 'white')
image.plot(x,y, bud1[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = '',
zlim = c(0, 30000))
plot(world,  add = T, border = 'white')
image.plot(x, y, im, col=magma(30),zlim=c(0,3),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
}
"


############################################################
### ENTIRE WORLD
############################################################

using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using MyUnitful
using AxisArrays
using ClimatePref

eco = readTOML("setupB.toml")
runTOML("runB.toml", eco)
eco = readTOML("examples/TOML/NH_sw_tp_3/setup.toml")

abun = mapslices(sum, eco.abundances.grid, 1)[1, :, :]
y = ustrip.(-0.375°:0.375°:90°)
x = ustrip.(-20°:0.375°:180°)
@rput abun; @rput x; @rput y
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
abun[is.na(abun)] = 0
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
jpeg(paste('examples/TOML/NH_sw_tp_3/plots/simulation_start.jpg'), quality=1000,
    width =1500, height =1000)
image.plot(x, y, log(abun+1), col=magma(30),zlim=c(0,15),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
"

R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(raster)
abun[is.na(abun)] = 0
world = shapefile('../../gbif/ne_10m_land/ne_10m_land.shp')
jpeg(paste('plots/simulation_start.jpg'), quality=1000,
    width =1500, height =1000)
image.plot(x, y, log(abun+1), col=magma(30),zlim=c(0,15),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
"

runTOML("examples/TOML/NH_sw_tp_3/run.toml", eco)
abun = readoutput("examples/TOML/NH_sw_tp_3/data/", "Run")
abun = reshape(abun, 534, 242, 1, 1285, 1)[:,:,1,:,1]


hab1 = ustrip.(eco.abenv.habitat.h1.matrix)
@rput hab1
hab2 = ustrip.(eco.abenv.habitat.h2.matrix)
@rput hab2
bud1 = ustrip.(eco.abenv.budget.b1.matrix)
@rput bud1
bud2 = ustrip.(eco.abenv.budget.b2.matrix)
@rput bud2
x = -20.0:0.375:180
y = -0.375:0.375:90
@rput x; @rput y
@rput abun; @rput hab1; @rput hab2; @rput bud1; @rput bud2
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
for (i in seq(0, 100, by = 10)){
jpeg(paste('examples/TOML/NH_sw_tp_5/plots/simulation_world', i, '.jpg'), quality=1000,
    width =1500, height =1000)
im = abun[ , , i]
im[is.na(im)] = 0
par(mfrow=c(2,2))
image.plot(x,y, hab1[,,i],  col=viridis(50), main = months[i],
zlim = c(-80, 45))
plot(world,  add = T, border = 'white')
image.plot(x,y, hab2[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = years[i%/%13 +1],
zlim = c(0, 400))
plot(world,  add = T, border = 'white')
image.plot(x,y, bud1[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = '',
zlim = c(0, 50000))
plot(world,  add = T, border = 'white')
image.plot(x, y, im, col=magma(30),zlim=c(0,5),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
}
"

abun = readoutput("examples/TOML/NH_sw_tp_5/data/", "Run")
abun = reshape(abun, 534, 242, 1, 1284, 1)[:,:,1,:,1]
@rput x; @rput y
@rput abun
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
casp = readOGR('/Users/claireh/Documents/worldglwd1/worldglwd1.shp', layer='worldglwd1')
for (i in seq(1, 1284, by = 1)){
png(paste0('examples/TOML/NH_sw_tp_5/plots/simulation_world', sprintf('%03d', i), '.png'),
    width =1069, height =484)
im = abun[ , , i]
im[is.na(im)] = 0
image.plot(x, y, im, col=magma(30),zlim=c(0,5),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
plot(casp, add = 2, col = 1, border = 'white')
dev.off()
}
"

abun2, beta = Simulation.readbeta("examples/TOML/NH_sw_tp_5/data/beta/", "Run")
abun2 = reshape(abun, 534, 242, 1, 169, 1)[:,:,1,:,1]
@rput abun2; @rput beta
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
casp = readOGR('/Users/claireh/Documents/worldglwd1/worldglwd1.shp', layer='worldglwd1')
for (i in seq(1, 169, by = 12)){
png(paste0('examples/TOML/NH_sw_tp_5/plots/simulation_beta', i, '.png'),
    width =1069, height =968)
par(mfrow = c(2, 1))
im = abun2[ , , i]
im[is.na(im)] = 0
image.plot(x, y, im, col=magma(30),zlim=c(0,20),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
plot(beta[1:i], xlim=c(0, 157),ylim=c(0, 4), type = 'l', xlab = 'Month',
    ylab = 'Number of distinct ecosystems')
points(i, beta[i], pch = 20, col = 'red')
dev.off()
}
"

abun = readoutput("examples/TOML/NH_sw_tp_5/data/", "Run")
abun = reshape(abun, 534, 242, 1, 1309, 1)[:,:,1,:,1]
@rput x; @rput y
@rput abun
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
for (i in seq(1, 1300, by = 10)){
jpeg(paste('examples/TOML/NH_sw_tp_5/plots/simulation_world', i, '.jpg'), quality=1000,
    width =1500, height =1000)
im = abun[ , , i]
im[is.na(im)] = 0
image.plot(x, y, im, col=magma(30),zlim=c(0,5),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
}
"

abun = readoutput("examples/TOML/NH_sw_tp_10/data/", "Run")
abun = reshape(abun, 534, 242, 1, 673, 1)[:,:,1,:,1]
@rput x; @rput y
@rput abun
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
for (i in seq(1, 600, by = 10)){
jpeg(paste('examples/TOML/NH_sw_tp_10/plots/simulation_world', i, '.jpg'), quality=1000,
    width =1500, height =1000)
im = abun[ , , i]
im[is.na(im)] = 0
image.plot(x, y, im, col=magma(30),zlim=c(0,10),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
}
"


abun = readoutput("examples/TOML/NH_sw_tp_50/data/", "Run")
abun = reshape(abun, 534, 242, 1, 649, 1)[:,:,1,:,1]
@rput x; @rput y
@rput abun
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
for (i in seq(1, 649, by = 12)){
png(paste0('examples/TOML/NH_sw_tp_50/plots/simulation_world', i, '.png'),
    width = 1068, height = 484)
im = abun[ , , i]
im[is.na(im)] = 0
image.plot(x, y, im, col=magma(30),zlim=c(0,20),
xlab='', ylab='', main ='')
plot(world,  add = T, border = 'white')
dev.off()
}
"
abun = readoutput("examples/TOML/NH_sw_tp_100/data/", "Run")
abun = reshape(abun, 534, 242, 1, 1081, 1)[:,:,1,:,1]
@rput x; @rput y
@rput abun
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
casp = readOGR('/Users/claireh/Documents/worldglwd1/worldglwd1.shp', layer='worldglwd1')
for (i in seq(1, 1081, by = 12)){
png(paste0('examples/TOML/NH_sw_tp_100/plots/simulation_world', i, '.png'),
    width = 1068, height = 484)
    im = abun[ , , i]
    im[is.na(im)] = 0
    image.plot(x, y, im, col=c('black', viridis(30)),zlim=c(0,100),
    xlab='', ylab='', main ='')
    plot(world,  add = T, border = 'white')
    plot(casp, add = T, col = 1, border = 'white')
dev.off()
}
"
abun, abun2 = Simulation.readall("examples/TOML/NH_sw_tp_100/data/all/", "Run")
abun = reshape(abun, 534, 242, 3, 1201, 1)[:,:,:,:,1]
@rput abun; @rput abun2
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
mains = c(' ','Number of distinct ecosystems', 'Number of remaining species')
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 100)
years = rep(c(1901:2010), each =12)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
casp = readOGR('/Users/claireh/Documents/worldglwd1/worldglwd1.shp', layer='worldglwd1')
for (i in seq(1200, 1200, by = 1)){
png(paste0('examples/TOML/NH_sw_tp_100/plots/simulation_alpha', sprintf('%04d',i+1), '.png'),
    width =2000, height =1000)
par(mfrow = c(2, 2), mar = c(5,4,4, 6))
image.plot(x,y, hab1[,,i],  col=magma(50), main = 'Temperature',
zlim = c(-80, 45), xlab='', ylab='', axes=F, cex.main = 4)
plot(world,  add = T, border = 'white')
image.plot(x,y, log(1 + hab2[,,rep(c(1:12), 100)[i]]),  col=viridis(50), main = 'Precipitation',
zlim = c(0, 10), xlab='', ylab='', axes=F, cex.main = 4)
plot(world,  add = T, border = 'white')
image.plot(x,y, bud1[,,rep(c(1:12), 100)[i]],  col=viridis(50), main = 'Solar radiation',
zlim = c(0, 50000), xlab='', ylab='', axes=F, cex.main = 4)
plot(world,  add = T, border = 'white')
im = abun[ , , 1, i]
im[is.na(im)] = 0
image.plot(x, y, im, col=c('black', viridis(100)[10:100]),zlim=c(0,100),
xlab='', ylab='',  axes=F, main = paste('December',  2010), cex.main = 4)
plot(world,  add = T, border = 'white')
plot(casp, add = T, col = 1, border = 'white')
dev.off()
}
"

R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
mains = c(' ','Number of distinct ecosystems', 'Number of remaining species')
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
casp = readOGR('/Users/claireh/Documents/worldglwd1/worldglwd1.shp', layer='worldglwd1')
for (i in seq(1082, 1201, by = 1)){
png(paste0('examples/TOML/NH_sw_tp_100/plots/simulation_beta', sprintf('%04d',i), '.png'),
    width =1069, height =968)
par(mfrow = c(2, 1), mar = c(5,6,4,4))
im = abun[ , , 2, i]
im[is.na(im)] = 0
image.plot(x, y, log(1+im), col=c('black',viridis(100)[10:100]),zlim=c(0,10),
xlab='', ylab='', main ='', axes =F)
plot(world,  add = T, border = 'white')
plot(casp, add = T, col = 1, border = 'white')
plot(abun2[2, 1:i], xlim = c(1, 1308), ylim = c(1, 20), type = 'l', xlab = 'Year',
    ylab = list(mains[2], cex=2), xaxt='n', cex.axis=1.5, cex.lab = 2)
axis(1, at=seq(1, 1308, 60), labels=seq(0, 109, 5), cex.axis = 1.5)
points(i, abun2[2, i], pch = 20, col = 'red')
dev.off()
}
"
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
mains = c(' ','Number of distinct ecosystems', 'Number of remaining species')
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
casp = readOGR('/Users/claireh/Documents/worldglwd1/worldglwd1.shp', layer='worldglwd1')
for (i in seq(1082, 1201, by = 1)){
png(paste0('examples/TOML/NH_sw_tp_100/plots/simulation_gamma', sprintf('%04d',i), '.png'),
    width =1069, height =968)
par(mfrow = c(2, 1), mar = c(5,6,4,4))
im = abun[ , , 3, i]
im[is.na(im)] = 0
image.plot(x, y, im, col=c('black', viridis(100)[10:100]), zlim=c(0,1000),
xlab='', ylab='', main ='', axes =F)
plot(world,  add = T, border = 'white')
plot(casp, add = T, col = 1, border = 'white')
plot(abun2[3, 1:i], xlim=c(0, 1308), ylim=c(0, 100), type = 'l', xlab = 'Year',
    ylab = list(mains[3], cex=2), xaxt='n', cex.axis=1.5, cex.lab = 2)
axis(1, at=seq(1, 1308, 60), labels=seq(0, 109, 5), cex.axis=1.5)
points(i, abun2[3, i], pch = 20, col = 'red')
dev.off()
}
"
hab1 = ustrip.(eco.abenv.habitat.h1.matrix)
@rput hab1
hab2 = ustrip.(eco.abenv.habitat.h2.matrix)
@rput hab2
bud1 = ustrip.(eco.abenv.budget.b1.matrix)
@rput bud1
bud2 = ustrip.(eco.abenv.budget.b2.matrix)
@rput bud2
x = -20.0:0.375:180
y = -0.375:0.375:90
abun = readoutput("examples/TOML/NH_sw_tp_1000/data/", "Run")
abun = reshape(abun, 534, 242, 1, 109, 1)[:,:,1,:,1]
@rput x; @rput y
@rput abun
R"library(fields);par(mfrow=c(1,1))
library(viridis)
library(rgdal)
months = rep(c('January', 'February', 'March', 'April', 'May', 'June', 'July',
'August', 'September', 'October', 'November', 'December'), 10)
years = c(1980:1989)
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
casp = readOGR('/Users/claireh/Documents/worldglwd1/worldglwd1.shp', layer='worldglwd1')
for (i in seq(1, 109, by = 12)){
png(paste0('examples/TOML/NH_sw_tp_1000/plots/simulation_world', i, '.png'),
    width = 1068, height = 484)
    par(mfrow=c(2,2))
    image.plot(x,y, hab1[,,i],  col=magma(50), main = 'Temperature',
    zlim = c(-80, 45))
    plot(world,  add = T, border = 'white')
    image.plot(x,y, log(1 + hab2[,,rep(c(1:12), 10)[i]]),  col=viridis(50), main = 'Precipitation',
    zlim = c(0, 10))
    plot(world,  add = T, border = 'white')
    image.plot(x,y, bud1[,,rep(c(1:12), 10)[i]],  col=viridis(50), main = 'Solar radiation',
    zlim = c(0, 50000))
im = abun[ , , i]
im[is.na(im)] = 0
image.plot(x, y, im, col=c('black', viridis(50)),zlim=c(0,1000),
xlab='', ylab='', main ='Alpha diversity')
plot(world,  add = T, border = 'white')
plot(casp, add = 2, col = 1, border = 'white')
dev.off()
}
"
