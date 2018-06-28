using Simulation
using JLD
using Diversity
using RCall
using Distributions
using Unitful
using Unitful.DefaultSymbols
using myunitful
using AxisArrays
using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/Temperature.jld",
 "Temperature")
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
file = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/World.tif"
world = extractfile(file)
europe = world[-10° .. 60°, 35° .. 80°]
eu = ustrip.(europe)
@rput eu
R"image(t(eu))"
dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_1990"
tempax2 = extractERA(dir2, "t2m", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_2000"
tempax3 = extractERA(dir3, "t2m", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_2010"
tempax4 = extractERA(dir4, "t2m", collect(361month:1month:38year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3}, srad.array[-10° .. 60°, 35° .. 80°,:])
srad = SolarBudget(srad, 1)
active = Array{Bool, 2}(.!isnan.(eu))

testtemp = tempax1
testtemp.array = tempax1.array[-10° .. 60°, 35° .. 80°, :]
# Create ecosystem
kernel = GaussianKernel(10.0km, numSpecies, 10e-4)
movement = BirthOnlyMovement(kernel)

native = repmat([true], numSpecies)
traits = Array(transpose(Temp[1,:]))
traits = TempBin{typeof(1.0°C)}(traits)
abun = Multinomial(individuals, numSpecies)
sppl = SpeciesList(numSpecies, traits, abun, energy_vec,
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
image.plot(t(abunmat), main = 'Abundances', col=magma(50))
   "
  R"pdf(file='plots/gbif_sim_era_sol.pdf', paper = 'a4r', height= 8.27, width=11.69 )
  library(rgdal)
  library(fields)
  months = c('January', 'February', 'March', 'April', 'May', 'June', 'July',
  'August', 'September', 'October', 'November', 'December')
  world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
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
 world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
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
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
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
using myunitful
using AxisArrays
using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/Temperature.jld",
 "Temperature")
 Rain = load("/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/Rainfall.jld",
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
file = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/World.tif"
world = extractfile(file)
europe = world[-10° .. 60°, 35° .. 80°]
eu = ustrip.(europe)
@rput eu
R"image(t(eu))"
dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_1990"
tempax2 = extractERA(dir2, "t2m", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_2000"
tempax3 = extractERA(dir3, "t2m", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_2010"
tempax4 = extractERA(dir4, "t2m", collect(361month:1month:38year))

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_1980"
precax1 = extractERA(dir1, "tp", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_1990"
precax2 = extractERA(dir2, "tp", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_2000"
precax3 = extractERA(dir3, "tp", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_2010"
precax4 = extractERA(dir4, "tp", collect(361month:1month:38year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/wc"
prec = extractworldclim(joinpath(dir, "wc2.0_5m_prec"))
prec.array = prec.array[-10° .. 60°, 35° .. 80°,:]
x = prec.array.axes[1]
y = prec.array.axes[2]
t = prec.array.axes[3]
prec.array = AxisArray(1.0.*(prec.array),
    Axis{:longitude}(x), Axis{:latitude}(y), Axis{:time}(t))


dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3},
            srad.array[-10° .. 60°, 35° .. 80°,:])
srad = SolarBudget(srad, 1)
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
  world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
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
 world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
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
using myunitful
using AxisArrays
using ClimatePref

# Load in temperature profiles
Temp = load("/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/Temperature.jld",
 "Temperature")
 Rain = load("/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/Rainfall.jld",
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
file = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/World.tif"
world = extractfile(file)
europe = world[-10° .. 60°, 35° .. 80°]
eu = ustrip.(europe)

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_1980"
tempax1 = extractERA(dir1, "t2m", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_1990"
tempax2 = extractERA(dir2, "t2m", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_2000"
tempax3 = extractERA(dir3, "t2m", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_moda_2010"
tempax4 = extractERA(dir4, "t2m", collect(361month:1month:38year))

dir1 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_1980"
precax1 = extractERA(dir1, "tp", collect(1.0month:1month:10year))
dir2 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_1990"
precax2 = extractERA(dir2, "tp", collect(121month:1month:20year))
dir3 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_2000"
precax3 = extractERA(dir3, "tp", collect(241month:1month:30year))
dir4 = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/era_interim_mdfa_prec_2010"
precax4 = extractERA(dir4, "tp", collect(361month:1month:38year))

dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/wc"
prec = extractworldclim(joinpath(dir, "wc2.0_5m_prec"))
prec.array = prec.array[-10° .. 60°, 35° .. 80°,:]
x = prec.array.axes[1]
y = prec.array.axes[2]
t = prec.array.axes[3]
prec.array = AxisArray(1.0.*(prec.array),
    Axis{:longitude}(x), Axis{:latitude}(y), Axis{:time}(t))


dir = "/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/wc"
srad = extractworldclim(joinpath(dir, "wc2.0_5m_srad"))
srad = convert(Array{typeof(2.0*day^-1*kJ*m^-2),3},
            srad.array[-10° .. 60°, 35° .. 80°,:])
srad = SolarBudget(srad, 1)
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
world = readOGR('/Users/claireh/Documents/PhD/GIT/ClimatePref.jl/data/ne_10m_land/ne_10m_land.shp', layer='ne_10m_land')
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
