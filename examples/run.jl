# Include simulation functions
using Simulation
using Distributions
#include("Load.jl")

# Create a partition
part=MatrixLandscape(reshape([1, 2, 3, 4, 5, 6, 7, 8], (2, 2, 2)), Habitats([1 2; 3 4]))

# Create an ecosystem
eco=Ecosystem(part, Species(), StringTraits(["A", "B"]))

# Calculate ordinariness of ecosystem
getordinariness!(eco)

b = β(eco)

sb = subdiv(β(eco), 3)
m = metadiv(β(eco), 3)

# Create a simple habitat matrix
mat=ones(10, 10)
# Populate habitat with 10,000 individuals from 50 species
LS=populate(50, 10000, Habitats(mat))
# Create ecosystem from habitat, making every species distinct and
# having the same trait
eco=Ecosystem(LS,Species(), StringTraits(repmat(["A"],50)))
# Create grid of species richness
sr=SR(eco, 10000)

# Plot
heatmap(LS.abundances[1, :, :])
heatmap(sr)

# Investigate alpha Diversity
alpha_bar=subdiv(ᾱ(eco),0)
heatmap(alpha_bar)

## Investigate multiple runs for single population with two species
using RCall
numSpecies=2
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = [2, 2]

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.01
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, timestep, l, s]

# Create ecosystem
sppl = SpeciesList(numSpecies, numTraits,
                  energy_vec)
abenv = MatrixAbioticEnv(numNiches, (1,1), 1000)
eco = Ecosystem(sppl, abenv)

# Run simulations 100 times
ab = run_sim(numSpecies, numTraits, Multinomial(50, numSpecies), energy_vec,
  numNiches, (1,1), 500,
  param, 1000, 1)

# Calculate mean and confidence limits
mean = ab[2][:, :, 1]
uc = ab[2][:, :, 1] + ab[3][:, :, 1]
lc = ab[2][:, :, 1] - ab[3][:, :, 1]

meanE = ab[4][:, :, 1]
ucE = ab[4][:, :, 1] + ab[5][:, :, 1]
lcE = ab[4][:, :, 1] - ab[5][:, :, 1]

# Input results to R
@rput mean;@rput meanE
@rput uc;@rput ucE
@rput lc;@rput lcE
@rput energy_vec; @rput numSpecies

# Plot in R
R"library(scales); par(mfcol=c(1,numSpecies), xaxs='i', yaxs='i');
row=1:nrow(meanE)
  for (j in 1:ncol(mean)){
    if (j==1) plot_fun=plot else plot_fun=lines
      plot_fun(row,mean[,j],
           ylab = list('Abundance',cex=1.4), xlab = list('Time',cex=1.4), type='l',
           col=j,
           ylim=c(0, max(uc)+50), cex.axis=1.2, cex.main=1.4);
      polygon(c(row,rev(row)),
              c(lc[,j],rev(uc[,j])),
              col=alpha(j, 0.3), border =F);
      legend('topright', legend= paste('Species', 1:numSpecies, 'energy=', energy_vec),
      pch=20, col=1:2)
    }
    plot(row,meanE,
                 ylab = list('Energy', cex=1.4), xlab = list('Time', cex=1.4), type='l',
                 col='black',
                 ylim=c(0, max(ucE)), cex.axis=1.2);
    polygon(c(row,rev(row)),
                    c(lcE,rev(ucE)),
                    col=alpha('black', 0.3), border =F)"

R"library(scales); par(mfcol=c(1,numSpecies), xaxs='i', yaxs='i');
row=1:nrow(mean)
for (j in 1:ncol(mean)){
  if (j==1) plot_fun=plot else plot_fun=lines
plot_fun(row, mean[,j],
     ylab = list('Abundance',cex=1.4), xlab = list('Time',cex=1.4), type='l',
     col=j,
     ylim=c(0, max(mean)+50), cex.axis=1.2, cex.main=1.4)
     }
     plot(row,meanE,
                 ylab = list('Energy', cex=1.4), xlab = list('Time', cex=1.4), type='l',
                 col='black',
                 ylim=c(0, max(meanE)), cex.axis=1.2)"

R"library(scales); par(mfcol=c(1,numSpecies+1), xaxs='i', yaxs='i');
row=1:nrow(meanE)
  for (j in 1:ncol(mean)){
      plot(row,mean[,j],
           ylab = list('Abundance',cex=1.4), xlab = list('Time',cex=1.4), type='l',
           col=j, main=paste('Species', j, ', energy=', energy_vec[j]),
           ylim=c(0, max(uc[,j])), cex.axis=1.2, cex.main=1.4);
      polygon(c(row,rev(row)),
              c(lc[,j],rev(uc[,j])),
              col=alpha(j, 0.3), border =F);
    }
    plot(row,meanE,
                 ylab = list('Energy', cex=1.4), xlab = list('Time', cex=1.4), type='l',
                 col='black',
                 ylim=c(0, max(ucE)), cex.axis=1.2);
    polygon(c(row,rev(row)),
                    c(lcE,rev(ucE)),
                    col=alpha('black', 0.3), border =F)"


## Run simulation on larger grid : 2 by 2
using RCall
numSpecies=2
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = [2, 4]

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, timestep, l, s]

# Create ecosystem
sppl = SpeciesList(numSpecies, numTraits, Multinomial(50, numSpecies),
                  energy_vec)
abenv = MatrixAbioticEnv(numNiches, (2,2), 500)

eco = Ecosystem(sppl, abenv)

# Run simulations 100 times
ab = run_sim(eco, param, 1000, 200)

# Calculate mean and confidence limits
mean = ab[2][:, :, :, 1]
uc = ab[2][:, :, :, 1] + ab[3][:, :, :, 1]
lc = ab[2][:, :, :, 1] - ab[3][:, :, :, 1]

meanE = ab[4][:, :, 1]
ucE = ab[4][:, :, 1] + ab[5][:, :, 1]
lcE = ab[4][:, :, 1] - ab[5][:, :, 1]

gridSize = 4

# Input results to R
@rput mean;@rput meanE
@rput uc;@rput ucE
@rput lc;@rput lcE
@rput energy_vec; @rput gridSize

# Plot in R
R"library(scales); par(mfcol=c(2,gridSize), xaxs='i', yaxs='i');
rows = 1:dim(mean)[1]
cols = 1:dim(mean)[2]
for (i in 1:dim(mean)[3]){
  for (j in cols){
    if(j ==1) plot_fun=plot else plot_fun=lines
      plot_fun(rows,mean[,j,i],
           ylab = list('Abundance', cex=1.4), xlab = list('Time', 1.4), type='l',
           col=j, main=paste('Population', i),
           ylim=c(0, max(uc[,j,])+50), cex.axis=1.2, cex.main=1.4);
      polygon(c(rows,rev(rows)),
              c(lc[,j,i],rev(uc[,j,i])),
              col=alpha(j, 0.3), border =F)
      legend('topright', legend= paste('Species', cols, ', energy=', energy_vec[cols]), col=cols, pch=20)
              }
              }
for (k in 1:dim(mean)[3]){
    plot(rows,meanE[,k],
                 ylab = list('Energy', cex=1.4), xlab = list('Time', cex=1.4), type='l',
                 col='black',
                 ylim=c(0, max(ucE[,k])), main = paste('Population', k), cex.axis=1.2, cex.main=1.4);
    polygon(c(rows,rev(rows)),
                    c(lcE[,k],rev(ucE[,k])),
                    col=alpha('black', 0.3), border =F)}
"

## Run simulation on larger grid : 3 by 3
using RCall
numSpecies=2
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = [2, 3]

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.1
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, timestep, l, s]

# Create ecosystem
sppl = SpeciesList(numSpecies, numTraits, Multinomial(500, numSpecies),
                   energy_vec)
abenv = MatrixAbioticEnv(numNiches, (3,3), 100000)
eco = Ecosystem(sppl, abenv, false)

# Run simulations 100 times
ab = run_sim(eco, param, 1000, 10)

# Calculate mean and confidence limits
# For abundance:
mean = ab[2][:, :, :, 1]
uc = ab[2][:, :, :, 1] + ab[3][:, :, :, 1]
lc = ab[2][:, :, :, 1] - ab[3][:, :, :, 1]

# For energy:
meanE = ab[4][:, :, 1]
ucE = ab[4][:, :, 1] + ab[5][:, :, 1]
lcE = ab[4][:, :, 1] - ab[5][:, :, 1]

gridSize = [size(eco.abenv.habitat.matrix,1), size(eco.abenv.habitat.matrix,2)]

# Input results to R
@rput mean;@rput meanE
@rput uc;@rput ucE
@rput lc;@rput lcE
@rput energy_vec; @rput gridSize

# Plot in R
R"library(scales); par(mfcol=c(gridSize[1], gridSize[1]*2), xaxs='i', yaxs='i');
rows = 1:dim(mean)[1]
cols = 1:dim(mean)[2]
for (i in 1:dim(mean)[3]){
  for (j in cols){
    if(j ==1) plot_fun=plot else plot_fun=lines
      plot_fun(rows,mean[,j,i],
           ylab = 'Abundance', xlab = 'Time', type='l',
           col=j, main=paste('Population', i),
           ylim=c(0, max(uc[,j,])));
      polygon(c(rows,rev(rows)),
              c(lc[,j,i],rev(uc[,j,i])),
              col=alpha(j, 0.3), border =F)
      legend('topright', legend= paste('Species', cols, ', energy=', energy_vec[cols]), col=cols, pch=20)
              }
              }
for (k in 1:dim(mean)[3]){
    plot(rows,meanE[,k],
                 ylab = 'Energy', xlab = 'Time', type='l',
                 col='black',
                 ylim=c(0, max(ucE[,k])), main = paste('Population', k));
    polygon(c(rows,rev(rows)),
                    c(lcE[,k],rev(ucE[,k])),
                    col=alpha('black', 0.3), border =F)}
"

## Run simulation over a grid and plot
numSpecies=4
numTraits=2
numNiches=2

# Set up how much energy each species consumes
energy_vec = RealEnergy(repmat([2], numSpecies))

# Set probabilities
birth = 0.6
death = 0.6
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, timestep, l, s]

grid = (5,5)
gridSize = 1
totalK = 1000
individuals=100

# Create ecosystem
movement = GaussianMovement(0.2, numSpecies, 10e-4)
spplist = SpeciesList(numSpecies, numTraits, Multinomial(individuals, numSpecies),
                   energy_vec, movement)
abenv = MatrixAbioticEnv(numNiches, grid, totalK, gridSize)
eco = Ecosystem(spplist, abenv, false)
plot_move(eco, 2, 2, 1)
# Run simulation grid
abun = run_sim_spatial(eco, param, 100, 1, 1, 1, false)

plot_abun(abun, numSpecies, grid[1])
