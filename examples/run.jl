#Load packages
using Diversity
using Diversity.ShortNames
using Simulation
using Plots
using PhyloTrees
using Distributions
#Pkg.update()
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

randtree=jtree(17, Exponential(0.1))
Plots.plot(randtree)

coaltree=jcoal(14, 5)
Plots.plot(coaltree)

tree=jcoal(17, 100)
Plots.plot(tree)
switch_rate=0.5
traits=["A","B", "C"]
trait_tree=assign_trait(tree,switch_rate, traits)
get_traits(trait_tree)
get_times(trait_tree)

Plots.plot(tree,markershape=:circle,
markercolor= [:blue,false],
markerstrokecolor=[:black,false, :red])


using StatsBase
using RCall
# Set up Habitat
species=50; individuals=10000
mat=create_habitat((10,10), ["A","B"], [0.4,0.6])

# Set up tree
tree=jcoal(50, 100)
assign_traits!(tree, 0.5, ["A","B"])
sp_trt=get_traits(tree, true)
# Create ecosystem
pop=populate(species, individuals, Niches(mat), sp_trt)
eco=Ecosystem(pop,Species(), StringTraits(sp_trt))

# Check species richness before
before=SR(eco, 100000)
@rlibrary("fields")
@rput before
R"image.plot(before)"
# Update by one timestep
update!(eco, 0.9, 0.6,0.5, 0.1)
# Check species richness after
after=SR(eco, 100000)
@rput after
R"image.plot(after)"


using StatsBase
using RCall
# p= amount of fragmentation, A = expected proportion of habitat
mat=random_habitat((100,100), ["A","B"], 0.5, [0.5,0.5])
hab= Array{Int64}(100,100)
hab[mat.=="A"]=1
hab[mat.=="B"]=2
@rput hab
R"
image(hab, legend = F, axes=F);
"

species=10; individuals=50000
# Set up tree
tree=jcoal(species, 100)
assign_traits!(tree, 0.2, ["A","B"])
sp_trt=get_traits(tree, true)
budg= Array{Float64}(50,50)
fill!(budg, 100)
energy=repmat([1], species)
# Try with a skewed distribution
pop=populate(species, individuals, Niches(mat), sp_trt, Budget(budg),
  Multinomial(individuals, rand(Dirichlet(species,1))))
eco=Ecosystem(pop, Species(), StringTraits(sp_trt), RealEnergy(energy))
maximum(mapslices(sum,eco.partition.abundances,1))

@rlibrary("fields")
a = subdiv(ᾱ(eco), 2)
@rput a
R"par(mfrow=c(1,2));image.plot(a,col=heat.colors(100), breaks=seq(0,max(a), length.out=101));image(hab, legend = F)"

p = subdiv(ρ̄(eco), 1)
@rput p
R"par(mfrow=c(1,2));image.plot(p,col=heat.colors(100), breaks=seq(0,max(p), length.out=101));image(hab, legend = F)"

b = subdiv(β(eco), 1)
@rput b
R"par(mfrow=c(1,2));image.plot(b,col=heat.colors(100), breaks=seq(0,max(b), length.out=101));image(hab, legend = F)"

g = subdiv(γ(eco), 1)
@rput g
R"par(mfrow=c(1,2));image.plot(b,col=heat.colors(100), breaks=seq(0,max(g), length.out=101));image(hab, legend = F)"

# Set up initial conditions
birth = 0.4
death = 0.4
move = 0.5
timestep = 1
# Check species richness before
before=subdiv(ᾱ(eco), 2)
@rlibrary("fields")
@rlibrary("grDevices")
@rput before
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(50)[1:20], breaks=seq(0,20,1));image(hab, legend = F)"
R"pdf(file='Before_steady.pdf', paper='a4r',onefile=T,width=11.69,height=6)"
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(20), breaks=seq(0,20,1));image(hab, legend = F)"
R"dev.off()"
for i in 1:10

# Update by one timestep
update!(eco, birth, death, move, timestep)
# Check species richness after
after=subdiv(ᾱ(eco), 2)
@rput after
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(50)[1:20], breaks=seq(0,20,1));image(hab, legend = F)"
maximum(mapslices(sum,eco.partition.abundances,1))
@rput i
#R"par(mfrow=c(1,2));image.plot(after,col=rainbow(20), breaks=seq(0,50,1));image(hab, legend = F)"
R"pdf(file=paste('After_steady', i, '.pdf', sep=''), paper='a4r',onefile=T,width=11.69,height=6)"
R"par(mfrow=c(1,2));print(image.plot(after,col=rainbow(20), breaks=seq(0,20,1)));image(hab, legend = F)"
R"dev.off()"
end
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(20), breaks=seq(0,20,1));image(hab, legend = F)"

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
move = 0.0
l = 1.0
s = 0.01
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, move, timestep, l, s]

# Create ecosystem
#sppl = SpeciesList(numSpecies, numTraits,
#                   energy_vec)
#abenv = MatrixAbioticEnv(numNiches, (1,1), 1000)
#eco = Ecosystem(sppl, abenv)

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
move = 0.1
l = 1.0
s = 0.0
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, move, timestep, l, s]

# Create ecosystem
#sppl = SpeciesList(numSpecies, numTraits, Multinomial(50, numSpecies),
#                   energy_vec)
#abenv = MatrixAbioticEnv(numNiches, (2,2), 500)
#start=zeros(100, 2)
#for j in 1:100
#  eco = Ecosystem(sppl, abenv)
#  start[j, :]=eco.spplist.abun
#end

#eco = Ecosystem(sppl, abenv)

# Run simulations 100 times
#ab = run_sim(eco, param, 1000, 100)
ab = run_sim(numSpecies, numTraits, Multinomial(50, numSpecies), energy_vec,
  numNiches, (2,2), 500,
  param, 1000, 200)

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
move = 0.1
l = 1.0
s = 0.1
timestep = 1

# Collect model parameters together (in this order!!)
param = [birth, death, move, timestep, l, s]

# Create ecosystem
sppl = SpeciesList(numSpecies, numTraits, Multinomial(500, numSpecies),
                   energy_vec)
abenv = MatrixAbioticEnv(numNiches, (3,3), 100000)
eco = Ecosystem(sppl, abenv)

# Run simulations 100 times
ab = run_sim(eco, param, 1000, 10)

# Calculate mean and confidence limits
mean = ab[2][:, :, :, 1]
uc = ab[2][:, :, :, 1] + ab[3][:, :, :, 1]
lc = ab[2][:, :, :, 1] - ab[3][:, :, :, 1]

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
