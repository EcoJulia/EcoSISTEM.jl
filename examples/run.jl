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



# p= amount of fragmentation, A = expected proportion of habitat
mat=random_habitat((50,50), ["A","B"], 0.5, [0.5,0.5])
hab= Array{Int64}(50,50)
hab[mat.=="A"]=1
hab[mat.=="B"]=2
@rput hab
R"image(hab, legend = F)"

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
eco=Ecosystem(pop,Species(), StringTraits(sp_trt), RealEnergy(energy))
maximum(mapslices(sum,eco.partition.abundances,1))

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
before=SR(eco)
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
after=SR(eco)
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
