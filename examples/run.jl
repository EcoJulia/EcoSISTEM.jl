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
# Set up Habitat
species=50; individuals=10000; mat=ones(10, 10)
mat=create_habitat((10,10), ["A","B"], [0.2,0.8])
# Set up tree
tree=jcoal(50, 100)
assign_traits!(tree, 0.5, ["A","B"])
sp_trt=get_traits(tree, true)
# Create ecosystem
pop=populate(species, individuals, Niches(mat), sp_trt)
eco=Ecosystem(pop,Species(), StringTraits(sp_trt))
# Check species richness before
before=SR(eco, 100000)
@rput before
R"image(before)"
# Update by one timestep
update!(eco, 0.9, 0.6,0.5, 1)
# Check species richness after
after=SR(eco, 100000)
@rput after
R"image(after)"
