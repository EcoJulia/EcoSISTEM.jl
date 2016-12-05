using Simulation
using Diversity
using Diversity.ShortNames
using Plots
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
