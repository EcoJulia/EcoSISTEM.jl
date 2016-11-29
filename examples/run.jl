using Simulation
using Diversity
using Diversity.β

# Create a partition
part=MatrixLandscape(reshape([1, 2, 3, 4], (1, 2, 2)), Habitats([1 2; 3 4]))

# Create an ecosystem
eco=Ecosystem(part, Species(), StringTraits(["A"]))

# Calculate ordinariness of ecosystem
getordinariness!(eco)

b = β(eco)

sb = subdiv(β(eco), 3)
