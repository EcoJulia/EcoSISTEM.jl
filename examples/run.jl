using Simulation
using Diversity
using Diversity.ShortNames
# Create a partition
part=MatrixLandscape(reshape([1, 2, 3, 4, 5, 6, 7, 8], (2, 2, 2)), Habitats([1 2; 3 4]))

# Create an ecosystem
eco=Ecosystem(part, Species(), StringTraits(["A", "B"]))

# Calculate ordinariness of ecosystem
getordinariness!(eco)

b = β(eco)

sb = subdiv(β(eco), 3)
m = metadiv(β(eco), 3)


mat=ones(10, 10)
LS=populate(50, 10000, Habitats(mat))
eco=Ecosystem(LS,Species(), StringTraits(repmat(["A"],50)))

