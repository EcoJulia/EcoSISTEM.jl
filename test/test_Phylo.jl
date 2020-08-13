using Phylo
using DataFrames
using Simulation
using Test

tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
reroot!(tree, "tip 2")
@test getroot(tree) == "NewRoot"

# Assign continuous trait
tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
trts = DataFrame([[1.0], [0.5]], [:start, :σ²])
assign_traits!(tree, trts)
@test nrow(get_traits(tree)) == 10
resettraits!(tree)
@test_throws ErrorException get_traits(tree)

# Assign discrete trait
tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
switch = [0.5, 0.5]
trts = DataFrame(trait1 = collect(1:2))
assign_traits!(tree, switch, trts)
@test nrow(get_traits(tree)) == 10
resettraits!(tree)
@test_throws ErrorException get_traits(tree)
