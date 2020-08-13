using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units
using Phylo

import Simulation: DiscreteTrait

# Gaussian trait
opts = fill(5.0°C, numSpecies)
vars = rand(Uniform(0, 25/9), numSpecies)  * °C
@test_nowarn GaussTrait(opts, vars)
@test Simulation.iscontinuous(GaussTrait(opts, vars)) == true
@test eltype(GaussTrait(opts, vars)) <: Unitful.Temperature

# Discrete trait
@test_nowarn DiscreteTrait(fill(1, 10))
@test Simulation.iscontinuous(DiscreteTrait(fill(1, 10))) == false
@test eltype(DiscreteTrait(fill(1, 10))) <: Int

# Temperature bin
@test_nowarn TempBin(fill(1, 10, 2))
@test Simulation.iscontinuous(TempBin(fill(1, 10, 2))) == true
@test eltype(TempBin(fill(1, 10, 2))) <: Unitful.Temperature

# Rainfall bin
@test_nowarn RainBin(fill(1, 10, 2))
@test Simulation.iscontinuous(RainBin(fill(1, 10, 2))) == true
@test eltype(RainBin(fill(1, 10, 2))) <: Unitful.Length
@test_nowarn TraitCollection2(TempBin(fill(1, 10, 2)), RainBin(fill(1, 10, 2)))

# Multiple traits
tr2 = TraitCollection2(TempBin(fill(1, 10, 2)), RainBin(fill(1, 10, 2)))
@test Simulation.iscontinuous(tr2) == [true, true]
@test eltype(tr2) == [typeof(1.0K), typeof(1.0mm)]
@test_nowarn TraitCollection3(GaussTrait(opts, vars), TempBin(fill(1, 10, 2)), RainBin(fill(1, 10, 2)))
tr3 = TraitCollection3(GaussTrait(opts, vars), TempBin(fill(1, 10, 2)), RainBin(fill(1, 10, 2)))
@test Simulation.iscontinuous(tr3) == [true, true, true]
@test eltype(tr3) == [typeof(1.0K), typeof(1.0K), typeof(1.0mm)]


# Test evolution of traits
tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
@test_nowarn DiscreteEvolve(2, tree, 0.5)
@test typeof(get_traits(tree)) <: DataFrame
@test maximum(get_traits(tree)[!, :trait1]) == 2

tree = rand(Ultrametric{BinaryTree{OneRoot, DataFrame, DataFrame}}(10))
@test_nowarn ContinuousEvolve(1.0, 0.1, tree)
@test typeof(get_traits(tree)) <: DataFrame
@test all(get_traits(tree)[!, :σ²] .== 0.1)
@test length(unique(get_traits(tree)[!, :start])) > 1
