using Diversity
using Cubature

import Diversity.getabundance, Diversity.getpartition, Diversity.gettypes
import Diversity.getordinariness!, Diversity.countsubcommunities
import Diversity.getnames
import Diversity.counttypes, Diversity.getnames
import Diversity.getordinariness, Diversity.getsimilarity

## Habitat types
abstract AbstractHabitat


# Habitats : matrix of float values
type Habitats <: AbstractHabitat
  matrix::Matrix{Float64}
  size::Real
end

# Niches : matrix of string values
type Niches <: AbstractHabitat
  matrix::Matrix{String}
  size::Real
end

# Env budget types
abstract AbstractBudget

type Budget <: AbstractBudget
  matrix::Matrix{Float64}
end

# Species trait types
abstract AbstractTraits

type StringTraits <: AbstractTraits
  traits::Vector{String}
end
type RealTraits <: AbstractTraits
  traits::Vector{Real}
end

type TraitRelationship
  matrix::Matrix{Real}
end

# Species energy types
abstract AbstractEnergy

type RealEnergy <: AbstractEnergy
  energy::Vector{Real}
end

abstract AbstractMovement

type GaussianMovement <: AbstractMovement
  var::Vector{Real}
  thresh::Float64
end

function GaussianMovement(move_var, numSpecies, pThresh)
  GaussianMovement(repmat([move_var], numSpecies), pThresh)
end

# Species list type - all info on species
type SpeciesList{FP <: AbstractFloat, M <: AbstractMatrix, T <: AbstractTraits,
                 E<: AbstractEnergy, TR <: Tree, MO <: AbstractMovement} <:
                             AbstractTypes
  similarity::M
  names::Vector{String}
  traits::T
  abun::Vector{Int64}
  energy::E
  phylo::TR
  movement::MO
end

function SpeciesList{FP <: AbstractFloat, T <: AbstractTraits,
    E<: AbstractEnergy, TR <: Tree, MO <: AbstractMovement}(similarity::AbstractMatrix{FP},
    names:: Vector{String}, traits::T, abun::Vector{Int64}, energy::E, phylo::TR,
    movement::MO)
    # Assign names
    names = map(x -> "$x", 1:size(similarity, 1))

    # Check similarity is square matrix
    size(similarity, 1) == size(similarity, 2) ||
    throw(DimensionMismatch("Similarity matrix is not square"))

    # Check dimensions of abundance and similarity match
    length(abun) == size(similarity, 1) ||
    throw(DimensionMismatch("Similarity matrix does not match abundances"))

    # Check similarity is bounded between 0 and 1
    minimum(similarity) ≥ 0 || throw(DomainError())
    maximum(similarity) ≤ 1 || warn("Similarity matrix has values above 1")
    SpeciesList{FP, typeof(similarity), T, E, TR, MO}(similarity, names, traits, abun, energy, phylo, movement)
end

function SpeciesList{E<: AbstractEnergy,
    MO <: AbstractMovement}(NumberSpecies::Int64,
    NumberTraits::Int64, abun_dist::Distribution, energy::E,
    movement::MO)

    # Create tree
    tree = jcoal(NumberSpecies, 100)
    # Create traits and assign to tips
    trts = map(string, 1:NumberTraits)
    assign_traits!(tree, 0.5, trts)
    # Get traits from tree
    sp_trt = StringTraits(get_traits(tree, true))
    # Create similarity matrix (for now identity)
    similarity = eye(NumberSpecies)
    names = map(x -> "$x", 1:size(similarity, 1))
    # Draw random set of abundances from distribution
    abun = rand(abun_dist)
    # error out when abun dist and NumberSpecies are not the same (same for energy dist)
    length(abun)==NumberSpecies || throw(DimensionMismatch("Abundance vector
                                          doesn't match number species"))
    length(energy.energy)==NumberSpecies || throw(DimensionMismatch("Energy vector
                                          doesn't match number species"))
    size(similarity)==(NumberSpecies,NumberSpecies) || throw(DimensionMismatch("
                                Similarity matrix doesn't match number species"))

  SpeciesList(similarity, names, sp_trt, abun, energy, tree, movement)
end

function counttypes(sl::SpeciesList)
    return size(sl.similarity, 1)
end

function getsimilarity(sl::SpeciesList)
    return sl.similarity
end

function floattypes{FP, M}(::SpeciesList{FP, M})
    return Set([FP])
end

function getnames(sl::SpeciesList)
    return sl.names
end

function getordinariness(sl::SpeciesList, a::AbstractArray)
    getsimilarity(sl) * a
end


#function SpeciesList(NumberSpecies::Int64, NumberTraits::Int64,
#                      abun_dist::Distribution, energy::AbstractVector,
#                      movement::GaussianMovement)
#  # Create tree
#  tree = jcoal(NumberSpecies, 100)
#  # Create traits and assign to tips
#  trts = map(string, 1:NumberTraits)
#  assign_traits!(tree, 0.5, trts)
#  # Get traits from tree
#  sp_trt = get_traits(tree, true)
#  # Create similarity matrix (for now identity)
#  similarity = eye(NumberSpecies)
#  # Draw random set of abundances from distribution
#  abun = rand(abun_dist)
#  # error out when abun dist and NumberSpecies are not the same (same for energy dist)
#  length(abun)==NumberSpecies || throw(DimensionMismatch("Abundance vector
#                                        doesn't match number species"))
#  length(energy)==NumberSpecies || throw(DimensionMismatch("Energy vector
#                                        doesn't match number species"))
#  size(similarity)==(NumberSpecies,NumberSpecies) || throw(DimensionMismatch("
#                              Similarity matrix doesn't match number species"))
#
#  SpeciesList(similarity, StringTraits(sp_trt), abun, RealEnergy(energy), tree, movement)
#end


# Abiotic environment types- all info about habitat and relationship to species
# traits
abstract AbstractAbiotic{H <: AbstractHabitat, B <:AbstractBudget} <: AbstractPartition

type MatrixAbioticEnv{H, B} <: AbstractAbiotic{H, B}
  habitat::H
  budget::B
  names::Vector{String}
end

function MatrixAbioticEnv(NumberNiches::Int64, dimension::Tuple, maxBud::Real, size::Real)
  # Create niches
  niches = map(string, 1:NumberNiches)
  # Create niche-like environment
  hab = random_habitat(dimension, niches, 0.5, repmat([0.5], NumberNiches))
  # Create empty budget and for now fill with one value
  bud = zeros(dimension)
  fill!(bud, maxBud)
  names = map(x -> "$x", 1:(dimension[1]*dimension[2]))
  MatrixAbioticEnv(Niches(hab, size), Budget(bud), names)
end

function countsubcommunities(mae::MatrixAbioticEnv)
  return length(mae.habitat.matrix)
end

function getnames(mae::MatrixAbioticEnv)
    return mae.names
end


# Matrix Landscape types - houses abundances (initially empty)
type MatrixLandscape{A <: AbstractArray} <: AbstractArray
  abundances::A
end

function MatrixLandscape(abenv::AbstractAbiotic, spplist::SpeciesList)
  # Create an array of zero abundances according to the size of the habitat
  # and the number of species
  abundances=zeros(length(spplist.abun),size(abenv.habitat.matrix,1),
                   size(abenv.habitat.matrix,2))
  MatrixLandscape(abundances)
end

# Ecosystem type - holds all information and populates ML
abstract AbstractEcosystem{A, Part <: AbstractStructuredPartition,
          S <: SpeciesList, AB <: AbstractAbiotic, R<: TraitRelationship,
          L <: AbstractArray} <: AbstractMetacommunity{Float64, A, Part}

type Ecosystem{A, Part, S, AB, R, L} <: AbstractEcosystem{A, Part, S, AB, R, L}
  partition::Part
  ordinariness::Nullable{A}
  spplist::S
  abenv::AB
  relationship::R
  lookup::L
end

function Ecosystem(spplist::SpeciesList, abenv::AbstractAbiotic, traits::Bool)
  sum(spplist.abun.*spplist.energy.energy)<= sum(abenv.budget.matrix) ||
  error("Environment does not have enough energy to support species")
  # Create matrix landscape of zero abundances
  ml = MatrixLandscape(abenv, spplist)
  # Populate this matrix with species abundances
  species = length(spplist.abun)
  populate!(ml, spplist, abenv, traits)
  # For now create an identity matrix for the species relationships
  A = typeof(ml.abundances)
  rel = eye(length(spplist.traits.traits), size(abenv.habitat.matrix,3))
  lookup_tab = map(x -> lookup(abenv.habitat.size, maximum(size(abenv.habitat.matrix)),
   x, sppl.movement.thresh), sppl.movement.var)
  Ecosystem(ml, Nullable{A}(), spplist, abenv, TraitRelationship(rel), lookup_tab)
end

# Function to copy an Ecosystem type for modification - similar to array copy
function copy_eco(eco::Ecosystem)
  # Collects species list and abiotic environment
  spplist = eco.spplist
  abenv = eco.abenv
  # Create new ML with abundances from ecosystem
  ml= MatrixLandscape(abenv, spplist)
  ml.abundances = copy(eco.partition.abundances)
  # Create relationships same as before
  A = typeof(ml.abundances)
  rel = eye(length(spplist.traits.traits), size(abenv.habitat.matrix,3))
  lookup_tab = eco.lookup
  # Create a new ecosystem with the same components
  Ecosystem(ml, Nullable{A}(), spplist, abenv, TraitRelationship(rel), lookup)
end

# Function like expand.grid in R - given variables it will create every combination in a grid
function expandgrid(args...)
    if length(args) == 0
        error("Nothing to expand")
    elseif length(args) == 1
        return args[1]
    else
        rest = expandgrid(args[2:end]...)
        ret  = Any[]
        for i in args[1]
            for r in rest
                push!(ret, vcat(i,r))
            end
        end
        return ret
    end
end

# Get grid of possible moves for abiotic environment
function get_grid(abenv)
  dims = size(abenv.habitat.matrix)
  combs = collect(-dims[1]:dims[2])
  eg = expandgrid(combs, combs)
  eg = hcat(eg...)'
  #sums = mapslices(sum, abs(eg), 2)
  #order = sortperm(sums[:,1])
  #eg[order, :]
end
function symmetric_grid(grid::Array{Real, 2})
   for x in 1:size(grid,1)
     if grid[x, 1] != grid[x, 2]
       grid = vcat(grid, hcat(grid[x, 2], grid[x, 1] , grid[x, 3]))
     end
   end
   for x in 1:size(grid,1)
     if (grid[x, 1] > 0)
       grid = vcat(grid, hcat(-grid[x, 1], grid[x, 2] , grid[x, 3]))
     end
     if (grid[x, 2] > 0)
       grid = vcat(grid, hcat(grid[x, 1], -grid[x, 2] , grid[x, 3]))
     end
     if (grid[x, 1] > 0 && grid[x, 2] > 0)
       grid = vcat(grid, hcat(-grid[x, 1], -grid[x, 2] , grid[x, 3]))
     end
   end
   grid
 end


function lookup(squareSize::Real, maxGridSize::Int64, dispersalSD::Real,
                pThresh::Float64)
  lookup_tab = Array{Real, 2}(1, 3)
   function disperse(r::AbstractArray)
     1/(π * dispersalSD^2)*exp(-((r[3]-r[1])^2+(r[4]-r[2])^2)/(dispersalSD^2))
   end
    k=0; m=0; count = 0
    while (k<= maxGridSize && m <=maxGridSize)
      count= count + 1
      calc_prob = pcubature(disperse, [0, 0, k*squareSize, m*squareSize],
      [squareSize, squareSize, (k+1)*squareSize, (m+1)*squareSize])[1] /
       squareSize^2
      if m == 0 && calc_prob < pThresh
        break
      end
        if count == 1
           lookup_tab[1, :] = hcat(Int64(k), Int64(m) , calc_prob)
           k = k + 1
        else
          if (calc_prob > pThresh && m <= k)
            lookup_tab = vcat(lookup_tab, hcat(Int64(k), Int64(m) , calc_prob))
            m = m + 1
          else m = 0; k = k + 1
          end
        end
    end

    isdefined(lookup_tab) || error("probability threshold too high")
    lookup_tab = symmetric_grid(lookup_tab)
    lookup_tab[:,1:2] = map(Int64,lookup_tab[:,1:2])
    info(sum(lookup_tab[:, 3]))
    lookup_tab[:, 3] = lookup_tab[:, 3]/sum(lookup_tab[:, 3])
    lookup_tab
end

#eg = expandgrid(1:dims[1], 1:dims[2])
#eg = hcat(eg...)'
#A = eg.== x & eg.== y
#keep = find(!map( (a,b)-> a==x && b==y, eg[:,1], eg[:,2]))
#eg = eg[keep,:]
#eg

#a = eco.spplist.movement.move_var[spp]
# Calculate pdf at points in lookup table

#coords = map((x, y) -> [(x-0.5, y-0.5), (x+0.5, y+0.5)], eco.lookup[:, 2], eco.lookup[:, 3])
#coords = hcat(coords...)
#probs = map((x, y) -> hcubature(disperse, x, y)[1], coords[1,:], coords[2,:])

# Normalise distribution
#eco.lookup[:, 1] = probs/sum(probs)

#points = rand(Uniform(-1,1), 100)
#n = map(x -> Normal(x, 0.5), points)
#m = map(y -> pdf(y, collect(minimum(grd):.1:maximum(grd))), n)
#dist = mapslices(sum, m , 1)
#len = size(grd, 1)
#num = reverse(collect(1:round(len/2, RoundUp)))
#append!(num, collect(1:round(len/2, RoundDown)))
#dist = pdf(Normal(1, 0.8), num)
#@rput dist
#R"plot(1:49,dist, type='l')"
