using StatsBase
using ProgressMeter

"""
    get_neighbours(mat::Matrix, x_coord::Int64, y_coord::Int64, chess::Int64=4)

Function to get the neighbours of a grid square in a matrix in 4 or 8 directions
"""
function get_neighbours(mat::Matrix, x_coord::Int64, y_coord::Int64, chess::Int64=4)
  # Calculate dimensions
  dims=size(mat)
  x_coord <= dims[1] && y_coord <= dims[2] || error("Coordinates outside grid")
  # Include 4 directions
  if chess==4
    neighbour_vec=[x_coord y_coord-1; x_coord y_coord+1; x_coord-1 y_coord;
     x_coord+1 y_coord]
  # Include 8 directions
  elseif chess==8
    neighbour_vec=[x_coord y_coord-1; x_coord y_coord+1; x_coord-1 y_coord;
     x_coord+1 y_coord; x_coord-1 y_coord-1; x_coord-1 y_coord+1;
      x_coord+1 y_coord-1; x_coord+1 y_coord+1]
  else
    # Give error if other number chosen than 4 or 8
    error("Can only calculate neighbours in 4 or 8 directions")
  end
  # Remove answers outside of the dimensions of the matrix
  remove=vcat(mapslices(all, [neighbour_vec.>=1 neighbour_vec[:,1].<=
    dims[1] neighbour_vec[:,2].<=dims[2]], 2)...)
  neighbour_vec=neighbour_vec[remove,:]
  neighbour_vec
end
"""
    update!(eco::Ecosystem,  birth::Float64, death::Float64,
       l::Float64, s::Float64, timestep::Real)
Function to update a Ecosystem after one timestep. It takes in parameters of
birth, death rates and longevity of species (l & s) to generate the abundances
of the species stochastically. Movement takes place across the landscape via
movement rates defined in the ecosystem.
"""
function update!(eco::Ecosystem, timestep::Float64)

  # Calculate abundance in overall grid (to be implemented later?)
  #abun=map(i->sum(eco.abundances[i,:,:]), 1:size(eco.abundances,1))

  # Calculate dimenions of habitat and number of species
  dims = length(eco.abenv.habitat.matrix)
  spp = size(eco.abundances.grid,1)
  net_migration = zeros(size(eco.abundances.matrix))
  params = eco.spplist.params
  # Loop through grid squares
  for i in 1:dims

      # Get the overall energy budget of that square
      K = eco.abenv.budget.matrix[i]
      randomise=collect(1:spp)
      randomise=randomise[randperm(length(randomise))]
      # Loop through species in chosen square
      for j in randomise

        # Get abundances of square we are interested in
        square = eco.abundances.matrix[:, i]

        if square[j] <= 0
          eco.abundances.matrix[j, i] = 0
        else
        # Get energy budgets of species in square
        ϵ̄ = eco.spplist.requirement.energy
        E = sum(square .* ϵ̄)

        # Alter rates by energy available in current pop & own requirements
        birth_energy = (ϵ̄[j])^(-params.l - params.s) * K / E
        death_energy = (ϵ̄[j])^(-params.l + params.s) * E / K

        # Calculate effective rates
        birthprob = params.birth[j] * timestep * birth_energy
        deathprob = params.death[j] * timestep * death_energy

        # If traits are same as habitat type then give birth "boost"
        #if eco.spplist.traits.trait[j] != eco.abenv.habitat.matrix[x, y]
        #  birthrate = birthrate * 0.8
        #end

        # If zero abundance then go extinct
        if square[j] == 0
          birthprob = 0
          deathprob = 0
        end

        # Throw error if rates exceed 1
        #birthprob <= 1 && deathprob <= 1 && moveprob <= 1 ||
        #  error("rates larger than one in binomial draw")


        # Put probabilities into 0 - 1
        probs = map(prob -> 1 - exp(-prob), [birthprob, deathprob])

        # Calculate how many births and deaths
        births = jbinom(1, Int(square[j]), probs[1])[1]
        deaths = jbinom(1, Int(square[j]), probs[2])[1]

        # Update population
        eco.abundances.matrix[j, i] = eco.abundances.matrix[j, i] +
          births - deaths

        # Perform gaussian movement
        move!(eco, eco.spplist.movement, i, j, net_migration, births)
      end
    end
  end
  eco.abundances.matrix .= eco.abundances.matrix .+ net_migration
end

function convert_coords(i::Int64, width::Int64)
  i = i - 1
  x = (i % width) + 1
  y = div(i, width)  + 1
  return (x, y)
end
function convert_coords(x::Int64, y::Int64, width::Int64)
  x = x - 1 ; y = y - 1
  i = x + width * y
  return i + 1
end

function calc_lookup_moves(i::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  width = size(eco.abenv.habitat.matrix, 1)
  (x, y) = convert_coords(i, width)
  # Draw moves from Multinomial dist
  table = copy(eco.lookup[spp])
  table[:X] += x; table[:Y] += y
  maxX = size(eco.abenv.habitat.matrix, 1)
  maxY = size(eco.abenv.habitat.matrix, 2)
  # Can't go over maximum dimension
  lower  = find(mapslices(x -> all(x .> 0), Array(table), 2))
  upperX = find(map(x -> (x .<= maxX), table[:,1]))
  upperY = find(map(x -> (x .<= maxY), table[:,2]))
  valid = intersect(lower, upperX, upperY)
  table = table[valid, :]
  table[:Prob] = table[:Prob]/sum(table[:Prob])

  moves = rand(Multinomial(abun,
          Vector{Float64}(table[:Prob])))
  # Add moves to lookup table
  table[:Moves] = moves
  table
end
"""
    move!(i::Int64, spp::Int64, eco::Ecosystem, grd::Array{Float64, 2})

Function to calculate the movement of species `spp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(eco::Ecosystem, ::AlwaysMovement, i::Int64, spp::Int64, grd::Array{Float64, 2}, ::Int64)
  width = size(eco.abenv.habitat.matrix, 1)
  full_abun = Int64(eco.abundances.matrix[spp, i])
  table = calc_lookup_moves(i, spp, eco, full_abun)
  # Lose moves from current grid square
  grd[spp, i] = grd[spp, i] - sum(table[:Moves])
  # Map moves to location in grid
  map(x -> grd[spp, convert_coords(table[x, :X], table[x, :Y], width)] = grd[spp,
  convert_coords(table[x, :X], table[x, :Y], width)] + table[x, :Moves], 1:nrow(table))
end


function move!(eco::Ecosystem, ::BirthOnlyMovement, i::Int64, spp::Int64, grd::Array{Float64, 2},
                births::Int64)
  width = size(eco.abenv.habitat.matrix, 1)
  table = calc_lookup_moves(i, spp, eco, births)
  # Lose moves from current grid square
  grd[spp, i] = grd[spp, i] - sum(table[:Moves])
  # Map moves to location in grid
  map(x -> grd[spp, convert_coords(table[x, :X], table[x, :Y], width)] = grd[spp,
  convert_coords(table[x, :X], table[x, :Y], width)] + table[x, :Moves], 1:nrow(table))
end



"""
    populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)

Function to populate a grid landscape given the abundances found in species list
and whether or not to include traits.
"""
function populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)
  # Calculate size of habitat
  dim = size(abenv.habitat.matrix)
  len = dim[1]*dim[2]
  grid = collect(1:len)
  # Set up copy of budget
  b=copy(abenv.budget.matrix)
  # Loop through species
  for i in eachindex(spplist.abun)
    # Get abundance of species
    abun=spplist.abun[i]
    # Get species trait
    pref=spplist.traits.trait[i]
    # Calculate weighting, giving preference to squares that match with trait
    wv= Vector{Float64}(grid)
    wv[find(reshape(abenv.habitat.matrix, (len,1))[grid].==pref)]= 0.5
    wv[find(reshape(abenv.habitat.matrix, (len,1))[grid].!=pref)]= 0.5
    # Loop through individuals
      while abun>0
        zs=findin(b[grid], 0)
        deleteat!(grid, zs)
        deleteat!(wv, zs)
        # Randomly choose position on grid (weighted)
        if traits pos=sample(grid, weights(wv)) else pos=sample(grid) end
      # Add individual to this location
      ml.grid[i,pos]=ml.grid[i,pos]+1
      abun=abun-1
      b[pos]=b[pos]-1
    end
  end
end
"""
    repopulate!(eco::Ecosystem, traits::Bool)
Function to repopulate an ecosystem `eco`, with option for including trait
preferences.
"""
function repopulate!(eco::Ecosystem, traits::Bool)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun), length(eco.spplist.abun)))
  populate!(eco.abundances, eco.spplist, eco.abenv, traits)
end
