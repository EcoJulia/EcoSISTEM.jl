using StatsBase
using ProgressMeter
using Query

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
function update!(eco::Ecosystem, timestep::Unitful.Time)
  # Calculate dimenions of habitat and number of species
  dims = length(eco.abenv.habitat.matrix)
  spp = size(eco.abundances.grid,1)
  net_migration = zeros(Int64, size(eco.abundances.matrix))
  params = eco.spplist.params
  # Loop through grid squares
  for i in 1:dims

      # Get the overall energy budget of that square
      width = size(eco.abenv.habitat.matrix, 1)
      (x, y) = convert_coords(i, width)
      K = eco.abenv.budget.matrix[x, y]
      randomise=collect(1:spp)
      randomise=randomise[randperm(length(randomise))]
      # Get abundances of square we are interested in
      currentabun = view(eco.abundances.matrix, :, i)

      # Get energy budgets of species in square
      ϵ̄ = eco.spplist.requirement.energy
      E = sum(convert(Vector{Float64}, currentabun) .* ϵ̄)
      # Traits
      ϵ̄real = map(ϵ̄, 1:spp) do epsilon, k
        epsilon/TraitFun(eco, i, k)
      end

      # Loop through species in chosen square
      for j in randomise
        if currentabun[j] <= 0
          eco.abundances.matrix[j, i] = 0
        else

        #params.s = 1 - exp(-params.s / TraitFun(eco, i, j))
        # Alter rates by energy available in current pop & own requirements
        birth_energy = ϵ̄[j]^-params.l * ϵ̄real[j]^-params.s * min(K/E, params.boost)
        death_energy = ϵ̄[j]^-params.l * ϵ̄real[j]^params.s * (E / K)

        # Calculate effective rates
        birth = uconvert(unit(timestep)^-1, params.birth[j])
        birthprob = birth * timestep * birth_energy
        death = uconvert(unit(timestep)^-1, params.death[j])
        deathprob = death * timestep * death_energy
      # If zero abundance then go extinct
        if currentabun[j] == 0
          birthprob = 0
          deathprob = 0
        end
        # Put probabilities into 0 - 1
        newbirthprob, newdeathprob = 1.0 - exp.(-[birthprob, deathprob])

        # Calculate how many births and deaths
        births = jbinom(1, currentabun[j], newbirthprob)[1]
        deaths = jbinom(1, currentabun[j], newdeathprob)[1]

        # Update population
        eco.abundances.matrix[j, i] += (births - deaths)

        # Perform gaussian movement
        move!(eco, eco.spplist.movement, i, j, net_migration, births)
      end
    end
  end
  eco.abundances.matrix .= eco.abundances.matrix .+ net_migration
  # Update environment
  getchangefun(eco)(eco, timestep)
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

  lookup = eco.lookup[spp]
  maxX = size(eco.abenv.habitat.matrix, 1) - x
  maxY = size(eco.abenv.habitat.matrix, 2) - y
  # Can't go over maximum dimension
  valid = find((lookup.x .> -x) .& (lookup.y .> -y) .&
   (lookup.x .<= maxX) .& (lookup.y .<= maxY))
  probs = lookup.p[valid]
  probs ./= sum(probs)
  moves = rand(Multinomial(abun, probs))
  return hcat(((lookup.x[valid] .+ x), (lookup.y[valid] .+ y), moves)...)
end
"""
    move!(i::Int64, spp::Int64, eco::Ecosystem, grd::Array{Int64, 2})

Function to calculate the movement of species `spp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(eco::Ecosystem, ::AlwaysMovement, i::Int64, spp::Int64,
  grd::Array{Int64, 2}, ::Int64)

  width = size(eco.abenv.habitat.matrix, 1)
  full_abun = Int64(eco.abundances.matrix[spp, i])
  table = calc_lookup_moves(i, spp, eco, full_abun)
  # Lose moves from current grid square
  grd[spp, i] = grd[spp, i] - sum(table[:, 3])
  # Map moves to location in grid
  locs = map(row -> convert_coords(table[row, 1], table[row, 2], width),
   1:size(table, 1))
  grd[spp, locs] .= grd[spp, locs] .+ table[:, 3]
  return eco
end
"""
    move!(i::Int64, spp::Int64, eco::Ecosystem, grd::Array{Int64, 2})

Function to calculate the movement of species `spp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(eco::Ecosystem, ::NoMovement, i::Int64, spp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  return eco
end


function move!(eco::Ecosystem, ::BirthOnlyMovement, i::Int64, spp::Int64, grd::Array{Int64, 2},
                births::Int64)
  width = size(eco.abenv.habitat.matrix, 1)
  table = calc_lookup_moves(i, spp, eco, births)
  # Lose moves from current grid square
  grd[spp, i] = grd[spp, i] - sum(table[:, 3])
  # Map moves to location in grid
  locs = map(row -> convert_coords(table[row, 1], table[row, 2], width),
   1:size(table, 1))
  grd[spp, locs] .= grd[spp, locs] .+ table[:, 3]
  return eco
end



"""
    populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)

Function to populate a grid landscape given the abundances found in species list
and whether or not to include traits.
"""
function populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic)
  # Calculate size of habitat
  dim = size(abenv.habitat.matrix)
  len = dim[1] * dim[2]
  grid = collect(1:len)
  # Set up copy of budget
  b = reshape(copy(abenv.budget.matrix), size(grid))
  # Loop through species
  for i in eachindex(spplist.abun)
    # Get abundance of species
    abun = spplist.abun[i]
    # Loop through individuals
      while abun>0
        # Randomly choose position on grid (weighted)
        pos = sample(grid[b .> 0])
      # Add individual to this location
      ml.matrix[i, pos] = ml.matrix[i, pos] .+ 1
      abun = abun .- 1
      b[pos] = b[pos] .- spplist.requirement.energy[i]
    end
  end
end
"""
    repopulate!(eco::Ecosystem, traits::Bool)
Function to repopulate an ecosystem `eco`, with option for including trait
preferences.
"""
function repopulate!(eco::Ecosystem)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun), length(eco.spplist.abun)))
  populate!(eco.abundances, eco.spplist, eco.abenv)
end
