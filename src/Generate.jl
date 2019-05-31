using StatsBase
using Compat
using LinearAlgebra
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
    dims[1] neighbour_vec[:,2].<=dims[2]], dims=2)...)
  neighbour_vec=neighbour_vec[remove,:]
  neighbour_vec
end
function get_neighbours(mat::Matrix, x_coord::Array{Int64,1},
     y_coord::Array{Int64,1}, chess::Int64=4)
     neighbours  =map(n -> get_neighbours(mat, x_coord[n], y_coord[n], chess),
      eachindex(x_coord))
      return vcat(neighbours...)
end

update!(eco::Ecosystem, timestep::Unitful.Time) =
    update!(eco, timestep, Val{Threads.nthreads()}())

"""
    update!(eco::Ecosystem, time::Unitful.Time)
Function to update a ecosystem abundances and environment for one timestep.
"""
function update!(eco::Ecosystem, timestep::Unitful.Time, ::Val{N}) where N

    # Calculate dimenions of habitat and number of species
    dims = _countsubcommunities(eco.abenv.habitat)
    spp = size(eco.abundances.grid,1)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall energy budget of that square
    update_energy_usage!(eco, Val{N}())

    # Loop through species in chosen square
    Threads.@threads for j in 1:spp
        rng = eco.abundances.seed[Threads.threadid()]
        # Loop through grid squares
        for i in 1:dims
            # Calculate how much birth and death should be adjusted
            adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, i, j)

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(i, width)
            # Check if grid cell currently active
            if eco.abenv.active[x, y]

                currentabun = @view eco.abundances.matrix[:, i]

                # Calculate effective rates
                birthprob = params.birth[j] * timestep * adjusted_birth
                deathprob = params.death[j] * timestep * adjusted_death

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                # Calculate how many births and deaths
                births = rand(rng, Poisson(currentabun[j] * newbirthprob))
                deaths = rand(rng, Binomial(currentabun[j], newdeathprob))

                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Calculate moves and write to cache
                move!(eco, eco.spplist.movement, i, j, eco.cache.netmigration, births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.matrix .+= eco.cache.netmigration

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
end

function update!(eco::Ecosystem, timestep::Unitful.Time, ::Val{1})

    # Calculate dimenions of habitat and number of species
    dims = _countsubcommunities(eco.abenv.habitat)
    spp = size(eco.abundances.grid,1)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall energy budget of that square
    update_energy_usage!(eco, Val{1}())

    # Loop through species in chosen square
    for j in 1:spp
        rng = eco.abundances.seed[Threads.threadid()]
        # Loop through grid squares
        for i in 1:dims
            # Calculate how much birth and death should be adjusted
            adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, i, j)

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(i, width)
            # Check if grid cell currently active
            if eco.abenv.active[x, y]

                currentabun = @view eco.abundances.matrix[:, i]

                # Calculate effective rates
                birthprob = params.birth[j] * timestep * adjusted_birth
                deathprob = params.death[j] * timestep * adjusted_death

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                # Calculate how many births and deaths
                births = rand(rng, Poisson(currentabun[j] * newbirthprob))
                #births = rand(rng, NegativeBinomial(currentabun[j]*(1 - newbirthprob)/newbirthprob))
                deaths = rand(rng, Binomial(currentabun[j], newdeathprob))
                #deaths = rand(rng, NegativeBinomial(currentabun[j]*(1 − newdeathprob)/newdeathprob))
                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Calculate moves and write to cache
                move!(eco, eco.spplist.movement, i, j, eco.cache.netmigration, births)
            end
        end
    end

    # Update abundances with all movements
    eco.abundances.matrix .+= eco.cache.netmigration

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
end

GLOBAL_funcdict["update!"] = update!

"""
    update_energy_usage!(eco::Ecosystem)
Function to calculate how much energy has been used up by the current species in each grid square in the ecosystem, `eco`. This function is parameterised on whether the species have one type of energy requirement or two.
"""
function update_energy_usage!(eco::Ecosystem{A, SpeciesList{Tr,  Req, B, C, D}, E}, ::Val{N}) where {N, A, B, C, D, E, Tr, Req <: Abstract1Requirement}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.requirement.energy

    # Loop through grid squares
    Threads.@threads for i in 1:size(eco.abundances.matrix, 2)
        eco.cache.totalE[i, 1] = ((@view eco.abundances.matrix[:, i]) ⋅ ϵ̄) * eco.spplist.requirement.exchange_rate
    end
    eco.cache.valid = true
end

function update_energy_usage!(eco::Ecosystem{A, SpeciesList{Tr,  Req, B, C, D}, E}, ::Val{N}) where {N, A, B, C, D, E, Tr, Req <: Abstract2Requirements}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.r1.energy
    ϵ̄2 = eco.spplist.requirement.r2.energy

    # Loop through grid squares
    Threads.@threads for i in 1:size(eco.abundances.matrix, 2)
        currentabun = @view eco.abundances.matrix[:, i]
        eco.cache.totalE[i, 1] = (currentabun ⋅ ϵ̄1) * eco.spplist.requirement.r1.exchange_rate
        eco.cache.totalE[i, 2] = (currentabun ⋅ ϵ̄2) * eco.spplist.requirement.r2.exchange_rate
    end
    eco.cache.valid = true
end

function update_energy_usage!(eco::Ecosystem{A, SpeciesList{Tr,  Req, B, C, D}, E}, ::Val{1}) where {A, B, C, D, E, Tr, Req <: Abstract1Requirement}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.requirement.energy

    # Loop through grid squares
    for i in 1:size(eco.abundances.matrix, 2)
        eco.cache.totalE[i, 1] = ((@view eco.abundances.matrix[:, i]) ⋅ ϵ̄) * eco.spplist.requirement.exchange_rate
    end
    eco.cache.valid = true
end

function update_energy_usage!(eco::Ecosystem{A, SpeciesList{Tr,  Req, B, C, D}, E}, ::Val{1}) where {A, B, C, D, E, Tr, Req <: Abstract2Requirements}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.r1.energy
    ϵ̄2 = eco.spplist.requirement.r2.energy

    # Loop through grid squares
    for i in 1:size(eco.abundances.matrix, 2)
        currentabun = @view eco.abundances.matrix[:, i]
        eco.cache.totalE[i, 1] = (currentabun ⋅ ϵ̄1) * eco.spplist.requirement.r1.exchange_rate
        eco.cache.totalE[i, 2] = (currentabun ⋅ ϵ̄2) * eco.spplist.requirement.r2.exchange_rate
    end
    eco.cache.valid = true
end

"""
    energy_adjustment(eco::Ecosystem, bud::AbstractBudget, i::Int64, spp::Int64)
Function to calculate how much birth and death rates should be adjusted by, according to how much energy is available, `bud`, in the grid square, `i`, and how much energy the species, `spp`, requires.
"""
function energy_adjustment(eco::Ecosystem, bud::AbstractBudget, i::Int64, spp::Int64)
    if typeof(eco.spplist.params) <: NoGrowth
        return 0.0, 0.0
    else
        width = getdimension(eco)[1]
        (x, y) = convert_coords(i, width)
        params = eco.spplist.params
        K = getbudget(eco)[x, y] * eco.spplist.requirement.exchange_rate
        # Get energy budgets of species in square
        ϵ̄ = eco.spplist.requirement.energy[spp] * eco.spplist.requirement.exchange_rate
        E = eco.cache.totalE[i, 1]
        # Traits
        ϵ̄real = ϵ̄/traitfun(eco, i, spp)
        # Alter rates by energy available in current pop & own requirements
        birth_energy = ϵ̄^-params.longevity * ϵ̄real^-params.survival * min(K/E, params.boost)
        death_energy = ϵ̄^-params.longevity * ϵ̄real^params.survival * (E / K)
    end
    return birth_energy, death_energy
end

function energy_adjustment(eco::Ecosystem, bud::BudgetCollection2, i::Int64, spp::Int64)
     width = getdimension(eco)[1]
     (x, y) = convert_coords(i, width)
     params = eco.spplist.params
    K1 = _getbudget(eco.abenv.budget, :b1)[x, y] * eco.spplist.requirement.r1.exchange_rate
    K2 = _getbudget(eco.abenv.budget, :b2)[x, y] * eco.spplist.requirement.r2.exchange_rate
    # Get abundances of square we are interested in
    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.r1.energy[spp] * eco.spplist.requirement.r1.exchange_rate
    E1 =  eco.cache.totalE[i, 1]
    ϵ̄2 = eco.spplist.requirement.r2.energy[spp] * eco.spplist.requirement.r2.exchange_rate
    E2 =  eco.cache.totalE[i, 2]
    ϵ̄real1 = ϵ̄1/traitfun(eco, i, spp)
    ϵ̄real2 = ϵ̄2/traitfun(eco, i, spp)
    # Alter rates by energy available in current pop & own requirements
    birth_energy = (ϵ̄1 * ϵ̄2)^-params.longevity * (ϵ̄real1 * ϵ̄real2)^-params.survival * min(K1/E1, K2/E2, params.boost)
    death_energy = (ϵ̄1 * ϵ̄2)^-params.longevity * (ϵ̄real1 * ϵ̄real2)^params.survival * max(E1/K1, E2/K2)
    return birth_energy, death_energy
end

"""
    convert_coords(i::Int64, width::Int64)
    convert_coords(x::Int64, y::Int64, width::Int64)
Function to convert coordinates from two-dimensional (`x`,`y`) format to one dimension (`i`), or vice versa, using the `width` of the grid. This function can also be applied to arrays of coordinates.
"""
function convert_coords(i::Int64, width::Int64)
  x = ((i - 1) % width) + 1
  y = div((i - 1), width)  + 1
  return (x, y)
end
function convert_coords(i::Array{Int64, 1}, width::Int64)
  x = ((i .- 1) .% width) .+ 1
  y = div.((i .- 1), width)  .+ 1
  return (x, y)
end
function convert_coords(x::Int64, y::Int64, width::Int64)
  i = x + width * (y - 1)
  return i
end

function convert_coords(x::Array{Int64, 1}, y::Array{Int64, 1}, width::Int64)
  i = x .+ (width .* (y .- 1))
  return i
end

"""
    calc_lookup_moves!(bound, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)

Function to calculate the number of moves taken by a species, `spp`, from a specific grid square location (`x`, `y`). There is a boundary condition, `bound`, which determines how the species can move across space (see AbstractBoundary). The total abundance of individuals is given in `abun`, which may be the number of births in the timestep, or total indiviuals.
"""
function calc_lookup_moves!(bound::NoBoundary, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
    lookup = eco.lookup[spp]
    maxX = Simulation.getdimension(eco)[1] - x
    maxY = Simulation.getdimension(eco)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        valid =  (-x < lookup.x[i] <= maxX) && (-y < lookup.y[i] <= maxY) && (eco.abenv.active[lookup.x[i], lookup.y[i]])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(eco.abundances.seed[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Cylinder, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
    lookup = eco.lookup[spp]
    maxX = Simulation.getdimension(eco)[1] - x
    maxY = Simulation.getdimension(eco)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, Simulation.getdimension(eco)[1]) + 1

        valid =  (-y < lookup.y[i] <= maxY) && (eco.abenv.active[newx, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(eco.abundances.seed[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Torus, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  lookup = eco.lookup[spp]
  maxX = Simulation.getdimension(eco)[1] - x
  maxY = Simulation.getdimension(eco)[2] - y
  # Can't go over maximum dimension
  for i in eachindex(lookup.x)
      newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, Simulation.getdimension(eco)[1]) + 1
      newy =  -y < lookup.y[i] <= maxY ? lookup.y[i] + y : mod(lookup.y[i] + y - 1, Simulation.getdimension(eco)[2]) + 1
      valid = eco.abenv.active[newx, newy]

      lookup.pnew[i] = valid ? lookup.p[i] : 0.0
  end
  lookup.pnew ./= sum(lookup.pnew)
  dist = Multinomial(abun, lookup.pnew)
  rand!(eco.abundances.seed[Threads.threadid()], dist, lookup.moves)
end

"""
    move!(eco::Ecosystem, ::AbstractMovement, i::Int64, spp::Int64, grd::Array{Int64, 2}, abun::Int64)

Function to calculate the movement of species `spp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a cached grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(eco::Ecosystem, ::AlwaysMovement, i::Int64, spp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  width, height = getdimension(eco)
  (x, y) = convert_coords(i, width)
  lookup = eco.lookup[spp]
  full_abun = eco.abundances.matrix[spp, i]
  calc_lookup_moves!(getboundary(eco.spplist.movement), x, y, spp, eco, full_abun)
  # Lose moves from current grid square
  grd[spp, i] -= full_abun
  # Map moves to location in grid
  mov = lookup.moves
  for i in eachindex(eco.lookup[spp].x)
      newx = mod(lookup.x[i] + x - 1, width) + 1
      newy = mod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(newx, newy, width)
      grd[spp, loc] += mov[i]
  end
  return eco
end

function move!(eco::Ecosystem, ::NoMovement, i::Int64, spp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  return eco
end

function move!(eco::Ecosystem, ::BirthOnlyMovement, i::Int64, spp::Int64,
    grd::Array{Int64, 2}, births::Int64)
  width, height = getdimension(eco)
  (x, y) = convert_coords(i, width)
  lookup = eco.lookup[spp]
  calc_lookup_moves!(getboundary(eco.spplist.movement), x, y, spp, eco, births)
  # Lose moves from current grid square
  grd[spp, i] -= births
  # Map moves to location in grid
  mov = lookup.moves
  for i in eachindex(lookup.x)
      newx = mod(lookup.x[i] + x - 1, width) + 1
      newy = mod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(newx, newy, width)
      grd[spp, loc] += mov[i]
  end
  return eco
end


function populate!(ml::GridLandscape, spplist::SpeciesList, abenv::AbstractAbiotic)
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b = reshape(ustrip.(_getbudget(abenv.budget)), size(grid))
    activity = reshape(copy(abenv.active), size(grid))
    b[.!activity] .= 0.0
    B = b./sum(b)
    # Loop through species
    for i in eachindex(spplist.abun)
        rand!(Multinomial(spplist.abun[i], B), (@view ml.matrix[i, :]))
    end
end

"""
    populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)

Function to populate a grid landscape given the abundances found in species list
and whether or not to include traits.
"""
function _populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic)
  # Calculate size of habitat
 dim = _getdimension(abenv.habitat)
 len = dim[1] * dim[2]
  grid = collect(1:len)
  # Set up copy of budget
  b = reshape(copy(_getbudget(abenv.budget)), size(grid))
  units = unit(b[1])
  activity = reshape(copy(abenv.active), size(grid))
  b[.!activity] .= 0.0 * units
  # Loop through species
  for i in eachindex(spplist.abun)
      # Get abundance of species
      abun = spplist.abun[i]
      # Loop through individuals
      while abun>0
          # Randomly choose position on grid (weighted)
          pos = sample(grid[b .> (0 * units)])
          # Add individual to this location
          ml.matrix[i, pos] = ml.matrix[i, pos] .+ 1
          abun = abun .- 1
          b[pos] = b[pos] .- spplist.requirement.energy[i]
      end
  end
end
function populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::GridAbioticEnv{H, BudgetCollection2{B1, B2}} where {H <:
                    AbstractHabitat, B1 <: AbstractBudget, B2 <: AbstractBudget})
    # Calculate size of habitat
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b1 = reshape(copy(_getbudget(abenv.budget, :b1)), size(grid))
    b2 = reshape(copy(_getbudget(abenv.budget, :b2)), size(grid))
    units1 = unit(b1[1])
    units2 = unit(b2[1])
    activity = reshape(copy(abenv.active), size(grid))
    b1[.!activity] .= 0.0 * units1
    b2[.!activity] .= 0.0 * units2
    # Loop through species
    for i in eachindex(spplist.abun)
        # Get abundance of species
        abun = spplist.abun[i]
        # Loop through individuals
        while abun>0
            # Randomly choose position on grid (weighted)
            pos = sample(grid[(b1 .> (0 * units1)) .& (b2 .> (0 * units2))])
            # Add individual to this location
            ml.matrix[i, pos] = ml.matrix[i, pos] .+ 1
            abun = abun .- 1
            b1[pos] = b1[pos] .- spplist.requirement.r1.energy[i]
            b2[pos] = b2[pos] .- spplist.requirement.r2.energy[i]
        end
    end
end
GLOBAL_funcdict["populate!"] = populate!


"""
    populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)

Function to populate a grid landscape given the abundances found in species list
and whether or not to include traits.
"""
function simplepopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AB, rel::R) where {AB <: AbstractAbiotic, R <: AbstractTraitRelationship}
  # Calculate size of habitat
  dim = _getdimension(abenv.habitat)
  numsquares = dim[1] * dim[2]
  numspp = length(spplist.names)
  maxrng = spplist.traits.mean .+ spplist.traits.var
  minrng = spplist.traits.mean .- spplist.traits.var
  hab = reshape(abenv.habitat.matrix, numsquares)
  probabilities = [_traitfun(abenv.habitat, spplist.traits, rel, i, spp, iscontinuous(abenv.habitat)) for i in 1:numsquares, spp in 1:numspp]
  # Loop through species
  for i in eachindex(spplist.abun)
      if spplist.native[i]
        # Get abundance of species
        abun = rand(Multinomial(spplist.abun[i], probabilities[:, i]./sum(probabilities[:, i])))
        # Add individual to this location
        ml.matrix[i, :] .+= abun
     end
   end
end
function simplerepopulate!(eco::Ecosystem)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun), length(eco.spplist.abun)))
  simplepopulate!(eco.abundances, eco.spplist, eco.abenv, eco.relationship)
end

GLOBAL_funcdict["simplepopulate!"] = simplepopulate!

"""
    repopulate!(eco::Ecosystem, abun::Int64)
Function to repopulate an ecosystem `eco`, with option for including trait
preferences. An additional `abun` parameter can be included, in order to repopulate the ecosystem with a specified number of individuals.
"""
function repopulate!(eco::Ecosystem)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun), length(eco.spplist.abun)))
  populate!(eco.abundances, eco.spplist, eco.abenv)
end
function repopulate!(eco::Ecosystem, abun::Int64)
    dim = _getdimension(eco.abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b = reshape(copy(_getbudget(eco.abenv.budget)), size(grid))
    units = unit(b[1])
    activity = reshape(copy(eco.abenv.active), size(grid))
    b[.!activity] .= 0.0 * units
    pos = sample(grid[b .> (0 * units)], abun)
    # Add individual to this location
    map(pos) do p
        eco.abundances.matrix[end, p] += 1
    end
end
GLOBAL_funcdict["repopulate!"] = repopulate!

"""
    reenergise!(eco::Ecosystem, budget::Union{Float64, Unitful.Quantity{Float64}}, grid::Tuple{Int64, Int64})
Function to refill an ecosystem `eco`, with energy from a budget value, `budget` and a grid size.
"""
function reenergise!(eco::Ecosystem, budget::Union{Float64, Unitful.Quantity{Float64}}, grid::Tuple{Int64, Int64})
    fill!(eco.abenv.budget.matrix, budget/(grid[1]*grid[2]))
end
function reenergise!(eco::Ecosystem, budget::Tuple{Unitful.Quantity{Float64}, Unitful.Quantity{Float64}}, grid::Tuple{Int64, Int64})
    fill!(eco.abenv.budget.b1.matrix, budget[1]/(grid[1]*grid[2]))
    fill!(eco.abenv.budget.b2.matrix, budget[2]/(grid[1]*grid[2]))
end
GLOBAL_funcdict["reenergise!"] = reenergise!
