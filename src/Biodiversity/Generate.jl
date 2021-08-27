using StatsBase
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

"""
    biodiversity_update!(eco::Ecosystem, time::Unitful.Time)
Function to update a ecosystem abundances and environment for one timestep.
"""
function biodiversity_update!(eco::Ecosystem, timestep::Unitful.Time)

    # Calculate dimenions of habitat and number of species
    dims = _countsubcommunities(eco.abenv.habitat)
    spp = size(eco.abundances.grid,1)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall energy budget of that square
    update_energy_usage!(eco)

    # Loop through species in chosen square
    Threads.@threads for j in 1:spp
        rng = eco.abundances.rngs[Threads.threadid()]
        # Loop through grid squares
        for i in 1:dims
            # Calculate how much birth and death should be adjusted
            adjusted_birth, adjusted_death = energy_adjustment(eco, eco.abenv.budget, i, j)

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(eco, i, width)
            # Check if grid cell currently active
            if eco.abenv.active[x, y] && (eco.cache.totalE[i, 1] > 0)
                # Calculate effective rates
                birthrate = params.birth[j] * timestep * adjusted_birth
                deathparam = params.death[j] * timestep * adjusted_death

                # Turn deathparam into probability and cancel units of birthrate
                birthrate += 0.0
                deathprob = 1.0 - exp(-deathparam)

                (birthrate >= 0) & (deathprob >= 0) || error("Birth: $birthrate \n Death: $deathprob \n \n i: $i \n j: $j")
                # Calculate how many births and deaths
                births = rand(rng, Poisson(eco.abundances.matrix[j, i] * birthrate))
                deaths = rand(rng, Binomial(eco.abundances.matrix[j, i], deathprob))

                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Calculate moves and write to cache
                move!(eco, eco.spplist.species.movement, i, j, eco.cache.netmigration, births)
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


"""
    update_energy_usage!(eco::Ecosystem)
Function to calculate how much energy has been used up by the current species in each grid square in the ecosystem, `eco`. This function is parameterised on whether the species have one type of energy requirement or two.
"""
function update_energy_usage!(eco::AbstractEcosystem{A, B, SpeciesList{SpeciesTypes{C, Req, D, E},  F, G}, H, I, J}) where {A, B, C, D, E, F, G, H, I, J, Req <: Abstract1Requirement}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.species.requirement.energy

    # Loop through grid squares
    Threads.@threads for i in 1:size(eco.abundances.matrix, 2)
        eco.cache.totalE[i, 1] = ((@view eco.abundances.matrix[:, i]) ⋅ ϵ̄) * eco.spplist.species.requirement.exchange_rate
    end
    eco.cache.valid = true
end

function update_energy_usage!(eco::AbstractEcosystem{A, B, SpeciesList{SpeciesTypes{C, Req, D, E},  F, G}, H, I, J}) where {A, B, C, D, E, F, G, H, I, J, Req <: Abstract2Requirements}
    !eco.cache.valid || return true

    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.species.requirement.r1.energy
    ϵ̄2 = eco.spplist.species.requirement.r2.energy

    # Loop through grid squares
    Threads.@threads for i in 1:size(eco.abundances.matrix, 2)
        currentabun = @view eco.abundances.matrix[:, i]
        eco.cache.totalE[i, 1] = (currentabun ⋅ ϵ̄1) * eco.spplist.species.requirement.r1.exchange_rate
        eco.cache.totalE[i, 2] = (currentabun ⋅ ϵ̄2) * eco.spplist.species.requirement.r2.exchange_rate
    end
    eco.cache.valid = true
end

"""
    energy_adjustment(eco::Ecosystem, bud::AbstractBudget, i::Int64, sp::Int64)
Function to calculate how much birth and death rates should be adjusted by, according to how much energy is available, `bud`, in the grid square, `i`, and how much energy the species, `sp`, requires.
"""
function energy_adjustment(eco::AbstractEcosystem, bud::AbstractBudget, i::Int64, sp::Int64)
    if typeof(eco.spplist.params) <: NoGrowth
        return 0.0, 0.0
    else
        width = getdimension(eco)[1]
        (x, y) = convert_coords(eco, i, width)
        params = eco.spplist.params
        K = getbudget(eco)[x, y] * eco.spplist.species.requirement.exchange_rate
        # Get energy budgets of species in square
        ϵ̄ = eco.spplist.species.requirement.energy[sp] * eco.spplist.species.requirement.exchange_rate
        E = eco.cache.totalE[i, 1]
        # Traits
        ϵ̄real = 1/traitfun(eco, i, sp, eco.spplist.species)
        # Alter rates by energy available in current pop & own requirements
        birth_energy = ϵ̄^-params.longevity * ϵ̄real^-params.survival * min(K/E, params.boost)
        death_energy = ϵ̄^-params.longevity * ϵ̄real^params.survival * (E / K)
    end
    return birth_energy, death_energy
end

function energy_adjustment(eco::AbstractEcosystem, bud::BudgetCollection2, i::Int64, sp::Int64)
    width = getdimension(eco)[1]
    (x, y) = convert_coords(eco, i, width)
    params = eco.spplist.params
    K1 = _getbudget(eco.abenv.budget.b1)[x, y] * eco.spplist.species.requirement.r1.exchange_rate
    K2 = _getbudget(eco.abenv.budget.b2)[x, y] * eco.spplist.species.requirement.r2.exchange_rate
    # Get abundances of square we are interested in
    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.species.requirement.r1.energy[sp] * eco.spplist.species.requirement.r1.exchange_rate
    ϵ̄2 = eco.spplist.species.requirement.r2.energy[sp] * eco.spplist.species.requirement.r2.exchange_rate
    E1 = eco.cache.totalE[i, 1]
    E2 = eco.cache.totalE[i, 2]
    ϵ̄real1 = 1/traitfun(eco, i, sp, eco.spplist.species)
    ϵ̄real2 = 1/traitfun(eco, i, sp, eco.spplist.species)
    # Alter rates by energy available in current pop & own requirements
    birth_energy = (ϵ̄1 * ϵ̄2)^-params.longevity * (ϵ̄real1 * ϵ̄real2)^-params.survival * min(K1/E1, K2/E2, params.boost)
    death_energy = (ϵ̄1 * ϵ̄2)^-params.longevity * (ϵ̄real1 * ϵ̄real2)^params.survival * max(E1/K1, E2/K2)
    return birth_energy, death_energy
end

"""
    convert_coords(eco, i::Int64, width::Int64)
    convert_coords(eco, x::Int64, y::Int64, width::Int64)
Function to convert coordinates from two-dimensional (`x`,`y`) format to one dimension (`i`), or vice versa, using the `width` of the grid. This function can also be applied to arrays of coordinates.
"""
function convert_coords(eco::AbstractEcosystem, i::Int64, width::Int64 = getdimension(eco)[1])
    x = ((i - 1) % width) + 1
    y = div((i - 1), width)  + 1
    return (x, y)
end
function convert_coords(eco::AbstractEcosystem, pos::Tuple{Int64, Int64}, width::Int64 = getdimension(eco)[1])
    i = pos[1] + width * (pos[2] - 1)
    return i
end

function convert_coords(i::Int64, width::Int64)
  x = ((i - 1) % width) + 1
  y = div((i - 1), width)  + 1
  return (x, y)
end
function convert_coords(x::Int64, y::Int64, width::Int64)
  i = x + width * (y - 1)
  return i
end
"""
    calc_lookup_moves!(bound, x::Int64, y::Int64, sp::Int64, eco::Ecosystem, abun::Int64)

Function to calculate the number of moves taken by a species, `sp`, from a specific grid square location (`x`, `y`). There is a boundary condition, `bound`, which determines how the species can move across space (see AbstractBoundary). The total abundance of individuals is given in `abun`, which may be the number of births in the timestep, or total indiviuals.
"""
function calc_lookup_moves!(bound::NoBoundary, x::Int64, y::Int64, sp::Int64, eco::AbstractEcosystem, abun::Int64)
    lookup = getlookup(eco, sp)
    maxX = getdimension(eco)[1] - x
    maxY = getdimension(eco)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        valid =  (-x < lookup.x[i] <= maxX) && (-y < lookup.y[i] <= maxY) && (eco.abenv.active[lookup.x[i] + x, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(eco.abundances.rngs[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Cylinder, x::Int64, y::Int64, sp::Int64, eco::AbstractEcosystem, abun::Int64)
    lookup = getlookup(eco, sp)
    maxX = getdimension(eco)[1] - x
    maxY = getdimension(eco)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, getdimension(eco)[1]) + 1

        valid =  (-y < lookup.y[i] <= maxY) && (eco.abenv.active[newx, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(eco.abundances.rngs[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Torus, x::Int64, y::Int64, sp::Int64, eco::AbstractEcosystem, abun::Int64)
  lookup = getlookup(eco, sp)
  maxX = getdimension(eco)[1] - x
  maxY = getdimension(eco)[2] - y
  # Can't go over maximum dimension
  for i in eachindex(lookup.x)
      newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, getdimension(eco)[1]) + 1
      newy =  -y < lookup.y[i] <= maxY ? lookup.y[i] + y : mod(lookup.y[i] + y - 1, getdimension(eco)[2]) + 1
      valid = eco.abenv.active[newx, newy]

      lookup.pnew[i] = valid ? lookup.p[i] : 0.0
  end
  lookup.pnew ./= sum(lookup.pnew)
  dist = Multinomial(abun, lookup.pnew)
  rand!(eco.abundances.rngs[Threads.threadid()], dist, lookup.moves)
end

"""
    move!(eco::Ecosystem, ::AbstractMovement, i::Int64, sp::Int64, grd::Array{Int64, 2}, abun::Int64)

Function to calculate the movement of species `sp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a cached grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(eco::AbstractEcosystem, ::AlwaysMovement, i::Int64, sp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  width, height = getdimension(eco)
  (x, y) = convert_coords(eco, i, width)
  lookup = getlookup(eco, sp)
  full_abun = eco.abundances.matrix[sp, i]
  calc_lookup_moves!(getboundary(eco.spplist.species.movement), x, y, sp, eco, full_abun)
  # Lose moves from current grid square
  grd[sp, i] -= full_abun
  # Map moves to location in grid
  mov = lookup.moves
  for i in eachindex(lookup.x)
      newx = mod(lookup.x[i] + x - 1, width) + 1
      newy = mod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(eco, (newx, newy), width)
      grd[sp, loc] += mov[i]
  end
  return eco
end

function move!(eco::AbstractEcosystem, ::NoMovement, i::Int64, sp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  return eco
end

function move!(eco::AbstractEcosystem, ::BirthOnlyMovement, i::Int64, sp::Int64,
    grd::Array{Int64, 2}, births::Int64)
  width, height = getdimension(eco)
  (x, y) = convert_coords(eco, i, width)
   lookup = getlookup(eco, sp)
  calc_lookup_moves!(getboundary(eco.spplist.species.movement), x, y, sp, eco, births)
  # Lose moves from current grid square
  grd[sp, i] -= births
  # Map moves to location in grid
  mov = lookup.moves
  for i in eachindex(lookup.x)
      newx = mod(lookup.x[i] + x - 1, width) + 1
      newy = mod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(eco, (newx, newy), width)
      grd[sp, loc] += mov[i]
  end
  return eco
end

"""
    populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)

Function to populate a grid landscape given the abundances found in species list according to availability of resources.
"""

function populate!(ml::GridLandscape, spplist::SpeciesList, abenv::AB, rel::R) where {AB <: AbstractAbiotic, R <: AbstractTraitRelationship}
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b = reshape(ustrip.(_getbudget(abenv.budget)), size(grid))
    activity = reshape(copy(abenv.active), size(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    B = b./sum(b)
    # Loop through species
    for i in eachindex(spplist.species.abun)
        rand!(Multinomial(spplist.species.abun[i], B), (@view ml.matrix[i, :]))
    end
end

function populate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::GridAbioticEnv{H, BudgetCollection2{B1, B2}}, rel::R) where {H <: AbstractHabitat, B1 <: AbstractBudget, B2 <: AbstractBudget, R <: AbstractTraitRelationship}
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
    B = (b1./sum(b1)) .* (b2./sum(b2))
    # Loop through species
    for i in eachindex(spplist.species.abun)
        rand!(Multinomial(spplist.species.abun[i], B./sum(B)), (@view ml.matrix[i, :]))
    end
end

"""
    repopulate!(eco::Ecosystem, abun::Int64)
Function to repopulate an ecosystem `eco`, with option for including trait preferences. An additional `abun` parameter can be included, in order to repopulate the ecosystem with a specified number of individuals.
"""
function repopulate!(eco::Ecosystem)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.species.abun = rand(Multinomial(sum(eco.spplist.species.abun), length(eco.spplist.species.abun)))
  populate!(eco.abundances, eco.spplist, eco.abenv, eco.relationship)
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


"""
    traitpopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic)

Function to populate a grid landscape given the abundances found in species list based upon how well the species traits match their environment.
"""
function traitpopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AB, rel::R) where {AB <: AbstractAbiotic, R <: AbstractTraitRelationship}
  # Calculate size of habitat
  dim = _getdimension(abenv.habitat)
  numsquares = dim[1] * dim[2]
  numspp = length(spplist.species.names)
  maxrng = spplist.species.traits.mean .+ spplist.species.traits.var
  minrng = spplist.species.traits.mean .- spplist.species.traits.var
  hab = reshape(abenv.habitat.matrix, numsquares)
  probabilities = [_traitfun(abenv.habitat, spplist.species.traits, rel, i, sp) for i in 1:numsquares, sp in 1:numspp]
  # Loop through species
  for i in eachindex(spplist.species.abun)
      if spplist.species.native[i]
        # Get abundance of species
        probs = probabilities[:, i]./sum(probabilities[:, i])
        probs[isnan.(probs)] .= 1/numsquares
        abun = rand(Multinomial(spplist.species.abun[i], probs))
        # Add individual to this location
        ml.matrix[i, :] .+= abun
     end
   end
end

"""
    repopulate!(eco::Ecosystem, abun::Int64)
Function to repopulate an ecosystem `eco`, with option for including trait preferences. An additional `abun` parameter can be included, in order to repopulate the ecosystem with a specified number of individuals.
"""
function traitrepopulate!(eco::Ecosystem)
  eco.abundances = emptygridlandscape(eco.abenv, eco.spplist)
  eco.spplist.species.abun = rand(Multinomial(sum(eco.spplist.species.abun), length(eco.spplist.species.abun)))
  traitpopulate!(eco.abundances, eco.spplist, eco.abenv, eco.relationship)
end

"""
    emptypopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AB, rel::R) where {AB <: EcoSISTEM.AbstractAbiotic, R <: EcoSISTEM.AbstractTraitRelationship}
"""
function emptypopulate!(ml::GridLandscape, spplist::SpeciesList,
                   abenv::AB, rel::R) where {AB <: EcoSISTEM.AbstractAbiotic, R <: EcoSISTEM.AbstractTraitRelationship}
@warn "Ecosystem not populated!"
end
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
