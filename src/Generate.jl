using StatsBase
using Query
using Compat
using LinearAlgebra

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
    dims = _countsubcommunities(eco.abenv.habitat)
    spp = size(eco.abundances.grid,1)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    # Set the overall energy budget of that square
    update_energy_usage!(eco)

    # Loop through species in chosen square
    Threads.@threads for j in 1:spp
    # Loop through grid squares
        for i in 1:dims
            (x, y) = convert_coords(i, width)
            birth_energy, death_energy = energy(eco, eco.abenv.budget, i,
                spp)
            if eco.abenv.active[x, y] && any((@view eco.abundances.matrix[:, i]) .≠ 0)
                currentabun = @view eco.abundances.matrix[:, i]
                # Calculate effective rates
                birthprob = params.birth[j] * timestep * birth_energy
                deathprob = params.death[j] * timestep * death_energy

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                # Calculate how many births and deaths
                births = rand(Binomial(currentabun[j], newbirthprob))
                deaths = rand(Binomial(currentabun[j], newdeathprob))
                # Update population
                eco.abundances.matrix[j, i] += (births - deaths)

                # Perform gaussian movement
                move!(eco, eco.spplist.movement, i, j, eco.cache.netmigration, births)
            end
        end
    end
    eco.abundances.matrix .+= eco.cache.netmigration
    invalidatecaches!(eco)
    # Update environment
    habitatupdate!(eco, timestep)
    budgetupdate!(eco, timestep)
end
GLOBAL_funcdict["update!"] = update!

function update_energy_usage!(eco::Ecosystem{A, SpeciesList{Tr,  Req, B, C, D}, E}) where {A,B, C, D, E, Tr, Req <: Abstract1Requirement}
    !eco.cache.valid || return true
    # Get energy budgets of species in square
    ϵ̄ = eco.spplist.requirement.energy
    # Loop through grid squares
    Threads.@threads for i in 1:size(eco.abundances.matrix, 2)
        eco.cache.totalE[i, 1] = ((@view eco.abundances.matrix[:, i]) ⋅ ϵ̄) * eco.spplist.requirement.exchange_rate
    end
    eco.cache.valid = true
end
function update_energy_usage!(eco::Ecosystem{A, SpeciesList{Tr,  Req}}) where {A, Tr, Req <: Abstract2Requirements}
    !eco.cache.valid || return true
    # Get energy budgets of species in square
    ϵ̄1 = eco.spplist.requirement.r1.energy
    ϵ̄2 = eco.spplist.requirement.r2.energy

    Threads.@threads for i in 1:size(eco.abundances.matrix, 2)
        currentabun = @view eco.abundances.matrix[:, i]
        eco.cache.totalE[i, 1] = (currentabun ⋅ ϵ̄1) * eco.spplist.requirement.r1.exchange_rate
        eco.cache.totalE[i, 2] = (currentabun ⋅ ϵ̄2) * eco.spplist.requirement.r2.exchange_rate
    end
    eco.cache.valid = true
end

function energy(eco::Ecosystem, bud::AbstractBudget, i::Int64, spp::Int64)
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

function energy(eco::Ecosystem, bud::BudgetCollection2, i::Int64, spp::Int64)
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
    birth_energy = (ϵ̄1^-params.longevity + ϵ̄2^-params.longevity) * (ϵ̄real1^-params.survival +
     ϵ̄real2^-params.survival) * min(K1/E1, K2/E2, params.boost)
    death_energy = (ϵ̄1^-params.longevity + ϵ̄2^-params.longevity) * (ϵ̄real1^params.survival +
    ϵ̄real2^params.survival) * max(E1 / K1, E2/K2)
    return birth_energy, death_energy
end


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
function calc_lookup_moves(bound::NoBoundary, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  lookup = eco.lookup[spp]
  maxX = getdimension(eco)[1] - x
  maxY = getdimension(eco)[2] - y
  # Can't go over maximum dimension
  valid = (lookup.x .> -x) .& (lookup.y .> -y) .&
   (lookup.x .<= maxX) .& (lookup.y .<= maxY)
  for i in eachindex(valid)
      if valid[i]
          valid[i] = valid[i] & (eco.abenv.active[lookup.x[i] .+ x,
          lookup.y[i] .+ y])
      end
  end
  lookup.pnew[.!valid] .= 0.0
  lookup.pnew[valid] = lookup.p[valid]
  lookup.pnew ./= sum(lookup.pnew)
  multinom = Multinomial(abun, lookup.pnew)
  draw = rand(multinom)
  lookup.moves .= draw
end

function calc_lookup_moves(bound::Cylinder, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  lookup = eco.lookup[spp]
  maxX = getdimension(eco)[1] - x
  maxY = getdimension(eco)[2] - y
  lookup.x[lookup.x .<= -x] .+= Simulation.getdimension(eco)[1]
  lookup.x[lookup.x .> maxX] .-= Simulation.getdimension(eco)[1]

  valid = (lookup.y .> -y) .& (lookup.y .<= maxY)
  for i in eachindex(lookup.x)
      if valid[i]
          valid[i] = valid[i] & (eco.abenv.active[lookup.x[i] .+ x,
            lookup.y[i] .+ y])
      end
  end

  lookup.pnew[.!valid] .= 0.0
  lookup.pnew[valid] = lookup.p[valid]
  lookup.pnew ./= sum(lookup.pnew)
  multinom = Multinomial(abun, lookup.pnew)
  draw = rand(multinom)
  lookup.moves .= draw
end

function calc_lookup_moves(bound::Torus, x::Int64, y::Int64, spp::Int64, eco::Ecosystem, abun::Int64)
  lookup = eco.lookup[spp]
  maxX = Simulation.getdimension(eco)[1] - x
  maxY = Simulation.getdimension(eco)[2] - y
  # Can't go over maximum dimension
  for i in eachindex(x)
      newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, Simulation.getdimension(eco)[1]) + 1
      newy =  -y < lookup.y[i] <= maxY ? lookup.y[i] + y : mod(lookup.y[i] + y - 1, Simulation.getdimension(eco)[2]) + 1

      valid = eco.abenv.active[newx, newy]

      lookup.pnew[i] = valid ? lookup.p[i] : 0.0
  end
  lookup.pnew ./= sum(lookup.pnew)
  multinom = Multinomial(abun, lookup.pnew)
  draw = rand(multinom)
  lookup.moves .= draw
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
  width = getdimension(eco)[1]
  (x, y) = convert_coords(i, width)
  full_abun = eco.abundances.matrix[spp, i]
  calc_lookup_moves(getboundary(eco.spplist.movement), x, y, spp, eco, full_abun)
  # Lose moves from current grid square
  grd[spp, i] -= full_abun
  # Map moves to location in grid
  mov = eco.lookup[spp].moves
  for i in eachindex(eco.lookup[spp].x)
      loc = convert_coords(eco.lookup[spp].x[i] + x, eco.lookup[spp].y[i] + y, width)
      grd[spp, loc] += mov[i]
  end
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


function move!(eco::Ecosystem, ::BirthOnlyMovement, i::Int64, spp::Int64,
    grd::Array{Int64, 2}, births::Int64)
  width = getdimension(eco)[1]
  (x, y) = convert_coords(i, width)
  table = calc_lookup_moves(getboundary(eco.spplist.movement), x, y, spp, eco, births)
  # Lose moves from current grid square
  grd[spp, i] -= births
  # Map moves to location in grid
  mov = eco.lookup[spp].moves
  for i in eachindex(eco.lookup[spp].x)
      loc = convert_coords(eco.lookup[spp].x[i] + x, eco.lookup[spp].y[i] + y, width)
      grd[spp, loc] += mov[i]
  end
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
    b1[.!activity] = 0.0 * units1
    b2[.!activity] = 0.0 * units2
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

function simplepopulate!(ml::GridLandscape, spplist::SpeciesList,
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
      if spplist.native[i]
        # Get abundance of species
        abun = spplist.abun[i]
        pos = sample(grid[b .> (0 * units)], abun)
        # Add individual to this location
        ml.matrix[i, pos] = ml.matrix[i, pos] .+ 1
     end
   end
end
GLOBAL_funcdict["simplepopulate!"] = simplepopulate!

"""
    repopulate!(eco::Ecosystem)
Function to repopulate an ecosystem `eco`, with option for including trait
preferences.
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
    reenergise!(eco::Ecosystem, budget::)
Function to repopulate an ecosystem `eco`, with option for including trait
preferences.
"""
function reenergise!(eco::Ecosystem, budget::Union{Float64, Unitful.Quantity{Float64}}, grid::Tuple{Int64, Int64})
    fill!(eco.abenv.budget.matrix, budget/(grid[1]*grid[2]))
end

GLOBAL_funcdict["reenergise!"] = reenergise!
