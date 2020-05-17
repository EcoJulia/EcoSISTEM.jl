using StatsBase
using Compat
using LinearAlgebra

"""
    update!(epi::EpiSystem, time::Unitful.Time)
Function to update disease class abundances and environment for one timestep.
"""
function update!(epi::EpiSystem, timestep::Unitful.Time)

    # Human movement loop - ignore for now while we focus on lockdown
    # humanmove!(epi, timestep)

    # Virus movement loop
    virusupdate!(epi, timestep)

    # Birth/death/infection/recovery loop of each class
    classupdate!(epi, timestep)

    # Invalidate all caches for next update
    invalidatecaches!(epi)

    # Update environment - habitat and energy budgets
    habitatupdate!(epi, timestep)
    applycontrols!(epi, timestep)
end

function humanmove!(epi::EpiSystem, timestep::Unitful.Time)
    dims = _countsubcommunities(epi.epienv.habitat)
    width = getdimension(epi)[1]
    classes = size(abundances(epi).matrix, 1)
    Threads.@threads for j in 2:classes
        rng = abundances(epi).seed[Threads.threadid()]
        # Loop through grid squares
        for i in 1:dims
            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(epi, i, width)
            # Check if grid cell currently active
            if epi.epienv.active[x, y]
                # Calculate moves and write to cache
                move!(epi, epi.epilist.movement, i, j, epi.cache.netmigration, abundances(epi).matrix[j, i])
            end
        end
    end
    # Update abundances with all movements
    abundances(epi).matrix .+= epi.cache.netmigration
end

function virusupdate!(epi::EpiSystem, timestep::Unitful.Time)
    dims = _countsubcommunities(epi.epienv.habitat)
    width = getdimension(epi)[1]
    params = epi.epilist.params
    id = Threads.threadid()
    rng = abundances(epi).seed[id]
    classes = findall((params.virus_growth .* timestep) .> 0)
    # Loop through grid squares
    Threads.@threads for j in classes
        for i in 1:dims
            # Calculate how much birth and death should be adjusted

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(epi, i, width)
            # Check if grid cell currently active
            if epi.epienv.active[x, y]
                traitmatch = traitfun(epi, i, 1)
                # Calculate effective rates
                birthrate = params.virus_growth[j] * timestep * abundances(epi).matrix[j, i]
                deathrate = params.virus_decay[j] * timestep * traitmatch^-1

                # Convert death rate into 0 - 1 probability
                deathprob = 1.0 - exp(-deathrate)

                (birthrate >= 0) & (deathprob >= 0) || error("Birth: $birthprob \n Death: $newdeathprob \n \n i: $i")
                # Calculate how many births and deaths
                births = rand(rng, Poisson(birthrate))
                deaths = rand(rng, Binomial(abundances(epi).matrix[1, i], deathprob))

                # Update population
                epi.cache.virusmigration[j, i] += births
                epi.cache.virusdecay[j, i] -= deaths
                virusmove!(epi, i, j, epi.cache.virusmigration, births)
            end
        end
    end
    vm = sum(epi.cache.virusmigration, dims = 1)[1, :]
    nm = sum(epi.cache.virusdecay, dims = 1)[1, :]
    abundances(epi).matrix[1, :] .+= (nm .+ vm)
    epi.cache.virusmigration[1, :] .+= vm
end

function sum_pop(M::Matrix{Int64}, start::Int64, i::Int64)
    N = 0
    for j in start:size(M, 1)
        N += M[j, i]
    end
    return N
end
"""
    classupdate!(epi::EpiSystem, timestep::Unitful.Time)
Function to update disease class abundances for one timestep. Dispatches differently depending on the class of model stored in the EpiSystem.
"""
function classupdate!(epi::EpiSystem, timestep::Unitful.Time)
    # Calculate dimenions of habitat and number of classes
    dims = _countsubcommunities(epi.epienv.habitat)
    params = epi.epilist.params
    width = getdimension(epi)[1]
    classes = size(abundances(epi).matrix, 1)
    # Loop through grid squares
    Threads.@threads for i in 1:dims
        rng = abundances(epi).seed[Threads.threadid()]
        susclass = findfirst(epi.epilist.names .== "Susceptible")
        N = sum_pop(abundances(epi).matrix, susclass, i)
        # Loop through classes in chosen square
        for j in susclass:classes
            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(epi, i, width)
            # Check if grid cell currently active
            if epi.epienv.active[x, y]
                # Births
                births = rand(rng, Binomial(abundances(epi).matrix[j, i],  params.births[j] * timestep))
                abundances(epi).matrix[susclass, i] += births

                # Calculate force of inf and env inf
                env_inf = (params.transition_virus[j, :] .* timestep .* abundances(epi).matrix[1, i]) ./ N

                force_inf = (params.transition_force[j, :] .* timestep .* epi.cache.virusmigration[1, i]) ./ N

                # Add to transitional probabilities
                trans_prob = (params.transition[j, :] .* timestep) .+ env_inf .+  force_inf
                trans_prob = 1.0 .- exp.(-1 .* trans_prob)

                # Make transitions
                trans = rand.(fill(rng, length(trans_prob)), Binomial.(abundances(epi).matrix[:, i],  trans_prob))
                abundances(epi).matrix[j, i] += sum(trans)
                abundances(epi).matrix[:, i] .-= trans
            end
        end
    end
end

"""
    populate!(ml::EpiLandscape, epilist::EpiList, epienv::EE, rel::R)
Function to populate an EpiLandscape with information on each disease class in the EpiList.
"""
function populate!(ml::EpiLandscape, epilist::EpiList, epienv::EE, rel::R) where {EE <: AbstractEpiEnv, R <: AbstractTraitRelationship}
    dim = _getdimension(epienv.habitat)
    len = dim[1] * dim[2]
    # Loop through classes
    for i in eachindex(epilist.abun)
        rand!(Multinomial(epilist.abun[i], len), (@view ml.matrix[i, :]))
    end
end

"""
    applycontrols!(epi::EpiSystem, timestep::Unitful.Time)
Function to apply control strategies to an EpiSystem for one timestep.
"""
function applycontrols!(epi::AbstractEpiSystem, timestep::Unitful.Time)
    _applycontrols!(epi, epi.epienv.control, timestep)
end

function _applycontrols!(epi::AbstractEpiSystem, controls::NoControl, timestep::Unitful.Time)
    return controls
end

function convert_coords(epi::AbstractEpiSystem, i::Int64, width::Int64 = getdimension(epi)[1])
    x = ((i - 1) % width) + 1
    y = div((i - 1), width)  + 1
    return (x, y)
end
function convert_coords(epi::AbstractEpiSystem, pos::Tuple{Int64, Int64}, width::Int64 = getdimension(epi)[1])
    i = pos[1] + width * (pos[2] - 1)
    return i
end


function calc_lookup_moves!(bound::NoBoundary, x::Int64, y::Int64, id::Int64, epi::AbstractEpiSystem, abun::Int64)
    lookup = getlookup(epi, id)
    maxX = getdimension(epi)[1] - x
    maxY = getdimension(epi)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        valid =  (-x < lookup.x[i] <= maxX) && (-y < lookup.y[i] <= maxY) && (epi.epienv.active[lookup.x[i] + x, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(abundances(epi).seed[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Cylinder, x::Int64, y::Int64, id::Int64, epi::AbstractEpiSystem, abun::Int64)
    lookup = getlookup(epi, id)
    maxX = getdimension(epi)[1] - x
    maxY = getdimension(epi)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, getdimension(epi)[1]) + 1

        valid =  (-y < lookup.y[i] <= maxY) && (epi.epienv.active[newx, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(abundances(epi).seed[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Torus, x::Int64, y::Int64, id::Int64, epi::AbstractEpiSystem, abun::Int64)
    lookup = getlookup(epi, id)
    maxX = getdimension(epi)[1] - x
    maxY = getdimension(epi)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, getdimension(epi)[1]) + 1
        newy =  -y < lookup.y[i] <= maxY ? lookup.y[i] + y : mod(lookup.y[i] + y - 1, getdimension(epi)[2]) + 1
        valid = epi.epienv.active[newx, newy]

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(abundances(epi).seed[Threads.threadid()], dist, lookup.moves)
end

function virusmove!(epi::AbstractEpiSystem, pos::Int64, id::Int64, grd::Array{Int64, 2}, births::Int64)
  width, height = getdimension(epi)
  (x, y) = convert_coords(epi, pos, width)
  lookup = getlookup(epi, id)
  calc_lookup_moves!(getboundary(epi.epilist.movement), x, y, id, epi, births)
  # Lose moves from current grid square
  grd[id, pos] -= births
  # Map moves to location in grid
  mov = lookup.moves
  for i in eachindex(lookup.x)
      newx = mod(lookup.x[i] + x - 1, width) + 1
      newy = mod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(epi, (newx, newy), width)
      grd[id, loc] += mov[i]
  end
  return epi
end

"""
    move!(epi::AbstractEpiSystem, ::AlwaysMovement, i::Int64, sp::Int64, grd::Array{Int64, 2}, ::Int64)

Function to calculate the movement of a disease class `sp` from a given position in the landscape `i`, using the lookup table found in the EpiSystem and updating the movement patterns on a cached grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead of the entire population.
"""
function move!(epi::AbstractEpiSystem, ::AlwaysMovement, i::Int64, sp::Int64, grd::Array{Int64, 2}, ::Int64)
  width, height = getdimension(epi)
  (x, y) = convert_coords(epi, i, width)
  lookup = getlookup(epi, sp)
  full_abun = abundances(epi).matrix[sp, i]
  calc_lookup_moves!(getboundary(epi.epilist.movement), x, y, sp, epi, full_abun)
  # Lose moves from current grid square
  grd[sp, i] -= full_abun
  # Map moves to location in grid
  mov = @view lookup.moves[:, Threads.threadid()]
  for i in eachindex(epi.lookup[sp].x)
      newx = mod(lookup.x[i] + x - 1, width) + 1
      newy = mod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(epi, (newx, newy), width)
      grd[sp, loc] += mov[i]
  end
  return epi
end

function move!(epi::AbstractEpiSystem, ::NoMovement, i::Int64, sp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  return epi
end

function move!(epi::AbstractEpiSystem, ::BirthOnlyMovement, i::Int64, sp::Int64,
    grd::Array{Int64, 2}, births::Int64)
  width, height = getdimension(epi)
  (x, y) = convert_coords(epi, i, width)
   lookup = getlookup(epi, sp)
  calc_lookup_moves!(getboundary(epi.epilist.movement), x, y, sp, epi, births)
  # Lose moves from current grid square
  grd[sp, i] -= births
  # Map moves to location in grid
  mov = @view lookup.moves[:, Threads.threadid()]
  for i in eachindex(lookup.x)
      newx = mod(lookup.x[i] + x - 1, width) + 1
      newy = mod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(epi, (newx, newy), width)
      grd[sp, loc] += mov[i]
  end
  return epi
end

function habitatupdate!(epi::AbstractEpiSystem, timestep::Unitful.Time)
  _habitatupdate!(epi, epi.epienv.habitat, timestep)
end
function _habitatupdate!(epi::AbstractEpiSystem, hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab}, timestep::Unitful.Time)
    hab.change.changefun(epi, hab, timestep)
end

function _habitatupdate!(epi::AbstractEpiSystem, hab::HabitatCollection2, timestep::Unitful.Time)
    _habitatupdate!(epi, hab.h1, timestep)
    _habitatupdate!(epi, hab.h2, timestep)
end

function TempChange(epi::AbstractEpiSystem, hab::ContinuousHab, timestep::Unitful.Time)
  v = uconvert(K/unit(timestep), hab.change.rate)
  hab.matrix .+= (v * timestep)
end
