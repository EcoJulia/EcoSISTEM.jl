using StatsBase
using Compat
using LinearAlgebra

"""
    update!(eco::Ecosystem, time::Unitful.Time)
Function to update a ecosystem abundances and environment for one timestep.
"""
function update!(epi::EpiSystem, timestep::Unitful.Time)

    # Birth/death loop of each class including virus
    classupdate!(epi, timestep)

    # Update abundances with all movements
    epi.abundances.matrix .+= epi.cache.netmigration

    # Calculate new infections based on virus spread
    newinfections!(epi, timestep)
    # And new recoveries based upon recovery rate
    newrecoveries!(epi, timestep)

    # Invalidate all caches for next update
    invalidatecaches!(epi)

    # Update environment - habitat and energy budgets
    habitatupdate!(epi, timestep)
    applycontrols!(epi, timestep)
end

# function virusupdate!(epi::EpiSystem, timestep::Unitful.Time)
#     # Calculate dimenions of habitat and number of species
#     dims = _countsubcommunities(epi.epienv.habitat)
#     params = epi.epilist.params
#     width = getdimension(epi)[1]
#     rng = epi.abundances.seed[Threads.threadid()]
#
#     Threads.@threads for i in 1:dims
#         # Calculate how much birth and death should be adjusted
#         adjusted_birth, adjusted_death = adjustment(epi, i, j)
#
#         # Convert 1D dimension to 2D coordinates
#         (x, y) = convert_coords(epi, i, width)
#         # Check if grid cell currently active
#         if epi.epienv.active[x, y]
#             # Calculate effective rates
#             birthprob = params.beta * timestep * adjusted_birth
#             deathprob = params.death[1] * timestep * adjusted_death
#
#             # Put probabilities into 0 - 1
#             newbirthprob = 1.0 - exp(-birthprob)
#             newdeathprob = 1.0 - exp(-deathprob)
#
#             (newbirthprob >= 0) & (newdeathprob >= 0) || error("Birth: $newbirthprob \n Death: $newdeathprob \n \n i: $i \n j: $j")
#             # Calculate how many births and deaths
#             births = rand(rng, Poisson(eco.abundances.matrix[1, i] * newbirthprob))
#             deaths = rand(rng, Binomial(eco.abundances.matrix[1, i], newdeathprob))
#
#             # Update population
#             epi.abundances.matrix[1, i] += (births - deaths)
#
#             # Calculate moves and write to cache
#             move!(eco, eco.spplist.movement, i, 1, eco.cache.netmigration, births)
#         end
#     end
# end

function classupdate!(epi::EpiSystem, timestep::Unitful.Time)
    # Calculate dimenions of habitat and number of classes
    dims = _countsubcommunities(epi.epienv.habitat)
    classes = epi.epilist.names
    params = epi.epilist.params
    width = getdimension(epi)[1]
    # Loop through classes in chosen square
    Threads.@threads for j in 1:length(classes)
        class = getclass(classes[j])
        rng = epi.abundances.seed[Threads.threadid()]
        # Loop through grid squares
        for i in 1:dims
            # Calculate how much birth and death should be adjusted
            adjusted_birth, adjusted_death = adjustment(epi, class, i, j)

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(epi, i, width)
            # Check if grid cell currently active
            if epi.epienv.active[x, y]
                # Calculate effective rates
                birthprob = params.birth[j] * timestep * adjusted_birth
                deathprob = params.death[j] * timestep * adjusted_death

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                (newbirthprob >= 0) & (newdeathprob >= 0) || error("Birth: $newbirthprob \n Death: $newdeathprob \n \n i: $i \n j: $j")
                # Calculate how many births and deaths
                births = rand(rng, Binomial(epi.abundances.matrix[j, i],  newbirthprob))
                deaths = rand(rng, Binomial(epi.abundances.matrix[j, i], newdeathprob))

                # Update population
                epi.abundances.matrix[j, i] += (births - deaths)

                # Calculate moves and write to cache
                move!(epi, epi.epilist.movement, i, j, epi.cache.netmigration, births)
            end
        end
    end
end

function newinfections!(epi::EpiSystem, timestep::Unitful.Time)
    rng = epi.abundances.seed[Threads.threadid()]
    gridsize = epi.epienv.habitat.size
    virus = epi.abundances.matrix[1, :]
    virusparam = epi.epilist.params.viralload * gridsize / unit(gridsize) * timestep
    infprob = virus ./ (virusparam .+ virus)
    infections = rand.(fill(rng, length(infprob)), Binomial.(epi.abundances.matrix[2, :], infprob))
    epi.abundances.matrix[2, :] .-= infections
    epi.abundances.matrix[3, :] .+= infections
end

function newrecoveries!(epi::EpiSystem, timestep::Unitful.Time)
    rng = epi.abundances.seed[Threads.threadid()]
    infecteds = epi.abundances.matrix[3, :]
    recoveries = rand.(fill(rng, length(infecteds)), Binomial.(infecteds, uconvert(NoUnits, epi.epilist.params.sigma * timestep)))
    epi.abundances.matrix[3, :] .-= recoveries
    epi.abundances.matrix[4, :] .+= recoveries
end


function adjustment(epi::AbstractEpiSystem, class::Virus, i::Int64, j::Int64)
    susceptibles = epi.abundances.matrix[2, i]
    infecteds = epi.abundances.matrix[3, i]
    total = sum(epi.abundances.matrix[:, i])
    traitmatch = 1/traitfun(epi, i, j)
    birth_boost = (susceptibles * infecteds/ total) * traitmatch^-1
    death_boost = (susceptibles * infecteds/ total)^-1 * traitmatch
    return birth_boost, death_boost
end
function adjustment(epi::AbstractEpiSystem, class::Susceptible, i::Int64, j::Int64)
    return 1, 1
end
function adjustment(epi::AbstractEpiSystem, class::Infected, i::Int64, j::Int64)
    return 1, 1
end
function adjustment(epi::AbstractEpiSystem, class::Recovered, i::Int64, j::Int64)
    return 1, 1
end

function populate!(ml::EpiLandscape, epilist::EpiList, epienv::EE, rel::R) where {EE <: AbstractEpiEnv, R <: AbstractTraitRelationship}
    dim = _getdimension(epienv.habitat)
    len = dim[1] * dim[2]
    # Loop through classes
    for i in eachindex(epilist.abun)
        rand!(Multinomial(epilist.abun[i], len), (@view ml.matrix[i, :]))
    end
end


function applycontrols!(epi::EpiSystem, timestep::Unitful.Time)
    _applycontrols!(epi, epi.epienv.control, timestep)
end

function _applycontrols!(epi::EpiSystem, controls::NoControl, timestep::Unitful.Time)
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


function calc_lookup_moves!(bound::NoBoundary, x::Int64, y::Int64, sp::Int64, epi::AbstractEpiSystem, abun::Int64)
    lookup = getlookup(epi, sp)
    maxX = getdimension(epi)[1] - x
    maxY = getdimension(epi)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        valid =  (-x < lookup.x[i] <= maxX) && (-y < lookup.y[i] <= maxY) && (epi.epienv.active[lookup.x[i] + x, lookup.y[i] + y])

        lookup.pnew[i] = valid ? lookup.p[i] : 0.0
    end
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(epi.abundances.seed[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Cylinder, x::Int64, y::Int64, sp::Int64, epi::AbstractEpiSystem, abun::Int64)
    lookup = getlookup(epi, sp)
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
    rand!(epi.abundances.seed[Threads.threadid()], dist, lookup.moves)
end

function calc_lookup_moves!(bound::Torus, x::Int64, y::Int64, sp::Int64, epi::AbstractEpiSystem, abun::Int64)
  lookup = getlookup(epi, sp)
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
  rand!(epi.abundances.seed[Threads.threadid()], dist, lookup.moves)
end

"""
    move!(eco::Ecosystem, ::AbstractMovement, i::Int64, sp::Int64, grd::Array{Int64, 2}, abun::Int64)

Function to calculate the movement of species `sp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a cached grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(epi::AbstractEpiSystem, ::AlwaysMovement, i::Int64, sp::Int64,
  grd::Array{Int64, 2}, ::Int64)
  width, height = getdimension(epi)
  (x, y) = convert_coords(epi, i, width)
  lookup = getlookup(epi, sp)
  full_abun = epi.abundances.matrix[sp, i]
  calc_lookup_moves!(getboundary(epi.epilist.movement), x, y, sp, epi, full_abun)
  # Lose moves from current grid square
  grd[sp, i] -= full_abun
  # Map moves to location in grid
  mov = lookup.moves
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
  mov = lookup.moves
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
  val = hab.change.rate
  v = uconvert(K/unit(timestep), val)
  hab.matrix .+= v * timestep
end
