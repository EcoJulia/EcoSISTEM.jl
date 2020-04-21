using StatsBase
using Compat
using LinearAlgebra

"""
    update!(epi::EpiSystem, time::Unitful.Time)
Function to update disease class abundances and environment for one timestep.
"""
function update!(epi::EpiSystem, timestep::Unitful.Time)

    # Human movement loop
    humanmove!(epi, timestep, getclass(epi))

    # Virus movement loop
    virusupdate!(epi, timestep, getclass(epi))

    # Birth/death/infection/recovery loop of each class
    classupdate!(epi, timestep, getclass(epi))

    # Invalidate all caches for next update
    invalidatecaches!(epi)

    # Update environment - habitat and energy budgets
    habitatupdate!(epi, timestep)
    applycontrols!(epi, timestep)
end

function humanmove!(epi::EpiSystem, timestep::Unitful.Time, model::MC) where MC <: ModelClass
    dims = _countsubcommunities(epi.epienv.habitat)
    width = getdimension(epi)[1]
    classes = size(epi.abundances.matrix, 1)
    Threads.@threads for j in 2:classes
        rng = epi.abundances.seed[Threads.threadid()]
        # Loop through grid squares
        for i in 1:dims
            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(epi, i, width)
            # Check if grid cell currently active
            if epi.epienv.active[x, y]
                # Calculate moves and write to cache
                move!(epi, epi.epilist.movement, i, j, epi.cache.netmigration, epi.abundances.matrix[j, i])
            end
        end
    end
    # Update abundances with all movements
    epi.abundances.matrix .+= epi.cache.netmigration
end

function virusupdate!(epi::EpiSystem, timestep::Unitful.Time, model::MC) where MC <: ModelClass
    dims = _countsubcommunities(epi.epienv.habitat)
    width = getdimension(epi)[1]
    params = epi.epilist.params
    id = Threads.threadid()
    rng = epi.abundances.seed[id]
    # Loop through grid squares
    Threads.@threads for i in 1:dims
        # Calculate how much birth and death should be adjusted
        birth_adjust, death_adjust = adjustment(epi, i, model)

        # Convert 1D dimension to 2D coordinates
        (x, y) = convert_coords(epi, i, width)
        # Check if grid cell currently active
        if epi.epienv.active[x, y]
            # Calculate effective rates
            birthprob = params.birth[1] * timestep * birth_adjust
            deathprob = params.death[1] * timestep * death_adjust

            # Put probabilities into 0 - 1
            newbirthprob = 1.0 - exp(-birthprob)
            newdeathprob = 1.0 - exp(-deathprob)

            (newbirthprob >= 0) & (newdeathprob >= 0) || error("Birth: $newbirthprob \n Death: $newdeathprob \n \n i: $i")
            # Calculate how many births and deaths
            births = rand(rng, Binomial(epi.abundances.matrix[1, i],  newbirthprob))
            deaths = rand(rng, Binomial(epi.abundances.matrix[1, i], newdeathprob))

            # Update population
            epi.abundances.matrix[1, i] += (births - deaths)

            virusmove!(epi, epi.epilist.movement, i, id, epi.cache.virusmigration, births)

        end
    end
    epi.abundances.matrix[1, :] .+= sum(epi.cache.virusmigration, dims = 1)[1, :]
end

"""
    classupdate!(epi::EpiSystem, timestep::Unitful.Time, model::MC) where MC <: ModelClass
Function to update disease class abundances for one timestep. Dispatches differently depending on the class of model stored in the EpiSystem.
"""
function classupdate!(epi::EpiSystem, timestep::Unitful.Time, model::SIR)
    # Calculate dimenions of habitat and number of classes
    dims = _countsubcommunities(epi.epienv.habitat)
    params = epi.epilist.params
    width = getdimension(epi)[1]
    classes = size(epi.abundances.matrix, 1)
    dict = epi.epilist.model.dict
    # Loop through grid squares
    Threads.@threads for i in 1:dims
        rng = epi.abundances.seed[Threads.threadid()]
        # Loop through classes in chosen square
        for j in 2:classes
            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(epi, i, width)
            # Check if grid cell currently active
            if epi.epienv.active[x, y]
                # Calculate effective rates
                birthprob = params.birth[j] * timestep
                deathprob = params.death[j] * timestep

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                (newbirthprob >= 0) & (newdeathprob >= 0) || error("Birth: $newbirthprob \n Death: $newdeathprob \n \n i: $i \n j: $j")
                # Calculate how many births and deaths
                births = rand(rng, Binomial(epi.abundances.matrix[j, i],  newbirthprob))
                deaths = rand(rng, Binomial(epi.abundances.matrix[j, i], newdeathprob))

                # Update population
                epi.abundances.matrix[dict["Susceptible"], i] += births
                epi.abundances.matrix[j, i] -= deaths
            end
        end
        # Infections
        newinfections!(epi, timestep, i, model)

        # Recoveries
        newrecoveries!(epi, timestep, i, model)
    end
end

function classupdate!(epi::EpiSystem, timestep::Unitful.Time, model::SEI2HRD)
    # Calculate dimenions of habitat and number of classes
    dims = _countsubcommunities(epi.epienv.habitat)
    params = epi.epilist.params
    width = getdimension(epi)[1]
    classes = size(epi.abundances.matrix, 1)
    dict = epi.epilist.model.dict
    # Loop through grid squares
    Threads.@threads for i in 1:dims
        # Loop through classes in chosen square
        for j in 2:(classes - 1)
            rng = epi.abundances.seed[Threads.threadid()]

            # Convert 1D dimension to 2D coordinates
            (x, y) = convert_coords(epi, i, width)
            # Check if grid cell currently active
            if epi.epienv.active[x, y]
                # Calculate effective rates
                birthprob = params.birth[j] * timestep
                deathprob = params.death[j] * timestep

                # Put probabilities into 0 - 1
                newbirthprob = 1.0 - exp(-birthprob)
                newdeathprob = 1.0 - exp(-deathprob)

                (newbirthprob >= 0) & (newdeathprob >= 0) || error("Birth: $newbirthprob \n Death: $newdeathprob \n \n i: $i \n j: $j")
                # Calculate how many births and deaths
                births = rand(rng, Binomial(epi.abundances.matrix[j, i],  newbirthprob))
                deaths = rand(rng, Binomial(epi.abundances.matrix[j, i], newdeathprob))

                # Update population
                epi.abundances.matrix[dict["Susceptible"], i] += births
                epi.abundances.matrix[j, i] -= deaths
                epi.abundances.matrix[dict["Dead"], i] += deaths
            end
        end
        # Infections
        newinfections!(epi, timestep, i, model)
        # Recoveries
        newrecoveries!(epi, timestep, i, model)
        # Deaths
        newdeaths!(epi, timestep, i, model)
    end
end

"""
    newinfections!(epi::EpiSystem, timestep::Unitful.Time)
Function to generate new infections based on viral load in each grid square.
"""
function newinfections!(epi::EpiSystem, timestep::Unitful.Time, pos::Int64, model::SIR)
    dict = epi.epilist.model.dict
    susclass = dict["Susceptible"]
    infclass = dict["Infected"]
    rng = epi.abundances.seed[Threads.threadid()]
    infprob = epi.abundances.matrix[1, pos] * (epi.epilist.params.beta * timestep)
    infprob <= 1 || error("Infection probability greater than 1.")
    infections = rand(rng, Binomial(epi.abundances.matrix[susclass, pos], infprob))
    epi.abundances.matrix[susclass, pos] -= infections
    epi.abundances.matrix[infclass, pos] += infections
end

function newinfections!(epi::EpiSystem, timestep::Unitful.Time, pos::Int64, model::SEI2HRD)
    rng = epi.abundances.seed[Threads.threadid()]
    # Get disease classes
    dict = epi.epilist.model.dict
    virclass = dict["Virus"]
    susclass = dict["Susceptible"]
    expclass = dict["Exposed"]
    inf1class = dict["AsymptomaticInfected"]
    inf2class = dict["SymptomaticInfected"]
    hospclass = dict["Hospitalised"]
    # Calc prob of infection
    expprob = epi.abundances.matrix[virclass, pos] * (epi.epilist.params.beta * timestep)
    expprob <= 1 || error("Infection probability greater than 1.")
    # New exposures
    exposures = rand(rng, Binomial(epi.abundances.matrix[susclass, pos], expprob))
    epi.abundances.matrix[susclass, pos] -= exposures
    epi.abundances.matrix[expclass, pos] += exposures

    # New asymptomatic infections
    infections = rand(rng, Binomial(epi.abundances.matrix[expclass, pos], epi.epilist.params.mu_1 * timestep))
    epi.abundances.matrix[expclass, pos] -= infections
    epi.abundances.matrix[inf1class, pos] += infections

    # New symptomatic infections
    infections2 = rand(rng, Binomial(epi.abundances.matrix[inf1class, pos], epi.epilist.params.mu_2 * timestep))
    epi.abundances.matrix[inf1class, pos] -= infections2
    epi.abundances.matrix[inf2class, pos] += infections2

    # New hospitalisations
    hosp = rand(rng, Binomial(epi.abundances.matrix[inf2class, pos], epi.epilist.params.hospitalisation * timestep))
    epi.abundances.matrix[inf2class, pos] -= hosp
    epi.abundances.matrix[hospclass, pos] += hosp
end


"""
    newrecoveries!(epi::EpiSystem, timestep::Unitful.Time)
Function to generate new recoveries based on a set recovery rate in each grid square.
"""
function newrecoveries!(epi::EpiSystem, timestep::Unitful.Time, pos::Int64, model::SIR)
    dict = epi.epilist.model.dict
    infclass = dict["Infected"]
    recclass = dict["Recovered"]
    rng = epi.abundances.seed[Threads.threadid()]
    recoveries = rand(rng, Binomial(epi.abundances.matrix[infclass, pos], uconvert(NoUnits, epi.epilist.params.sigma * timestep)))
    epi.abundances.matrix[infclass, pos] -= recoveries
    epi.abundances.matrix[recclass, pos] += recoveries
end

function newrecoveries!(epi::EpiSystem, timestep::Unitful.Time, pos::Int64, model::SEI2HRD)
    dict = epi.epilist.model.dict
    inf1class = dict["AsymptomaticInfected"]
    inf2class = dict["SymptomaticInfected"]
    hospclass = dict["Hospitalised"]
    recclass = dict["Recovered"]
    rng = epi.abundances.seed[Threads.threadid()]

    # Asymptomatic recoveries
    recoveries_1 = rand(rng, Binomial(epi.abundances.matrix[inf1class, pos], uconvert(NoUnits, epi.epilist.params.sigma_1 * timestep)))
    epi.abundances.matrix[inf1class, pos] -= recoveries_1
    epi.abundances.matrix[recclass, pos] += recoveries_1

    # Symptomatic recoveries
    recoveries_2 = rand(rng, Binomial(epi.abundances.matrix[inf2class, pos], uconvert(NoUnits, epi.epilist.params.sigma_2 * timestep)))
    epi.abundances.matrix[inf2class, pos] -= recoveries_2
    epi.abundances.matrix[recclass, pos] += recoveries_2

    # Hospital recoveries
    recoveries_hosp = rand(rng, Binomial(epi.abundances.matrix[hospclass, pos], uconvert(NoUnits, epi.epilist.params.sigma_hospital * timestep)))
    epi.abundances.matrix[hospclass, pos] -= recoveries_hosp
    epi.abundances.matrix[recclass, pos] += recoveries_hosp

end

"""
    newdeaths!(epi::EpiSystem, timestep::Unitful.Time)
Function to generate new deaths based on a death rate in each grid square.
"""
function newdeaths!(epi::EpiSystem, timestep::Unitful.Time, pos::Int64, model::SEI2HRD)
    rng = epi.abundances.seed[Threads.threadid()]
    dict = epi.epilist.model.dict
    inf2class = dict["SymptomaticInfected"]
    hospclass = dict["Hospitalised"]
    deathclass = dict["Dead"]

    # Deaths at home
    deaths_home = rand(rng, Binomial(epi.abundances.matrix[inf2class, pos], uconvert(NoUnits, epi.epilist.params.death_home * timestep)))
    epi.abundances.matrix[inf2class, pos] -= deaths_home
    epi.abundances.matrix[deathclass, pos] += deaths_home

    # Deaths at hospital
    deaths_hosp = rand(rng, Binomial(epi.abundances.matrix[hospclass, pos], uconvert(NoUnits, epi.epilist.params.death_hospital * timestep)))
    epi.abundances.matrix[hospclass, pos] -= deaths_hosp
    epi.abundances.matrix[deathclass, pos] += deaths_hosp
end


"""
    adjustment(epi::AbstractEpiSystem, class::DiseaseClass, i::Int64, j::Int64)
Function to calculate match between environment and birth/death for each disease class. This is assumed to be 1 for all susceptible, infected, recovered and dependent on a trait function for the virus.
"""
function adjustment(epi::AbstractEpiSystem, pos::Int64, model::SIR)
    dict = epi.epilist.model.dict
    infclass = dict["Infected"]
    virclass = dict["Virus"]
    traitmatch = traitfun(epi, pos, 1)
    birth_boost = epi.abundances.matrix[infclass, pos] * traitmatch
    death_boost = epi.abundances.matrix[virclass, pos] * traitmatch^-1
    return birth_boost, death_boost
end

function adjustment(epi::AbstractEpiSystem, pos::Int64, model::SEI2HRD)
    dict = epi.epilist.model.dict
    inf1class = dict["AsymptomaticInfected"]
    inf2class = dict["SymptomaticInfected"]
    virclass = dict["Virus"]
    traitmatch = traitfun(epi, pos, 1)
    birth_boost = (epi.abundances.matrix[inf1class, pos] + epi.abundances.matrix[inf2class, pos]) * traitmatch
    death_boost = epi.abundances.matrix[virclass, pos] * traitmatch^-1
    return birth_boost, death_boost
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


function calc_lookup_moves!(bound::NoBoundary, x::Int64, y::Int64, sp::Int64, epi::AbstractEpiSystem, abun::Int64)
    id = Threads.threadid()
    lookup = getlookup(epi, sp)
    maxX = getdimension(epi)[1] - x
    maxY = getdimension(epi)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        valid =  (-x < lookup.x[i] <= maxX) && (-y < lookup.y[i] <= maxY) && (epi.epienv.active[lookup.x[i] + x, lookup.y[i] + y])

        lookup.pnew[i, id] = valid ? lookup.p[i] : 0.0
    end
    pnew = @view lookup.pnew[:, id]
    pnew ./= sum(pnew)
    dist = Multinomial(abun, pnew)
    m = @view lookup.moves[:, id]
    rand!(epi.abundances.seed[Threads.threadid()], dist, m)
end

function calc_lookup_moves!(bound::Cylinder, x::Int64, y::Int64, sp::Int64, epi::AbstractEpiSystem, abun::Int64)
    id = Threads.threadid()
    lookup = getlookup(epi, sp)
    maxX = getdimension(epi)[1] - x
    maxY = getdimension(epi)[2] - y
    # Can't go over maximum dimension
    for i in eachindex(lookup.x)
        newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, getdimension(epi)[1]) + 1

        valid =  (-y < lookup.y[i] <= maxY) && (epi.epienv.active[newx, lookup.y[i] + y])

        lookup.pnew[i, id] = valid ? lookup.p[i] : 0.0
    end
    pnew = @view lookup.pnew[:, id]
    pnew ./= sum(pnew)
    dist = Multinomial(abun, pnew)
    m = @view lookup.moves[:, id]
    rand!(epi.abundances.seed[Threads.threadid()], dist, m)
end

function calc_lookup_moves!(bound::Torus, x::Int64, y::Int64, sp::Int64, epi::AbstractEpiSystem, abun::Int64)
    id = Threads.threadid()
  lookup = getlookup(epi, sp)
  maxX = getdimension(epi)[1] - x
  maxY = getdimension(epi)[2] - y
  # Can't go over maximum dimension
  for i in eachindex(lookup.x)
      newx = -x < lookup.x[i] <= maxX ? lookup.x[i] + x : mod(lookup.x[i] + x - 1, getdimension(epi)[1]) + 1
      newy =  -y < lookup.y[i] <= maxY ? lookup.y[i] + y : mod(lookup.y[i] + y - 1, getdimension(epi)[2]) + 1
      valid = epi.epienv.active[newx, newy]

      lookup.pnew[i, id] = valid ? lookup.p[i] : 0.0
  end
  pnew = @view lookup.pnew[:, id]
  pnew ./= sum(pnew)
  dist = Multinomial(abun, pnew)
  m = @view lookup.moves[:, id]
  rand!(epi.abundances.seed[Threads.threadid()], dist, m)
end

function virusmove!(epi::AbstractEpiSystem, ::AlwaysMovement, pos::Int64, id::Int64, grd::Array{Int64, 2}, ::Int64)
  width, height = getdimension(epi)
  (x, y) = convert_coords(epi, pos, width)
  lookup = getlookup(epi, 1)
  full_abun = epi.abundances.matrix[1, pos]
  calc_lookup_moves!(getboundary(epi.epilist.movement), x, y, 1, epi, full_abun)
  # Lose moves from current grid square
  grd[id, pos] -= full_abun
  # Map moves to location in grid
  mov = @view lookup.moves[:, id]
  for i in eachindex(epi.lookup[1].x)
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
  full_abun = epi.abundances.matrix[sp, i]
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
