using StatsBase
using Compat
using LinearAlgebra

"""
    update!(epi::EpiSystem, time::Unitful.Time)
Function to update disease and virus class abundances and environment for one timestep.
"""
function update!(epi::EpiSystem, timestep::Unitful.Time)

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

"""
    virusupdate!(epi::EpiSystem, time::Unitful.Time)
Function to update virus abundances and disperse for one timestep.
"""
function virusupdate!(epi::EpiSystem, timestep::Unitful.Time)
    dims = _countsubcommunities(epi.epienv.habitat)
    width = getdimension(epi)[1]
    params = epi.epilist.params
    id = Threads.threadid()
    rng = epi.abundances.seed[id]
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
                birthrate = params.virus_growth[j] * timestep * human(epi.abundances)[j, i]
                deathrate = params.virus_decay[j] * timestep * traitmatch^-1

                # Convert death rate into 0 - 1 probability
                deathprob = 1.0 - exp(-deathrate)

                (birthrate >= 0) & (deathprob >= 0) || error("Birth: $birthprob \n Death: $deathprob \n \n i: $i")
                # Calculate how many births and deaths
                births = rand(rng, Poisson(birthrate))
                deaths = rand(rng, Binomial(virus(epi.abundances)[1, i], deathprob))

                # Update population
                epi.cache.virusmigration[j, i] += births
                epi.cache.virusdecay[j, i] -= deaths
                virusmove!(epi, i, j, epi.cache.virusmigration, births)
            end
        end
    end
    vm = sum(epi.cache.virusmigration, dims = 1)[1, :]
    nm = sum(epi.cache.virusdecay, dims = 1)[1, :]
    virus(epi.abundances)[1, :] .+= (nm .+ vm)
    virus(epi.abundances)[2, :] = vm
end

"""
    sum_pop(M::Matrix{Int64}, i::Int64)
Function to sum a population matrix, `M`, without memory allocation, at a grid location `i`.
"""
function sum_pop(M::Matrix{Int64}, i::Int64)
    N = 0
    for j in 1:size(M, 1)
        N += M[j, i]
    end
    return N
end


"""
    classupdate!(epi::EpiSystem, timestep::Unitful.Time)
Function to update disease class abundances for one timestep.
"""
function classupdate!(epi::EpiSystem, timestep::Unitful.Time)
    # Calculate dimenions of habitat and number of classes
    dims = _countsubcommunities(epi.epienv.habitat)
    params = epi.epilist.params
    width = getdimension(epi)[1]
    classes = size(human(epi.abundances), 1)

    # Loop through grid squares
    Threads.@threads for i in 1:dims
        # will result +=/-= 0 at end of inner loop, so safe to skip
        all(iszero, human(epi.abundances)[:, i]) && continue

        rng = epi.abundances.seed[Threads.threadid()]
        N = sum_pop(epi.abundances.matrix, i)
        # Convert 1D dimension to 2D coordinates
        (x, y) = convert_coords(epi, i, width)
        # Check if grid cell currently active. If inactive skip the inner loop
        epi.epienv.active[x, y] || continue
        # Loop through classes in chosen square
        for j in 1:classes
            # Births
            births = rand(rng, Binomial(human(epi.abundances)[j, i],  params.births[j] * timestep))
            human(epi.abundances)[1, i] += births

            # Note transposition of transition matrices to make iteration over k faster
            # Calculate force of inf and env inf
            for k in 1:size(params.transition_virus, 1)
                iszero(human(epi.abundances)[k, i]) && continue # will result +=/-= 0 at end of loop
                env_inf = (params.transition_virus[j, k] * timestep * virus(epi.abundances)[1, i]) / (N^params.freq_vs_density_env)

                force_inf = (params.transition_force[j, k] * timestep * virus(epi.abundances)[2, i]) / (N^params.freq_vs_density_force)

                # Add to transitional probabilities
                trans_val = (params.transition[j, k] * timestep) + ((1 - params.force_vs_env) * env_inf) + (params.force_vs_env * force_inf)
                trans_prob = 1 - exp(-trans_val)
                iszero(trans_prob) && continue # will result +=/-= 0 at end of loop

                # Make transitions
                trans = rand(rng, Binomial(human(epi.abundances)[k, i], trans_prob))
                human(epi.abundances)[j, i] += trans
                human(epi.abundances)[k, i] -= trans
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
    for i in eachindex(epilist.human.abun)
        rand!(Multinomial(epilist.human.abun[i], len), (@view ml.matrix[i, :]))
    end
    for i in eachindex(epilist.virus.abun)
        rand!(Multinomial(epilist.virus.abun[i], len), (@view ml.matrix_v[i, :]))
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

    # lookup.pnew[i] would be set to zero at the top of the loop, but do it here instead
    # incase we return early when we discover there's no hard work to do
    fill!(lookup.pnew, zero(eltype(lookup.pnew)))

    # If all lookup.p are zero or abun is zero then there's no need to do anything.
    # This is because the insertion of non-zero rands into lookup.moves depends on dist
    # producing non rands, which is only true if both abun and lookup.p are non-zero
    if iszero(abun) || all(iszero, lookup.p)
      fill!(lookup.moves, zero(eltype(lookup.moves)))
      return nothing
    end

    maxX = getdimension(epi)[1]
    maxY = getdimension(epi)[2]
    # the for loop looks a bit grotty because it's been optimised
    # make references to containers on structs to avoid getproperty overheads
    lookuppnew, lookupp, lookupx, lookupy = lookup.pnew, lookup.p, lookup.x, lookup.y
    epiepienvactive = epi.epienv.active
    @inbounds for i in eachindex(lookupx)
        lookuppi = lookupp[i]
        iszero(lookuppi) && continue # then lookuppnew is already correct
        lookupxi = lookupx[i] + x
        (0 < lookupxi <= maxX) || continue # Can't go over maximum dimension
        lookupyi = lookupy[i] + y
        (0 < lookupyi <= maxY) || continue # Can't go over maximum dimension
        epiepienvactive[lookupxi, lookupyi] || continue # skip if inactive
        lookuppnew[i] = lookuppi
    end
    # Case that produces NaN (if sum(lookup.pnew) also pnew always .>= 0) dealt with above
    lookup.pnew ./= sum(lookup.pnew)
    dist = Multinomial(abun, lookup.pnew)
    rand!(epi.abundances.seed[Threads.threadid()], dist, lookup.moves)
    return nothing
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
    rand!(epi.abundances.seed[Threads.threadid()], dist, lookup.moves)
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
    rand!(epi.abundances.seed[Threads.threadid()], dist, lookup.moves)
end


"""
    virusmove!(epi::AbstractEpiSystem, pos::Int64, id::Int64, grd::Array{Int64, 2}, newvirus::Int64)

Function to calculate the movement of force of infection `id` from a given position in the landscape `pos`, using the lookup table found in the EpiSystem and updating the movement patterns on a cached grid, `grd`. The number of new virus is provided, so that movement only takes place as part of the generation process.
"""
function virusmove!(epi::AbstractEpiSystem, pos::Int64, id::Int64, grd::Array{Int64, 2}, newvirus::Int64)
  width, height = getdimension(epi)
  (x, y) = convert_coords(epi, pos, width)
  lookup = getlookup(epi, id)
  calc_lookup_moves!(getboundary(epi.epilist.human.movement), x, y, id, epi, newvirus)
  # Lose moves from current grid square
  grd[id, pos] -= newvirus
  # Map moves to location in grid
  @inbounds for (i, move) in enumerate(lookup.moves)
      newx = quickmod(lookup.x[i] + x - 1, width) + 1
      newy = quickmod(lookup.y[i] + y - 1, height) + 1
      loc = convert_coords(epi, (newx, newy), width)
      grd[id, loc] += move
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

function quickmod(x::T, h::T) where {T<:Integer}
    while !(-1 < x < h)
        if x <= -1
            x += h
        elseif x >= h
            x -= h
        end
    end
    return x
end
