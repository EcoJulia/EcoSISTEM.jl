using JLSO
using SparseArrays

"""
    AbstractEpiSystem

Abstract supertype for all disease system types.
"""
abstract type AbstractEpiSystem{Part <: AbstractEpiEnv, EL <: EpiList, TR <: AbstractTraitRelationship} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                    Matrix{Float64}, EL, Part} end


mutable struct EpiCache <: AbstractCache
  virusmigration::Array{Float64, 2}
  valid::Bool
end

@enum MovementType homeMovement workMovement

struct EpiLookup <: AbstractLookup
  homelookup::SparseMatrixCSC{Float64, Int32}
  worklookup::SparseMatrixCSC{Float64, Int32}
  function EpiLookup(homelookup::SparseMatrixCSC{Float64, Int32}, worklookup::SparseMatrixCSC{Float64, Int32})
      all(0 .<= homelookup.nzval .<= 1) || error("Home lookup values must be between 0 and 1")
      all(0 .<= worklookup.nzval .<= 1) || error("Work lookup values must be between 0 and 1")
      return new(homelookup, worklookup)
  end
end

"""
    EpiSystem{EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractRelationship} <: AbstractEpiSystem{EE, EL, ER}

EpiSystem houses information on different disease classes, `epilist`, the environment, `epienv`, and their relationship to one another, `relationship`.

See `help?>plot_epidynamics` and `help?>plot_epiheatmaps` for relevant plotting functions.
"""
mutable struct EpiSystem{U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG}, EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractTraitRelationship} <: AbstractEpiSystem{EE, EL, ER}
  abundances::EpiLandscape{U, VecRNGType}
  epilist::EL
  epienv::EE
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::ER
  lookup::EpiLookup
  cache::EpiCache
  initial_infected::Int64
  ordered_active::Vector{Int64}
  transitions::Union{Missing, TransitionList}
end
function EpiSystem(abundances::EpiLandscape{U, VecRNGType}, epilist::EL, epienv::EE,
    ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::EpiLookup,
    cache::EpiCache, initial_infected::Int64, transitions::Union{Missing, TransitionList}
    ) where {U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG},
    EE <: AbstractEpiEnv, EL <: EpiList, ER <: AbstractTraitRelationship}
  total_pop = sum(abundances.matrix, dims = 1)[1, :]
  sorted_grid_ids = sortperm(total_pop, rev = true)
  sorted_grid_ids = sorted_grid_ids[total_pop[sorted_grid_ids] .> 0]
  return EpiSystem(abundances, epilist, epienv, ordinariness, relationship, lookup, cache, initial_infected, sorted_grid_ids, transitions)
end

function EpiSystem(popfun::F, epilist::EpiList, epienv::GridEpiEnv,
      rel::AbstractTraitRelationship, intnum::U; initial_infected = 0,
      rngtype::Type{R} = Random.MersenneTwister,
      transitions = missing) where {F<:Function, U <: Integer, R <: Random.AbstractRNG}

  # Create matrix landscape of zero abundances
  ml = emptyepilandscape(epienv, epilist, intnum, rngtype)

  # Populate this matrix with species abundances
  popfun(ml, epilist, epienv, rel)
  initial_pop = sum(ml.matrix, dims = 1)
  # Create lookup table of all moves and their probabilities
  home_lookup = genlookups(epienv, epilist.human.movement.home)
  work_lookup = genlookups(epienv, epilist.human.movement.work, initial_pop)
  lookup = EpiLookup(home_lookup, work_lookup)
  vm = zeros(Float64, size(ml.matrix))
  return EpiSystem(ml, epilist, epienv, missing, rel, lookup, EpiCache(vm, false), initial_infected, transitions)
end

function EpiSystem(epilist::EpiList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
        intnum::U = Int64(1); initial_infected = 0, rngtype::Type{R} = Random.MersenneTwister,
        transitions = missing
        ) where {U <: Integer, R <: Random.AbstractRNG}
    return EpiSystem(populate!, epilist, epienv, rel, intnum, initial_infected = initial_infected, rngtype = rngtype, transitions = transitions)
end

function EpiSystem(epilist::EpiList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
        initial_population::A, intnum::U = Int64(1); initial_infected = 0,
        rngtype::Type{R} = Random.MersenneTwister,
        transitions = missing) where {U <: Integer, A <: AbstractArray, R <: Random.AbstractRNG}
    if size(initial_population) != size(epienv.active)
        msg = "size(initial_population)==$(size(initial_population)) != " *
            "size(epienv.active)==$(size(epienv.active))"
        throw(DimensionMismatch(msg))
    end
    epienv.active .&= .!_inactive.(initial_population)

    # Create matrix landscape of zero abundances
    ml = emptyepilandscape(epienv, epilist, intnum, rngtype)

    # Create lookup table of all moves and their probabilities
    home_lookup = genlookups(epienv, epilist.human.movement.home)
    work_lookup = genlookups(epienv, epilist.human.movement.work, initial_population[1:end])
    lookup = EpiLookup(home_lookup, work_lookup)

    vm = zeros(Float64, size(ml.matrix))

    epi = EpiSystem(ml, epilist, epienv, missing, rel, lookup, EpiCache(vm, false), initial_infected, transitions)

    # Add in the initial susceptible population
    # TODO Need to fix code so it doesn't rely on name of susceptible class
    idx = findfirst(occursin.("Susceptible", epilist.human.names))
    if idx == nothing
        msg = "epilist has no Susceptible category. epilist.names = $(epilist.human.names)"
        throw(ArgumentError(msg))
    end
    # Modify active cells based on new population
    initial_population = convert_population(initial_population, intnum)
    epi.abundances.grid[idx, :, :] .+= initial_population

    total_pop = sum(human(epi.abundances), dims = 1)[1, :]
    sorted_grid_ids = sortperm(total_pop, rev = true)
    sorted_grid_ids = sorted_grid_ids[total_pop[sorted_grid_ids] .> 0]
    epi.ordered_active = sorted_grid_ids
    return epi
end

"""
    isapprox(epi_1::AbstractEpiSystem, epi_2::AbstractEpiSystem; kwargs...)

Compare two `EpiSystem`s for approximate equality. Specifically, compares the
`EpiLandscape`s of the two systems.

## Keyword arguments
- Anything to pass to `Base.isapprox`.

!!! note
    You may want to pass in `atol` or `rtol` to loosen the equality tolerance.
"""
function Base.isapprox(epi_1::AbstractEpiSystem, epi_2::AbstractEpiSystem; kwargs...)
    return isapprox(epi_1.abundances, epi_2.abundances; kwargs...)
end

save(path::String, system::EpiSystem) = JLSO.save(path, :episystem => system)
load(path::String, obj_type::Type{EpiSystem}) = JLSO.load(path)[:episystem]

function getsize(epi::AbstractEpiSystem)
  return _getsize(epi.epienv.habitat)
end

function getgridsize(epi::AbstractEpiSystem)
  return _getgridsize(epi.epienv.habitat)
end

function getdimension(epi::AbstractEpiSystem)
    return _getdimension(epi.epienv.habitat)
end


function gettraitrel(epi::AbstractEpiSystem)
  return epi.relationship
end

function gethabitat(epi::AbstractEpiSystem)
  return epi.epienv.habitat
end

import Diversity.API: _getabundance
function _getabundance(epi::AbstractEpiSystem, input::Bool)
    if input
        return epi.abundances.matrix
    else
        return _calcabundance(_gettypes(epi), epi.abundances.matrix / sum(epi.abundances.matrix))[1]
    end
end

import Diversity.API: _getmetaabundance
function _getmetaabundance(epi::AbstractEpiSystem)
  return sumoversubcommunities(epi, _getabundance(epi))
end


import Diversity.API: _getpartition
function _getpartition(epi::AbstractEpiSystem)
  return epi.epienv
end
import Diversity.API: _gettypes
function _gettypes(epi::AbstractEpiSystem)
    return epi.epilist
end
import Diversity.API: _getordinariness!
function _getordinariness!(epi::AbstractEpiSystem)
    if ismissing(epi.ordinariness)
        relab = getabundance(epi, false)
        epi.ordinariness = _calcordinariness(epi.epilist, relab)
    end
    return epi.ordinariness
end

import Diversity.API._getscale
function _getscale(epi::AbstractEpiSystem)
    return _calcabundance(_gettypes(epi), getabundance(epi, false))[2]
end

function invalidatecaches!(epi::AbstractEpiSystem)
    epi.ordinariness = missing
    epi.cache.virusmigration .= 0
    epi.cache.valid = false
end

function getdispersaldist(epi::AbstractEpiSystem, id::Int64)
  dist = epi.epilist.human.movement.home.kernels[id].dist
  return dist
end
function getdispersaldist(epi::AbstractEpiSystem, id::String)
  num = findfirst(epi.epilist.human.names .== id)
  getdispersaldist(epi, num)
end

function getdispersalvar(epi::AbstractEpiSystem, id::Int64)
    var = (epi.epilist.human.movement.home.kernels[id].dist)^2 * pi / 4
    return var
end
function getdispersalvar(epi::AbstractEpiSystem, id::String)
    num = findfirst(epi.epilist.human.names .== id)
    getdispersalvar(epi, num)
end

function getlookup(epi::AbstractEpiSystem, id::Int64, movetype::MovementType)
    if movetype == homeMovement
        return epi.lookup.homelookup[id, :]
    elseif movetype == workMovement
        return epi.lookup.worklookup[id, :]
    else
        return error("No other movement types currently implemented")
    end
end

function getlookup(epi::AbstractEpiSystem, id::Int64)
    return epi.lookup.homelookup[id, :], epi.lookup.worklookup[id, :]
end

function genlookups(epienv::AbstractEpiEnv, mov::Commuting, pop_size)
    total_size = (size(epienv.active, 1) * size(epienv.active, 2))
    # Column access so Js should be source grid cells
    Js = Int32.(mov.home_to_work[!, :from])
    # Is should be destination grid cells
    Is = Int32.(mov.home_to_work[!, :to])
    Vs = mov.home_to_work[!, :count]
    work = sparse(Is, Js, Vs, total_size, total_size)
    # Divide through by total population size
    work.nzval ./= pop_size[Is]
    work.nzval[isnan.(work.nzval)] .= 0
    # Make sure each row adds to one (probability of movement)
    summed = map(j -> sum(work[:, j]), unique(Js))
    summed[summed .== 0] .= 1.0
    work.nzval ./= summed
    return work
end
function genlookups(epienv::GridEpiEnv, mov::AlwaysMovement)
    total_size = (size(epienv.active, 1) * size(epienv.active, 2))
    # Generate grid ids and x,y coords for active cells only
    grid_locs = 1:total_size
    activity = epienv.active[1:end]
    grid_locs = grid_locs[activity]
    xys = convert_coords.(grid_locs, size(epienv.active, 2))

    # Collate all movement related parameters
    grid_size = _getgridsize(epienv.habitat)
    sd = [(2 .* k.dist) ./ sqrt(pi) for k in mov.kernels][activity]
    relsize =  grid_size ./ sd
    thresh = [k.thresh for k in mov.kernels][activity]
    grid_size /= unit(grid_size)

    # Calculate lookup probabilities for each grid location
    res = map((i, r, t) -> EcoSISTEM.genlookups(i, grid_locs, xys, grid_size, r, t, epienv), grid_locs, relsize, thresh)

    # Column vectors are source grid cells (repeated for each destination calculated)
    Js = vcat([fill(grid_locs[r], length(res[r][1])) for r in eachindex(res)]...)
    # Row vectors are destination grid cells
    Is = vcat([r[1] for r in res]...)
    Vs = vcat([r[2] for r in res]...)
    return sparse(Int32.(Is), Int32.(Js), Vs, total_size, total_size)
end

function genlookups(from::Int64, to::Vector{Int64}, xys::Array{Tuple{Int64,Int64},1}, grid_size::Float64, relsize::Float64, thresh::Float64, epienv::GridEpiEnv)
    x, y = xys[to .== from][1]
    maxX = ceil(Int64, x + 1/relsize); minX = ceil(Int64, x - 1/relsize)
    maxY = floor(Int64, y + 1/relsize); minY = floor(Int64, y - 1/relsize)
    keep = [(i[1] <= maxX) & (i[2] <= maxY) & (i[1] >= minX) & (i[2] >= minY) for i in xys]
    to = to[keep]
    probs = [_lookup((x = x, y = y), (x = i[1], y = i[2]), relsize, _gaussian_disperse) for i in xys[keep]]
    keep = probs .> thresh
    probs = probs[keep]
    probs ./= sum(probs)
    return to[keep], probs
end

function _lookup(from::NamedTuple, to::NamedTuple, relSquareSize::Float64, dispersalfn::F) where {F<:Function}
    return calc_prob = hcubature(r -> dispersalfn(r),
      [from.y *relSquareSize - relSquareSize, from.x * relSquareSize - relSquareSize, to.y * relSquareSize - relSquareSize, to.x * relSquareSize - relSquareSize],
      [from.y * relSquareSize, from.x * relSquareSize, to.y * relSquareSize, to.x * relSquareSize],
      maxevals= 100, rtol = 0.01)[1] / relSquareSize^2
end
