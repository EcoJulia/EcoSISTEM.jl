using Diversity
if VERSION > v"1.0.0"
    using HCubature
else
    using Cubature
end
using DataFrames
using Unitful
using EcoSISTEM.Units
using Missings
using Compat
using RecipesBase


import Diversity: _calcabundance

abstract type AbstractLookup end
"""
    Lookup

Lookup houses information on `x`, `y` grid locations and the probability of
occurrence at the location for the species in question `p`. `pnew` and `moves`
are initially empty storage and written over by the movement step in update!().
`pnew` is the recalculated probability based on which directions are available
and `moves` is the number of moves to that grid location in that step.
"""
mutable struct Lookup
  x::Vector{Int64}
  y::Vector{Int64}
  p::Vector{Float64}
  pnew::Vector{Float64}
  moves::Vector{Int64}
end

mutable struct SpeciesLookup  <: AbstractLookup
    species::Vector{Lookup}
end

abstract type AbstractCache end
"""
    Cache

Cache houses an integer array of moves made by all species in a timestep for the
update! function, `netmigration`.
"""
mutable struct Cache <: AbstractCache
  netmigration::Array{Int64, 2}
  totalE::Matrix{Float64}
  valid::Bool
end


Lookup(df::DataFrame) = Lookup(df[!, :X], df[!, :Y], df[!, :Prob],
zeros(Float64, nrow(df)),zeros(Int64, nrow(df)))


"""
    tematch(sppl::SpeciesList, abenv::AbstractAbiotic)

Function to check that the types of a trait list and habitat list are
the same for a species list (`sppl`) and abiotic environment (`abenv`).
"""
function tematch(sppl::SpeciesList, abenv::AbstractAbiotic)
    (eltype(sppl.species.traits) == eltype(abenv.habitat)) &&
    (iscontinuous(sppl.species.traits) == iscontinuous(abenv.habitat))
end
"""
    trmatch(sppl::SpeciesList, traitrel::AbstractTraitRelationship)

Function to check that the types of a trait list and trait relationship list are
the same for a species list (`sppl`) and trait relationship (`traitrel`).
"""
function trmatch(sppl::SpeciesList, traitrel::AbstractTraitRelationship)
    eltype(sppl.species.traits) == eltype(traitrel) &&
    (iscontinuous(sppl.species.traits) == iscontinuous(traitrel))
end



function create_cache(sppl::SpeciesList, ml::GridLandscape)
  nm = zeros(Int64, size(ml.matrix))
  totalE = zeros(Float64, (size(ml.matrix, 2), numrequirements(typeof(sppl.species.requirement))))
  return Cache(nm, totalE, false)
end

abstract type
    AbstractEcosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
        TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <: AbstractMetacommunity{Float64, Matrix{Int64},
                                        Matrix{Float64}, SL, Part}
end


mutable struct Ecosystem{L <: AbstractLandscape, Part <: AbstractPartition, SL <: AbstractSpeciesList,
    TR <: AbstractTraitRelationship, LU <: AbstractLookup, C <: AbstractCache} <: AbstractEcosystem{L, Part, SL, TR, LU, C}
  abundances::L
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::LU
  cache::C
  transitions::Union{Missing, TransitionList}
end

function Ecosystem(popfun::F, spplist::SpeciesList{T, Req}, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship) where {F<:Function, T, Req}

    # Create matrix landscape of zero abundances
  ml = emptygridlandscape(abenv, spplist)
  # Populate this matrix with species abundances
  popfun(ml, spplist, abenv, rel)
  # Create lookup table of all moves and their probabilities
  lookup = SpeciesLookup(collect(map(k -> genlookups(abenv.habitat, k), getkernels(spplist.species.movement))))
  cache = create_cache(spplist, ml)
  transitions = create_transitions(spplist, abenv)
  return Ecosystem{typeof(ml), typeof(abenv), typeof(spplist), typeof(rel), typeof(lookup), typeof(cache)}(ml, spplist, abenv,
  missing, rel, lookup, cache, transitions)
end

function Ecosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
   rel::AbstractTraitRelationship)
   return Ecosystem(populate!, spplist, abenv, rel)
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
function Base.isapprox(eco_1::AbstractEcosystem, eco_2::AbstractEcosystem; kwargs...)
    return isapprox(eco_1.abundances, eco_2.abundances; kwargs...)
end

save(path::String, system::Ecosystem) = JLSO.save(path, :ecosystem => system)
load(path::String, obj_type::Type{Ecosystem}) = JLSO.load(path)[:ecosystem]

function addspecies!(eco::Ecosystem, abun::Int64)
    eco.abundances.matrix = vcat(eco.abundances.matrix, zeros(1, size(eco.abundances.matrix, 2)))
    eco.abundances.grid = reshape(eco.abundances.matrix, (counttypes(eco.spplist, true)+1, _getdimension(eco.abenv.habitat)...))
    repopulate!(eco, abun)
    push!(eco.spplist.species.names, string.(counttypes(eco.spplist, true)+1))
    append!(eco.spplist.species.abun, abun)
    append!(eco.spplist.species.native, true)
    addtraits!(eco.spplist.species.traits)
    addmovement!(eco.spplist.species.movement)
    addparams!(eco.spplist.params)
    addrequirement!(eco.spplist.species.requirement)
    addtypes!(eco.spplist.species.types)
end
function addtraits!(tr::GaussTrait)
    append!(tr.mean, tr.mean[end])
    append!(tr.var, tr.var[end])
end

function addtraits!(tr::DiscreteTrait)
    append!(tr.val, rand(tr.val))
end

addmovement!(mv::AbstractMovement) = push!(mv.kernels, mv.kernels[end])

function addparams!(pr::AbstractParams)
    append!(pr.birth, pr.birth[end])
    append!(pr.death, pr.death[end])
end

addrequirement!(rq::AbstractRequirement) = append!(rq.energy, rq.energy[end])

function addtypes!(ut::UniqueTypes)
    ut = UniqueTypes(ut.num+1)
end

import Diversity.API: _getabundance
function _getabundance(eco::AbstractEcosystem, input::Bool)
    if input
        return eco.abundances.matrix
    else
        return _calcabundance(_gettypes(eco), eco.abundances.matrix / sum(eco.abundances.matrix))[1]
    end
end

import Diversity.API: _getmetaabundance
function _getmetaabundance(eco::AbstractEcosystem)
  return sumoversubcommunities(eco, _getabundance(eco))
end

import Diversity.API: _getpartition
function _getpartition(eco::AbstractEcosystem)
  return eco.abenv
end
import Diversity.API: _gettypes
function _gettypes(eco::AbstractEcosystem)
    return eco.spplist
end
import Diversity.API: _getordinariness!
function _getordinariness!(eco::AbstractEcosystem)
    if ismissing(eco.ordinariness)
        relab = getabundance(eco, false)
        eco.ordinariness = _calcordinariness(eco.spplist, relab)
    end
    return eco.ordinariness
end

import Diversity.API._getscale
function _getscale(eco::AbstractEcosystem)
    return _calcabundance(_gettypes(eco), getabundance(eco, false))[2]
end

function invalidatecaches!(eco::AbstractEcosystem)
    _invalidatecaches!(eco, eco.cache)
end
function _invalidatecaches!(eco::A, cache::Cache) where A <: AbstractEcosystem
    eco.ordinariness = missing
    eco.cache.netmigration .= 0
    eco.cache.valid = false
end

"""
    gettraitrel(eco::Ecosystem)

Function to extract trait relationships.
"""
function gettraitrel(eco::AbstractEcosystem)
  return eco.relationship
end

"""
    gethabitat(eco::Ecosystem)

Function to extract habitat from Ecosystem object.
"""
function gethabitat(eco::AbstractEcosystem)
  return eco.abenv.habitat
end
"""
    getbudget(eco::Ecosystem)

Function to extract budget from Ecosystem object.
"""
function getbudget(eco::AbstractEcosystem)
    return _getbudget(eco.abenv.budget)
end

"""
    getsize(eco::Ecosystem)

Function to extract size of habitat from Ecosystem object.
"""
function getsize(eco::AbstractEcosystem)
  return _getsize(eco.abenv.habitat)
end

"""
    getgridsize(eco::Ecosystem)

Function to extract grid cell size of habitat from Ecosystem object.
"""
function getgridsize(eco::AbstractEcosystem)
  return _getgridsize(eco.abenv.habitat)
end

"""
    getdimension(eco::Ecosystem)

Function to extract dimension of habitat from Ecosystem object.
"""
function getdimension(eco::AbstractEcosystem)
    return _getdimension(eco.abenv.habitat)
end

"""
    getdispersaldist(eco::Ecosystem)

Function to extract average dispersal distance of species from Ecosystem object.
Returns a vector of distances, unless a specific species is provided as a String
or Integer.
"""
function getdispersaldist(eco::A, sp::Int64) where A <: AbstractEcosystem
    dist = eco.spplist.species.movement.kernels[sp].dist
    return dist
end
function getdispersaldist(eco::A, sp::String) where A <: AbstractEcosystem
    sp = findfirst(eco.spplist.species.names .== sp)
    dist = eco.spplist.species.movement.kernels[sp].dist
    return dist
end

"""
    getdispersalvar(eco::Ecosystem)

Function to extract dispersal varaince of species from Ecosystem object.
Returns a vector of distances, unless a specific species is provided as a String
or Integer.
"""
function getdispersalvar(eco::A, sp::Int64) where A <: AbstractEcosystem
    var = (eco.spplist.species.movement.kernels[sp].dist)^2 * pi / 4
    return var
end
function getdispersalvar(eco::A, sp::String) where A <: AbstractEcosystem
    sp = findfirst(eco.spplist.species.names .== sp)
    var = (eco.spplist.species.movement.kernels[sp].dist)^2 * pi / 4
    return var
end

"""
    getlookup(eco::Ecosystem)

Function to extract movement lookup table of species from Ecosystem object.
"""
function getlookup(eco::A, sp::Int64) where A <: AbstractEcosystem
    return _getlookup(eco.lookup, sp)
end

function _getlookup(lookup::SpeciesLookup, sp::Int64)
    return lookup.species[sp]
end

"""
    resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝐓^-1)})

Function to reset the rate of habitat change for a species.
"""
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, typeof(𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, Unitful.Dimensions{()})
end
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, typeof(𝚯*𝐓^-1)})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, typeof(dimension(1K)))
end
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, 𝐓^-1})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, Unitful.Dimensions{()})
end
function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, 𝚯*𝐓^-1})
    eco.abenv.habitat.change = HabitatUpdate(
    eco.abenv.habitat.change.changefun, rate, typeof(dimension(1K)))
end

function resettime!(eco::AbstractEcosystem)
    _resettime!(eco.abenv.habitat)
end


function _symmetric_grid(grid::DataFrame)
   for x in 1:nrow(grid)
     if grid[x, 1] != grid[x, 2]
       push!(grid, hcat(grid[x, 2], grid[x, 1] , grid[x, 3]))
     end
   end
   for x in 1:nrow(grid)
     if (grid[x, 1] > 0)
       push!(grid, hcat(-grid[x, 1], grid[x, 2] , grid[x, 3]))
     end
     if (grid[x, 2] > 0)
       push!(grid, hcat(grid[x, 1], -grid[x, 2] , grid[x, 3]))
     end
     if (grid[x, 1] > 0 && grid[x, 2] > 0)
       push!(grid, hcat(-grid[x, 1], -grid[x, 2] , grid[x, 3]))
     end
   end
   grid
 end

 # Define gaussian kernel function
function _gaussian_disperse(r)
  exp(-((r[3]-r[1])^2+(r[4]-r[2])^2)) / π
end

function _2Dt_disperse(r, b)
    return((b - 1)/(π)) * (1 + ((r[3]-r[1])^2+(r[4]-r[2])^2))^-b
end

"""
    genlookups(hab::AbstractHabitat, mov::GaussianMovement)

Function to generate lookup tables, which hold information on the probability
of moving to neighbouring squares.
"""
function genlookups(hab::AbstractHabitat, mov::GaussianKernel)
  sd = (2 * mov.dist) / sqrt(pi)
  relsize =  _getgridsize(hab) ./ sd
  m = maximum(_getdimension(hab))
  p = mov.thresh
  return Lookup(_lookup(relsize, m, p, _gaussian_disperse))
end
function genlookups(hab::AbstractHabitat, mov::LongTailKernel)
    sd = (2 * mov.dist) / sqrt(pi)
    relsize =  _getgridsize(hab) ./ sd
    m = maximum(_getdimension(hab))
    p = mov.thresh
    b = mov.shape
    return EcoSISTEM.Lookup(EcoSISTEM._lookup(relsize, m, p, b, EcoSISTEM._2Dt_disperse))
end

function _lookup(relSquareSize::Float64, maxGridSize::Int64,
                pThresh::Float64, dispersalfn::F) where {F<:Function}
  # Create empty array
  lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

  # Loop through directions until probability is below threshold
  k = 0
  m = 0
  count = 0
  while (k <= maxGridSize && m <= maxGridSize)
    count = count + 1
    calc_prob = hcubature(r -> dispersalfn(r),
      [0, 0, k*relSquareSize, m*relSquareSize],
      [relSquareSize, relSquareSize, (k+1)*relSquareSize, (m+1)*relSquareSize],
      maxevals= 10000)[1] / relSquareSize^2
    if m == 0 && calc_prob < pThresh
      break
    end
    if count == 1
      push!(lookup_tab, [k m calc_prob])
      k = k + 1
    elseif (calc_prob > pThresh && m <= k)
      push!(lookup_tab, [k m calc_prob])
      m = m + 1
    else
      m = 0
      k = k + 1
    end
  end
  # If no probabilities can be calculated, threshold is too high
  nrow(lookup_tab) != 0 || error("probability threshold too high")
  # Find all other directions
  lookup_tab = _symmetric_grid(lookup_tab)
  #info(sum(lookup_tab[:, 3]))
  # Normalise
  lookup_tab[!, :Prob] = lookup_tab[!, :Prob]/sum(lookup_tab[!, :Prob])
  lookup_tab
end


function _lookup(relSquareSize::Float64, maxGridSize::Int64,
                pThresh::Float64, b::Float64, dispersalfn::F
                ) where {F<:Function}
  # Create empty array
  lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

  # Loop through directions until probability is below threshold
  k = 0
  m = 0
  count = 0
  while (k <= maxGridSize && m <= maxGridSize)
    count = count + 1
    calc_prob = hcubature(r -> dispersalfn(r, b),
      [0, 0, k*relSquareSize, m*relSquareSize],
      [relSquareSize, relSquareSize, (k+1)*relSquareSize, (m+1)*relSquareSize],
      maxevals=10000)[1] / relSquareSize^2
    if m == 0 && calc_prob < pThresh
      break
    end
    if count == 1
      push!(lookup_tab, [k m calc_prob])
       k = k + 1
    elseif (calc_prob > pThresh && m <= k)
      push!(lookup_tab, [k m calc_prob])
      m = m + 1
    else
      m = 0
      k = k + 1
    end
  end
  # If no probabilities can be calculated, threshold is too high
  nrow(lookup_tab) != 0 || error("probability threshold too high")
  # Find all other directions
  lookup_tab = _symmetric_grid(lookup_tab)
  #info(sum(lookup_tab[:, 3]))
  # Normalise
  lookup_tab[!, :Prob] = lookup_tab[!, :Prob]/sum(lookup_tab[!, :Prob])
  lookup_tab
end


"""
    CachedEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
        TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}

CachedEcosystem houses the same information as Ecosystem (see ?Ecosystem), but
holds the time period abundances as a CachedGridLandscape, so that they may
be present or missing.
"""
mutable struct CachedEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
    TR <: AbstractTraitRelationship} <: AbstractEcosystem{CachedGridLandscape, Part, SL, TR, SpeciesLookup, Cache}
  abundances::CachedGridLandscape
  spplist::SL
  abenv::Part
  ordinariness::Union{Matrix{Float64}, Missing}
  relationship::TR
  lookup::SpeciesLookup
  cache::Cache
end

"""
    CachedEcosystem(eco::Ecosystem, outputfile::String, rng::StepRangeLen)

Function to create a CachedEcosystem given an existing ecosystem, `eco`,
output folder to which the simulations are saved, `outputfile`, and a range of
times over which to simulate, `rng`.
"""
function CachedEcosystem(eco::Ecosystem, outputfile::String, rng::StepRangeLen)
    if size(eco.abenv.habitat, 3) > 1
        size(eco.abenv.habitat, 3) == length(rng) || error("Time range does not match habitat")
    end
    abundances = CachedGridLandscape(outputfile, rng)
    abundances.matrix[1] = eco.abundances
  CachedEcosystem{typeof(eco.abenv), typeof(eco.spplist), typeof(eco.relationship)}(abundances,
  eco.spplist, eco.abenv, eco.ordinariness, eco.relationship, eco.lookup, eco.cache)
end



using JLSO
using SparseArrays

function create_cache(sppl::SpeciesList, ml::EpiLandscape)
  vm = zeros(Float64, size(ml.matrix))
  return EpiCache(vm, false)
end

mutable struct EpiCache <: AbstractCache
  virusmigration::Array{Float64, 2}
  initial_infected::Int64
  ordered_active::Vector{Int64}
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

function Ecosystem(abundances::EpiLandscape{U, VecRNGType}, epilist::EL, epienv::EE,
    ordinariness::Union{Matrix{Float64}, Missing}, relationship::ER, lookup::EpiLookup,
    vm::Array{Float64, 2}, initial_infected::Int64, valid::Bool, transitions::Union{Missing, TransitionList}
    ) where {U <: Integer, VecRNGType <: AbstractVector{<:Random.AbstractRNG},
    EE <: AbstractEpiEnv, EL <: SpeciesList, ER <: AbstractTraitRelationship}
  total_pop = sum(abundances.matrix, dims = 1)[1, :]
  sorted_grid_ids = sortperm(total_pop, rev = true)
  sorted_grid_ids = sorted_grid_ids[total_pop[sorted_grid_ids] .> 0]
  cache = EpiCache(vm, initial_infected, sorted_grid_ids, valid)
  return Ecosystem(abundances, epilist, epienv, ordinariness, relationship, lookup, cache, transitions)
end

function Ecosystem(popfun::F, epilist::SpeciesList, epienv::GridEpiEnv,
      rel::AbstractTraitRelationship, intnum::U; initial_infected = 0,
      rngtype::Type{R} = Random.MersenneTwister,
      transitions = missing) where {F<:Function, U <: Integer, R <: Random.AbstractRNG}

  # Create matrix landscape of zero abundances
  ml = emptyepilandscape(epienv, epilist, intnum, rngtype)

  # Populate this matrix with species abundances
  popfun(ml, epilist, epienv, rel)
  initial_pop = sum(ml.matrix, dims = 1)
  # Create lookup table of all moves and their probabilities
  home_lookup = genlookups(epienv, epilist.species.movement.home)
  work_lookup = genlookups(epienv, epilist.species.movement.work, initial_pop)
  lookup = EpiLookup(home_lookup, work_lookup)
  vm = zeros(Float64, size(ml.matrix))
  return Ecosystem(ml, epilist, epienv, missing, rel, lookup, vm, initial_infected, false, transitions)
end

function Ecosystem(epilist::SpeciesList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
        intnum::U = Int64(1); initial_infected = 0, rngtype::Type{R} = Random.MersenneTwister,
        transitions = missing
        ) where {U <: Integer, R <: Random.AbstractRNG}
    return Ecosystem(populate!, epilist, epienv, rel, intnum, initial_infected = initial_infected, rngtype = rngtype, transitions = transitions)
end

function Ecosystem(epilist::SpeciesList, epienv::GridEpiEnv, rel::AbstractTraitRelationship,
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
    home_lookup = genlookups(epienv, epilist.species.movement.home)
    work_lookup = genlookups(epienv, epilist.species.movement.work, initial_population[1:end])
    lookup = EpiLookup(home_lookup, work_lookup)

    vm = zeros(Float64, size(ml.matrix))

    epi = Ecosystem(ml, epilist, epienv, missing, rel, lookup, vm, initial_infected, false, transitions)

    # Add in the initial susceptible population
    # TODO Need to fix code so it doesn't rely on name of susceptible class
    idx = findfirst(occursin.("Susceptible", epilist.species.names))
    if idx == nothing
        msg = "epilist has no Susceptible category. epilist.names = $(epilist.species.names)"
        throw(ArgumentError(msg))
    end
    # Modify active cells based on new population
    initial_population = convert_population(initial_population, intnum)
    epi.abundances.grid[idx, :, :] .+= initial_population

    total_pop = sum(human(epi.abundances), dims = 1)[1, :]
    sorted_grid_ids = sortperm(total_pop, rev = true)
    sorted_grid_ids = sorted_grid_ids[total_pop[sorted_grid_ids] .> 0]
    epi.cache.ordered_active = sorted_grid_ids
    return epi
end

function _invalidatecaches!(eco::Ecosystem, cache::EpiCache)
    eco.ordinariness = missing
    eco.cache.virusmigration .= 0
    eco.cache.valid = false
end

function getlookup(epi::Ecosystem, id::Int64, movetype::MovementType)
    if movetype == homeMovement
        return epi.lookup.homelookup[id, :]
    elseif movetype == workMovement
        return epi.lookup.worklookup[id, :]
    else
        return error("No other movement types currently implemented")
    end
end

function _getlookup(lookup::EpiLookup, id::Int64)
    return lookup.homelookup[id, :], lookup.worklookup[id, :]
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
