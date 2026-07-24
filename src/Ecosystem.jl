# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity
using HCubature
using DataFrames
using Unitful
using EcoSISTEM.Units
using Missings
using RecipesBase
using Random

using Diversity.API: _calcabundance

"""
    makerngs(seed::Integer, n::Integer)

Build a vector of `n` independent, deterministically-seeded random number
generators, one per species. Species `j` is seeded as `Xoshiro(hash((seed, j)))`
so its random stream is a pure function of `(seed, j)` — independent of how
species are distributed across threads or MPI processes. This is what makes
simulation results reproducible across different thread and process counts (each
species is always processed by exactly one task on one rank, drawing in a fixed
cell order). See [`getrng`](@ref).

Note: this per-species scheme is sufficient only because no single species' draws
are ever split across ranks/tasks. If a species' cells were ever partitioned
across ranks, a per-`(species, cell)` counter-based generator would be needed
instead.
"""
function makerngs(seed::Integer, n::Integer)
    return [Random.Xoshiro(hash((seed, j))) for j in 1:n]
end

"""
    Lookup

Lookup houses information on `x`, `y` grid locations and the probability of
occurrence at the location for the species in question `p`. `pnew` and `moves`
are initially empty storage and written over by the movement step in
[`update!`](@ref)(). `pnew` is the recalculated probability based on which
directions are available and `moves` is the number of moves to that grid
location in that step.
"""
struct Lookup
    x::Vector{Int64}
    y::Vector{Int64}
    p::Vector{Float64}
    pnew::Vector{Float64}
    moves::Vector{Int64}
end

"""
    Cache

Cache houses an integer array of moves made by all species in a timestep for the
[`update!`](@ref) function, `netmigration`.
"""
mutable struct Cache
    netmigration::Matrix{Int64}
    totalE::Matrix{Float64}
    valid::Bool
end

function Lookup(df::DataFrame)
    return Lookup(df[!, :X],
                  df[!, :Y],
                  df[!, :Prob],
                  zeros(Float64, nrow(df)),
                  zeros(Int64, nrow(df)))
end

# Check the abundance matrix `m` (species × subcommunities) has one row per species in `sppl` and one
# column per subcommunity in `part`. A pure dimension check.
function _mcmatch(m::AbstractMatrix, sppl::SpeciesList, part::AbstractHabitat)
    return _counttypes(sppl, true) == size(m, 1) &&
           _countsubcommunities(part) == size(m, 2)
end

"""
    tematch(sppl::SpeciesList, habitat::AbstractHabitat)

Check that the types of a tolerance list and regime list are the same for a species
list (`sppl`) and abiotic environment (`habitat`).
"""
function tematch(sppl::SpeciesList, habitat::AbstractHabitat)
    return (eltype(sppl.tolerance) == eltype(habitat.regime)) &&
           (iscontinuous(sppl.tolerance) == iscontinuous(habitat.regime))
end

"""
    nfmatch(sppl::SpeciesList, traitrel::AbstractNicheFit)

Check that the types of a tolerance list and trait nichefit list are the same
for a species list (`sppl`) and trait nichefit (`traitrel`).
"""
function nfmatch(sppl::SpeciesList, traitrel::AbstractNicheFit)
    return eltype(sppl.tolerance) == eltype(traitrel) &&
           (iscontinuous(sppl.tolerance) == iscontinuous(traitrel))
end

# Format an `iscontinuous(...)` result — a `Bool` for a single trait/regime/nichefit, or a
# `Vector{Bool}` for a collection — as "continuous"/"discrete" label(s).
_kindlabel(b::Bool) = b ? "continuous" : "discrete"
_kindlabel(bs::AbstractVector{Bool}) = "[" * join(_kindlabel.(bs), ", ") * "]"

"""
    AbstractEcosystem{Part <: AbstractHabitat, SL <: SpeciesList,
        NF <: AbstractNicheFit} <: AbstractMetacommunity{Float64,
            Matrix{Int64}, Matrix{Float64}, SL, Part}

Abstract supertype for all ecosystem types and a subtype of
AbstractMetacommunity.
"""
abstract type AbstractEcosystem{Part <: AbstractHabitat,
                                SL <: SpeciesList,
                                NF <: AbstractNicheFit} <:
              AbstractMetacommunity{Float64, Matrix{Int64}, Matrix{Float64}, SL,
                                    Part} end

"""
    Ecosystem{Part <: AbstractHabitat} <:
       AbstractEcosystem{Part, SL, NF}

Ecosystem houses information on species and their interaction with their
environment. For species, it holds abundances and locations, as well as
properties such as trait information, `spplist`, and movement types, `lookup`.
For environments, it provides information on environmental conditions and
available resources,`habitat`. Finally, there is a slot for the nichefit
between the environment and the characteristics of the species, `nichefit`.
"""
mutable struct Ecosystem{Part <: AbstractHabitat,
                         SL <: SpeciesList,
                         NF <: AbstractNicheFit} <:
               AbstractEcosystem{Part, SL, NF}
    abundances::GridLandscape
    spplist::SL
    habitat::Part
    ordinariness::Union{Matrix{Float64}, Missing}
    nichefit::NF
    lookup::Vector{Lookup}
    cache::Cache
    rngs::Vector{Random.Xoshiro}

    function Ecosystem{Part, SL, NF}(abundances::GridLandscape,
                                     spplist::SL,
                                     habitat::Part,
                                     ordinariness::Union{Matrix{Float64},
                                                         Missing},
                                     nichefit::NF,
                                     lookup::Vector{Lookup},
                                     cache::Cache,
                                     rngs::Vector{Random.Xoshiro}) where {Part <:
                                                                          AbstractHabitat,
                                                                          SL <:
                                                                          SpeciesList,
                                                                          NF <:
                                                                          AbstractNicheFit}
        tematch(spplist, habitat) ||
            error("Species tolerances and regime are incompatible: tolerances are " *
                  "$(_kindlabel(iscontinuous(spplist.tolerance))) " *
                  "$(eltype(spplist.tolerance)), the regime is " *
                  "$(_kindlabel(iscontinuous(habitat.regime))) " *
                  "$(eltype(habitat.regime)). Pair a continuous trait (e.g. " *
                  "NicheTolerance) with a continuous regime (simplehabitat / " *
                  "tempgradhabitat), or a discrete trait (DiscreteTolerance) with a " *
                  "discrete regime (simplenichehabitat).")
        nfmatch(spplist, nichefit) ||
            error("Species tolerances and the trait nichefit are incompatible: " *
                  "tolerances are $(_kindlabel(iscontinuous(spplist.tolerance))) " *
                  "$(eltype(spplist.tolerance)), the nichefit is " *
                  "$(_kindlabel(iscontinuous(nichefit))) " *
                  "$(eltype(nichefit)). Use NicheSuitability with continuous tolerances, " *
                  "or MatchSuitability / LandCoverSuitability with discrete tolerances.")
        _mcmatch(abundances.matrix, spplist, habitat) ||
            error("Dimension mismatch: the abundance matrix " *
                  "($(size(abundances.matrix, 1)) × $(size(abundances.matrix, 2))) " *
                  "does not match the $(_counttypes(spplist, true)) species and " *
                  "$(_countsubcommunities(habitat)) subcommunities of the species " *
                  "list and environment.")
        return new{Part, SL, NF}(abundances,
                                 spplist,
                                 habitat,
                                 ordinariness,
                                 nichefit,
                                 lookup,
                                 cache,
                                 rngs)
    end
end

@recipe function f(::AbstractMovement, eco::AbstractEcosystem, sp::Int64)
    l = eco.lookup[sp]
    maxX = maximum(l.x)
    maxY = maximum(l.y)
    x, y = round(Int64, maxX / 2), round(Int64, maxY / 2)
    # Can't go over maximum dimension
    valid = findall((l.x .> -x) .& (l.y .> -y) .& (l.x .<= (maxX - x)) .&
                    (l.y .<= (maxY - y)))
    probs = l.p[valid]
    probs ./= sum(probs)
    xs = (l.x[valid] .+ x)
    ys = (l.y[valid] .+ y)
    A = zeros(maxX, maxY)
    for i in eachindex(xs)
        A[xs[i], ys[i]] .= probs[i]
    end
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> "Movement kernel (km)"
    return xrange(getregime(eco)), yrange(getregime(eco)), A
end

"""
    Ecosystem(spplist::SpeciesList, habitat::GridHabitat,
        nichefit::AbstractNicheFit)

Create an `Ecosystem` given a species list, an abiotic environment and trait
nichefit. An optional population function can be added, `popfun`, which
defaults to generic random filling of the ecosystem. A `seed` may be supplied to
make the run reproducible: it deterministically seeds one random number
generator per species (see [`makerngs`](@ref)), so results are identical
regardless of the number of threads used. If no `seed` is given, one is drawn at
random.
"""
function Ecosystem(popfun::F,
                   spplist::SpeciesList{T, DM},
                   habitat::GridHabitat,
                   nichefit::AbstractNicheFit;
                   seed::Integer = rand(UInt64)) where {F <: Function, T, DM}

    # Check there is enough resource to support number of individuals at set up
    #all(getdemand(spplist) .<= getavailablesupply(habitat)) ||
    #error("Environment does not have enough resource to support species")
    # Create matrix landscape of zero abundances
    ml = emptygridlandscape(habitat, spplist)
    # One deterministically-seeded RNG per species, so births/deaths/dispersal
    # and the initial population draw are reproducible across thread counts
    rngs = makerngs(seed, size(ml.matrix, 1))
    # Populate this matrix with species abundances
    popfun(ml, spplist, habitat, nichefit, rngs)
    # Create lookup table of all moves and their probabilities
    lookup_tab = collect(map(k -> genlookups(habitat.regime, k),
                             getkernels(spplist.movement)))
    nm = zeros(Int64, size(ml.matrix))
    totalE = zeros(Float64, (size(ml.matrix, 2), numdemands(DM)))
    return Ecosystem{typeof(habitat), typeof(spplist), typeof(nichefit)}(ml,
                                                                         spplist,
                                                                         habitat,
                                                                         missing,
                                                                         nichefit,
                                                                         lookup_tab,
                                                                         Cache(nm,
                                                                               totalE,
                                                                               false),
                                                                         rngs)
end

function Ecosystem(spplist::SpeciesList,
                   habitat::GridHabitat,
                   nichefit::AbstractNicheFit;
                   seed::Integer = rand(UInt64))
    return Ecosystem(populate!, spplist, habitat, nichefit; seed = seed)
end
@doc (@doc Ecosystem) Ecosystem(::SpeciesList,
                                ::GridHabitat,
                                ::AbstractNicheFit)

"""
    addspecies!(eco::Ecosystem, abun::Int64)

Add a new species to an existing [`Ecosystem`](@ref) with initial abundance
`abun`, copying trait, movement, parameter, demand, and type information
from the last existing species.
"""
function addspecies!(eco::Ecosystem, abun::Int64)
    eco.abundances.matrix = vcat(eco.abundances.matrix,
                                 zeros(1, size(eco.abundances.matrix, 2)))
    eco.abundances.grid = reshape(eco.abundances.matrix,
                                  (counttypes(eco.spplist, true) + 1,
                                   _getdimension(eco.habitat.regime)...))
    # Give the new species its own RNG stream, derived from the previous last one
    push!(eco.rngs, Random.Xoshiro(rand(eco.rngs[end], UInt64)))
    repopulate!(eco, abun)
    push!(eco.spplist.names, string.(counttypes(eco.spplist, true) + 1))
    append!(eco.spplist.abun, abun)
    append!(eco.spplist.native, true)
    addtolerance!(eco.spplist.tolerance)
    addmovement!(eco.spplist.movement)
    addparams!(eco.spplist.params)
    adddemand!(eco.spplist.demand)
    return addtypes!(eco.spplist.types)
end

function addtolerance!(tolerance::NicheTolerance)
    return push!(tolerance.dists, tolerance.dists[end])
end

function addtolerance!(tolerance::DiscreteTolerance)
    return append!(tolerance.val, rand(tolerance.val))
end

addmovement!(mv::AbstractMovement) = push!(mv.kernels, mv.kernels[end])

function addparams!(params::AbstractParams)
    append!(params.birth, params.birth[end])
    return append!(params.death, params.death[end])
end

adddemand!(d::AbstractDemand) = append!(d.resource, d.resource[end])

function addtypes!(ut::UniqueTypes)
    return ut = UniqueTypes(ut.num + 1)
end

"""
    CachedEcosystem{Part <: AbstractHabitat, SL <: SpeciesList,
        NF <: AbstractNicheFit} <: AbstractEcosystem{Part, SL, NF}

CachedEcosystem houses the same information as [`Ecosystem`](@ref) (see
?Ecosystem), but holds the time period abundances as a
[`CachedGridLandscape`](@ref), so that they may be present or missing.
"""
mutable struct CachedEcosystem{Part <: AbstractHabitat,
                               SL <: SpeciesList,
                               NF <: AbstractNicheFit} <:
               AbstractEcosystem{Part, SL, NF}
    abundances::CachedGridLandscape
    spplist::SL
    habitat::Part
    ordinariness::Union{Matrix{Float64}, Missing}
    nichefit::NF
    lookup::Vector{Lookup}
    cache::Cache
    rngs::Vector{Random.Xoshiro}
end

"""
    CachedEcosystem(eco::Ecosystem, outputfile::String, times::StepRangeLen;
                    saveinterval::Unitful.Time = step(times))

Create a CachedEcosystem given an existing [`Ecosystem`](@ref), `eco`, an output
folder to which the simulations are saved, `outputfile`, and a range of times
over which to simulate, `times`. The step of `times` is the simulation timestep;
`saveinterval` controls how often checkpoints are written to disk (a multiple of
the timestep, defaulting to every step). Because the simulation always advances
by the timestep, results are independent of `saveinterval`.
"""
function CachedEcosystem(eco::Ecosystem, outputfile::String,
                         times::StepRangeLen;
                         saveinterval::Unitful.Time = step(times))
    if size(eco.habitat.regime, 3) > 1
        size(eco.habitat.regime, 3) == length(times) ||
            error("Time range does not match regime")
    end
    abundances = CachedGridLandscape(outputfile, times;
                                     saveinterval = saveinterval)
    abundances.matrix[1] = eco.abundances
    return CachedEcosystem{typeof(eco.habitat), typeof(eco.spplist),
                           typeof(eco.nichefit)}(abundances,
                                                 eco.spplist,
                                                 eco.habitat,
                                                 eco.ordinariness,
                                                 eco.nichefit,
                                                 eco.lookup,
                                                 eco.cache,
                                                 eco.rngs)
end

import Diversity.API: _getabundance
function _getabundance(eco::AbstractEcosystem, raw::Bool)
    if raw
        return eco.abundances.matrix
    else
        return _calcabundance(_gettypes(eco),
                              eco.abundances.matrix /
                              sum(eco.abundances.matrix))[1]
    end
end

function _getabundance(cache::CachedEcosystem, raw::Bool)
    if all(ismissing.(cache.abundances.matrix))
        error("Abundances are missing")
    else
        id = findall(.!ismissing.(cache.abundances.matrix))[end]
        abun = cache.abundances.matrix[id]
    end

    if raw
        return abun.matrix
    else
        return abun.matrix / sum(abun.matrix)
    end
end

import Diversity.API: _getpartition
function _getpartition(eco::AbstractEcosystem)
    return eco.habitat
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
    eco.ordinariness = missing
    eco.cache.netmigration .= 0
    return eco.cache.valid = false
end

"""
    getnichefit(eco::Ecosystem)

Extract niche fits.
"""
function getnichefit(eco::AbstractEcosystem)
    return eco.nichefit
end

"""
    getregime(eco::Ecosystem)

Extract regime from [`Ecosystem`](@ref) object.
"""
function getregime(eco::AbstractEcosystem)
    return eco.habitat.regime
end

"""
    getsupply(eco::Ecosystem)

Extract supply from [`Ecosystem`](@ref) object.
"""
function getsupply(eco::AbstractEcosystem)
    return _getsupply(eco.habitat.supply)
end

"""
    getsize(eco::Ecosystem)

Extract size of regime from [`Ecosystem`](@ref) object.
"""
function getsize(eco::AbstractEcosystem)
    return _getsize(eco.habitat.regime)
end

"""
    getgridsize(eco::Ecosystem)

Extract grid cell size of regime from [`Ecosystem`](@ref) object.
"""
function getgridsize(eco::AbstractEcosystem)
    return _getgridsize(eco.habitat.regime)
end

"""
    getdimension(eco::Ecosystem)

Extract dimension of regime from [`Ecosystem`](@ref) object.
"""
function getdimension(eco::AbstractEcosystem)
    return _getdimension(eco.habitat.regime)
end

"""
    getdispersaldist(eco::Ecosystem)

Extract average dispersal distance of species from [`Ecosystem`](@ref) object.
Returns a vector of distances, unless a specific species is provided as a String
or Integer.
"""
function getdispersaldist(eco::AbstractEcosystem, sp::Int64)
    dist = eco.spplist.movement.kernels[sp].dist
    return dist
end

function getdispersaldist(eco::AbstractEcosystem, sp::String)
    num = findall(eco.spplist.names .== sp)[1]
    return getdispersaldist(eco, num)
end
@doc (@doc getdispersaldist) getdispersaldist(::AbstractEcosystem, ::String)

"""
    getdispersalvar(eco::Ecosystem)

Extract dispersal varaince of species from [`Ecosystem`](@ref) object. Returns a
vector of distances, unless a specific species is provided as a String or
Integer.
"""
function getdispersalvar(eco::AbstractEcosystem, sp::Int64)
    var = (eco.spplist.movement.kernels[sp].dist)^2 * pi / 4
    return var
end

function getdispersalvar(eco::AbstractEcosystem, sp::String)
    num = findall(eco.spplist.names .== sp)[1]
    return getdispersalvar(eco, num)
end
@doc (@doc getdispersalvar) getdispersalvar(::AbstractEcosystem, ::String)
"""
    getlookup(eco::Ecosystem)

Extract movement lookup table of species from [`Ecosystem`](@ref) object.
"""
function getlookup(eco::AbstractEcosystem, sp::Int64)
    return eco.lookup[sp]
end

function getlookup(eco::AbstractEcosystem, sp::String)
    num = findall(eco.spplist.names .== sp)[1]
    return getlookup(eco, num)
end
@doc (@doc getlookup) getlookup(::AbstractEcosystem, ::String)

"""
    getrng(eco::AbstractEcosystem, sp::Int64)

Return the per-species random number generator for global species index `sp`.
Because each species has its own stream and is processed by exactly one task per
timestep, random draws are both thread-safe and reproducible independent of the
number of threads or MPI processes (see [`makerngs`](@ref)).
"""
function getrng(eco::AbstractEcosystem, sp::Int64)
    return eco.rngs[sp]
end

"""
    resetrate!(eco::Ecosystem, rate::Quantity{Float64, typeof(𝐓^-1)})

Reset the rate of regime change for a species.
"""
function resetrate!(eco::AbstractEcosystem,
                    rate::Quantity{Float64, typeof(𝐓^-1)})
    return eco.habitat.regime.dynamics = LayerUpdate(eco.habitat.regime.dynamics.changefun,
                                                     rate,
                                                     Unitful.Dimensions{()})
end

function resetrate!(eco::AbstractEcosystem,
                    rate::Quantity{Float64, typeof(𝚯 * 𝐓^-1)})
    return eco.habitat.regime.dynamics = LayerUpdate(eco.habitat.regime.dynamics.changefun,
                                                     rate,
                                                     typeof(dimension(1K)))
end

function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, 𝐓^-1})
    return eco.habitat.regime.dynamics = LayerUpdate(eco.habitat.regime.dynamics.changefun,
                                                     rate,
                                                     Unitful.Dimensions{()})
end

function resetrate!(eco::AbstractEcosystem, rate::Quantity{Float64, 𝚯 * 𝐓^-1})
    return eco.habitat.regime.dynamics = LayerUpdate(eco.habitat.regime.dynamics.changefun,
                                                     rate,
                                                     typeof(dimension(1K)))
end

function resettime!(eco::AbstractEcosystem)
    return _resettime!(eco.habitat.regime)
end

function _symmetric_grid(grid::DataFrame)
    for x in 1:nrow(grid)
        if grid[x, 1] != grid[x, 2]
            push!(grid, hcat(grid[x, 2], grid[x, 1], grid[x, 3]))
        end
    end
    for x in 1:nrow(grid)
        if (grid[x, 1] > 0)
            push!(grid, hcat(-grid[x, 1], grid[x, 2], grid[x, 3]))
        end
        if (grid[x, 2] > 0)
            push!(grid, hcat(grid[x, 1], -grid[x, 2], grid[x, 3]))
        end
        if (grid[x, 1] > 0 && grid[x, 2] > 0)
            push!(grid, hcat(-grid[x, 1], -grid[x, 2], grid[x, 3]))
        end
    end
    return grid
end

# Define gaussian kernel function
function _gaussian_disperse(r)
    return exp(-((r[3] - r[1])^2 + (r[4] - r[2])^2)) / π
end

function _2Dt_disperse(r, b)
    return ((b - 1) / (π)) * (1 + ((r[3] - r[1])^2 + (r[4] - r[2])^2))^-b
end

"""
    genlookups(regime::AbstractRegime, kernel::GaussianMovement)

Generate lookup tables, which hold information on the probability of moving to
neighbouring squares.
"""
function genlookups(regime::AbstractRegime, kernel::GaussianKernel)
    sd = (2 * kernel.dist) / sqrt(pi)
    relsize = uconvert(NoUnits, _getgridsize(regime) / sd)
    m = maximum(_getdimension(regime))
    p = kernel.thresh
    return Lookup(_lookup(relsize, m, p, _gaussian_disperse))
end

function genlookups(regime::AbstractRegime, kernel::LongTailKernel)
    sd = (2 * kernel.dist) / sqrt(pi)
    relsize = uconvert(NoUnits, _getgridsize(regime) / sd)
    m = maximum(_getdimension(regime))
    p = kernel.thresh
    b = kernel.shape
    return EcoSISTEM.Lookup(EcoSISTEM._lookup(relsize, m, p, b,
                                              EcoSISTEM._2Dt_disperse))
end

function _lookup(relSquareSize::Float64,
                 maxGridSize::Int64,
                 pThresh::Float64,
                 dispersalfn::F) where {F <: Function}
    # Create empty array
    lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

    # Loop through directions until probability is below threshold
    k = 0
    m = 0
    count = 0
    while (k <= maxGridSize && m <= maxGridSize)
        count = count + 1
        calc_prob = hcubature(r -> dispersalfn(r),
                              [0, 0, k * relSquareSize, m * relSquareSize],
                              [
                                  relSquareSize,
                                  relSquareSize,
                                  (k + 1) * relSquareSize,
                                  (m + 1) * relSquareSize
                              ],
                              maxevals = 10000)[1] / relSquareSize^2
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
    lookup_tab[!, :Prob] = lookup_tab[!, :Prob] / sum(lookup_tab[!, :Prob])
    return lookup_tab
end

function _lookup(relSquareSize::Float64,
                 maxGridSize::Int64,
                 pThresh::Float64,
                 b::Float64,
                 dispersalfn::F) where {F <: Function}
    # Create empty array
    lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

    # Loop through directions until probability is below threshold
    k = 0
    m = 0
    count = 0
    while (k <= maxGridSize && m <= maxGridSize)
        count = count + 1
        calc_prob = hcubature(r -> dispersalfn(r, b),
                              [0, 0, k * relSquareSize, m * relSquareSize],
                              [
                                  relSquareSize,
                                  relSquareSize,
                                  (k + 1) * relSquareSize,
                                  (m + 1) * relSquareSize
                              ],
                              maxevals = 10000)[1] / relSquareSize^2
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
    lookup_tab[!, :Prob] = lookup_tab[!, :Prob] / sum(lookup_tab[!, :Prob])
    return lookup_tab
end
