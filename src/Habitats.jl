# SPDX-License-Identifier: LGPL-3.0-or-later

using StatsBase
using Diversity
using Diversity.API
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using RecipesBase
using EcoBase
using Measures: AbsoluteLength

import Base.eltype
import EcoBase: xmin, ymin, xcellsize, ycellsize, xcells, ycells, cellsize,
                cells, xrange, yrange, xmax, ymax, indices, coordinates
import Diversity.API._countsubcommunities
import Diversity.countsubcommunities

const px = AbsoluteLength(0.254)

"""
    AbstractHabitat

Abstract supertype for all habitat types
"""
abstract type AbstractHabitat{H} <: EcoBase.AbstractGrid end

function Diversity.countsubcommunities(ah::AbstractHabitat)
    return _countsubcommunities(ah)
end
xmin(ah::AbstractHabitat) = 0
ymin(ah::AbstractHabitat) = 0
xcellsize(ah::AbstractHabitat) = Float64(ah.size / km)
ycellsize(ah::AbstractHabitat) = Float64(ah.size / km)
xcells(ah::AbstractHabitat) = size(ah.matrix, 1)
ycells(ah::AbstractHabitat) = size(ah.matrix, 2)
function indices(ah::AbstractHabitat)
    return hcat(collect.(convert_coords.(eachindex(ah.matrix), xcells(ah)))...)'
end
indices(ah::AbstractHabitat, idx) = indices(ah)[:, idx]
coordinates(ah::AbstractHabitat) = indices(ah)

"""
    HabitatUpdate{F <: Function, DT}

Stores the update rule for a habitat. `changefun` is the function applied to the
habitat matrix at each timestep, and `rate` is the rate of change with
appropriate units.
"""
struct HabitatUpdate{F <: Function, DT}
    changefun::F
    rate::DT
end

"""
    HabitatUpdate(changefun::F, rate::DT, ::Type{D})

Construct a `HabitatUpdate` with a type-checked rate. Errors if `rate * 1month`
does not have dimensions `D`.
"""
function HabitatUpdate(changefun::F,
                       rate::DT,
                       ::Type{D}) where {F <: Function, DT,
                                         D <: Unitful.Dimensions}
    dimension(rate * 1month) isa D ||
        error("Failed to match types $(rate * 1month) vs $D")
    return HabitatUpdate(changefun, rate)
end

"""
    ContinuousHab{C <: Number} <: AbstractHabitat{C}

This habitat subtype houses a habitat matrix `matrix` of any units, a grid
square size `size` and [`HabitatUpdate`](@ref) type `change`.
"""
mutable struct ContinuousHab{C <: Number} <: AbstractHabitat{C}
    matrix::Matrix{C}
    size::Unitful.Length
    change::HabitatUpdate
end

iscontinuous(::ContinuousHab) = true

function Base.eltype(hab::ContinuousHab{C}) where {C}
    return C
end
@recipe function f(H::ContinuousHab{C}) where {C}
    unitdict = Dict(K => "Temperature (K)",
                    °C => "Temperature (°C)",
                    mm => "Rainfall (mm)",
                    kJ => "Solar Radiation (kJ)")
    h = ustrip.(H.matrix)
    seriestype := :heatmap
    grid --> false
    right_margin --> 0.1px
    margin --> 10px
    aspect_ratio --> 1
    title --> unitdict[unit(C)]
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    return xrange(H), yrange(H), h
end
@recipe function f(H::ContinuousHab{C}) where {C <: Unitful.Temperature}
    unitdict = Dict(K => "Temperature (K)",
                    °C => "Temperature (°C)",
                    mm => "Rainfall (mm)",
                    kJ => "Solar Radiation (kJ)")
    h = ustrip.(uconvert.(°C, H.matrix))
    seriestype := :heatmap
    grid --> false
    right_margin --> 0.1px
    margin --> 10px
    aspect_ratio --> 1
    title --> unitdict[unit(C)]
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    return xrange(H), yrange(H), h
end

"""
    ContinuousTimeHab{C <: Number, M <: AbstractArray{C, 3}} <: AbstractHabitat{C}

This habitat subtype houses a habitat matrix `matrix` of any units, the time
slice of the habitat matrix currently being operated on `time`, a grid square
size `size` and [`HabitatUpdate`](@ref) type `change`.
"""
mutable struct ContinuousTimeHab{C <: Number, M <: AbstractArray{C, 3}} <:
               AbstractHabitat{C}
    matrix::M
    time::Int64
    size::Unitful.Length
    change::HabitatUpdate
end
@recipe function f(H::ContinuousTimeHab{C, M},
                   time::Int64) where {C, M <: AbstractArray{C, 3}}
    h = ustrip.(H.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(C)]
    clims --> (minimum(h[:, :, time]) * 0.99, maximum(h[:, :, time]) * 1.01)
    return xrange(H), yrange(H), h[:, :, time]
end

function iscontinuous(hab::ContinuousTimeHab{C, M}) where
    {C, M <: AbstractArray{C, 3}}
    return true
end
function Base.eltype(hab::ContinuousTimeHab{C, M}) where
    {C, M <: AbstractArray{C, 3}}
    return C
end
function _resettime!(hab::ContinuousTimeHab)
    return hab.time = 1
end

function Diversity.API._countsubcommunities(hab::ContinuousHab)
    return length(hab.matrix)
end

function Diversity.API._countsubcommunities(hab::ContinuousTimeHab)
    return length(hab.matrix[:, :, 1])
end

"""
    DiscreteHab{D} <: AbstractHabitat{D}

Habitat subtype with a discrete `matrix` of element type `D`, a grid cell
`size`, and a [`HabitatUpdate`](@ref) rule `change` applied at each timestep.
"""
mutable struct DiscreteHab{D} <: AbstractHabitat{D}
    matrix::Matrix{D}
    size::Unitful.Length
    change::HabitatUpdate
end
@recipe function f(H::DiscreteHab{D}) where {D}
    h = ustrip.(H.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(D)]
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    return xrange(H), yrange(H), h
end

iscontinuous(hab::DiscreteHab) = false
function Base.eltype(hab::DiscreteHab{D}) where {D}
    return D
end
function Diversity.API._countsubcommunities(hab::DiscreteHab)
    return length(hab.matrix)
end

function _getdimension(hab::Union{DiscreteHab, ContinuousHab,
                                  ContinuousTimeHab})
    return (size(hab.matrix, 1), size(hab.matrix, 2))
end
function _getsize(hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab})
    x = hab.size * size(hab.matrix, 1)
    y = hab.size * size(hab.matrix, 2)
    return x * y
end
import Base.size
function Base.size(hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab}, d)
    return size(hab.matrix, d)
end

function _getgridsize(hab::Union{DiscreteHab, ContinuousHab, ContinuousTimeHab})
    return hab.size
end

"""
    HabitatCollection2{H1, H2} <: AbstractHabitat{Tuple{H1, H2}}

Composite habitat pairing two sub-habitats `h1` and `h2`, allowing
multi-variable abiotic environments (e.g. temperature and rainfall together).
"""
struct HabitatCollection2{H1, H2} <: AbstractHabitat{Tuple{H1, H2}}
    h1::H1
    h2::H2
end
function iscontinuous(hab::HabitatCollection2{H1, H2}) where {H1, H2}
    return [iscontinuous(hab.h1), iscontinuous(hab.h2)]
end
function Base.eltype(hab::HabitatCollection2)
    return [eltype(hab.h1), eltype(hab.h2)]
end
@recipe function f(H::HabitatCollection2{H1, H2}) where {H1, H2}
    x, y = H.h1, H.h2
    layout := 2
    @series begin
        subplot := 1
        x
    end
    @series begin
        subplot := 2
        y
    end
end

function _resettime!(hab::HabitatCollection2)
    _resettime!(hab.h1)
    return _resettime!(hab.h2)
end

"""
    HabitatCollection3{H1, H2, H3} <: AbstractHabitat{Tuple{H1, H2, H3}}

Composite habitat combining three sub-habitats `h1`, `h2`, and `h3`.
"""
struct HabitatCollection3{H1, H2, H3} <:
       AbstractHabitat{Tuple{H1, H2, H3}}
    h1::H1
    h2::H2
    h3::H3
end
function iscontinuous(hab::HabitatCollection3)
    return [iscontinuous(hab.h1), iscontinuous(hab.h2), iscontinuous(hab.h3)]
end
function Base.eltype(hab::HabitatCollection3)
    return [eltype(hab.h1), eltype(hab.h2), eltype(hab.h3)]
end

@recipe function f(H::HabitatCollection3{H1, H2, H3}) where {H1, H2, H3}
    x, y, z = H.h1, H.h2, H.h3
    layout := 3
    @series begin
        subplot := 1
        x
    end
    @series begin
        subplot := 2
        y
    end
    @series begin
        subplot := 3
        z
    end
end

function _resettime!(hab::HabitatCollection3)
    _resettime!(hab.h1)
    _resettime!(hab.h2)
    return _resettime!(hab.h3)
end

function _getdimension(hab::Union{HabitatCollection2, HabitatCollection3})
    return (size(hab.h1.matrix, 1), size(hab.h2.matrix, 2))
end
function _getsize(hab::Union{HabitatCollection2, HabitatCollection3})
    return _getsize(hab.h1)
end
function Base.size(hab::Union{HabitatCollection2, HabitatCollection3}, d)
    return size(hab.h1, d)
end

function _getgridsize(hab::Union{HabitatCollection2, HabitatCollection3})
    return _getgridsize(hab.h1)
end

"""
    gethabitat(hab::H, pos::Int64) where H <: AbstractHabitat

Return the habitat value at 1D grid position `pos`, converting to 2D coordinates
internally.
"""
function gethabitat(hab::H, pos::Int64) where {H <: AbstractHabitat}
    x, y = convert_coords(pos, size(hab.matrix, 1))
    return hab.matrix[x, y]
end
@doc (@doc gethabitat) gethabitat(::H, ::Symbol) where {H <: AbstractHabitat}
function gethabitat(hab::H, field::Symbol) where {H <: AbstractHabitat}
    return getfield(hab, field)
end
"""
    gethabitat(hab::ContinuousTimeHab, pos::Int64)

Return the habitat value at position `pos` for the current time slice of a
[`ContinuousTimeHab`](@ref).
"""
function gethabitat(hab::ContinuousTimeHab, pos::Int64)
    x, y = convert_coords(pos, size(hab.matrix, 1))
    return hab.matrix[x, y, hab.time]
end

function Diversity.API._countsubcommunities(hab::HabitatCollection2)
    return _countsubcommunities(hab.h1)
end
function Diversity.API._countsubcommunities(hab::HabitatCollection3)
    return _countsubcommunities(hab.h1)
end

# Function to create a habitat from a discrete set of types according to the
# Saura-Martinez-Millan algorithm (2000)
function _percolate!(M::AbstractMatrix, clumpiness::Real)
    for i in eachindex(M)
        if rand(Uniform(0, 1)) < clumpiness
            M[i] = 1
        end
    end
end

# Function to create clusters from percolated grid
function _identify_clusters!(M::AbstractMatrix)
    # Begin cluster count
    count = 1
    # Loop through each grid square in M
    for x in Base.axes(M, 1)
        for y in Base.axes(M, 2)

            # If square is marked as 1, then apply cluster finding algorithm
            if M[x, y] == 1.0
                # Find neighbours of M at this location
                neighbours = get_neighbours(M, y, x)
                # Find out if any of the neighbours also have a value of 1, thus, have
                # not been assigned a cluster yet
                cluster = vcat(mapslices(x -> M[x[1], x[2]] .== 1, neighbours,
                                         dims = 2)...)
                # Find out if any of the neighbours have a value > 1, thus, have already
                # been assigned a cluster
                already = vcat(mapslices(x -> M[x[1], x[2]] .> 1, neighbours,
                                         dims = 2)...)
                # If any already assigned neighbours, then assign the grid square to this
                # same type
                if any(already)
                    neighbours = neighbours[already, :]
                    M[x, y] = M[neighbours[1, 1], neighbours[1, 2]]
                    # If none are assigned yet, then create a new cluster
                else
                    count = count + 1
                    neighbours = neighbours[cluster, :]
                    M[x, y] = count
                    map(i -> M[neighbours[i, 1], neighbours[i, 2]] = count,
                        Base.axes(neighbours, 1))
                end
            end
        end
    end
end

function _fill_in!(T, M, types, wv)
    # Loop through grid of clusters
    for x in Base.axes(M, 1)
        for y in Base.axes(M, 2)
            # If square is zero then it is yet to be assigned
            if M[x, y] == 0
                # Find neighbours of square on string grid
                neighbours = get_neighbours(T, y, x, 8)
                # Check if they have already been assigned
                already = vcat(mapslices(x -> isassigned(T, x[1], x[2]),
                                         neighbours, dims = 2)...)
                # If any already assigned then sample from most frequent neighbour traits
                if any(already)
                    neighbours = neighbours[already, :]
                    # Find all neighbour traits
                    neighbour_traits = map(i -> T[neighbours[i, 1],
                                                  neighbours[i, 2]],
                                           Base.axes(neighbours, 1))
                    # Find which one is counted most often
                    ind = argmax(map(x -> sum(neighbour_traits .== x), types))
                    # Assign this type to the grid square in T
                    T[x, y] = types[ind]
                    # If none are assigned in entire grid already,
                    # sample randomly from traits
                elseif all(M .<= 1)
                    T[x, y] = sample(types, wv)
                    # If some are assigned in grid, sample from these
                else
                    T[x, y] = sample(T[M .> 1])
                end
            end
        end
    end
end
"""
    randomniches(dimension::Tuple, types::Vector{Int64}, clumpiness::Float64,
        weights::Vector, gridsquaresize::Unitful.Length)

Create a [`DiscreteHab`](@ref) habitat of dimension `dimension`, made up of
integer niche types `types` with relative weightings `weights` and spatial
clumpiness controlled by `clumpiness`. Cell size is set by `gridsquaresize`.
"""
function randomniches(dimension::Tuple,
                      types::Vector{Int64},
                      clumpiness::Float64,
                      weights::Vector,
                      gridsquaresize::Unitful.Length)
    # Check that the proportion of coverage for each type matches the number
    # of types and that they add up to 1
    length(weights) == length(types) ||
        error("There must be an area proportion for each type")
    sum(weights) == 1 || error("Proportion of habitats must add up to 1")
    # Create weighting from proportion habitats
    wv = Weights(weights)

    # Create an empty grid of the right dimension
    M = zeros(dimension)

    # If the dimensions are too small for the algorithm, just use a weighted sample
    if dimension[1] <= 2 || dimension[2] <= 2
        T = sample(types, Weights(weights), dimension)
    else
        # Percolation step
        _percolate!(M, clumpiness)
        # Select clusters and assign types
        _identify_clusters!(M)
        # Create a string grid of the same dimensions
        T = Array{Int64}(undef, dimension)
        # Fill in T with clusters already created
        map(x -> T[M .== x] .= sample(types, wv), 1:maximum(M))
        # Fill in undefined squares with most frequent neighbour
        _fill_in!(T, M, types, wv)
    end
    habitatupdate = HabitatUpdate(NoChange, 0.0 / s, Unitful.Dimensions{()})
    return DiscreteHab(T, gridsquaresize, habitatupdate)
end

"""
    simplehabitat(val::Unitful.Quantity, size::Unitful.Length,
    dim::Tuple{Int64, Int64})

Create a [`ContinuousHab`](@ref) habitat of dimension `dim`, with cell `size`
and filled value, `val`.
"""
function simplehabitat(val::Unitful.Quantity, size::Unitful.Length,
                       dim::Tuple{Int64, Int64})
    M = fill(val, dim)
    func = ChangeLookup[unit(val)]
    rate = 0.0 * unit(val) / s
    habitatupdate = HabitatUpdate(func, rate, typeof(dimension(val)))
    return ContinuousHab(M, size, habitatupdate)
end

"""
    simplehabitat(val::Float64, size::Unitful.Length, dim::Tuple{Int64, Int64})

Create a dimensionless [`ContinuousHab`](@ref) filled with `val`. Uses
[`NoChange`](@ref) as the update rule.
"""
function simplehabitat(val::Float64, size::Unitful.Length,
                       dim::Tuple{Int64, Int64})
    M = fill(val, dim)

    habitatupdate = HabitatUpdate(NoChange, 0.0 / s, Unitful.Dimensions{()})
    return ContinuousHab(M, size, habitatupdate)
end

"""
    tempgrad(minT::Unitful.Temperature{Float64}, maxT::Unitful.Temperature{Float64},
      size::Unitful.Length{Float64},
      dim::Tuple{Int64, Int64}, rate::Quantity{Float64, 𝚯*𝐓^-1})

Create a [`ContinuousHab`](@ref) habitat with a temperature gradient.
"""
function tempgrad(minT::Unitful.Temperature{Float64},
                  maxT::Unitful.Temperature{Float64},
                  size::Unitful.Length{Float64},
                  dim::Tuple{Int64, Int64},
                  rate::Quantity{Float64, 𝚯 * 𝐓^-1})
    dim[1] > 1 ||
        error("First dimension should be greater than 1 for temperature gradient")
    M = Array{typeof(minT)}(undef, dim)
    total = dim[1]
    temp_range = collect(range(minT, stop = maxT, length = total))
    map(1:total) do seq
        return M[seq, :] .= temp_range[seq]
    end
    habitatupdate = HabitatUpdate(TempChange, rate, typeof(dimension(minT)))
    return ContinuousHab(M, size, habitatupdate)
end

"""
    raingrad(minR::Unitful.Length{Float64}, maxR::Unitful.Length{Float64},
      size::Unitful.Length{Float64},
      dim::Tuple{Int64, Int64}, rate::Quantity{Float64, 𝐋*𝐓^-1})

Create a [`ContinuousHab`](@ref) habitat with a rainfall gradient.
"""
function raingrad(minR::Unitful.Length{Float64},
                  maxR::Unitful.Length{Float64},
                  size::Unitful.Length{Float64},
                  dim::Tuple{Int64, Int64},
                  rate::Quantity{Float64, 𝐋 * 𝐓^-1})
    dim[1] > 1 ||
        error("First dimension should be greater than 1 for temperature gradient")
    M = Array{typeof(minR)}(undef, dim)
    total = dim[1]
    rain_range = collect(range(minR, stop = maxR, length = total))
    map(1:total) do seq
        return M[seq, :] .= rain_range[seq]
    end
    habitatupdate = HabitatUpdate(RainfallChange, rate, typeof(dimension(minR)))
    return ContinuousHab(M, size, habitatupdate)
end
