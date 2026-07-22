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

# `AbstractRegime`, the `HabitatUpdate` dynamics, and the `*Regime`/`RegimeCollection` condition-layer
# types now live in `Layer.jl` (as `AbstractLayer{Habitat}` + `ContinuousLayer`/…). The methods below
# dispatch on those; the old `*Hab`/`HabitatCollection` names are deprecated aliases (deprecations.jl).

function Diversity.countsubcommunities(regime::AbstractRegime)
    return _countsubcommunities(regime)
end
xmin(regime::AbstractRegime) = 0
ymin(regime::AbstractRegime) = 0
xcellsize(regime::AbstractRegime) = Float64(regime.size / km)
ycellsize(regime::AbstractRegime) = Float64(regime.size / km)
xcells(regime::AbstractRegime) = size(regime.matrix, 1)
ycells(regime::AbstractRegime) = size(regime.matrix, 2)
function indices(regime::AbstractRegime)
    return hcat(collect.(convert_coords.(eachindex(regime.matrix),
                                         xcells(regime)))...)'
end
indices(regime::AbstractRegime, idx) = indices(regime)[:, idx]
coordinates(regime::AbstractRegime) = indices(regime)

iscontinuous(::ContinuousRegime) = true

function Base.eltype(regime::ContinuousRegime{C}) where {C}
    return C
end
@recipe function f(H::ContinuousRegime{C}) where {C}
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
@recipe function f(H::ContinuousRegime{C}) where {C <: Unitful.Temperature}
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

@recipe function f(H::ContinuousTimeRegime{C, M},
                   time::Int64) where {C, M <: AbstractArray{C, 3}}
    h = ustrip.(H.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(C)]
    clims --> (minimum(h[:, :, time]) * 0.99, maximum(h[:, :, time]) * 1.01)
    return xrange(H), yrange(H), h[:, :, time]
end

function iscontinuous(regime::ContinuousTimeRegime{C, M}) where
    {C, M <: AbstractArray{C, 3}}
    return true
end
function Base.eltype(regime::ContinuousTimeRegime{C, M}) where
    {C, M <: AbstractArray{C, 3}}
    return C
end
function _resettime!(regime::ContinuousTimeRegime)
    return regime.time = 1
end

function Diversity.API._countsubcommunities(regime::ContinuousRegime)
    return length(regime.matrix)
end

function Diversity.API._countsubcommunities(regime::ContinuousTimeRegime)
    return length(regime.matrix[:, :, 1])
end

@recipe function f(H::DiscreteRegime{D}) where {D}
    h = ustrip.(H.matrix)
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> unitdict[unit(D)]
    clims --> (minimum(h) * 0.99, maximum(h) * 1.01)
    return xrange(H), yrange(H), h
end

iscontinuous(regime::DiscreteRegime) = false
function Base.eltype(regime::DiscreteRegime{D}) where {D}
    return D
end
function Diversity.API._countsubcommunities(regime::DiscreteRegime)
    return length(regime.matrix)
end

function _getdimension(regime::Union{DiscreteRegime, ContinuousRegime,
                                     ContinuousTimeRegime})
    return (size(regime.matrix, 1), size(regime.matrix, 2))
end
function _getsize(regime::Union{DiscreteRegime, ContinuousRegime,
                                ContinuousTimeRegime})
    x = regime.size * size(regime.matrix, 1)
    y = regime.size * size(regime.matrix, 2)
    return x * y
end
import Base.size
function Base.size(regime::Union{DiscreteRegime, ContinuousRegime,
                                 ContinuousTimeRegime}, d)
    return size(regime.matrix, d)
end

function _getgridsize(regime::Union{DiscreteRegime, ContinuousRegime,
                                    ContinuousTimeRegime})
    return regime.size
end

function iscontinuous(regime::RegimeCollection2{H1, H2}) where {H1, H2}
    return [iscontinuous(regime.one), iscontinuous(regime.two)]
end
function Base.eltype(regime::RegimeCollection2)
    return [eltype(regime.one), eltype(regime.two)]
end
@recipe function f(H::RegimeCollection2{H1, H2}) where {H1, H2}
    x, y = H.one, H.two
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

function _resettime!(regime::RegimeCollection2)
    _resettime!(regime.one)
    return _resettime!(regime.two)
end

function iscontinuous(regime::RegimeCollection3)
    return [
        iscontinuous(regime.one),
        iscontinuous(regime.two),
        iscontinuous(regime.three)
    ]
end
function Base.eltype(regime::RegimeCollection3)
    return [eltype(regime.one), eltype(regime.two), eltype(regime.three)]
end

@recipe function f(H::RegimeCollection3{H1, H2, H3}) where {H1, H2, H3}
    x, y, z = H.one, H.two, H.three
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

function _resettime!(regime::RegimeCollection3)
    _resettime!(regime.one)
    _resettime!(regime.two)
    return _resettime!(regime.three)
end

function _getdimension(regime::Union{RegimeCollection2, RegimeCollection3})
    return (size(regime.one.matrix, 1), size(regime.two.matrix, 2))
end
function _getsize(regime::Union{RegimeCollection2, RegimeCollection3})
    return _getsize(regime.one)
end
function Base.size(regime::Union{RegimeCollection2, RegimeCollection3}, d)
    return size(regime.one, d)
end

function _getgridsize(regime::Union{RegimeCollection2, RegimeCollection3})
    return _getgridsize(regime.one)
end

"""
    getregime(regime::H, pos::Int64) where H <: AbstractRegime

Return the regime value at 1D grid position `pos`, converting to 2D coordinates
internally.
"""
function getregime(regime::H, pos::Int64) where {H <: AbstractRegime}
    x, y = convert_coords(pos, size(regime.matrix, 1))
    return regime.matrix[x, y]
end
@doc (@doc getregime) getregime(::H, ::Symbol) where {H <: AbstractRegime}
function getregime(regime::H, field::Symbol) where {H <: AbstractRegime}
    return getfield(regime, field)
end
"""
    getregime(regime::ContinuousTimeRegime, pos::Int64)

Return the regime value at position `pos` for the current time slice of a
[`ContinuousTimeRegime`](@ref).
"""
function getregime(regime::ContinuousTimeRegime, pos::Int64)
    x, y = convert_coords(pos, size(regime.matrix, 1))
    return regime.matrix[x, y, regime.time]
end

function Diversity.API._countsubcommunities(regime::RegimeCollection2)
    return _countsubcommunities(regime.one)
end
function Diversity.API._countsubcommunities(regime::RegimeCollection3)
    return _countsubcommunities(regime.one)
end

# Function to create a regime from a discrete set of types according to the
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

Create a [`DiscreteRegime`](@ref) regime of dimension `dimension`, made up of
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
        error("There must be a weight for each type")
    sum(weights) == 1 || error("Weights of regimes must sum to 1")
    # Create weighting from proportion regimes
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
    regimeupdate = HabitatUpdate(NoChange, 0.0 / s, Unitful.Dimensions{()})
    return DiscreteRegime(T, gridsquaresize, regimeupdate)
end

"""
    simpleregime(val::Unitful.Quantity, size::Unitful.Length,
    dim::Tuple{Int64, Int64}, axis::Type{<:NicheAxis} = Unclassified)

Create a [`ContinuousRegime`](@ref) regime of dimension `dim`, with cell `size`
and filled value, `val`, on niche axis `axis`.
"""
function simpleregime(val::Unitful.Quantity, size::Unitful.Length,
                      dim::Tuple{Int64, Int64},
                      axis::Type{<:NicheAxis} = Unclassified)
    M = fill(val, dim)
    # A static regime (rate 0); its change function comes from the layer's niche `axis` via
    # `dynamics(axis)`.
    rate = 0.0 * unit(val) / s
    regimeupdate = HabitatUpdate(dynamics(axis()), rate,
                                 typeof(dimension(val)))
    return ContinuousRegime(M, size, regimeupdate)
end

"""
    simpleregime(val::Float64, size::Unitful.Length, dim::Tuple{Int64, Int64},
    axis::Type{<:NicheAxis} = Unclassified)

Create a dimensionless [`ContinuousRegime`](@ref) filled with `val`. Its update
rule comes from `axis` via `dynamics(axis)`.
"""
function simpleregime(val::Float64, size::Unitful.Length,
                      dim::Tuple{Int64, Int64},
                      axis::Type{<:NicheAxis} = Unclassified)
    M = fill(val, dim)
    regimeupdate = HabitatUpdate(dynamics(axis()), 0.0 / s,
                                 Unitful.Dimensions{()})
    return ContinuousRegime(M, size, regimeupdate)
end

"""
    tempgrad(minT::Unitful.Temperature{Float64}, maxT::Unitful.Temperature{Float64},
      size::Unitful.Length{Float64},
      dim::Tuple{Int64, Int64}, rate::Quantity{Float64, 𝚯*𝐓^-1})

Create a [`ContinuousRegime`](@ref) regime with a temperature gradient.
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
    regimeupdate = HabitatUpdate(TempChange, rate, typeof(dimension(minT)))
    return ContinuousRegime(M, size, regimeupdate)
end

"""
    raingrad(minR::Unitful.Length{Float64}, maxR::Unitful.Length{Float64},
      size::Unitful.Length{Float64},
      dim::Tuple{Int64, Int64}, rate::Quantity{Float64, 𝐋*𝐓^-1})

Create a [`ContinuousRegime`](@ref) regime with a rainfall gradient.
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
    regimeupdate = HabitatUpdate(RainfallChange, rate, typeof(dimension(minR)))
    return ContinuousRegime(M, size, regimeupdate)
end
