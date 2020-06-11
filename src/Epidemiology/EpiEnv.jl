using Diversity
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref

using Diversity.API

"""
    AbstractEpiEnv{H <: AbstractHabitat, C <: AbstractControl} <:
    AbstractPartition{H}

Abstract supertype for all epi environment types and a subtype of AbstractPartition.
"""
abstract type AbstractEpiEnv{H <: AbstractHabitat, C <: AbstractControl} <:
   AbstractPartition{H} end

"""
    GridEpiEnv{H, C} <: AbstractAbiotic{H, C}

This epi environment type holds a habitat and control strategy, as well as a string of
subcommunity names, and initial susceptible population.
"""
mutable struct GridEpiEnv{H, C, A, P} <: AbstractEpiEnv{H, C}
    habitat::H
    active::A
    control::C
    initial_population::P
    names::Vector{String}
    function (::Type{GridEpiEnv{H, C}})(
        habitat::H,
        active::A,
        control::C,
        initial_population::P=zeros(Int, _getdimension(habitat)),
        names::Vector{String}=map(x -> "$x", 1:countsubcommunities(habitat))
    ) where {H, C, A <: AbstractMatrix{Bool}, P <: AbstractMatrix{Int}}
        countsubcommunities(habitat) == length(names) ||
            error("Number of subcommunities must match subcommunity names")
        if size(habitat.matrix) != size(active)
            throw(DimensionMismatch(
                "size(habitat.matrix)=$(size(habitat.matrix)) != " *
                "size(active)=$(size(active))"
            ))
        end
        if size(initial_population) != size(active)
            throw(DimensionMismatch(
                "size(initial_population)=$(size(initial_population)) != " *
                "size(active)=$(size(active))"
            ))
        end
        return new{H, C, A, P}(habitat, active, control, initial_population, names)
    end
end

import Diversity.API: _countsubcommunities
function _countsubcommunities(epienv::GridEpiEnv)
  return countsubcommunities(epienv.habitat)
end
import Diversity.API: _getsubcommunitynames
function _getsubcommunitynames(epienv::GridEpiEnv)
    return epienv.names
end

"""
    _shrink_to_active(M::AbstractMatrix, active::AbstractMatrix{<:Bool})

Shrink the matrix `M` to the minimum rectangular region which contains all active cells, as
defined by `active`. Returns the shrunk matrix.
"""
function _shrink_to_active(M::AM, active::A) where {AM <: AbstractMatrix, A <: AbstractMatrix{<: Bool}}
    if size(M) != size(active)
        throw(DimensionMismatch("size(M)=$(size(M)) != size(active)=$(size(active))"))
    end
    # Find indices of non-missing values
    idx = Tuple.(findall(active))
    # Separate into row and column indices
    row_idx = first.(idx)
    col_idx = last.(idx)
    # Return the shrunk region
    shrunk_rows = minimum(row_idx):maximum(row_idx)
    shrunk_cols = minimum(col_idx):maximum(col_idx)
    #return M[shrunk_rows, shrunk_cols]
    return _construct_shrunk_matrix(M, shrunk_rows, shrunk_cols)
end

"""
    _construct_shrunk_matrix

Construct a shrunk matrix by selecting certain rows and columns specified by `row_idxs` and
`col_idxs` from AbstractMatrix `M`.

Return an AxisArray{T, 2}. The axes will be the selected subset of the original axes if `M`
is an AxisArray. If `M` is a normal matrix, the axes of the returned AxisArray are the
selected coordinates.
"""
function _construct_shrunk_matrix(M::Matrix, row_idxs, col_idxs)::AxisArray
    return AxisArray(
        M[row_idxs, col_idxs];
        row_idxs=row_idxs,
        col_idxs=col_idxs,
    )
end

function _construct_shrunk_matrix(M::AxisArray, row_idxs, col_idxs)::AxisArray
    return M[row_idxs, col_idxs]
end

"""
    function simplehabitatAE(
        val::Union{Float64, Unitful.Quantity{Float64}},
        dimension::Tuple{Int64, Int64},
        area::Unitful.Area{Float64},
        active::AbstractMatrix{Bool},
        control::C,
        initial_population::AbstractMatrix{<:Integer}=zeros(Int, dimension),
    )

Function to create a simple `ContinuousHab` type epi environment. It creates a
`ContinuousHab` filled with a given value `val`, of dimensions `dimension` and specified
area `area`. If a Bool matrix `active` of active grid squares is included, this is used,
else one is created with all grid cells active.

!!! note
    The simulation grid will be shrunk so that it tightly wraps the active values
"""
function simplehabitatAE(
    val::Union{Float64, Unitful.Quantity{Float64}},
    dimension::Tuple{Int64, Int64},
    area::Unitful.Area{Float64},
    active::M,
    control::C,
    initial_population::P=zeros(Int, dimension),
) where {C <: AbstractControl, M <: AbstractMatrix{Bool}, P <: AbstractMatrix{<:Integer}}
    if typeof(val) <: Unitful.Temperature
        val = uconvert(K, val)
    end
    area = uconvert(km^2, area)
    gridsquaresize = sqrt(area / (dimension[1] * dimension[2]))

    # Shrink to active region
    # This doesn't change the gridsquaresize
    initial_population = _shrink_to_active(initial_population, active)
    active = _shrink_to_active(active, active)
    dimension = size(active)

    hab = simplehabitat(val, gridsquaresize, dimension)
    return GridEpiEnv{typeof(hab), typeof(control)}(hab, active, control, initial_population)
end

function simplehabitatAE(
    val::Union{Float64, Unitful.Quantity{Float64}},
    dimension::Tuple{Int64, Int64},
    area::Unitful.Area{Float64},
    control::C,
) where C <: AbstractControl
    active = fill(true, dimension)
    return simplehabitatAE(val, dimension, area, active, control)
end

"""
    simplehabitatAE(
        val::Union{Float64, Unitful.Quantity{Float64}},
        area::Unitful.Area{Float64},
        control::C,
        initial_population::AbstractMatrix{<:Real},
    )

Create a simple `ContinuousHab` type epi environment from a specified `initial_population`
matrix.

## Inputs
- `val`: Fill the habitat with this value
- `initial_population`: Used to derive the dimensions of the habitat, and the initial
    susceptible population. Values in `initial_population` which are `NaN` or `Missing` are
    used to mask off inactive areas. `initial_population` will be rounded to integers.
- `area`: The area of the habitat
- `control`: The control to apply

!!! note
    The simulation grid will be shrunk so that it tightly wraps the active values in
    `initial_population`.
"""
function simplehabitatAE(
    val::Union{Float64, Unitful.Quantity{Float64}},
    area::Unitful.Area{Float64},
    control::C,
    initial_population::M,
) where {C <: AbstractControl, M <: AbstractMatrix{<:Real}}
    inactive(x) = isnan(x) || ismissing(x)
    if all(inactive.(initial_population))
        throw(ArgumentError("initial_population is all NaN / missing"))
    end
    dimension = size(initial_population)
    active = Matrix{Bool}(.!inactive.(initial_population))
    initial_population = _convert_population(initial_population, active)
    return simplehabitatAE(val, dimension, area, active, control, initial_population)
end

"""
    _convert_population

Convert populatioin matrix to Int matrix by filling in the inactive area with 0 population
and rounding the active area.
"""
function _convert_population(
    initial_population::Matrix{<:Real},
    active::AbstractMatrix{Bool}
)::Matrix{<:Int}
    initial_population[.!active] .= 0
    initial_population = Int.(round.(initial_population))
    return initial_population
end

function _convert_population(
    initial_population::AxisArray{<:Real, 2},
    active::AbstractMatrix{Bool}
)::AxisArray{<:Int, 2}
    # NOTE: this is a workaround as logical indexing directly on AxisArray leads to
    #   stackoverflow. see issue: https://github.com/JuliaArrays/AxisArrays.jl/issues/179
    initial_population.data[.!active] .= 0
    return AxisArray(
            Int.(round.(initial_population.data)),
            initial_population.axes
        )
end
