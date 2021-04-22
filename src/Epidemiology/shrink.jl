
"""
    _shrink_to_active(M::AbstractMatrix, active::AbstractMatrix{<:Bool})

Shrink the matrix `M` to the minimum rectangular region which contains all active cells, as
defined by `active`. Returns the shrunk matrix.

If `active is not provided, automatically determines the active region by masking out
entries which are `NaN` or `missing`.
"""
function shrink_to_active(M::AM, active::A) where {AM <: AbstractMatrix, A <: AbstractMatrix{<: Bool}}
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

function _shrink_to_active(A::AbstractArray, active::AbstractMatrix{<:Bool})
    if size(A[:, :, 1]) != size(active)
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
    return _construct_shrunk_matrix(A, shrunk_rows, shrunk_cols)
end

function shrink_to_active(M::AM) where {AM <: AbstractMatrix}
    active = .!_inactive.(M)
    return shrink_to_active(M, active)
end

function shrink_to_active(A::AA) where {AA <: AbstractArray}
    M = dropdims(sum(Float64.(A), dims=3), dims=3)
    M[M .≈ 0.0] .= NaN
    active = .!_inactive.(M)
    return _shrink_to_active(A, active)
end

_inactive(x) = ismissing(x) || isnan(x)

function findactive(A::AbstractArray{T, 3}) where T
    M = dropdims(sum(Float64.(A), dims=3), dims=3)
    M[M .== 0] .= NaN
    return .!_inactive.(M)
end

function findactive(M::AbstractArray{T, 2}) where T
    return .!_inactive.(M)
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
function _construct_shrunk_matrix(M::AbstractMatrix, row_idxs, col_idxs)::AxisArray
     return AxisArray(
         M[row_idxs, col_idxs];
         row_idxs=row_idxs,
         col_idxs=col_idxs,
     )
 end

function _construct_shrunk_matrix(M::AxisArray{T, 3}, row_idxs, col_idxs)::AxisArray{T, 3} where T
    return M[row_idxs, col_idxs, :]
end

function _construct_shrunk_matrix(M::AxisArray, row_idxs, col_idxs)::AxisArray
    return M[row_idxs, col_idxs]
 end

"""
    function convert_population(
        initial_population,
        intnum::U = Int64(1)
    )

Convert population matrix to Int matrix by filling in the inactive area with 0 population
and rounding the active area.
"""
function convert_population(
    initial_population::Matrix,
    intnum::U = Int64(1)
) where U <: Integer
    # Don't modify the arg
    initial_population = copy(initial_population)
    active = .!_inactive.(initial_population)
    initial_population[.!active] .= 0
    initial_population = U.(round.(initial_population))
    return initial_population
end

function convert_population(
    initial_population::AxisArray,
    intnum::U = Int64(1)
) where U <: Integer
    # Don't modify the arg
    initial_population = copy(initial_population)
    active = .!_inactive.(initial_population)
    # NOTE: this is a workaround as logical indexing directly on AxisArray leads to
    #   stackoverflow. see issue: https://github.com/JuliaArrays/AxisArrays.jl/issues/179
    initial_population.data[.!active] .= 0
    return AxisArray(
        U.(round.(initial_population.data)),
        initial_population.axes
    )
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
