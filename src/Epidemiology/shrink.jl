
"""
    _shrink_to_active(M::AbstractMatrix, active::AbstractMatrix{<:Bool})

Shrink the matrix `M` to the minimum rectangular region which contains all active cells, as
defined by `active`. Returns the shrunk matrix.
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

function shrink_to_boundingbox(M::AM) where {AM <: AbstractMatrix}
    inactive(x) = isnan(x) || ismissing(x) || x==0
    return shrink_to_active(M, .!inactive.(M))
end
