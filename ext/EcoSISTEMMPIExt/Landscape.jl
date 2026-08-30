# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The distributed abundance matrix: `MPIGridLandscape` and the two collective operations that keep
# its two decompositions agreeing.
#
# A rank owns a block of **species rows** and a block of **grid columns**, and the hot loop needs
# both views: demographics run per species, dispersal per cell. `synchronise_from_rows!` and
# `synchronise_from_cols!` are how one view is rebuilt from the other, and they are collective — every
# rank must reach them, in the same order, or the run deadlocks.

import EcoSISTEM
using MPI
using Random
using BlockArrays: mortar
using DimensionalData: DimArray, Dim

"""
    MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple}

Ecosystem abundances distributed across MPI nodes. `rows_matrix` is the local
slice of the abundance matrix (species owned by this node × all grid cells).
`cols_vector` is the flattened column-partitioned view (all species × grid cells
owned by this node). `reshaped_cols` holds reshaped views per MPI block.
`rows_tuple` and `cols_tuple` are named tuples with `total`, `first`, `last`,
and `counts` fields describing the partitioning. `dimrows` and `dimcols` are
labelled views of the same memory, carrying the global species and cell indices
this rank holds. Random draws during simulation use Julia's task-local default
RNG, so no generator state is stored here.
"""
struct MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple, DR, DC} <:
       EcoSISTEM.MPIGridLandscape
    rows_matrix::Matrix{Int64}
    cols_vector::Vector{Int64}
    reshaped_cols::Vector{RA}
    rows_tuple::NT
    cols_tuple::NT
    dimrows::DR
    dimcols::DC

    function MPIGridLandscape(sppcounts::Vector{Int32},
                              sccounts::Vector{Int32},
                              rows_matrix::Matrix{Int64},
                              cols_vector::Vector{Int64})
        rank = MPI.Comm_rank(MPI.COMM_WORLD)

        totalspp = sum(sppcounts)
        totalsc = sum(sccounts)

        lastsp = sum(sppcounts[1:(rank + 1)])
        firstsp = lastsp - sppcounts[rank + 1] + 1

        lastsc = sum(sccounts[1:(rank + 1)])
        firstsc = lastsc - sccounts[rank + 1] + 1

        sppindices = [0; cumsum(sppcounts) .* sccounts[rank + 1]]
        scindices = [0; cumsum(sccounts) .* sppcounts[rank + 1]]

        reshaped_cols = map(eachindex(sccounts)) do i
            return reshape(view(cols_vector,
                                (sppindices[i] + 1):sppindices[i + 1]),
                           Int64(sppcounts[i]),
                           Int64(sccounts[rank + 1]))
        end
        rows = (total = totalspp,
                first = firstsp,
                last = lastsp,
                counts = sccounts .* sppcounts[rank + 1])
        cols = (total = totalsc,
                first = firstsc,
                last = lastsc,
                counts = sppcounts .* sccounts[rank + 1])
        (rows.last - rows.first + 1) * cols.total == length(rows_matrix) ||
            error("rows_matrix size mismatch: $(rows.last - rows.first + 1) * $(cols.total) !=$(length(rows_matrix))")
        (cols.last - cols.first + 1) * rows.total == length(cols_vector) ||
            error("cols_vector size mismatch: $(cols.last - cols.first + 1) * $(rows.total) !=$(length(cols_vector))")

        # The labelled views. Both are views rather than copies, and both say which slice of the
        # whole run this rank actually holds -- which matters more here than serially, where the
        # matrix simply is everything. The labels are global indices: the species names are not
        # available at this point, and adding them would change a released public signature.
        dimrows = DimArray(rows_matrix,
                           (Dim{:species}(firstsp:lastsp),
                            Dim{:location}(1:totalsc)))
        # `cols_vector` arrives block-stacked -- one block per contributing rank -- because that is
        # what `Alltoallv!` produces. `mortar` presents those blocks as one ordinary
        # `(all species x this rank's cells)` matrix without copying, so a caller need not know.
        # It is for inspection only: a scalar index has to find its block first, which measured
        # about eight times the cost of walking `reshaped_cols` directly, so the hot loop keeps
        # doing that.
        dimcols = DimArray(mortar(reshape(reshaped_cols, length(reshaped_cols),
                                          1)),
                           (Dim{:species}(1:totalspp),
                            Dim{:location}(firstsc:lastsc)))
        return new{typeof(reshaped_cols[1]), typeof(rows), typeof(dimrows),
                   typeof(dimcols)}(rows_matrix,
                                    cols_vector,
                                    reshaped_cols,
                                    rows,
                                    cols,
                                    dimrows,
                                    dimcols)
    end
end

EcoSISTEM.MPIGridLandscape(args...) = MPIGridLandscape(args...)

"""
    empty_mpi_gridlandscape(sppcounts::Vector{Int32}, sccounts::Vector{Int32})

Create an empty MPIGridLandscape given information about the MPI setup.
"""
function EcoSISTEM.empty_mpi_gridlandscape(sppcounts::Vector{Int32},
                                           sccounts::Vector{Int32})
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    rows_matrix = zeros(Int64, sppcounts[rank + 1], sum(sccounts))
    cols_vector = zeros(Int64, sum(sppcounts) * sccounts[rank + 1])

    return MPIGridLandscape(sppcounts, sccounts, rows_matrix, cols_vector)
end

"""
    synchronise_from_rows!(ml::MPIGridLandscape)

Synchronise abundance data from the row-partitioned `rows_matrix` to the
column-partitioned `cols_vector` across all MPI nodes using `Alltoallv!`.
"""
function EcoSISTEM.synchronise_from_rows!(ml::MPIGridLandscape)
    return MPI.Alltoallv!(MPI.VBuffer(ml.rows_matrix, ml.rows_tuple.counts),
                          MPI.VBuffer(ml.cols_vector, ml.cols_tuple.counts),
                          MPI.COMM_WORLD)
end

"""
    synchronise_from_cols!(ml::MPIGridLandscape)

Synchronise abundance data from the column-partitioned `cols_vector` back to the
row-partitioned `rows_matrix` across all MPI nodes using `Alltoallv!`.
"""
function EcoSISTEM.synchronise_from_cols!(ml::MPIGridLandscape)
    return MPI.Alltoallv!(MPI.VBuffer(ml.cols_vector, ml.cols_tuple.counts),
                          MPI.VBuffer(ml.rows_matrix, ml.rows_tuple.counts),
                          MPI.COMM_WORLD)
end
