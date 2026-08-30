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
using DimensionalData: DimArray, Dim, X, Y, dims
using DimensionalData.Lookups: Categorical, NoLookup

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
struct MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple, DR, DG,
                        DC} <:
       EcoSISTEM.MPIGridLandscape
    rows_matrix::Matrix{Int64}
    cols_vector::Vector{Int64}
    reshaped_cols::Vector{RA}
    rows_tuple::NT
    cols_tuple::NT
    dimrows::DR
    dimgrid::DG
    dimcols::DC

    function MPIGridLandscape(sppcounts::Vector{Int32},
                              sccounts::Vector{Int32},
                              rows_matrix::Matrix{Int64},
                              cols_vector::Vector{Int64},
                              names::Vector{String},
                              grid::EcoSISTEM.StudyGrid)
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
        # This rank's species, by name, over every cell -- the direct analogue of the serial
        # landscape's `dimmatrix`, and its `dimgrid`, which is why there are two of them: `rows`
        # covers all the cells, so the real Y and X coordinates are available and meaningful.
        # Cell descriptors, computed on demand from the grid: the lookup holds a grid reference
        # and nothing else, so labelling costs the same few bytes at any grid size. `dimcols` takes
        # a view of the same vector, which is how a rank's own cells get their real names.
        cells = EcoSISTEM.CellNames(grid)
        yx = (grid.y, grid.x)
        myspecies = Dim{:species}(Categorical(names[firstsp:lastsp]))
        dimrows = DimArray(rows_matrix, (myspecies, Dim{:location}(cells)))
        dimgrid = DimArray(reshape(rows_matrix,
                                   (lastsp - firstsp + 1, length.(yx)...)),
                           (myspecies, yx...))
        # `cols_vector` arrives block-stacked -- one block per contributing rank -- because that is
        # what `Alltoallv!` produces. `mortar` presents those blocks as one ordinary
        # `(all species x this rank's cells)` matrix without copying, so a caller need not know.
        # It is for inspection only: a scalar index has to find its block first, which measured
        # about eight times the cost of walking `reshaped_cols` directly, so the hot loop keeps
        # doing that.
        # Every species, by name, over this rank's cells -- which carry their real descriptors,
        # the same ones the serial landscape uses. They are a contiguous run of the flat ordering
        # rather than a rectangle, so there is no Y/X view of them to give.
        dimcols = DimArray(mortar(reshape(reshaped_cols, length(reshaped_cols),
                                          1)),
                           (Dim{:species}(Categorical(names)),
                            Dim{:location}(view(cells, firstsc:lastsc))))
        return new{typeof(reshaped_cols[1]), typeof(rows), typeof(dimrows),
                   typeof(dimgrid), typeof(dimcols)}(rows_matrix,
                                                     cols_vector,
                                                     reshaped_cols,
                                                     rows,
                                                     cols,
                                                     dimrows,
                                                     dimgrid,
                                                     dimcols)
    end
end

EcoSISTEM.MPIGridLandscape(args...) = MPIGridLandscape(args...)

# The distributed method of `EcoSISTEM.empty_landscape`: same habitat and species list as the
# serial one, plus the partition. The names and grid dimensions are derived here exactly as the
# serial method derives them, rather than being passed in already computed -- which is what stops
# the two sides drifting apart, and what keeps this signature stable when the derivation changes.
function EcoSISTEM.empty_landscape(habitat::EcoSISTEM.GridHabitat,
                                   spplist::EcoSISTEM.SpeciesList,
                                   sppcounts::Vector{Int32},
                                   sccounts::Vector{Int32})
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    rows_matrix = zeros(Int64, sppcounts[rank + 1], sum(sccounts))
    cols_vector = zeros(Int64, sum(sppcounts) * sccounts[rank + 1])

    return MPIGridLandscape(sppcounts, sccounts, rows_matrix, cols_vector,
                            spplist.names, EcoSISTEM.getcoords(habitat))
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
