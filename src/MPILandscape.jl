using MPI

"""
    MPIGridLandscape

MPIEcosystem abundances housed in the landscape, shared across multiple nodes.

"""
mutable struct MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple}
  rows_matrix::Matrix{Int64}
  cols_vector::Vector{Int64}
  reshaped_cols::Vector{RA}
  rows_tuple::NT
  cols_tuple::NT
  seed::Vector{MersenneTwister}

  function MPIGridLandscape(sppcounts::Vector{Int32}, sccounts::Vector{Int32}, rows_matrix::Matrix{Int64}, cols_vector::Vector{Int64})
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    totalspp = sum(sppcounts)
    totalsc = sum(sccounts)

    lastsp = sum(sppcounts[1:(rank + 1)])
    firstsp = lastsp - sppcounts[rank + 1] + 1

    lastsc = sum(sccounts[1:(rank + 1)])
    firstsc = lastsc - sccounts[rank + 1] + 1

    sppindices = [0; cumsum(sppcounts) .* sccounts[rank + 1]]
    scindices = [0; cumsum(sccounts) .* sppcounts[rank + 1]]

    reshaped_cols = map(1:length(sccounts)) do i
      reshape(view(cols_vector, (sppindices[i] + 1) : sppindices[i + 1]),
              Int64(sppcounts[i]), Int64(sccounts[rank + 1]))
    end
    rows = (total = totalspp, first = firstsp, last = lastsp,
            counts = sccounts .* sppcounts[rank + 1])
    cols = (total = totalsc, first = firstsc, last = lastsc,
            counts = sppcounts .* sccounts[rank + 1])
    (rows.last - rows.first + 1) * cols.total == length(rows_matrix) ||
      error("rows_matrix size mismatch: $(rows.last - rows.first + 1) * $(cols.total) !=$(length(rows_matrix))")
    (cols.last - cols.first + 1) * rows.total == length(cols_vector) ||
      error("cols_vector size mismatch: $(cols.last - cols.first + 1) * $(rows.total) !=$(length(cols_vector))")

    return new{typeof(reshaped_cols[1]), typeof(rows)}(rows_matrix,
      cols_vector, reshaped_cols, rows, cols,
      [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
end

"""
    emptyMPIgridlandscape(sppcounts::Vector{Int32}, sccounts::Vector{Int32})

Function to create an empty MPIGridLandscape given information about the MPI setup.
"""
function emptyMPIgridlandscape(sppcounts::Vector{Int32},
                               sccounts::Vector{Int32})
  rank = MPI.Comm_rank(MPI.COMM_WORLD)

  rows_matrix = zeros(Int64, sppcounts[rank + 1], sum(sccounts))
  cols_vector = zeros(Int64, sum(sppcounts) * sccounts[rank + 1])

  return MPIGridLandscape(sppcounts, sccounts, rows_matrix, cols_vector)
end

function synchronise_from_rows!(ml::MPIGridLandscape)
  MPI.Alltoallv!(MPI.VBuffer(ml.rows_matrix, ml.rows_tuple.counts),
                 MPI.VBuffer(ml.cols_vector, ml.cols_tuple.counts),
                 MPI.COMM_WORLD)
end

function synchronise_from_cols!(ml::MPIGridLandscape)
  MPI.Alltoallv!(MPI.VBuffer(ml.cols_vector, ml.cols_tuple.counts),
                 MPI.VBuffer(ml.rows_matrix, ml.rows_tuple.counts),
                 MPI.COMM_WORLD)
end
