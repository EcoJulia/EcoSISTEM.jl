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

  function MPIGridLandscape(speciescounts::Vector{Int32}, sccounts::Vector{Int32}, rows_matrix::Matrix{Int64}, cols_vector::Vector{Int64})
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    totalsc = sum(sccounts)
    totalspecies = sum(speciescounts)

    lastspecies = sum(speciescounts[1:(rank + 1)])
    firstspecies = lastspecies - speciescounts[rank + 1] + 1

    lastsc = sum(sccounts[1:(rank + 1)])
    firstsc = lastsc - sccounts[rank + 1] + 1

    spindices = [0; cumsum(speciescounts) .* sccounts[rank + 1]]
    scindices = [0; cumsum(sccounts) .* speciescounts[rank + 1]]

    reshaped_cols = map(1:length(sccounts)) do i
      reshape(view(cols_vector, (spindices[i] + 1) : spindices[i + 1]), Int64(speciescounts[i]), Int64(sccounts[rank + 1]))
    end
    rows = (total = totalsc, first = firstsc, last = lastsc, counts = sccounts .* speciescounts[rank + 1])
    cols = (total = totalspecies, first = firstspecies, last = lastspecies, counts = speciescounts .* sccounts[rank + 1])
    return new{typeof(reshaped_cols[1]), typeof(rows)}(rows_matrix,
    cols_vector, reshaped_cols, rows, cols,
    [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
end


"""
    emptyMPIgridlandscape(speciescounts::Vector{Int32},
    sccounts::Vector{Int32})

Function to create an empty MPIGridLandscape given information about the MPI setup.
"""
function emptyMPIgridlandscape(speciescounts::Vector{Int32},
  sccounts::Vector{Int32})
  rank = MPI.Comm_rank(MPI.COMM_WORLD)

  totalsc = sum(sccounts)
  totalspecies = sum(speciescounts)
  rows_matrix = zeros(Int64, speciescounts[rank + 1], totalsc)
  cols_vector = zeros(Int64, totalspecies * sccounts[rank + 1])

  return MPIGridLandscape(speciescounts, sccounts,
  rows_matrix, cols_vector)
end

function synchronise_from_rows!(ml::MPIGridLandscape)
  MPI.Alltoallv!(ml.rows_matrix, ml.cols_vector,
   ml.rows_tuple.counts, ml.cols_tuple.counts, MPI.COMM_WORLD)
end

function synchronise_from_cols!(ml::MPIGridLandscape)
  MPI.Alltoallv!(ml.cols_vector, ml.rows_matrix,
  ml.cols_tuple.counts, ml.rows_tuple.counts, MPI.COMM_WORLD)
end
