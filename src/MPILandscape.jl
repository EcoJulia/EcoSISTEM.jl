using MPI

"""
    MPIGridLandscape

MPIEcosystem abundances housed in the landscape, shared across multiple nodes.

"""
mutable struct MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple}
  horizontal_matrix::Matrix{Int64}
  vertical_vector::Vector{Int64}
  reshaped_vertical::Vector{RA}
  horizontal_tuple::NT
  vertical_tuple::NT
  seed::Vector{MersenneTwister}

  function MPIGridLandscape(speciescounts::Vector{Int32}, sccounts::Vector{Int32}, horizontal_matrix::Matrix{Int64}, vertical_vector::Vector{Int64})
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    totalsc = sum(sccounts)
    totalspecies = sum(speciescounts)

    lastspecies = sum(speciescounts[1:(rank + 1)])
    firstspecies = lastspecies - speciescounts[rank + 1] + 1

    lastsc = sum(sccounts[1:(rank + 1)])
    firstsc = lastsc - sccounts[rank + 1] + 1

    spindices = [0; cumsum(speciescounts) .* sccounts[rank + 1]]
    scindices = [0; cumsum(sccounts) .* speciescounts[rank + 1]]

    reshaped_vertical = map(1:length(sccounts)) do i
      reshape(view(vertical_vector, (spindices[i] + 1) : spindices[i + 1]), Int64(speciescounts[i]), Int64(sccounts[rank + 1]))
    end
    horizontal = (total = totalsc, first = firstsc, last = lastsc, counts = sccounts .* speciescounts[rank + 1])
    vertical  = (total = totalspecies, first = firstspecies, last = lastspecies, counts = speciescounts .* sccounts[rank + 1])
    return new{typeof(reshaped_vertical[1]), typeof(horizontal)}(horizontal_matrix,
    vertical_vector, reshaped_vertical, horizontal, vertical,
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
  horizontal_matrix = zeros(Int64, speciescounts[rank + 1], totalsc)
  vertical_vector = zeros(Int64, totalspecies * sccounts[rank + 1])

  return MPIGridLandscape(speciescounts, sccounts,
  horizontal_matrix, vertical_vector)
end

function synchronise_from_horizontal!(ml::MPIGridLandscape)
  MPI.Alltoallv!(ml.horizontal_matrix, ml.vertical_vector,
   ml.horizontal_tuple.counts, ml.vertical_tuple.counts, MPI.COMM_WORLD)
end

function synchronise_from_vertical!(ml::MPIGridLandscape)
  MPI.Alltoallv!(ml.vertical_vector, ml.horizontal_matrix,
  ml.vertical_tuple.counts, ml.horizontal_tuple.counts, MPI.COMM_WORLD)
end
