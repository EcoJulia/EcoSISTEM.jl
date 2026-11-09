# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A distributed ecosystem's state saved and restored. Every rank holds its own species' rows and
# streams, so saving gathers them on the first rank, which writes the file; starting from one, the
# first rank reads it and sends each rank its own rows.

using Random

# The whole abundance matrix on the first rank, empty on the others.
EcoSISTEM._wholeabundances(eco::MPIEcosystem) = EcoSISTEM.gatherabundance(eco)

# Every species' stream on the first rank, one column a species in species order, each rank sending
# those of its own species; empty on the others.
function EcoSISTEM._wholestreams(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    words = fieldcount(Random.Xoshiro)
    own = eco.firstsp:(eco.firstsp + eco.sppcounts[rank + 1] - 1)
    mine = EcoSISTEM._streamwords(eco.rngs[own])
    if rank == 0
        whole = Matrix{UInt64}(undef, words, sum(eco.sppcounts))
        MPI.Gatherv!(mine, VBuffer(whole, Int32.(words .* eco.sppcounts)), 0,
                     comm)
        return whole
    end
    MPI.Gatherv!(mine, VBuffer(nothing), 0, comm)
    return Matrix{UInt64}(undef, words, 0)
end

# Each rank's own rows of a saved state, sent by the first rank, which alone holds it: a row from its
# species' column of the state, a species the state does not hold left empty. The columns are then
# rebuilt from the rows.
function EcoSISTEM._placesaved!(eco::MPIEcosystem, saved, columns)
    comm = MPI.COMM_WORLD
    rows = eco.abundances.rows_matrix
    if MPI.Comm_rank(comm) == 0
        bounds = vcat(0, cumsum(eco.sppcounts))
        for rank in 0:(MPI.Comm_size(comm) - 1)
            block = zeros(Int64, eco.sppcounts[rank + 1], size(rows, 2))
            for (i, row) in enumerate((bounds[rank + 1] + 1):bounds[rank + 2])
                column = columns[row]
                iszero(column) || (block[i, :] .= view(saved, column, :))
            end
            rank == 0 ? (rows .= block) :
            MPI.Send(block, comm, dest = rank, tag = 0)
        end
    else
        MPI.Recv!(rows, comm, source = 0, tag = 0)
    end
    EcoSISTEM.synchronise_from_rows!(eco.abundances)
    return eco
end
