# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Conformance to Diversity.jl's generics for an `MPIEcosystem` - the distributed counterpart of
# `src/DiversityInterface.jl`, and here for the same reason it is a separate file there: these are
# methods **we** supply on **someone else's** generic, which is a different job from declaring our
# own types, so they do not sit beside them.
#
# **Every hook here has a serial twin in `src/DiversityInterface.jl`**, and the pair should be read
# together - a distributed answer that disagrees with the serial one is a defect whatever else is
# true of it. The two `gather*` functions at the foot are the exception and have no twin: they exist
# only because a distributed run has an answer no single rank holds.
#
# A call graph cannot see any caller of the hooks - their only other mention is the `import` line -
# so every one of them reads as dead code and is not. Deleting one silently breaks the interface for
# a distributed run.
#
# Two of them are COLLECTIVE - `_getmetaabundance` calls `Allgatherv` and `_getweight` calls
# `Allreduce` - so every rank must reach them, in the same order, or the run deadlocks. That
# constrains where a guard or an early return may be placed, and is invisible when reading the
# serial file alone.

# **The normalisation constant is the metacommunity's total, never a rank's own.** Every quantity
# below is a share of the whole, so each divides by this. Collective: an `Allreduce`, so every rank
# must reach it.
#
# Taking `sum(rows_matrix)` instead was the single mistake behind two measured defects at once - the
# weights summed to the **rank count** rather than to 1 (2.0 at two ranks, 4.0 at four), and every
# abundance differed from the serial answer in its last digits because each block was scaled by its
# own total.
function _globaltotal(eco::MPIEcosystem)
    return MPI.Allreduce(sum(eco.abundances.rows_matrix),
                         +, MPI.COMM_WORLD)
end

using Diversity.API: _calcabundance, _gettypes
# This rank's block of the abundances, normalised by the **metacommunity** total - the internal input
# every quantity here is built from, and deliberately not what `getabundance` hands a caller (see the
# refusal below). In the rows decomposition that block is this rank's species over every cell.
#
# The non-raw form goes through `_calcabundance`, so a phylogeny's species abundances are mapped onto
# its branches; skipping it would answer a different question for a `PhyloBranches` species list and
# the same one for `UniqueTypes`, which is what makes such an omission invisible until a tree is used.
function _localabundance(eco::MPIEcosystem, raw::Bool)
    raw && return eco.abundances.rows_matrix
    return _calcabundance(_gettypes(eco),
                          eco.abundances.rows_matrix / _globaltotal(eco))[1]
end

# A block is never an honest answer to `getabundance`, because Diversity's own consumers each
# reduce it over a **different** axis: `_getmetaabundance` sums over cells, `_getweight` over
# species, `_getordinariness!` multiplies by the full similarity matrix. Between them they need every
# element, so whichever axis is partitioned, a block silently answers at least one of them wrongly -
# which is exactly what happened while this returned `rows_matrix`.
#
# So it refuses, and names the function that does assemble the whole matrix honestly. `gatherabundance`
# is collective and costs a full copy on the root, which is why it is a separate call rather than the
# quiet default.
import Diversity.API: _getabundance
function _getabundance(::MPIEcosystem, ::Bool)
    return error("`getabundance` asks for the whole metacommunity's abundances, and no single rank " *
                 "holds them. Use `gatherabundance(eco)`, which assembles them on the root rank, " *
                 "or run serially with `build_ecosystem(...; distributed = false)`.")
end

# Each type's total across the **whole** metacommunity. In the rows decomposition a rank holds every
# cell for its own species, so summing over cells gives that rank's species complete, and the blocks
# concatenate. Collective: `Allgatherv`, so every rank must reach it.
#
# The arity matters and is easy to get wrong: Diversity calls `_getmetaabundance(m, raw::Bool)`.
# A one-argument method is never reached, so the gather silently does not happen and Diversity's own
# two-argument default answers from one rank's block instead - measured as a metaabundance vector of
# 4 species at two ranks and 2 at four, where the metacommunity has 7.
import Diversity.API: _getmetaabundance
#
# `vec` is load-bearing: `sum(..., dims = 2)` gives an n-by-1 **matrix**, and the buffer has to be a
# vector. The receive buffer is the whole metacommunity's length while the send buffer is this rank's
# block, which is what `Allgatherv!` needs and what the single-buffer in-place form cannot express.
function _getmetaabundance(eco::MPIEcosystem, raw::Bool)
    ab = vec(sum(_localabundance(eco, raw), dims = 2))
    out = similar(ab, Int64(sum(eco.sppcounts)))
    MPI.Allgatherv!(ab, MPI.VBuffer(out, eco.sppcounts), MPI.COMM_WORLD)
    return out
end

# Each subcommunity's share of the total abundance. A cell's species are split across ranks in this
# decomposition, so each rank contributes its own species' part and the parts are summed. Collective,
# and the result is the whole per-cell vector on every rank - it is only as long as the grid, so it
# costs nothing to hold everywhere.
import Diversity.API: _getweight
function _getweight(eco::MPIEcosystem)
    w = sum(_localabundance(eco, false), dims = 1)
    return MPI.Allreduce(w, +, MPI.COMM_WORLD)[1, :]
end

# The scale that pairs with those abundances, from the same `_calcabundance` call Diversity's own
# `Metacommunity` uses. Without a method here the `AbstractEcosystem` one is reached, and it asks for
# `eco.abundances.matrix` - a field an `MPIGridLandscape` has not got, so it threw a `FieldError`
# rather than answering.
import Diversity.API: _getscale
function _getscale(eco::MPIEcosystem)
    return _calcabundance(_gettypes(eco),
                          eco.abundances.rows_matrix / _globaltotal(eco))[2]
end

import Diversity.API: _getordinariness!
using Diversity.API: _calcsimilarity, _calcabundance
# **`[MPI-DUP]` - this had BOTH problems the serial path had, and a `grep` over `src/` finds
# neither.** It cached into `eco.ordinariness` (stale the moment anything mutated the ecosystem),
# and it hardcoded the scale as `one(eltype(relab))` - right only when the types carry no
# similarity, and so wrong for exactly the phylogenetic case a scale exists for.
# Now: no cache, and **one `_calcabundance` call giving both the processed abundances and the
# real scale**, as Diversity's own `Metacommunity` does.
# **Computed in the COLUMN partition, and that is what makes it right rather than merely faster.**
# Ordinariness multiplies the similarity matrix by the abundances, so a species' ordinariness in a
# cell depends on every *other* species in that cell. In the row partition a cell's species are split
# across ranks, so the similarity matrix had to be cut to this rank's species on both axes - which
# silently discards every similarity between a species here and a species elsewhere. That is exact
# for `UniqueTypes`, whose similarity is the identity, and wrong for any phylogeny, where branches
# are shared across the partition. Every MPI test in the suite uses `UniqueTypes`, so nothing caught
# it.
#
# `parent(dimcols)` is every species over this rank's own cells, so the full similarity matrix
# applies with no slice at all, and the result needs no communication: a cell's value is complete on
# the rank that owns it.
#
# The blocked array is read from directly rather than materialised: the broadcast into
# `cache.relativeabun` already copies it once, and materialising first would copy it twice.
# `Matrix{Float64}(undef, ...)` rather than `similar` because `similar` on a blocked array returns
# another blocked array, which the cache's declared `Matrix{Float64}` will not hold.
# The one-argument `_calcordinariness(eco)` method this replaces was not in Diversity's contract
# at all - that generic is `(types, abundances, scale)`.
# Fills `cache.relativeabun` in place, exactly as the serial path does - this is where the big
# runs happen, so the ~1.9 GiB-a-call normalisation matters most here.
function _getordinariness!(eco::MPIEcosystem)
    cols = parent(eco.abundances.dimcols)
    size(eco.cache.relativeabun) == size(cols) ||
        (eco.cache.relativeabun = Matrix{Float64}(undef, size(cols)))
    eco.cache.relativeabun .= cols ./ _globaltotal(eco)
    processed, scale = _calcabundance(_gettypes(eco), eco.cache.relativeabun)
    return _calcsimilarity(eco.spplist, scale) * processed
end

# Each type's total ordinariness across the **whole** metacommunity. Diversity's default sums
# `_getordinariness!` over the second dimension, which under the column partition is only this rank's
# cells, so the sum has to be completed across ranks. Collective.
#
# Without this every measure that divides by it - beta, rho and gamma - silently uses one rank's
# share of the metacommunity as though it were all of it.
import Diversity.API: _getmetaordinariness!
function _getmetaordinariness!(eco::MPIEcosystem)
    ord = vec(sum(_getordinariness!(eco), dims = 2))
    return MPI.Allreduce(ord, +, MPI.COMM_WORLD)
end

# The inputs every Diversity measure is built from, in the column partition: all species over this
# rank's own cells, normalised by the metacommunity total, and the matching slice of the global
# per-cell weights. The two must be cut the same way - a full-length weight vector against a block of
# columns is a shape error at best and a silently wrong power mean at worst.
#
# Materialised here, unlike in `_getordinariness!`: this array is *stored* in the measure rather than
# broadcast into an existing buffer, so it has to be a concrete matrix either way.
function _measureinputs(eco::MPIEcosystem)
    cells = (eco.abundances.cols_tuple.first):(eco.abundances.cols_tuple.last)
    ab = Matrix(parent(eco.abundances.dimcols)) ./ _globaltotal(eco)
    return (abundances = _calcabundance(_gettypes(eco), ab)[1],
            weights = _getweight(eco)[cells])
end

# Build one of Diversity's measures from those inputs. Every one of its own constructors opens with
# `getabundance(meta)`, which refuses for a distributed ecosystem - a block is not the metacommunity's
# abundances, and saying so is the whole point of the refusal - so the measures are built here
# instead, from the block explicitly.
#
# `value` is the only thing that differs between the seven, which is why it arrives as a function of
# the three quantities they are all written in terms of.
function _mpimeasure(::Type{M}, eco::MPIEcosystem, value) where {M}
    inputs = _measureinputs(eco)
    ab, ws = inputs.abundances, inputs.weights
    v = value(_getordinariness!(eco), _getmetaordinariness!(eco), ws)
    return M{eltype(ab), typeof(ab), typeof(v), typeof(eco)}(ab, ws, v, eco)
end

# A measure built on a distributed ecosystem covers **this rank's cells**, so that is what its rows
# must be labelled with. Diversity's default asks the ecosystem, which correctly answers with every
# cell in the grid - and `subdiv` then broadcasts a 39-long vector of values against 77 names and
# throws.
#
# Defined on the **measure** rather than on the ecosystem deliberately: `getsubcommunitynames(eco)`
# and `countsubcommunities(eco)` both describe the whole grid and must go on doing so - the hot loop
# reads the latter for exactly that. It is the measure, not the ecosystem, that is rank-local.
function Diversity.getsubcommunitynames(dm::Diversity.DiversityMeasure{FP, A,
                                                                       D,
                                                                       <:MPIEcosystem}) where {FP <:
                                                                                               AbstractFloat,
                                                                                               A <:
                                                                                               AbstractMatrix,
                                                                                               D <:
                                                                                               AbstractArray}
    eco = Diversity._getmeta(dm)
    cells = (eco.abundances.cols_tuple.first):(eco.abundances.cols_tuple.last)
    return Diversity.getsubcommunitynames(eco)[cells]
end

# The seven, each mirroring the body of Diversity's own constructor exactly. A subcommunity measure
# built this way is **complete for this rank's cells**: a cell's value is a reduction down its own
# column, and the column partition puts every species of that cell on one rank.
#
# A *metacommunity* value taken from one of these would be this rank's cells alone, which is what
# `gatherdiversity` exists to assemble correctly.
function Diversity.RawAlpha(eco::MPIEcosystem)
    return _mpimeasure(Diversity.RawAlpha, eco,
                       (o, m, w) -> o .^ -1)
end
function Diversity.NormalisedAlpha(eco::MPIEcosystem)
    return _mpimeasure(Diversity.NormalisedAlpha, eco, (o, m, w) -> w' ./ o)
end
function Diversity.RawBeta(eco::MPIEcosystem)
    return _mpimeasure(Diversity.RawBeta, eco,
                       (o, m, w) -> o ./ m)
end
function Diversity.NormalisedBeta(eco::MPIEcosystem)
    return _mpimeasure(Diversity.NormalisedBeta, eco,
                       (o, m, w) -> o ./ (m .* w'))
end
function Diversity.RawRho(eco::MPIEcosystem)
    return _mpimeasure(Diversity.RawRho, eco,
                       (o, m, w) -> m ./ o)
end
function Diversity.NormalisedRho(eco::MPIEcosystem)
    return _mpimeasure(Diversity.NormalisedRho, eco,
                       (o, m, w) -> (m .* w') ./ o)
end
function Diversity.Gamma(eco::MPIEcosystem)
    return _mpimeasure(Diversity.Gamma, eco,
                       (o, m, w) -> fill!(similar(w), 1)' ./ m)
end

# == Bringing a distributed answer back ============================================================
# `gatherabundance` and `gatherdiversity` are **ours**, not hooks on a Diversity generic, so strictly
# they are not conformance. They live here anyway because that is what they are *for*: both exist to
# hand Diversity's world a whole-metacommunity answer that no single rank holds, and reading them
# beside the hooks they assemble is what makes either legible.
#
# They are also the counterpart to the refusal above. `getabundance` declines to answer for a
# distributed ecosystem precisely because these two do it honestly and say what they cost.

using DataFrames: DataFrame, ncol

function EcoSISTEM.gatherabundance(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    true_abuns = zeros(Int64, counttypes(eco), countsubcommunities(eco))
    if rank == 0
        # **Counted over the COLUMN partition, because that is what is being sent.**
        # `reshaped_cols` holds *every* species for *this rank's* cells, so a rank contributes
        # `counttypes(eco) * sccounts[rank]` values. The previous
        # `sppcounts .* sum(sccounts)` - this rank's *species* by *all* cells - is the row
        # partition, and the two agree only when species and cells both divide evenly across the
        # ranks. They always did in the old test fixture (8 species on a 4 × 4 grid), so the
        # mismatch never showed; with 7 species on 77 cells rank 0 sends 273 values while the old
        # expression asked for 308, and `MPI_Gatherv` fails outright.
        output_vbuf = VBuffer(true_abuns,
                              Int32.(counttypes(eco) .* eco.sccounts))
    else
        output_vbuf = VBuffer(nothing)
    end
    MPI.Gatherv!(vcat(eco.abundances.reshaped_cols...)[1:end], output_vbuf, 0,
                 comm)
    return true_abuns
end

# **A concatenation, not a combination**, and that is the whole of the change. Each rank computes
# the measure for its own cells in the column partition, where a cell holds every species, so its
# values are already the final ones - there is nothing to combine across ranks, only to assemble.
#
# What this replaced took each rank's per-cell vector and merged them with a power mean weighted by
# `sum(rows_matrix)`, one scalar for the whole rank. That is wrong twice over: the correct weight for
# combining a *cell* varies by cell rather than by rank, and in the column partition no cell is split
# across ranks in the first place, so no weighting arises at all.
#
# `Allgatherv!` rather than a gather to the root: the result is one value per cell, so it is small
# enough to hold everywhere, and a diversity that existed only on rank 0 would make the answer depend
# on which rank asked for it.
#
# Rows are grouped by `q` and then by cell, so each order is gathered separately - concatenating the
# ranks' whole frames would interleave the orders instead of the cells.
function EcoSISTEM.gatherdiversity(eco::MPIEcosystem, divmeasure::F,
                                   q) where {F <: Function}
    comm = MPI.COMM_WORLD
    mine = divmeasure(eco, q)
    counts = Int32.(eco.sccounts)
    total = Int64(sum(counts))
    names = Diversity.getsubcommunitynames(eco)
    ncol(mine) == 8 ||
        error("`gatherdiversity` cannot assemble a measure carrying extra output columns " *
              "(this one has $(ncol(mine))); the added columns of a phylogeny are not per-cell " *
              "and have no defined way to concatenate.")
    frames = map(unique(mine[!, :q])) do order
        rows = mine[mine[!, :q] .== order, :]
        full = similar(rows[!, :diversity], total)
        MPI.Allgatherv!(rows[!, :diversity], MPI.VBuffer(full, counts), comm)
        first_row = first(eachrow(rows))
        return DataFrame(div_type = first_row.div_type,
                         measure = first_row.measure,
                         q = order,
                         type_level = first_row.type_level,
                         type_name = first_row.type_name,
                         partition_level = first_row.partition_level,
                         partition_name = names,
                         diversity = full)
    end
    return reduce(append!, frames)
end
