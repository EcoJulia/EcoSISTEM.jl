# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Conformance to Diversity.jl's generics for an `MPIEcosystem` - the distributed counterpart of
# `src/DiversityInterface.jl`, and here for the same reason it is a separate file there: these are
# methods **we** supply on **someone else's** generic, which is a different job from declaring our
# own types, so they do not sit beside them.
#
# **Every one of these has a serial twin in `src/DiversityInterface.jl`**, and the pair should be
# read together - a distributed answer that disagrees with the serial one is a defect whatever else
# is true of it.
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
