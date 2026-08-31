# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Conformance to Diversity.jl's generics for an `MPIEcosystem` — the distributed counterpart of
# `src/DiversityInterface.jl`, and here for the same reason it is a separate file there: these are
# methods **we** supply on **someone else's** generic, which is a different job from declaring our
# own types, so they do not sit beside them.
#
# **Every one of these has a serial twin in `src/DiversityInterface.jl`**, and the pair should be
# read together — a distributed answer that disagrees with the serial one is a defect whatever else
# is true of it.
#
# Two of them are COLLECTIVE — `_getmetaabundance` calls `Allgatherv` and `_getweight` calls
# `Allreduce` — so every rank must reach them, in the same order, or the run deadlocks. That
# constrains where a guard or an early return may be placed, and is invisible when reading the
# serial file alone.

import Diversity.API: _getabundance
using Diversity.API: _calcabundance, _gettypes
# A rank's own share of the abundances, in the **rows** decomposition — this rank's species over
# every cell, which is what a per-species reduction needs. The non-raw form goes through
# `_calcabundance`, so a phylogeny's species abundances are mapped onto its branches; skipping it
# would answer a different question for a `PhyloBranches` species list and the same one for
# `UniqueTypes`, which is what makes such an omission invisible until a tree is used.
function _getabundance(eco::MPIEcosystem, raw::Bool)
    if raw
        return eco.abundances.rows_matrix
    else
        return _calcabundance(_gettypes(eco),
                              eco.abundances.rows_matrix /
                              sum(eco.abundances.rows_matrix))[1]
    end
end

# Each type's total across the **whole** metacommunity, gathered from every rank. Collective: it
# calls `Allgatherv`, so every rank must reach it or the run deadlocks.
import Diversity.API: _getmetaabundance
function _getmetaabundance(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    ab = sum(_getabundance(eco), dims = 2)
    return MPI.Allgatherv(MPI.VBuffer(ab, eco.sppcounts), comm)
end

# Each subcommunity's share of the total abundance, summed across ranks. Collective, as
# `_getmetaabundance` is, and for the same reason: the weights are a property of the whole
# metacommunity, not of one rank's block.
import Diversity.API: _getweight
function _getweight(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    w = sum(_getabundance(eco, false), dims = 1)
    return MPI.Allreduce(w, +, comm)[1, :]
end

import Diversity.API: _getordinariness!
using Diversity.API: _calcsimilarity, _calcabundance
# **`[MPI-DUP]` — this had BOTH problems the serial path had, and a `grep` over `src/` finds
# neither.** It cached into `eco.ordinariness` (stale the moment anything mutated the ecosystem),
# and it hardcoded the scale as `one(eltype(relab))` — right only when the types carry no
# similarity, and so wrong for exactly the phylogenetic case a scale exists for.
# Now: no cache, and **one `_calcabundance` call giving both the processed abundances and the
# real scale**, as Diversity's own `Metacommunity` does.
# What stays MPI-specific is the slice: each rank owns a block of species rows, so the
# similarity matrix is cut to `sp_rng` on both axes before multiplying.
# The one-argument `_calcordinariness(eco)` method this replaces was not in Diversity's contract
# at all — that generic is `(types, abundances, scale)`.
# Fills `cache.relativeabun` in place, exactly as the serial path does — this is where the big
# runs happen, so the ~1.9 GiB-a-call normalisation matters most here.
function _getordinariness!(eco::MPIEcosystem)
    rows = eco.abundances.rows_matrix
    size(eco.cache.relativeabun) == size(rows) ||
        (eco.cache.relativeabun = similar(rows, Float64))
    eco.cache.relativeabun .= rows ./ sum(rows)
    processed, scale = _calcabundance(_gettypes(eco), eco.cache.relativeabun)
    sp_rng = (eco.abundances.rows_tuple.first):(eco.abundances.rows_tuple.last)
    return _calcsimilarity(eco.spplist, scale)[sp_rng, sp_rng] * processed
end
