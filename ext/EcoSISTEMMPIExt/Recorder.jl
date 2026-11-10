# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The recorders on a distributed ecosystem. Every rank runs the callback a recorder is called from, so
# each rank takes part in the gather a recorder needs, and none returns before it.

# The abundances gathered from every rank, written into this occurrence's slice by the root.
function (recorder::EcoSISTEM.RecordAbundance)(eco::MPIEcosystem,
                                               occurrence::NamedTuple)
    EcoSISTEM._checkslot(recorder.storage, occurrence.count)
    abuns = EcoSISTEM.gatherabundance(eco)
    recorder.record[] = EcoSISTEM.provenance(eco)
    MPI.Comm_rank(MPI.COMM_WORLD) == 0 || return recorder
    recorder.storage[1:size(abuns, 1), :, occurrence.count] = abuns
    return recorder
end

# Each rank's own species summed over the grid, gathered into their places, and each rank's share of
# every cell summed across the ranks, so every rank holds both totals whole.
function (recorder::EcoSISTEM.RecordTotals)(eco::MPIEcosystem,
                                            occurrence::NamedTuple)
    EcoSISTEM._checkslot(recorder.species, occurrence.count, 2)
    comm = MPI.COMM_WORLD
    matrix = eco.abundances.rows_matrix
    mine = vec(sum(matrix, dims = 2))
    species = similar(mine, Int64(sum(eco.sppcounts)))
    MPI.Allgatherv!(mine, MPI.VBuffer(species, eco.sppcounts), comm)
    recorder.species[:, occurrence.count] = species
    recorder.cells[:, occurrence.count] = MPI.Allreduce(vec(sum(matrix,
                                                                dims = 1)),
                                                        +, comm)
    recorder.record[] = EcoSISTEM.provenance(eco)
    return recorder
end

# Each rank's cells gathered into the whole metacommunity's answer, which every rank holds and writes.
function (recorder::EcoSISTEM.RecordDiversity)(eco::MPIEcosystem,
                                               occurrence::NamedTuple)
    _recorddiversity!(recorder, recorder.measure, eco, occurrence.count)
    recorder.record[] = EcoSISTEM.provenance(eco)
    return recorder
end

# The storage form: the measure function's values for each rank's cells, by `gatherdiversity`.
function _recorddiversity!(recorder::EcoSISTEM.RecordDiversity, measure,
                           eco::MPIEcosystem, count::Integer)
    EcoSISTEM._checkslot(recorder.storage, count)
    frame = EcoSISTEM.gatherdiversity(eco, measure, recorder.qs)
    return EcoSISTEM._writediversity!(recorder, frame, count)
end

using Tables: Tables
# The level form: every measure built by Diversity from this rank's snapshot, its subcommunity
# values at every order gathered in one `Allgatherv`, and each metacommunity value their power mean
# over the whole grid's weights, as Diversity's own metacommunity level is defined - so every rank
# holds the whole answer. Every refusal comes before the first collective, alike on every rank.
function _recorddiversity!(recorder::EcoSISTEM.RecordDiversity, ::Nothing,
                           eco::MPIEcosystem, count::Integer)
    EcoSISTEM._checkslots(recorder, count)
    eco.spplist.types isa UniqueTypes ||
        error("`RecordDiversity` records diversity under MPI only for species with no " *
              "similarity between them; this species list's types are " *
              "$(nameof(typeof(eco.spplist.types))). Record it serially, with " *
              "`build_ecosystem(...; distributed = false)`.")
    measures = EcoSISTEM._recordedmeasures(recorder)
    qs = recorder.qs
    made = _ranksnapshot(eco)
    ncells = length(made.snapshot.weights)
    mine = Array{Float64, 3}(undef, length(measures), length(qs), ncells)
    for (i, measure) in enumerate(measures)
        columns = subdiv(Tables.columntable, measure(made.snapshot), qs)
        length(columns.diversity) == ncells * length(qs) ||
            error("`subdiv` gave $(length(columns.diversity)) rows for `$measure`, where " *
                  "$(ncells * length(qs)) were expected.")
        for (j, q) in enumerate(qs)
            start = (j - 1) * ncells
            EcoSISTEM._checkblock(columns, start + 1, subcommunityDiversity, q,
                                  measure)
            for c in 1:ncells
                mine[i, j, c] = columns.diversity[start + c]
            end
        end
    end
    whole = Array{Float64, 3}(undef, length(measures), length(qs),
                              length(made.weights))
    MPI.Allgatherv!(mine,
                    MPI.VBuffer(whole,
                                Int32.(length(measures) * length(qs) .*
                                       eco.sccounts)), MPI.COMM_WORLD)
    for (k, measure) in enumerate(recorder.measures.subcommunity)
        i = findfirst(==(measure), measures)
        for j in eachindex(qs)
            recorder.subcommunity[k, :, j, count] = view(whole, i, j, :)
        end
    end
    for (k, measure) in enumerate(recorder.measures.metacommunity)
        i = findfirst(==(measure), measures)
        for (j, q) in enumerate(qs)
            recorder.metacommunity[k, j,
                                   count] = Diversity.powermean(view(whole,
                                                                     i, j,
                                                                     :),
                                                                one(q) -
                                                                q,
                                                                made.weights)
        end
    end
    return recorder
end
