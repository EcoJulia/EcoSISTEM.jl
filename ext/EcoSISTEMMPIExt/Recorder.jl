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

# Each rank's cells gathered into the whole metacommunity's answer, which every rank holds and writes.
function (recorder::EcoSISTEM.RecordDiversity)(eco::MPIEcosystem,
                                               occurrence::NamedTuple)
    EcoSISTEM._checkslot(recorder.storage, occurrence.count)
    frame = EcoSISTEM.gatherdiversity(eco, recorder.measure, recorder.qs)
    EcoSISTEM._writediversity!(recorder, frame, occurrence.count)
    recorder.record[] = EcoSISTEM.provenance(eco)
    return recorder
end

# The abundances gathered from every rank, saved by the root.
function (recorder::EcoSISTEM.SaveAbundance)(eco::MPIEcosystem,
                                             occurrence::NamedTuple)
    abuns = EcoSISTEM.gatherabundance(eco)
    recorder.record[] = EcoSISTEM.provenance(eco)
    MPI.Comm_rank(MPI.COMM_WORLD) == 0 || return recorder
    EcoSISTEM._saveabundance(recorder, abuns, occurrence.count)
    return recorder
end
