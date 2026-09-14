# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a run keeps as it goes: recorders, callable values handed to `simulate!` in place of a
# callback.

using Unitful
using JLD2: @save

"""
    AbstractRecorder

A value that keeps something from a run each time a [`simulate!`](@ref) schedule fires -
[`RecordAbundance`](@ref), [`RecordDiversity`](@ref) or [`SaveAbundance`](@ref). Pass one to
`simulate!` in place of a callback, or call `recorder(eco, occurrence)` from a callback of your own to
keep several things from one run. Each holds the [`Provenance`](@ref) of the run as it stood at its
last write, which [`provenance`](@ref) returns.
"""
abstract type AbstractRecorder end

"""
    RecordAbundance(storage::AbstractArray)
    RecordAbundance(eco::Ecosystem, ntimes::Integer; maxspecies = length(eco.spplist.abun))

Record every species' abundance in every cell into `storage`, one slice per occurrence:
`storage[:, :, k]` holds species by cells at occurrence `k`. Under MPI the ranks' abundances are
gathered and the root writes them, every rank taking part in the gather.

The second form allocates the storage with [`generate_storage`](@ref). A run reaching more
occurrences than the storage has room for, or gaining more species, is refused when it tries to
record them.

# Arguments

  - `storage`: an integer array of at least species × cells × occurrences, as from
    `generate_storage`.
  - `eco`: the ecosystem whose species and cells size the storage.
  - `ntimes`: how many occurrences to leave room for.
  - `maxspecies`: how many species to leave room for, counting any an intervention adds.
"""
struct RecordAbundance{S <: AbstractArray} <: AbstractRecorder
    storage::S
    record::Base.RefValue{Union{Nothing, Provenance}}

    function RecordAbundance(storage::AbstractArray)
        return new{typeof(storage)}(storage,
                                    Ref{Union{Nothing, Provenance}}(nothing))
    end
end

function RecordAbundance(eco::Ecosystem, ntimes::Integer;
                         maxspecies::Integer = length(eco.spplist.abun))
    return RecordAbundance(generate_storage(eco, Int64(ntimes), 1,
                                            maxspecies = Int64(maxspecies)))
end

# Record the abundances at this occurrence into its slice, and take the run's provenance as it now
# stands.
function (recorder::RecordAbundance)(eco::Ecosystem, occurrence::NamedTuple)
    _checkslot(recorder.storage, occurrence.count)
    _record!(recorder.storage, eco, occurrence.count)
    recorder.record[] = provenance(eco)
    return recorder
end

function Base.show(io::IO, recorder::RecordAbundance)
    return print(io, "RecordAbundance(", join(size(recorder.storage), " × "),
                 " storage)")
end

"""
    RecordDiversity(storage::AbstractArray, measure, qs::AbstractVector)

Record a subcommunity diversity `measure` - such as Diversity's `norm_sub_alpha` - at each of the
orders `qs` into `storage`, one slice per occurrence: `storage[:, :, k]` holds cells by orders at
occurrence `k`. Under MPI each rank's cells are gathered with [`gatherdiversity`](@ref), so every rank
holds and writes the whole answer.

# Arguments

  - `storage`: a float array of at least cells × orders × occurrences, as from
    `generate_storage(eco, length(qs), ntimes, reps)`.
  - `measure`: the subcommunity diversity function, returning a `DataFrame` with a `:diversity`
    column. A metacommunity or individual measure is refused.
  - `qs`: the orders to measure at.
"""
struct RecordDiversity{S <: AbstractArray, F, Q <: AbstractVector} <:
       AbstractRecorder
    storage::S
    measure::F
    qs::Q
    record::Base.RefValue{Union{Nothing, Provenance}}

    function RecordDiversity(storage::AbstractArray, measure,
                             qs::AbstractVector)
        return new{typeof(storage), typeof(measure), typeof(qs)}(storage,
                                                                 measure, qs,
                                                                 Ref{Union{Nothing,
                                                                           Provenance}}(nothing))
    end
end

# Record the measure at this occurrence into its slice, and take the run's provenance as it now
# stands.
function (recorder::RecordDiversity)(eco::Ecosystem, occurrence::NamedTuple)
    _checkslot(recorder.storage, occurrence.count)
    _writediversity!(recorder, recorder.measure(eco, recorder.qs),
                     occurrence.count)
    recorder.record[] = provenance(eco)
    return recorder
end

function Base.show(io::IO, recorder::RecordDiversity)
    return print(io, "RecordDiversity(", join(size(recorder.storage), " × "),
                 " storage, ", recorder.measure, ", ", recorder.qs, ")")
end

"""
    SaveAbundance(folder::AbstractString, name::AbstractString)

Save every species' abundance in every cell to a JLD2 file at each occurrence, in `folder`, named
`name` followed by the occurrence's number counted from `00` - so `<name>00.jld2` holds the first -
with the run's provenance beside each as `<name>00.provenance.toml`. Under MPI the ranks' abundances
are gathered and the root writes the files.

# Arguments

  - `folder`: the directory to write into, created if it does not exist.
  - `name`: the start of each file's name.
"""
struct SaveAbundance <: AbstractRecorder
    folder::String
    name::String
    record::Base.RefValue{Union{Nothing, Provenance}}

    function SaveAbundance(folder::AbstractString, name::AbstractString)
        return new(String(folder), String(name),
                   Ref{Union{Nothing, Provenance}}(nothing))
    end
end

# Save the abundances at this occurrence, and the run's provenance as it now stands beside them.
function (recorder::SaveAbundance)(eco::Ecosystem, occurrence::NamedTuple)
    recorder.record[] = provenance(eco)
    _saveabundance(recorder, eco.abundances.matrix, occurrence.count)
    return recorder
end

function Base.show(io::IO, recorder::SaveAbundance)
    return print(io, "SaveAbundance(", repr(recorder.folder), ", ",
                 repr(recorder.name), ")")
end

# == Functions ==================================================================================

"""
    provenance(recorder::AbstractRecorder)

Return the [`Provenance`](@ref) of the run a recorder has kept from, as it stood at the recorder's last
write - how far the run had gone, and every intervention that had acted by then - or `nothing` before
the recorder has written.
"""
provenance(recorder::AbstractRecorder) = recorder.record[]

# Refuse an occurrence the storage has no slice for, rather than let the write fail on an index.
function _checkslot(storage::AbstractArray, count::Integer)
    count <= size(storage, 3) && return nothing
    return error("this run has reached occurrence $count, but the recording has room for " *
                 "$(size(storage, 3)). A run starting from zero has the starting state as its " *
                 "first occurrence, so size the storage for every multiple of `every` from zero " *
                 "to `duration`; a run continuing an earlier one has one fewer.")
end

# Refuse a diversity result that is not one value per subcommunity at each order - a metacommunity or
# individual measure - since recording or assembling it cell by cell gives values that mean nothing.
function _checksubcommunitylevel(frame, measure, what::AbstractString)
    levels = unique(frame[!, :partition_level])
    typelevels = unique(frame[!, :type_level])
    levels == ["subcommunity"] && typelevels == ["types"] && return nothing
    given = typelevels == ["types"] ? join(levels, ", ") :
            join(typelevels, ", ") * " by " * join(levels, ", ")
    return error("$what takes a subcommunity diversity measure, one value per cell at each " *
                 "order, such as `norm_sub_alpha`; `$measure` gives $given diversity, which " *
                 "is not yet supported here.")
end

# Write one occurrence of a diversity measure into its slice: the measure's rows run through the cells
# for each order in turn.
function _writediversity!(recorder::RecordDiversity, frame, count::Integer)
    _checksubcommunitylevel(frame, recorder.measure, "`RecordDiversity`")
    diversity = frame[!, :diversity]
    orders = length(recorder.qs)
    recorder.storage[:, :, count] = reshape(diversity,
                                            length(diversity) ÷ orders, orders)
    return recorder
end

# Save one occurrence's abundances and the provenance held for it, numbered from `00`.
function _saveabundance(recorder::SaveAbundance, abun::AbstractMatrix,
                        count::Integer)
    mkpath(recorder.folder)
    stem = joinpath(recorder.folder, recorder.name * lpad(count - 1, 2, '0'))
    @save stem * ".jld2" abun
    write_provenance(stem * ".provenance.toml", recorder.record[])
    return recorder
end
