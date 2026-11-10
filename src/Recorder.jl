# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a run keeps as it goes: recorders, callable values handed to `simulate!` in place of a
# callback.

using Unitful
using Diversity: Metacommunity, diversity, subcommunityDiversity,
                 metacommunityDiversity, DiversityMeasure
using Tables: Tables

"""
    AbstractRecorder

A value that keeps something from a run each time a [`simulate!`](@ref) schedule fires -
[`RecordAbundance`](@ref), [`RecordTotals`](@ref), [`RecordDiversity`](@ref) or
[`SaveAbundance`](@ref). Pass one to
`simulate!` in place of a callback, or call `recorder(eco, occurrence)` from a callback of your own to
keep several things from one run. Each holds the [`Provenance`](@ref) of the run as it stood at its
last write, which [`provenance`](@ref) returns.
"""
abstract type AbstractRecorder end

"""
    RecordAbundance(storage::AbstractArray)
    RecordAbundance(eco::Ecosystem, ntimes::Integer)

Record every species' abundance in every cell into `storage`, one slice per occurrence:
`storage[:, :, k]` holds species by cells at occurrence `k`. Under MPI the ranks' abundances are
gathered and the root writes them, every rank taking part in the gather.

The second form allocates the storage with [`generate_storage`](@ref), a row for every species of
the ecosystem's pool. A run reaching more occurrences than the storage has room for is refused when
it tries to record them.

# Arguments

  - `storage`: an integer array of at least species × cells × occurrences, as from
    `generate_storage`.
  - `eco`: the ecosystem whose species and cells size the storage.
  - `ntimes`: how many occurrences to leave room for.
  - `maxspecies`: deprecated, and ignored: every species a run may hold has its row already.
"""
struct RecordAbundance{S <: AbstractArray} <: AbstractRecorder
    storage::S
    record::Base.RefValue{Union{Nothing, Provenance}}

    function RecordAbundance(storage::AbstractArray)
        return new{typeof(storage)}(storage,
                                    Ref{Union{Nothing, Provenance}}(nothing))
    end
end

function RecordAbundance(eco::Ecosystem, ntimes::Integer; maxspecies = nothing)
    isnothing(maxspecies) || _maxspeciesdeprecated(:RecordAbundance)
    return RecordAbundance(generate_storage(eco, Int64(ntimes), 1))
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
    RecordTotals(eco::AbstractEcosystem, ntimes::Integer)

Record every species' individuals summed over the grid, and every cell's summed over the species,
one column per occurrence: `recorder.species[:, k]` holds each species' total at occurrence `k`, in
the species list's order, and `recorder.cells[:, k]` each cell's, in the grid's order. Every species
of the pool has its row, an arrival's included; an occurrence not yet reached holds zeros. Under MPI
the ranks' totals are summed, every rank taking part, and every rank holds them.

They are small enough to keep at every occurrence where the abundances are not, and say which
species a run loses and when, and where its individuals are.

# Arguments

  - `eco`: the ecosystem whose species and cells size the storage.
  - `ntimes`: how many occurrences to leave room for.

# Fields

  - `species`: the species' totals, species × occurrences.
  - `cells`: the cells' totals, cells × occurrences.
"""
struct RecordTotals <: AbstractRecorder
    species::Matrix{Int64}
    cells::Matrix{Int64}
    record::Base.RefValue{Union{Nothing, Provenance}}

    function RecordTotals(eco::AbstractEcosystem, ntimes::Integer)
        ntimes >= 0 ||
            error("`RecordTotals` needs room for zero or more occurrences, not $ntimes.")
        return new(zeros(Int64, length(eco.spplist.names), ntimes),
                   zeros(Int64, countsubcommunities(eco), ntimes),
                   Ref{Union{Nothing, Provenance}}(nothing))
    end
end

# Sum the abundances both ways into this occurrence's columns, and take the run's provenance as it
# now stands.
function (recorder::RecordTotals)(eco::Ecosystem, occurrence::NamedTuple)
    _checkslot(recorder.species, occurrence.count, 2)
    matrix = eco.abundances.matrix
    sum!(view(recorder.species, :, occurrence.count), matrix)
    sum!(reshape(view(recorder.cells, :, occurrence.count), 1, :), matrix)
    recorder.record[] = provenance(eco)
    return recorder
end

function Base.show(io::IO, recorder::RecordTotals)
    return print(io, "RecordTotals(", size(recorder.species, 1), " species, ",
                 size(recorder.cells, 1), " cells, ",
                 size(recorder.species, 2), " occurrences)")
end

"""
    RecordDiversity(eco::AbstractEcosystem, ntimes::Integer; qs, subcommunity = [],
                    metacommunity = [])
    RecordDiversity(storage::AbstractArray, measure, qs::AbstractVector)

Record diversity at each of the orders `qs`, one slice per occurrence.

The first form records Diversity's measure types - `NormalisedAlpha`, `RawBeta`, `Gamma` and the
rest - at the subcommunity level, a value per cell, and at the metacommunity level, one value for
the whole grid, each measure at the levels it is listed under: `recorder.subcommunity[i, :, j, k]`
holds the `i`th subcommunity measure in every cell at order `qs[j]` and occurrence `k`, and
`recorder.metacommunity[i, j, k]` the `i`th metacommunity measure. An occurrence not yet reached
holds `NaN`. Every measure of one occurrence is taken from one snapshot of the abundances, and a
measure listed at both levels is computed once for both. Under MPI each rank measures its own cells
and the values are gathered, the metacommunity value taken from them as Diversity takes it, so every
rank holds the whole answer; a species list with similarity between its species, such as a
phylogeny, is refused there.

The second form records a subcommunity diversity function - such as Diversity's `norm_sub_alpha` -
into `storage`: `storage[:, :, k]` holds cells by orders at occurrence `k`. Under MPI each rank's
cells are gathered with [`gatherdiversity`](@ref), so every rank holds and writes the whole answer.

# Arguments

  - `eco`: the ecosystem whose cells size the storage.
  - `ntimes`: how many occurrences to leave room for.
  - `qs`: the orders to measure at.
  - `subcommunity`: the measure types to record in every cell.
  - `metacommunity`: the measure types to record for the whole grid. At least one measure is needed
    between the two.
  - `storage`: a float array of at least cells × orders × occurrences, as from
    `generate_storage(eco, length(qs), ntimes, reps)`.
  - `measure`: the subcommunity diversity function, returning a `DataFrame` with a `:diversity`
    column. A metacommunity or individual measure is refused.

# Fields

  - `subcommunity`: measures × cells × orders × occurrences, or `nothing` where the first form
    lists no subcommunity measure, and in the second form.
  - `metacommunity`: measures × orders × occurrences, or `nothing` likewise.
  - `measures`: the measure types at each level, a named tuple of `subcommunity` and
    `metacommunity`.
  - `storage`, `measure`: the second form's, as given; `nothing` in the first.
  - `qs`: the orders, as given.
"""
struct RecordDiversity{S, F, Q <: AbstractVector, B, M} <: AbstractRecorder
    storage::S
    measure::F
    qs::Q
    subcommunity::B
    metacommunity::M
    measures::@NamedTuple{subcommunity::Vector{Type},
                          metacommunity::Vector{Type}}
    record::Base.RefValue{Union{Nothing, Provenance}}

    function RecordDiversity(eco::AbstractEcosystem, ntimes::Integer;
                             qs::AbstractVector{<:Real},
                             subcommunity = Type[], metacommunity = Type[])
        measures = (subcommunity = _measuretypes(subcommunity, :subcommunity),
                    metacommunity = _measuretypes(metacommunity,
                                                  :metacommunity))
        isempty(measures.subcommunity) && isempty(measures.metacommunity) &&
            error("`RecordDiversity` needs a measure to record: list Diversity's measure types " *
                  "as `subcommunity`, for a value in every cell, or `metacommunity`, for one " *
                  "value for the whole grid.")
        isempty(qs) &&
            error("`RecordDiversity` needs at least one order in `qs`.")
        ntimes >= 0 ||
            error("`RecordDiversity` needs room for zero or more occurrences, not $ntimes.")
        cells = countsubcommunities(eco)
        sub = isempty(measures.subcommunity) ? nothing :
              fill(NaN, length(measures.subcommunity), cells, length(qs),
                   ntimes)
        meta = isempty(measures.metacommunity) ? nothing :
               fill(NaN, length(measures.metacommunity), length(qs), ntimes)
        return new{Nothing, Nothing, typeof(qs), typeof(sub), typeof(meta)}(nothing,
                                                                            nothing,
                                                                            qs,
                                                                            sub,
                                                                            meta,
                                                                            measures,
                                                                            Ref{Union{Nothing,
                                                                                      Provenance}}(nothing))
    end

    function RecordDiversity(storage::AbstractArray, measure,
                             qs::AbstractVector)
        return new{typeof(storage), typeof(measure), typeof(qs), Nothing,
                   Nothing}(storage, measure, qs, nothing, nothing,
                            (subcommunity = Type[], metacommunity = Type[]),
                            Ref{Union{Nothing, Provenance}}(nothing))
    end
end

# Record this occurrence into its slices, and take the run's provenance as it now stands.
function (recorder::RecordDiversity)(eco::Ecosystem, occurrence::NamedTuple)
    _recorddiversity!(recorder, recorder.measure, eco, occurrence.count)
    recorder.record[] = provenance(eco)
    return recorder
end

function Base.show(io::IO, recorder::RecordDiversity{Nothing})
    return print(io, "RecordDiversity(subcommunity = ",
                 _measurenames(recorder.measures.subcommunity),
                 ", metacommunity = ",
                 _measurenames(recorder.measures.metacommunity), ", qs = ",
                 recorder.qs, ")")
end

function Base.show(io::IO, recorder::RecordDiversity)
    return print(io, "RecordDiversity(", join(size(recorder.storage), " × "),
                 " storage, ", recorder.measure, ", ", recorder.qs, ")")
end

"""
    SaveAbundance(folder::AbstractString, name::AbstractString)

Save the run's state at each occurrence, as [`write_state`](@ref) does, to a JLD2 file in `folder`
named `name` followed by the occurrence's number counted from `00` - so `<name>00.jld2` holds the
first - with the run's provenance beside each as `<name>00.provenance.toml`. Each holds every
species' abundance in every cell, and is a state a later run can start from with
`build_ecosystem(species, habitat; from = path)`. Under MPI the ranks' states are gathered and the
root writes the files.

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

# Save the state at this occurrence, and the run's provenance as it now stands beside it; under MPI a
# collective call, every rank taking part.
function (recorder::SaveAbundance)(eco::AbstractEcosystem,
                                   occurrence::NamedTuple)
    recorder.record[] = provenance(eco)
    write_state(joinpath(recorder.folder,
                         recorder.name * lpad(occurrence.count - 1, 2, '0') *
                         ".jld2"), eco)
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

# Refuse an occurrence the storage has no slice for, rather than let the write fail on an index;
# `dim` is the storage's dimension of occurrences.
function _checkslot(storage::AbstractArray, count::Integer, dim::Integer = 3)
    count <= size(storage, dim) && return nothing
    return error("this run has reached occurrence $count, but the recording has room for " *
                 "$(size(storage, dim)). A run starting from zero has the starting state as its " *
                 "first occurrence, so size the storage for every multiple of `every` from zero " *
                 "to `duration`; a run continuing an earlier one has one fewer.")
end

# The same for each level a diversity recorder keeps, before any measure is taken.
_checkslot(::Nothing, ::Integer, ::Integer) = nothing

function _checkslots(recorder::RecordDiversity, count::Integer)
    _checkslot(recorder.subcommunity, count, 4)
    return _checkslot(recorder.metacommunity, count, 3)
end

# The measure types listed for one level, refusing anything that is not a concrete Diversity measure
# type - a measure function among them, which the storage form takes - and any listed twice.
function _measuretypes(given, level::Symbol)
    measures = Type[]
    for m in given
        if !(m isa Type && !isabstracttype(m) &&
             typeintersect(m, DiversityMeasure) !== Union{})
            hint = m isa Function ?
                   "; a function such as `norm_sub_alpha` is recorded by " *
                   "`RecordDiversity(storage, measure, qs)`" : ""
            error("`$level` lists Diversity's measure types - `NormalisedAlpha`, `RawBeta`, " *
                  "`Gamma` and the rest - and was given $(repr(m))$hint.")
        end
        m in measures && error("`$level` lists `$m` twice.")
        push!(measures, m)
    end
    return measures
end

# A level's measure types as a list of their names, for `show`.
function _measurenames(measures::Vector{Type})
    return "[" * join(nameof.(measures), ", ") * "]"
end

# Each measure recorded, once whatever its levels, in the order first listed.
function _recordedmeasures(recorder::RecordDiversity)
    return unique(vcat(recorder.measures.subcommunity,
                       recorder.measures.metacommunity))
end

# The levels one measure is recorded at, subcommunity first, as `diversity` takes them.
function _recordedlevels(recorder::RecordDiversity, measure::Type)
    levels = Diversity.DiversityLevel[]
    measure in recorder.measures.subcommunity &&
        push!(levels, subcommunityDiversity)
    measure in recorder.measures.metacommunity &&
        push!(levels, metacommunityDiversity)
    return levels
end

# The storage form: the measure function's values, cell by cell, into the caller's array.
function _recorddiversity!(recorder::RecordDiversity, measure, eco::Ecosystem,
                           count::Integer)
    _checkslot(recorder.storage, count)
    return _writediversity!(recorder, measure(eco, recorder.qs), count)
end

# The level form: every measure from one snapshot - the scratch filled once with the abundances'
# shares and wrapped as a `Metacommunity`, which caches the ordinariness every measure is built
# from - and one `diversity` call per measure, which computes its subcommunity values once for both
# levels.
# The snapshot aliases the scratch, so nothing between here and the last write may run the
# ecosystem's own Diversity hooks, which refill it.
function _recorddiversity!(recorder::RecordDiversity, ::Nothing, eco::Ecosystem,
                           count::Integer)
    _checkslots(recorder, count)
    snapshot = Metacommunity(_relativeabundance!(eco), eco.spplist, eco.habitat)
    for measure in _recordedmeasures(recorder)
        levels = _recordedlevels(recorder, measure)
        columns = diversity(Tables.columntable, levels, (measure,), snapshot,
                            recorder.qs)
        _writelevels!(recorder, measure, levels, columns, count)
    end
    return recorder
end

# Write one measure's rows into each level's slice, walking the order `diversity` gives them in - by
# level, then order, then cell - and checking each block's level and order as it starts, so that a
# change in that order fails rather than writing values into the wrong slots.
function _writelevels!(recorder::RecordDiversity, measure::Type, levels,
                       columns::NamedTuple, count::Integer)
    row = 0
    for level in levels
        if level == subcommunityDiversity
            i = findfirst(==(measure), recorder.measures.subcommunity)
            slots = view(recorder.subcommunity, i, :, :, count)
        else
            i = findfirst(==(measure), recorder.measures.metacommunity)
            slots = reshape(view(recorder.metacommunity, i, :, count), 1, :)
        end
        for (j, q) in enumerate(recorder.qs)
            _checkblock(columns, row + 1, level, q, measure)
            for c in axes(slots, 1)
                row += 1
                slots[c, j] = columns.diversity[row]
            end
        end
    end
    row == length(columns.diversity) ||
        error("`diversity` gave $(length(columns.diversity)) rows for `$measure`, where " *
              "$row were expected.")
    return recorder
end

# Refuse a block of `diversity`'s rows that does not start at the level and order expected.
function _checkblock(columns::NamedTuple, row::Integer, level, q, measure::Type)
    expected = level == subcommunityDiversity ? "subcommunity" : "metacommunity"
    (row <= length(columns.diversity) &&
     columns.partition_level[row] == expected && columns.q[row] == q) &&
        return nothing
    return error("`diversity`'s rows for `$measure` are not in the order recorded: row $row was " *
                 "expected to start the $expected values at order $q.")
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
                 "order, such as `norm_sub_alpha`; `$measure` gives $given diversity. " *
                 "`RecordDiversity(eco, ntimes; qs, metacommunity = [...])` records " *
                 "Diversity's measure types for the whole grid.")
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
