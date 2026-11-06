# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Filters for a table of records, such as occurrences, deciding which rows a run can use before
# anything is placed from them. A filter reports what it removes, and a placement matches its table's
# species names by the same matcher, so the two report alike.

using CSV: CSV
using DataFrames: AbstractDataFrame, DataFrame, keepat!
using Tables: Tables

"""
    AbstractRecordFilter

A filter keeping the rows of a table of records that a run can use, for `filter!` (in place, on a
`DataFrame`) and `filter` (returning a new table): [`IsActiveSpecies`](@ref). Each reports what it
removes.
"""
abstract type AbstractRecordFilter end

"""
    IsActiveSpecies(species; column = :species, level = ProblemNotice(), verbosity = :normal)

A filter keeping the rows of a table that name one of the active species of `species`, a
[`SpeciesPool`](@ref) or a [`SpeciesList`](@ref): `filter!(IsActiveSpecies(pool), records)` removes the
other rows of a `DataFrame` in place, and `filter(IsActiveSpecies(pool), records)` returns the rows kept
as a new table - a `DataFrame` for a `DataFrame`, a `CSV.File` or any other column table, a named
tuple of columns for one. A row's species is found as a placement finds it: by a species' own name or
a copy's label, by a name the pool's data used, by the name authority's identifier, or by any name the
authority resolves to one of the species. What the filter removes is reported at `level`, each kind
with how many rows it took and its names, the commonest first: names the authority does not recognise
(with no authority, names of no species), names it recognises that are of no species held, species
that are not active, and rows naming no species.

# Arguments

  - `species`: the pool or the species list whose active species are kept.
  - `column`: the column naming each row's species.
  - `level`: how the report is raised - a [`ProblemNotice`](@ref), the default, not at all at
    `verbosity = :silent`; a [`ProblemWarning`](@ref); or a [`ProblemError`](@ref), which stops the
    filter before the table is changed.
  - `verbosity`: `:silent`, `:normal` or `:verbose`; `:verbose` lists every name, the others the
    first ten.

# Fields

  - `species`, `column`, `level`, `verbosity`: as given.
"""
struct IsActiveSpecies{S <: Union{SpeciesPool, SpeciesList},
                       L <: AbstractProblemSeverity} <: AbstractRecordFilter
    species::S
    column::Symbol
    level::L
    verbosity::Symbol

    function IsActiveSpecies(species::Union{SpeciesPool, SpeciesList};
                             column::Symbol = :species,
                             level::AbstractProblemSeverity = ProblemNotice(),
                             verbosity::Symbol = :normal)
        _checkverbosity(verbosity)
        return new{typeof(species), typeof(level)}(species, column, level,
                                                   verbosity)
    end
end

function Base.show(io::IO, f::IsActiveSpecies)
    print(io, "IsActiveSpecies(", f.species)
    f.column === :species || print(io, ", column = ", repr(f.column))
    f.level isa ProblemNotice ||
        print(io, ", level = ", nameof(typeof(f.level)), "()")
    f.verbosity === :normal || print(io, ", verbosity = ", repr(f.verbosity))
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", f::IsActiveSpecies)
    println(io, "IsActiveSpecies")
    println(io, "  species   ", f.species)
    println(io, "  column    ", repr(f.column))
    println(io, "  level     ", nameof(typeof(f.level)), "()")
    print(io, "  verbosity ", repr(f.verbosity))
    return nothing
end

function Base.filter!(f::AbstractRecordFilter, table::AbstractDataFrame)
    return keepat!(table, _keptrows(f, table))
end

function Base.filter!(::AbstractRecordFilter, ::CSV.File)
    return error("a `CSV.File` cannot be changed in place; `filter` returns the rows it keeps as " *
                 "a new `DataFrame`.")
end

function Base.filter(f::AbstractRecordFilter, table::AbstractDataFrame)
    return table[_keptrows(f, table), :]
end

function Base.filter(f::AbstractRecordFilter, file::CSV.File)
    return DataFrame(file, copycols = false)[_keptrows(f, file), :]
end

function Base.filter(f::AbstractRecordFilter, table::Tables.AbstractColumns)
    return DataFrame(table, copycols = false)[_keptrows(f, table), :]
end

function Base.filter(f::AbstractRecordFilter, table::NamedTuple)
    keep = _keptrows(f, table)
    return map(column -> column[keep], table)
end

# Which rows of `table` the filter keeps, its report raised first, so that a filter stopped by it
# leaves the table as it was.
function _keptrows(f::IsActiveSpecies, table)
    columns = Tables.columns(table)
    names = Tables.columnnames(columns)
    f.column in names ||
        error("the table names no species in a `$(f.column)` column; it has " *
              "$(_columnsphrase(names)). Give the column naming them as `column`.")
    index = _matchindex(f.species)
    given = Tables.getcolumn(columns, f.column)
    matched = _matchspecies(given, Iterators.repeated(1, length(given)), index)
    _reportproblems(_matchproblems(matched.dropped, f.level, "record",
                                   "removed",
                                   index, f.verbosity), f.verbosity)
    return .!isnothing.(matched.rows)
end

# A table's column names, for a message.
function _columnsphrase(names)
    return isempty(names) ? "no columns" :
           join(("`$n`" for n in names), ", ")
end

# A row's species as a `String`, refused where it is not a name.
_rowspecies(name::AbstractString, ::Int) = String(name)

function _rowspecies(name, row::Int)
    return error("row $row of a table names its species as `$(repr(name))`; a species is named by " *
                 "its name or its name authority's identifier.")
end

# How each row's species name matches the species `index` describes: the row of the species it names,
# `nothing` where it reaches no active species, and the weight each dropped name took (a record's 1, a
# placement's individuals) by kind - names not recognised (by the authority, or with none by any
# species), names recognised but of no species held, species not active, and rows naming none. Each
# distinct name is looked up once.
function _matchspecies(names, weights, index)
    rows = Vector{Union{Nothing, Int}}(nothing, length(names))
    dropped = (unrecognised = Dict{String, Int}(), absent = Dict{String, Int}(),
               inactive = Dict{String, Int}(), unnamed = Ref(0))
    seen = Dict{String,
                NamedTuple{(:row, :recognised),
                           Tuple{Union{Nothing, Int}, Bool}}}()
    for (r, (name, w)) in enumerate(zip(names, weights))
        if ismissing(name)
            dropped.unnamed[] += w
            continue
        end
        given = _rowspecies(name, r)
        hit = get!(seen, given) do
            return _matchrow(index, given)
        end
        if isnothing(hit.row)
            tally = hit.recognised ? dropped.absent : dropped.unrecognised
            tally[given] = get(tally, given, 0) + w
        elseif !index.active[hit.row]
            species = index.speciesname[hit.row]
            dropped.inactive[species] = get(dropped.inactive, species, 0) + w
        else
            rows[r] = hit.row
        end
    end
    return (rows = rows, dropped = dropped)
end

# A name's species row for a filter or a placement, and whether the name is recognised: found
# directly (a species' own name, a name of the pool's data, an identifier), else through the
# authority, which recognises names of species not held as well as of those held.
function _matchrow(index, name::AbstractString)
    row = something(_directrow(index, name), get(index.id, name, nothing),
                    Some(nothing))
    isnothing(row) || return (row = row, recognised = true)
    isnothing(index.resolve) && return (row = nothing, recognised = false)
    resolved = index.resolve(name, nothing)
    return (row = _resolvedrow(index, resolved),
            recognised = !isnothing(resolved.name))
end

# The report of what matching names dropped, at `level`: each kind with the `noun`s it took ("record",
# "individual") and what became of them (`done`), and its names, the heaviest first.
function _matchproblems(dropped, level::AbstractProblemSeverity,
                        noun::AbstractString, done::AbstractString, index,
                        verbosity::Symbol)
    authority = index.authority
    unrecognised = isnothing(authority) ? "of no species in $(index.holder)" :
                   "not recognised by $(_authorityphrase(authority))"
    absent = isnothing(authority) ? "" :
             "recognised by $(_authorityphrase(authority)) but of no species in " *
             index.holder
    activate = index.holder == "the pool" ?
               "`activate_species!` activates them." :
               "Activate them in their pool with `activate_species!` and build the species " *
               "again to keep them."
    problems = Problem[]
    for (tally, code, kind, what, after) in ((dropped.unrecognised, :not_recognised, "name",
        unrecognised, ""),
        (dropped.absent, :not_in_pool, "name", absent, ""),
        (dropped.inactive, :inactive_species, "species",
        "not active", " " * activate))
        isempty(tally) && continue
        total = sum(values(tally))
        k = length(tally)
        heaviest = sort!(collect(tally), by = p -> (-last(p), first(p)))
        push!(problems,
              Problem(level, code,
                      "$(_countword(total, noun)) of $k " *
                      (kind == "name" && k != 1 ? "names" : kind) * " $what " *
                      "$(_were(total)) $done: " *
                      _namelist(["`$n` ($w)" for (n, w) in heaviest],
                                verbosity) *
                      "." * after))
    end
    iszero(dropped.unnamed[]) ||
        push!(problems,
              Problem(level, :unnamed_species,
                      "$(_countword(dropped.unnamed[], noun)) naming no species " *
                      "$(_were(dropped.unnamed[])) $done."))
    return problems
end

# A count of a noun, as a phrase.
function _countword(n::Integer, noun::AbstractString)
    return n == 1 ? "1 $noun" : "$n $(noun)s"
end

# The verb agreeing with a count.
_were(n::Integer) = n == 1 ? "was" : "were"
