# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Filters for a table of records, such as occurrences, deciding which rows a run can use before
# anything is placed from them. A filter reports what it removes, and a placement matches its table's
# species names by the same matcher, so the two report alike.

using CSV: CSV
using DataFrames: AbstractDataFrame, DataFrame, keepat!
using Tables: Tables
using Unitful
using Unitful.DefaultSymbols

"""
    AbstractRecordFilter

A filter keeping the rows of a table of records that a run can use, for `filter!` (in place, on a
`DataFrame`) and `filter` (returning a new table): [`IsActiveSpecies`](@ref) by species,
[`InActiveCell`](@ref) by position. Each reports what it removes. Under MPI
each is a collective call: the first rank decides and reports, and every rank keeps the rows of its
own table by that decision, so every rank filters the same table, in the same order.
"""
abstract type AbstractRecordFilter end

"""
    IsActiveSpecies(species; column = :species, level = ProblemNotice(), verbosity = :normal)

A filter keeping the rows of a table that name one of the active species of `species`, a
[`SpeciesPool`](@ref), a [`SpeciesList`](@ref) or an ecosystem, whose species list it reads:
`filter!(IsActiveSpecies(pool), records)` removes the other rows of a `DataFrame` in place, and
`filter(IsActiveSpecies(pool), records)` returns the rows kept as a new table - a `DataFrame` for a
`DataFrame`, a `CSV.File` or any other column table, a named tuple of columns for one. A row's species is found as a placement finds it: by a species' own name or
a copy's label, by a name the pool's data used, by the name authority's identifier, or by any name the
authority resolves to one of the species. What the filter removes is reported at `level`, each kind
with how many rows it took and its names, the commonest first: names the authority does not recognise
(with no authority, names of no species), names it recognises that are of no species held, species
that are not active, and rows naming no species. Under MPI it is a collective call: the first rank
matches the names and reports, reading the authority alone, and every rank keeps the rows of its own
table by those answers, so every rank filters the same table, in the same order.

# Arguments

  - `species`: the pool, the species list or the ecosystem whose active species are kept.
  - `column`: the column naming each row's species.
  - `level`: how the report is raised - a [`ProblemNotice`](@ref), the default, not at all at
    `verbosity = :silent`; a [`ProblemWarning`](@ref); or a [`ProblemError`](@ref), which stops the
    filter before the table is changed.
  - `verbosity`: `:silent`, `:normal` or `:verbose`; `:verbose` lists every name, the others the
    first ten.

# Fields

  - `species`, `column`, `level`, `verbosity`: as given.
"""
struct IsActiveSpecies{S <: Union{SpeciesPool, SpeciesList, AbstractEcosystem},
                       L <: AbstractProblemSeverity} <: AbstractRecordFilter
    species::S
    column::Symbol
    level::L
    verbosity::Symbol

    function IsActiveSpecies(species::Union{SpeciesPool, SpeciesList,
                                            AbstractEcosystem};
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

"""
    InActiveCell(cells; latitude = :decimalLatitude, longitude = :decimalLongitude, snap = 0.0km,
                 level = ProblemNotice(), verbosity = :normal)

A filter keeping the rows of a table whose position lies in one of the active cells of `cells`, for
`filter!` and `filter` as [`IsActiveSpecies`](@ref) is: of a [`StudyArea`](@ref) or its
[`StudyAreaReport`](@ref), the cells the area was decided to simulate, known before anything is built
on it; of a [`GridHabitat`](@ref), or of an ecosystem, its habitat's, its cells as they stand when the
filter is called, which may be fewer, since a habitat can be narrower than the area it was built on
and cells can be deactivated as a run proceeds. A position is a latitude and a longitude in degrees, placed in the grid as
[`getcellat`](@ref) places it. A row in an inactive cell is kept where the centre of an active cell
lies within `snap` of its own cell's centre - a ground distance, straight on a projected grid and
along a great circle on a geographic one - since a placement given the same `snap` moves it there; the
filter itself moves nothing. What it removes is reported at `level`: rows off the grid, rows in
inactive cells with no active cell that near, and rows with no position; the rows kept by `snap` are
noted.

# Arguments

  - `cells`: the study area, its report, the habitat or the ecosystem whose active cells are kept.
  - `latitude`, `longitude`: the columns giving each row's position, in degrees; Darwin Core's
    names by default, as an occurrence download from GBIF has them.
  - `snap`: how far from an active cell a row in an inactive one may be and still be kept, as a
    length; `0.0km`, the default, keeps no row in an inactive cell.
  - `level`, `verbosity`: as for [`IsActiveSpecies`](@ref).

# Fields

  - `cells`, `latitude`, `longitude`, `snap`, `level`, `verbosity`: as given.
"""
struct InActiveCell{C <: Union{StudyArea, StudyAreaReport, GridHabitat,
                          AbstractEcosystem},
                    S <: Unitful.Length, L <: AbstractProblemSeverity} <:
       AbstractRecordFilter
    cells::C
    latitude::Symbol
    longitude::Symbol
    snap::S
    level::L
    verbosity::Symbol

    function InActiveCell(cells::Union{StudyArea, StudyAreaReport, GridHabitat,
                                       AbstractEcosystem};
                          latitude::Symbol = :decimalLatitude,
                          longitude::Symbol = :decimalLongitude,
                          snap::Unitful.Length = 0.0km,
                          level::AbstractProblemSeverity = ProblemNotice(),
                          verbosity::Symbol = :normal)
        _checksnap(snap)
        _checkverbosity(verbosity)
        return new{typeof(cells), typeof(snap), typeof(level)}(cells, latitude,
                                                               longitude, snap,
                                                               level, verbosity)
    end
end

function Base.show(io::IO, f::InActiveCell)
    print(io, "InActiveCell(", f.cells)
    f.latitude === :decimalLatitude ||
        print(io, ", latitude = ", repr(f.latitude))
    f.longitude === :decimalLongitude ||
        print(io, ", longitude = ", repr(f.longitude))
    iszero(f.snap) || print(io, ", snap = ", f.snap)
    f.level isa ProblemNotice ||
        print(io, ", level = ", nameof(typeof(f.level)), "()")
    f.verbosity === :normal || print(io, ", verbosity = ", repr(f.verbosity))
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", f::InActiveCell)
    println(io, "InActiveCell")
    println(io, "  cells     ", f.cells)
    println(io, "  position  ", repr(f.latitude), ", ", repr(f.longitude))
    println(io, "  snap      ", f.snap)
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
# leaves the table as it was. Under MPI the first rank looks the names up and reports - so it alone
# reads a name authority - and every rank keeps its rows by the first rank's answers, so every rank
# must filter the same table, in the same order.
function _keptrows(f::IsActiveSpecies, table)
    columns = Tables.columns(table)
    names = Tables.columnnames(columns)
    f.column in names ||
        error("the table names no species in a `$(f.column)` column; it has " *
              "$(_columnsphrase(names)). Give the column naming them as `column`.")
    index = _matchindex(f.species)
    given = Tables.getcolumn(columns, f.column)
    weights = Iterators.repeated(1, length(given))
    answers = _madeonroot() do
        answers = _nameanswers(index, given)
        matched = _matchspecies(given, weights, index, answers)
        _reportproblems(_matchproblems(matched.dropped, f.level, "record",
                                       "removed", index, f.verbosity),
                        f.verbosity)
        return answers
    end
    return .!isnothing.(_matchspecies(given, weights, index, answers).rows)
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

# Each distinct name of `names` with its species row and whether it is recognised, each looked up once.
function _nameanswers(index, names)
    answers = Dict{String,
                   NamedTuple{(:row, :recognised),
                              Tuple{Union{Nothing, Int}, Bool}}}()
    for (r, name) in enumerate(names)
        ismissing(name) && continue
        given = _rowspecies(name, r)
        haskey(answers, given) || (answers[given] = _matchrow(index, given))
    end
    return answers
end

# How each row's species name matches the species `index` describes, by `answers` for its names: the
# row of the species it names, `nothing` where it reaches no active species, and the weight each
# dropped name took (a record's 1, a placement's individuals) by kind - names not recognised (by the
# authority, or with none by any species), names recognised but of no species held, species not
# active, and rows naming none.
function _matchspecies(names, weights, index, answers::AbstractDict)
    rows = Vector{Union{Nothing, Int}}(nothing, length(names))
    dropped = (unrecognised = Dict{String, Int}(), absent = Dict{String, Int}(),
               inactive = Dict{String, Int}(), unnamed = Ref(0))
    for (r, (name, w)) in enumerate(zip(names, weights))
        if ismissing(name)
            dropped.unnamed[] += w
            continue
        end
        given = _rowspecies(name, r)
        hit = get(answers, given, nothing)
        isnothing(hit) &&
            error("row $r names `$given`, which the first MPI rank did not match: under MPI the " *
                  "first rank matches the names for every rank, so every rank filters the same " *
                  "table.")
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

# Which rows of `table` a filter by position keeps, its report raised first, so that a filter stopped
# by it leaves the table as it was. Under MPI the first rank decides and reports, and every rank keeps
# its rows by the first rank's decision.
function _keptrows(f::InActiveCell, table)
    columns = Tables.columns(table)
    names = Tables.columnnames(columns)
    for column in (f.latitude, f.longitude)
        column in names ||
            error("the table gives no position in a `$column` column; it has " *
                  "$(_columnsphrase(names)). Give the columns holding each row's position as " *
                  "`latitude` and `longitude`.")
    end
    latitudes = Tables.getcolumn(columns, f.latitude)
    longitudes = Tables.getcolumn(columns, f.longitude)
    keep = _madeonroot() do
        placed = _snapcells(f.cells, latitudes, longitudes, f.snap)
        _reportproblems(_placingproblems(placed, f.level, f.snap, "record",
                                         "removed", "kept"), f.verbosity)
        return .!isnothing.(placed.cells)
    end
    length(keep) == length(latitudes) ||
        error("the table has $(length(latitudes)) rows and the first MPI rank's has " *
              "$(length(keep)): under MPI the first rank decides for every rank, so every rank " *
              "filters the same table.")
    return keep
end

# Refuse a `snap` no distance can be.
function _checksnap(snap::Unitful.Length)
    (isfinite(snap) && snap >= zero(snap)) ||
        error("`snap` is how far a position may be from an active cell, so a finite, non-negative " *
              "length; got $snap.")
    return nothing
end

# Each record's cell in the grid of `x` - a study area, its report or a habitat - with a record in an
# inactive cell moved to the nearest active cell within `snap`, or `nothing` for a record off the
# grid, with no active cell that near, or with no position; and how many records each of those took.
# A cell's nearest active cell is found once, however many records it holds.
function _snapcells(x, latitudes, longitudes, snap::Unitful.Length)
    rows = [r
            for (r, (lat, long)) in enumerate(zip(latitudes, longitudes))
            if !(ismissing(lat) || ismissing(long))]
    places = [LatLong(_rowdegrees(latitudes[r], r, "latitude") * °,
                      _rowdegrees(longitudes[r], r, "longitude") * °)
              for r in rows]
    found = getcellat(x, places)
    yx = _gridyx(x)
    active = _activemask(x)
    cells = Vector{Union{Nothing, CartesianIndex{2}}}(nothing,
                                                      length(latitudes))
    nearest = Dict{CartesianIndex{2}, Union{Nothing, CartesianIndex{2}}}()
    offgrid = inactive = snapped = 0
    for (r, cell) in zip(rows, found)
        if isnothing(cell)
            offgrid += 1
        elseif active[cell]
            cells[r] = cell
        else
            target = get!(nearest, cell) do
                return _nearestactive(yx, active, cell, snap)
            end
            if isnothing(target)
                inactive += 1
            else
                snapped += 1
                cells[r] = target
            end
        end
    end
    return (cells = cells, offgrid = offgrid, inactive = inactive,
            snapped = snapped, unplaced = length(latitudes) - length(rows))
end

# A row's latitude or longitude in degrees, refused where it is not a number.
function _rowdegrees(x::Union{Real, Unitful.Quantity}, ::Int, ::AbstractString)
    return _degrees(x)
end

function _rowdegrees(x, row::Int, what::AbstractString)
    return error("row $row of a table gives its $what as `$(repr(x))`; a position is a latitude " *
                 "and a longitude in degrees.")
end

# The active cell nearest to `cell` whose centre lies within `snap` of its own, ties to the lower
# `(y, x)`; `nothing` where none does. Along a row the distance grows with the columns between, so
# each row is searched outwards until it passes `snap`, and only the rows whose own distance from
# `cell`'s is within it.
function _nearestactive(yx, active::AbstractMatrix{Bool},
                        cell::CartesianIndex{2},
                        snap::Unitful.Length)
    iszero(snap) && return nothing
    best = nothing
    nearest = snap
    for row in axes(active, 1)
        _rowdistance(yx, cell[1], row) <= snap || continue
        for (start, step) in ((cell[2], -1), (cell[2] + 1, 1))
            col = start
            while col in axes(active, 2)
                candidate = CartesianIndex(row, col)
                d = _celldistance(yx, cell, candidate)
                d <= snap || break
                if active[candidate] &&
                   (isnothing(best) || d < nearest ||
                    (d == nearest && Tuple(candidate) < Tuple(best)))
                    best, nearest = candidate, d
                end
                col += step
            end
        end
    end
    return best
end

# The centre of cell `i` on axis `d`, from its own edges, so whichever point of a cell its coordinate
# labels.
function _cellcentre(d, i::Int)
    lo, hi = _celledges(d, i)
    return (lo + hi) / 2
end

# The ground distance between two cells' centres.
function _celldistance(yx, a::CartesianIndex{2}, b::CartesianIndex{2})
    return _grounddistance(_cellcentre(yx[1], a[1]), _cellcentre(yx[2], a[2]),
                           _cellcentre(yx[1], b[1]), _cellcentre(yx[2], b[2]))
end

# The least ground distance between a cell of row `r1` and one of row `r2`, which no cell of `r2`
# lies nearer than: their difference in latitude, or in northing.
function _rowdistance(yx, r1::Int, r2::Int)
    return _arc(abs(_cellcentre(yx[1], r1) - _cellcentre(yx[1], r2)))
end

# A distance on the ground: a length as it is, an angle along a great circle of the package's one
# Earth radius, `LONGITUDE_DEGREE_LENGTH` to a degree.
_arc(d::Unitful.Length) = d
function _arc(angle::Unitful.DimensionlessQuantity)
    return uconvert(km, angle * LONGITUDE_DEGREE_LENGTH)
end

# The ground distance between two points: straight on a projected grid, along a great circle on a
# geographic one.
function _grounddistance(y1::Unitful.Length, x1, y2, x2)
    return sqrt((y1 - y2)^2 + (x1 - x2)^2)
end

function _grounddistance(φ1::Unitful.DimensionlessQuantity, λ1, φ2, λ2)
    h = sin((φ2 - φ1) / 2)^2 + cos(φ1) * cos(φ2) * sin((λ2 - λ1) / 2)^2
    return _arc(2 * asin(sqrt(min(h, 1.0))) * u"rad")
end

# The report of what placing records on a grid dropped, at `level`: records off the grid, in inactive
# cells with no active cell within `snap`, and with no position, saying what became of them (`done`);
# and, as a notice, those an active cell within `snap` saved, saying what became of those (`snapped`:
# a filter keeps them, a placement moves them).
function _placingproblems(placed::NamedTuple, level::AbstractProblemSeverity,
                          snap::Unitful.Length, noun::AbstractString,
                          done::AbstractString, snapped::AbstractString)
    near = iszero(snap) ? "" : ", with no active cell within $snap,"
    problems = Problem[]
    for (n, code, where) in ((placed.offgrid, :off_grid, "off the grid"),
        (placed.inactive, :inactive_cell,
        _inactivecells(placed.inactive) * near),
        (placed.unplaced, :no_position, "with no position"))
        iszero(n) && continue
        push!(problems,
              Problem(level, code,
                      "$(_countword(n, noun)) $where $(_were(n)) $done."))
    end
    iszero(placed.snapped) ||
        push!(problems,
              Problem(ProblemNotice(), :snapped,
                      "$(_countword(placed.snapped, noun)) $(_inactivecells(placed.snapped)) " *
                      "within $snap of an active cell $(_were(placed.snapped)) $snapped."))
    return problems
end

# Where records in inactive cells are, agreeing with how many.
function _inactivecells(n::Integer)
    return n == 1 ? "in an inactive cell" : "in inactive cells"
end
