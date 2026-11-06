# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Where a run's individuals start. A placement says how many individuals of which species go into which
# cells: `build_ecosystem` calls it with the species, the habitat and a `PlacementContext`, checks the
# table it returns, and only then writes it into the abundances. The shipped placements weight the
# cells by supply, by each species' own suitability or by any layer, read a table, or place nothing;
# any function of the same three arguments is a placement too.

using Random
using Tables: Tables

"""
    AbstractPlacement

A rule saying how many individuals of each species start in each cell: [`BySupply`](@ref),
[`BySuitability`](@ref), [`ByLayer`](@ref), [`Weighted`](@ref), [`FromTable`](@ref),
[`FromOccurrences`](@ref), [`NoAbundance`](@ref), any of them filled to capacity by
[`ToCapacity`](@ref), and several added with `+` ([`SumOfPlacements`](@ref)). Called with
`(species, habitat, context)` - a [`SpeciesList`](@ref), a [`GridHabitat`](@ref) and a
[`PlacementContext`](@ref) - one returns a [`PlacedAbundances`](@ref). [`build_ecosystem`](@ref) takes
one as its first argument, as it takes `DefaultEcosystem()` or any function of the same three
arguments.
"""
abstract type AbstractPlacement end

# How many individuals a weighted placement places: a community `total`, split across the active
# species, `each` species' own, or neither, when it places the default.
struct _Total{T}
    value::T
end

struct _Each{T}
    value::T
end

struct _DefaultCount end

# How many individuals `FromOccurrences` places for each record it places.
struct _PerRecord
    value::Int
end

"""
    BySupply(axes...; total = nothing, each = nothing)

A placement weighting each cell by what the habitat supplies there: a species' individuals are
spread among the cells it may be placed in at random, each cell's chance in proportion to its
supply. With several supplies each is normalised over those cells and the results multiplied, so a
cell needs every resource to be worth placing in.

On its own it is a placement; inside [`Weighted`](@ref) it is one weight among several, and takes no
numbers there.

# Arguments

  - `axes`: the supplies weighting the cells, each by its axis (`SolarRadiation`) or, where two share
    an axis, by its own name (`:summer`); none for every supply.
  - `total`: the community's individuals, split across the active species at random: a whole
    number, or a density (`2.0/km^2`) of the area of the cells it may place in.
  - `each`: each active species' own individuals: a whole number or a density for all of them, a
    [`TraitColumn`](@ref) of the species' traits, a [`DrawFrom`](@ref) drawn from each species' own
    stream, or a function of a species' row of its traits returning one of those.

Give at most one of `total` and `each`. With neither it places 1000 individuals per active species,
as a `total`, and says so.

# Fields

  - `axes`: the supplies named, as given.
  - `numbers`: how many it places - the `total`, `each`, or the default.
"""
struct BySupply{A <: Tuple, N} <: AbstractPlacement
    axes::A
    numbers::N

    function BySupply(axes::Union{Type{<:NicheAxis}, Symbol}...;
                      total = nothing, each = nothing)
        numbers = _numbers(total, each)
        return new{typeof(axes), typeof(numbers)}(axes, numbers)
    end
end

"""
    BySuitability(axes...; total = nothing, each = nothing)

A placement weighting each cell by how well it suits each species: a species' individuals are spread
among the cells it may be placed in at random, each cell's chance in proportion to the species'
suitability there, scored by the run's nichefit as the dynamics score it. It is the one weight that
differs between species.

On its own it is a placement; inside [`Weighted`](@ref) it is one weight among several, and takes no
numbers there.

# Arguments

  - `axes`: the members of the regime scored, each by its axis (`Temperature`) or, where two share an
    axis, by its own name; none for the whole regime. Several are combined as the nichefit combines
    them.
  - `total`, `each`: how many individuals, as for [`BySupply`](@ref).

# Fields

  - `axes`: the members named, as given.
  - `numbers`: how many it places - the `total`, `each`, or the default.
"""
struct BySuitability{A <: Tuple, N} <: AbstractPlacement
    axes::A
    numbers::N

    function BySuitability(axes::Union{Type{<:NicheAxis}, Symbol}...;
                           total = nothing, each = nothing)
        numbers = _numbers(total, each)
        return new{typeof(axes), typeof(numbers)}(axes, numbers)
    end
end

"""
    ByLayer(layer; total = nothing, each = nothing, provenance = ProvenanceRecord[])

A placement weighting each cell by the values of a layer of your own - a vegetation index, a
density surface - which need be neither the habitat's regime nor its supply. A species'
individuals are spread among the cells it may be placed in at random, each cell's chance in
proportion to the layer's value there.

On its own it is a placement; inside [`Weighted`](@ref) it is one weight among several, and takes no
numbers there.

# Arguments

  - `layer`: a spec (any [`AbstractSpec`](@ref), built on the habitat's own grid when placing, as
    [`materialise`](@ref) builds one) or a layer already built on that grid. Its values must be finite
    and non-negative wherever the placement may place.
  - `total`, `each`: how many individuals, as for [`BySupply`](@ref).
  - `provenance`: [`ProvenanceRecord`](@ref)s of the data the layer was made from, where a spec does not
    record its own files.

# Fields

  - `layer`: the spec or layer, as given.
  - `numbers`: how many it places - the `total`, `each`, or the default.
  - `records`: the records given as `provenance`.
"""
struct ByLayer{L, N} <: AbstractPlacement
    layer::L
    numbers::N
    records::Vector{ProvenanceRecord}

    function ByLayer(layer::Union{AbstractSpec, AbstractLayer};
                     total = nothing, each = nothing,
                     provenance = ProvenanceRecord[])
        layer isa LayerCollection &&
            error("`ByLayer` weights by one layer, and was given a collection of " *
                  "$(length(layer)); give each its own `ByLayer` inside `Weighted`.")
        numbers = _numbers(total, each)
        return new{typeof(layer), typeof(numbers)}(layer, numbers,
                                                   _recordvector(provenance))
    end
end

"""
    Weighted(weights...; total = nothing, each = nothing)

A placement weighting each cell by several weights at once: each weight - a [`BySupply`](@ref),
[`BySuitability`](@ref) or [`ByLayer`](@ref) - is normalised over the cells a species may be placed
in, and the results multiplied, so a cell must score on every weight to be worth placing in. A
species' individuals are then spread among those cells at random in proportion to the product.

# Arguments

  - `weights`: two or more weights, none given numbers of its own.
  - `total`, `each`: how many individuals, as for [`BySupply`](@ref).

# Fields

  - `weights`: the weights, as given.
  - `numbers`: how many it places - the `total`, `each`, or the default.
"""
struct Weighted{W <: Tuple, N} <: AbstractPlacement
    weights::W
    numbers::N

    function Weighted(weights::AbstractPlacement...; total = nothing,
                      each = nothing)
        isempty(weights) &&
            error("`Weighted` needs the weights to multiply: `BySupply`, `BySuitability` or " *
                  "`ByLayer`.")
        for w in weights
            w isa _Weight ||
                error("`Weighted` multiplies weights - `BySupply`, `BySuitability` or `ByLayer` - " *
                      "and was given $(_placementphrase(w)). Add placements with `+` instead.")
            w.numbers isa _DefaultCount ||
                error("$(_placementphrase(w)) is a weight inside `Weighted`, so give its `total` " *
                      "or `each` to the `Weighted` instead.")
        end
        numbers = _numbers(total, each)
        return new{typeof(weights), typeof(numbers)}(weights, numbers)
    end
end

"""
    NoAbundance()

A placement placing nothing: the run starts with every species absent, to be added later.
"""
struct NoAbundance <: AbstractPlacement end

"""
    FromTable(data; provenance = ProvenanceRecord[])

A placement reading the individuals from a table, one row for a number of individuals of one
species in one cell.

A row names its species in a `species` column, by its name in the species list (or a name the data
of the pool it came from used, or a copy's label), by its name authority's identifier, or by any name
the authority resolves to one of the species, as [`IsActiveSpecies`](@ref) finds them. It names its
place by a location - `latitude` and `longitude` columns in degrees, or a `location` column of
[`LatLong`](@ref)s or [`SpatialLocation`](@ref)s - or, in a table in memory, by the grid's own key for
the cell, a `cell` column of what [`getcellat`](@ref) returns. An optional `count` column gives each
row's individuals (one each without it), and an optional `time` column an elapsed time or a date
(placed through the run's epoch and calendar); rows at build must be at or before its start. Rows of
one species and cell add.

# Arguments

  - `data`: a path to a CSV file, a [`ProvenanceRecord`](@ref) of one (the file it names is read, and
    the record kept), or any Tables.jl table in memory. A table that can be read only once - a
    `CSV.Rows` - places once, and [`repopulate!`](@ref) refuses it.
  - `provenance`: [`ProvenanceRecord`](@ref)s of the data the table was made from.

# Fields

  - `data`: the path, record or table, as given.
  - `records`: the records given as `provenance`.
"""
struct FromTable{D} <: AbstractPlacement
    data::D
    records::Vector{ProvenanceRecord}

    function FromTable(data; provenance = ProvenanceRecord[])
        _checktabledata(data, "FromTable")
        return new{typeof(data)}(data, _recordvector(provenance))
    end
end

"""
    FromOccurrences(data; species = :species, latitude = :decimalLatitude,
                    longitude = :decimalLongitude, date = nothing, time = nothing, snap = 0.0km,
                    total = nothing, per_record = nothing, level = ProblemWarning(),
                    provenance = ProvenanceRecord[])
    FromOccurrences(rule, data; ...)

A placement reading the individuals from a table of occurrence records, a row for each record of a
species at a place, as a download from GBIF has them. Each record's species is found as
[`IsActiveSpecies`](@ref) finds it and its place as [`InActiveCell`](@ref) places it, a record in an
inactive cell moved into the nearest active cell within `snap`; the records that cannot be placed are
reported at `level`. By default each record is one individual, so a cell holds as many of a species as
it has records of it, and the species are in the proportions their records give them. `per_record`
places that many individuals for each record instead; `total` spreads a number of individuals over all
the records in proportion to them, drawn from the run's community stream, so a species with no records
gets none. `rule`, a function usually written as a `do` block, makes any other table from the records.

Records are undated unless `date` or `time` names a column. At build every dated record must be at or
before the run's start; a record whose date cannot be read is dropped, counted in a notice.

# Arguments

  - `rule`: a function of `(rows, habitat, context)` returning what a placement returns - a table of
    species, place, count and, where dated, time, as [`FromTable`](@ref) reads one, or a
    [`PlacedAbundances`](@ref). `rows` is a `DataFrame` of the records placed: every column of the
    data, with each record's `species` as the species list names it, its `cell` (after `snap`) and,
    where dated, its `time` as elapsed time, in place of any columns of those names. Not with `total`
    or `per_record`.
  - `data`: a path to a CSV file, a [`ProvenanceRecord`](@ref) of one, or a Tables.jl table, as for
    [`FromTable`](@ref).
  - `species`, `latitude`, `longitude`: the columns naming each record's species and giving its
    position in degrees; Darwin Core's names by default.
  - `date`: the column of each record's date - dates, or Darwin Core's `eventDate` text, where a
    partial date or a range is taken at its start - or a tuple of columns, `(:year, :month)` or
    `(:year, :month, :day)`; placed through the run's epoch and calendar. Not with `time`.
  - `time`: the column of each record's elapsed time. Not with `date`.
  - `snap`: as for [`InActiveCell`](@ref).
  - `total`: a whole number of individuals for the community, spread over the records; not a
    density, since records have no area for one to multiply.
  - `per_record`: a whole number of individuals for each record. Not with `total`.
  - `level`: how the records that cannot be placed are reported, a [`ProblemWarning`](@ref) by
    default, as for [`IsActiveSpecies`](@ref).
  - `provenance`: [`ProvenanceRecord`](@ref)s of the data the table was made from.

# Fields

  - `rule`: the function given, or `nothing`.
  - `data`, `species`, `latitude`, `longitude`, `date`, `time`, `snap`, `level`: as given.
  - `numbers`: the `total` or `per_record` given, or neither.
  - `records`: the records given as `provenance`.
"""
struct FromOccurrences{R, D, T, S <: Unitful.Length, N,
                       L <: AbstractProblemSeverity} <: AbstractPlacement
    rule::R
    data::D
    species::Symbol
    latitude::Symbol
    longitude::Symbol
    date::T
    time::Union{Nothing, Symbol}
    snap::S
    numbers::N
    level::L
    records::Vector{ProvenanceRecord}

    function FromOccurrences(rule::Union{Nothing, Function}, data;
                             species::Symbol = :species,
                             latitude::Symbol = :decimalLatitude,
                             longitude::Symbol = :decimalLongitude,
                             date::Union{Nothing, Symbol,
                                         Tuple{Vararg{Symbol}}} = nothing,
                             time::Union{Nothing, Symbol} = nothing,
                             snap::Unitful.Length = 0.0km,
                             total = nothing, per_record = nothing,
                             level::AbstractProblemSeverity = ProblemWarning(),
                             provenance = ProvenanceRecord[])
        _checktabledata(data, "FromOccurrences")
        _checkdatecolumns(date, time)
        _checksnap(snap)
        isnothing(rule) || (isnothing(total) && isnothing(per_record)) ||
            error("`FromOccurrences` takes a rule making the table, or `total` or `per_record` " *
                  "giving its numbers, not both.")
        _refusedensity(total, "total")
        numbers = _occurrencenumbers(total, per_record)
        return new{typeof(rule), typeof(data), typeof(date), typeof(snap),
                   typeof(numbers), typeof(level)}(rule, data, species,
                                                   latitude, longitude, date,
                                                   time, snap, numbers, level,
                                                   _recordvector(provenance))
    end
end

FromOccurrences(data; kw...) = FromOccurrences(nothing, data; kw...)

"""
    ToCapacity(placement; fraction = 1)

A placement filling each cell to what its resources support, taking `placement`'s table as the shares
of the species there: a cell holds individuals until, in its scarcest resource, what they demand is
`fraction` of what the cell supplies - where, in the dynamics, births begin to be held back - each
species keeping its share, rounded to the nearest whole individual. A cell `placement` puts nothing
in stays empty. The shares may be whole counts or any non-negative weights, such as recency-weighted
records from a [`FromOccurrences`](@ref) rule; their places and species are found as any placement's
are, and a `time` column is ignored.

A filled cell is at its equilibrium only where suitability does not enter - `survival` zero, or a
suitability of one - and births and deaths have equal base rates: births then balance deaths exactly
when demand meets supply. Otherwise a cell settles above or below capacity, as its species' rates and
suitabilities set, and a run filled this way starts away from that.

# Arguments

  - `placement`: the placement whose table gives the shares - a shipped one, `DefaultEcosystem()` or
    a function of `(species, habitat, context)`.
  - `fraction`: the share of capacity each cell is filled to, a non-negative number; `1`, the
    default, fills it.

# Fields

  - `placement`, `fraction`: as given.
"""
struct ToCapacity{P} <: AbstractPlacement
    placement::P
    fraction::Float64

    function ToCapacity(placement; fraction::Real = 1)
        _checkplacement(placement)
        (isfinite(fraction) && fraction >= 0) ||
            error("`fraction` is the share of capacity each cell is filled to, a finite, " *
                  "non-negative number; got $fraction.")
        return new{typeof(placement)}(placement, Float64(fraction))
    end
end

"""
    SumOfPlacements(placements...)

Several placements added: each places its individuals as it would alone, and their rows add, a
species' individuals in a cell being the sum of every placement's. Written with `+`:
`BySupply(total = 5000) + FromTable("founders.csv")`.

# Fields

  - `placements`: the placements added, in order.
"""
struct SumOfPlacements{P <: Tuple} <: AbstractPlacement
    placements::P

    function SumOfPlacements(placements::AbstractPlacement...)
        length(placements) >= 2 ||
            error("a `SumOfPlacements` adds two or more placements; write them with `+`.")
        return new{typeof(placements)}(placements)
    end
end

function Base.:+(a::AbstractPlacement, b::AbstractPlacement)
    return SumOfPlacements(_summands(a)...,
                           _summands(b)...)
end

"""
    PlacementContext(; cells, streams, community, nichefit, epoch = nothing,
                     calendar = ExactDates(), time = 0.0s)

What a placement is given besides the species and the habitat: where it may place, the random
streams it draws from, and when it is.

# Arguments and fields

  - `cells`: the cells it may place in, as the grid's own keys - the active cells, at build.
  - `streams`: a random stream per species, by name, made from the run's seed and the species' own key
    in its pool, so a species' draws depend on no other species and leave the streams its dynamics
    draw from as they were.
  - `community`: the stream a community `total` is split across the species from.
  - `nichefit`: the run's nichefit, which scores a species' suitability for a cell.
  - `epoch`: the date the run starts on, or `nothing` for a run with no dates.
  - `calendar`: how dates become elapsed time.
  - `time`: the elapsed time the individuals are placed at, zero at build.
"""
struct PlacementContext{C <: AbstractVector, S <: AbstractDict,
                        R <: AbstractRNG, NF <: AbstractNicheFit}
    cells::C
    streams::S
    community::R
    nichefit::NF
    epoch::Union{Nothing, Dates.TimeType}
    calendar::AbstractRunCalendar
    time::typeof(1.0s)

    function PlacementContext(; cells::AbstractVector,
                              streams::AbstractDict, community::AbstractRNG,
                              nichefit::AbstractNicheFit,
                              epoch::Union{Nothing, Dates.TimeType} = nothing,
                              calendar::AbstractRunCalendar = ExactDates(),
                              time::Unitful.Time = 0.0s)
        return new{typeof(cells), typeof(streams), typeof(community),
                   typeof(nichefit)}(cells, streams, community, nichefit,
                                     epoch, calendar,
                                     float(uconvert(s, time)))
    end
end

function Base.show(io::IO, c::PlacementContext)
    return print(io, "PlacementContext(", length(c.cells), " cells, ",
                 length(c.streams), " species, t = ", c.time, ")")
end

"""
    PlacedAbundances(table; records = ProvenanceRecord[], problems = Problem[])

What a placement returns: the individuals it places, as a table, with the records of what it read and
anything worth saying about it. [`build_ecosystem`](@ref) checks the table - its species, its cells
and its counts - before writing it.

# Arguments and fields

  - `table`: any Tables.jl table, its rows as [`FromTable`](@ref) reads them: `species`, a place
    (`cell`, `location`, or `latitude` and `longitude`), and optionally `count` and `time`.
  - `records`: the [`ProvenanceRecord`](@ref)s of the data it read, kept in the ecosystem's provenance.
  - `problems`: [`Problem`](@ref)s to report, at the severity each carries.
"""
struct PlacedAbundances{T}
    table::T
    records::Vector{ProvenanceRecord}
    problems::Vector{Problem}

    function PlacedAbundances(table; records = ProvenanceRecord[],
                              problems = Problem[])
        Tables.istable(table) ||
            error("`PlacedAbundances` holds a table of the individuals placed - a named tuple of " *
                  "columns, a `DataFrame`, or any Tables.jl table; got a value of type $(typeof(table)).")
        return new{typeof(table)}(table, _recordvector(records),
                                  collect(Problem, problems))
    end
end

function Base.show(io::IO, p::PlacedAbundances)
    return print(io, "PlacedAbundances(", _rowsphrase(p.table), ", ",
                 length(p.records), " records, ", length(p.problems),
                 " problems)")
end

# --- Display -----------------------------------------------------------------
# Each placement prints as it would be written, so a run's provenance can say how it was placed.

function Base.show(io::IO, p::Union{BySupply, BySuitability})
    return print(io, nameof(typeof(p)), "(",
                 _argumentsphrase(map(_axisphrase, p.axes), p.numbers), ")")
end

function Base.show(io::IO, p::ByLayer)
    return print(io, "ByLayer(",
                 _argumentsphrase((sprint(show, p.layer),), p.numbers), ")")
end

function Base.show(io::IO, p::Weighted)
    return print(io, "Weighted(",
                 _argumentsphrase(map(w -> sprint(show, w), p.weights),
                                  p.numbers), ")")
end

# Written out, since a singleton's own display names its module wherever that is not in scope.
Base.show(io::IO, ::NoAbundance) = print(io, "NoAbundance()")

function Base.show(io::IO, p::FromTable)
    return print(io, "FromTable(", _dataphrase(p.data), ")")
end

function Base.show(io::IO, p::FromOccurrences)
    print(io, "FromOccurrences(")
    isnothing(p.rule) || print(io, _placementphrase(p.rule), ", ")
    print(io, _dataphrase(p.data))
    p.species === :species || print(io, ", species = ", repr(p.species))
    p.latitude === :decimalLatitude ||
        print(io, ", latitude = ", repr(p.latitude))
    p.longitude === :decimalLongitude ||
        print(io, ", longitude = ", repr(p.longitude))
    isnothing(p.date) || print(io, ", date = ", repr(p.date))
    isnothing(p.time) || print(io, ", time = ", repr(p.time))
    iszero(p.snap) || print(io, ", snap = ", p.snap)
    for (k, v) in pairs(_numberkeywords(p.numbers))
        print(io, ", ", k, " = ", v)
    end
    p.level isa ProblemWarning ||
        print(io, ", level = ", nameof(typeof(p.level)), "()")
    return print(io, ")")
end

function Base.show(io::IO, p::ToCapacity)
    print(io, "ToCapacity(", _placementphrase(p.placement))
    isone(p.fraction) || print(io, ", fraction = ", p.fraction)
    return print(io, ")")
end

function Base.show(io::IO, p::SumOfPlacements)
    return join(io, (sprint(show, q) for q in p.placements), " + ")
end

# == Functions ==================================================================================

# --- Calling a placement ----------------------------------------------------

# A placement's numbers as given: what it places, `_Total`, `_Each`, or the default.
function _numbers(::Nothing, ::Nothing)
    return _DefaultCount()
end

_numbers(total, ::Nothing) = _Total(_checkcount(total, "total"))

_numbers(::Nothing, each) = _Each(_checkeach(each))

function _numbers(total, each)
    return error("give a placement `total` or `each`, not both: `total` is the community's, " *
                 "split across the species, and `each` each species' own.")
end

# A number of individuals as given: a whole number, or a density per area.
function _checkcount(x::Real, what::AbstractString)
    (isinteger(x) && x >= 0) ||
        error("`$what` must be a whole number of individuals, or a density per area such as " *
              "`2.0/km^2`; got $x.")
    return x
end

function _checkcount(x::Unitful.AbstractQuantity, what::AbstractString)
    (dimension(x) == dimension(u"km^-2") && x >= zero(x)) ||
        error("`$what` must be a whole number of individuals, or a non-negative density per area " *
              "such as `2.0/km^2`; got $x.")
    return x
end

function _checkcount(x, what::AbstractString)
    return error("`$what` must be a whole number of individuals, or a density per area such as " *
                 "`2.0/km^2`; got a value of type $(typeof(x)).")
end

# Each species' own number as given: a number for all, or a form giving each its own.
_checkeach(x::Number) = _checkcount(x, "each")

_checkeach(x::Union{TraitColumn, DrawFrom, Function}) = x

function _checkeach(x)
    return error("`each` is a number of individuals or a density for every species, a " *
                 "`TraitColumn`, a `DrawFrom`, or a function of a species' row of its traits; got " *
                 "a $(typeof(x)).")
end

# A placement's numbers as the keywords that would give them, to rebuild it.
_numberkeywords(::_DefaultCount) = NamedTuple()

_numberkeywords(n::_Total) = (total = n.value,)

_numberkeywords(n::_Each) = (each = n.value,)

_numberkeywords(n::_PerRecord) = (per_record = n.value,)

# The three placements that are also weights, which `Weighted` multiplies.
const _Weight = Union{BySupply, BySuitability, ByLayer}

# The placements a sum adds: a sum's own, or the placement itself.
_summands(p::SumOfPlacements) = p.placements

_summands(p::AbstractPlacement) = (p,)

# What a placement's data is, as a phrase: a path or a record as itself, a table by its size.
_dataphrase(path::AbstractString) = repr(path)

_dataphrase(record::ProvenanceRecord) = sprint(show, record)

_dataphrase(table) = "a $(nameof(typeof(table))), " * _rowsphrase(table)

# How many rows a table has, where that is known without reading it.
function _rowsphrase(table)
    table isa AbstractVector && return _rowsword(length(table))
    Tables.columnaccess(table) || return "rows read once"
    return _rowsword(Tables.rowcount(Tables.columns(table)))
end

# A number of rows, as words.
_rowsword(n::Integer) = n == 1 ? "1 row" : "$n rows"

# A number of individuals, as words.
_individualsword(n::Integer) = n == 1 ? "1 individual" : "$n individuals"

# An axis or a member's name as it was written.
_axisphrase(axis::Type) = string(nameof(axis))

_axisphrase(name::Symbol) = repr(name)

# A placement's arguments as written: its positional ones, then its numbers.
function _argumentsphrase(positional, numbers)
    numberphrase = _numbersphrase(numbers)
    return join(isnothing(numberphrase) ? positional :
                (positional..., numberphrase), ", ")
end

_numbersphrase(::_DefaultCount) = nothing

_numbersphrase(n::_Total) = "total = $(n.value)"

_numbersphrase(n::_Each) = "each = $(n.value)"

# A placement as a reader would recognise it: as written, or a function by its name.
_placementphrase(p::AbstractPlacement) = sprint(show, p)

function _placementphrase(f::Function)
    name = string(nameof(f))
    return startswith(name, '#') ? "an anonymous function" : "`$name`"
end

_placementphrase(p) = "a $(nameof(typeof(p)))"

# Refuse what cannot be a placement before anything is built: a placement is a shipped one, the
# default, or anything callable.
function _checkplacement(placement)
    placement isa Union{AbstractPlacement, DefaultEcosystem, Function} &&
        return nothing
    isempty(methods(placement)) ||
        return nothing
    return error("a placement says how many individuals start where - `BySupply()`, " *
                 "`FromTable(...)`, `DefaultEcosystem()` or a function of `(species, habitat, " *
                 "context)`; got $(_placementphrase(placement)).")
end

# What a placement returns, as `PlacedAbundances`: its own, or the bare table it returned.
function _callplacement(placement, species, habitat, context)
    return _asplaced(placement(species, habitat, context), placement)
end

_asplaced(placed::PlacedAbundances, ::Any) = placed

function _asplaced(x, placement)
    Tables.istable(x) && return PlacedAbundances(x)
    return error("a placement returns a `PlacedAbundances`, or a table of the individuals it " *
                 "places; $(_placementphrase(placement)) returned a value of type $(typeof(x)).")
end

# A placement ready to call on every MPI rank's behalf: a layer weight's spec built on the habitat's
# grid, which under MPI every rank must do together, before the root alone calls it.
_prepared(p, ::GridHabitat) = p

function _prepared(p::ByLayer, habitat::GridHabitat)
    return ByLayer(materialise(p.layer, habitat.area);
                   _numberkeywords(p.numbers)...,
                   provenance = vcat(p.records, _layerrecords(p.layer)))
end

function _prepared(p::Weighted, habitat::GridHabitat)
    return Weighted(map(w -> _prepared(w, habitat), p.weights)...;
                    _numberkeywords(p.numbers)...)
end

function _prepared(p::ToCapacity, habitat::GridHabitat)
    return ToCapacity(_prepared(p.placement, habitat), fraction = p.fraction)
end

function _prepared(p::SumOfPlacements, habitat::GridHabitat)
    return SumOfPlacements(map(q -> _prepared(q, habitat), p.placements)...)
end

# The records of the files a layer spec read, under `:abundance`: a read's own, a combination's
# members', none for a synthetic layer.
_layerrecords(spec::RasterSpec) = _readinputs(spec, role = :abundance)

function _layerrecords(spec::ConstructedRasterSpec)
    return reduce(vcat, (_layerrecords(l) for l in spec.layers),
                  init = ProvenanceRecord[])
end

_layerrecords(spec::Varying) = _layerrecords(_unwrapspec(spec))

_layerrecords(::Any) = ProvenanceRecord[]

# Whether a placement can be called again - every one but a table that can be read only once.
_reusable(::Any) = true

_reusable(p::FromTable) = _rereadable(p.data)

_reusable(p::FromOccurrences) = _rereadable(p.data)

_reusable(p::ToCapacity) = _reusable(p.placement)

_reusable(p::SumOfPlacements) = all(_reusable, p.placements)

# Whether a table's data can be read again: a file, a vector of rows, or a table of columns.
_rereadable(::Union{AbstractString, ProvenanceRecord, AbstractVector}) = true

_rereadable(table) = Tables.columnaccess(table)

# A table's data checked as it is given to the placement `name`: a file that exists, a record naming
# one, or a table.
function _checktabledata(path::AbstractString, name::AbstractString)
    isfile(path) ||
        error("`$name` was given the path $(repr(path)), and there is no file there.")
    return nothing
end

function _checktabledata(record::ProvenanceRecord, name::AbstractString)
    isnothing(record.location) &&
        error("$record names no file: a record `Provenancer.record_provenance` returned, or " *
              "`provenance(path)` read, says where its file is, and `$name` takes one of those.")
    return nothing
end

function _checktabledata(table, name::AbstractString)
    Tables.istable(table) ||
        error("`$name` reads a path to a CSV file, a `ProvenanceRecord` of one, or a Tables.jl " *
              "table; got a value of type $(typeof(table)).")
    return nothing
end

# --- The shipped placements ---------------------------------------------------

function (p::BySupply)(species::SpeciesList, habitat::GridHabitat,
                       context::PlacementContext)
    return _placeweighted(p, (p,), species, habitat, context)
end

function (p::BySuitability)(species::SpeciesList, habitat::GridHabitat,
                            context::PlacementContext)
    return _placeweighted(p, (p,), species, habitat, context)
end

function (p::ByLayer)(species::SpeciesList, habitat::GridHabitat,
                      context::PlacementContext)
    built = _prepared(p, habitat)
    return _placeweighted(p, (built,), species, habitat, context,
                          built.records)
end

function (p::Weighted)(species::SpeciesList, habitat::GridHabitat,
                       context::PlacementContext)
    built = _prepared(p, habitat)
    records = reduce(vcat,
                     (w.records for w in built.weights if w isa ByLayer),
                     init = ProvenanceRecord[])
    return _placeweighted(p, built.weights, species, habitat, context,
                          records)
end

function (::NoAbundance)(::SpeciesList, ::GridHabitat, ::PlacementContext)
    return PlacedAbundances(_emptyplaced())
end

function (p::FromTable)(::SpeciesList, ::GridHabitat, ::PlacementContext)
    got = _readtable(p.data, :abundance)
    return PlacedAbundances(got.table, records = vcat(got.records, p.records))
end

function (p::FromOccurrences)(species::SpeciesList, habitat::GridHabitat,
                              context::PlacementContext)
    got = _readtable(p.data, :abundance, column = p.species)
    resolved = _occurrencerows(p, got.table, species, habitat, context)
    records = vcat(got.records, p.records)
    if !isnothing(p.rule)
        placed = _asplaced(p.rule(resolved.rows, habitat, context), p.rule)
        return PlacedAbundances(placed.table,
                                records = vcat(records, placed.records),
                                problems = vcat(resolved.problems,
                                                placed.problems))
    end
    dated = !(isnothing(p.date) && isnothing(p.time))
    dated && _checkstarted(resolved.rows, context)
    table = _occurrencecounts(p.numbers, resolved.rows, species, context,
                              resolved.problems, dated)
    return PlacedAbundances(table, records = records,
                            problems = resolved.problems)
end

function (p::ToCapacity)(species::SpeciesList, habitat::GridHabitat,
                         context::PlacementContext)
    placed = _callplacement(p.placement, species, habitat, context)
    shares = _sharesof(placed.table, species, habitat, context)
    filled = _filltocapacity(shares, species, habitat, p.fraction)
    return PlacedAbundances(filled.table, records = placed.records,
                            problems = vcat(placed.problems, shares.problems,
                                            filled.problems))
end

function (p::SumOfPlacements)(species::SpeciesList, habitat::GridHabitat,
                              context::PlacementContext)
    parts = [_callplacement(q, species, habitat, context)
             for q in p.placements]
    tables = [_canonicaltable(part.table, habitat, context) for part in parts]
    return PlacedAbundances(_vcattables(tables, context),
                            records = reduce(vcat, (q.records for q in parts)),
                            problems = reduce(vcat,
                                              (q.problems for q in parts)))
end

# A table of no individuals, in the form the weighted placements return.
function _emptyplaced()
    return (species = String[], cell = CartesianIndex{2}[], count = Int[])
end

# --- Weights ----------------------------------------------------------------

# A weighted placement's table: each species' individuals - from `numbers` - among the cells it may
# be placed in, each cell's chance the product of `weights`.
function _placeweighted(p, weights::Tuple, species, habitat, context,
                        records = ProvenanceRecord[])
    problems = Problem[]
    allowed = _allowedmask(habitat, context.cells)
    counts = _placecounts(p.numbers, _placementphrase(p), species, habitat,
                          context, allowed, problems)
    weight = _combinedweight(map(w -> _weightof(w, species, habitat, context,
                                                allowed), weights))
    table = _drawweighted(weight, counts, species, habitat, context, problems)
    return PlacedAbundances(table, records = records, problems = problems)
end

# The cells a placement may place in, as a mask over the grid's cells in storage order.
function _allowedmask(habitat::GridHabitat, cells::AbstractVector)
    shape = getgridshape(habitat)
    mask = falses(prod(shape))
    linear = LinearIndices(shape)
    for c in cells
        mask[linear[c]] = true
    end
    return mask
end

# A weight's values on the cells, for species `sp`: a function returning them normalised over the
# allowed cells, or `nothing` where the weight is zero on every one. Shared weights are worked out
# once, and only suitability differs between species.
function _weightof(w::BySupply, species, habitat, context, allowed)
    names = _membernames(habitat.supply, w.axes, "supply")
    shared = _supplyweight(habitat.supply, names, allowed)
    return _ -> shared
end

function _weightof(w::BySuitability, species, habitat, context, allowed)
    names = _membernames(habitat.regime, w.axes, "regime")
    whole = isempty(w.axes)
    return sp -> _suitabilityweight(habitat.regime, species.tolerance,
                                    context.nichefit, names, whole, allowed,
                                    sp)
end

function _weightof(w::ByLayer, species, habitat, context, allowed)
    shared = _layerweight(w.layer, allowed)
    return _ -> shared
end

# Several weights as one: the product of their values, normalised again; one weight as itself.
_combinedweight(weights::Tuple{Any}) = only(weights)

function _combinedweight(weights::Tuple)
    return function (sp)
        parts = map(w -> w(sp), weights)
        any(isnothing, parts) && return nothing
        product = _fold(parts) do a, b
            return a .* b
        end
        return _normalised(product)
    end
end

# `values` divided by their total, or `nothing` where that is zero or not a number.
function _normalised(values::AbstractVector)
    total = sum_kbn(values)
    (iszero(total) || !isfinite(total)) && return nothing
    return values ./ total
end

# One supply's values on the allowed cells, normalised; several supplies' each normalised, multiplied
# and normalised again.
function _supplyweight(supply::AbstractLayer, ::AbstractVector, allowed)
    b = reshape(parent(ustrip.(_getsupply(supply))), length(allowed))
    b[.!allowed] .= 0.0
    return _normalised(b)
end

function _supplyweight(supply::LayerCollection, names::AbstractVector, allowed)
    members = Tuple(supply[n] for n in names)
    fractions = _zipmap(members) do member
        b = reshape(parent(copy(_getsupply(member))), length(allowed))
        b[.!allowed] .= zero(eltype(b))
        return b ./ sum_kbn(b)
    end
    product = _fold(fractions) do f1, f2
        return f1 .* f2
    end
    return _normalised(product)
end

# One species' suitability for each allowed cell, normalised: the whole regime's as the dynamics
# score it, or the named members', combined as the nichefit combines them.
function _suitabilityweight(regime, tolerance, nichefit, names, whole::Bool,
                            allowed, sp::Integer)
    s = zeros(Float64, length(allowed))
    for pos in eachindex(allowed)
        allowed[pos] || continue
        s[pos] = whole ? _suitability(regime, tolerance, nichefit, pos, sp) :
                 _membersuitability(regime, tolerance, nichefit, names, pos,
                                    sp)
    end
    return _normalised(s)
end

function _membersuitability(regime, tolerance, nichefit, names, pos, sp)
    results = map(n -> _suitability(regime[n], tolerance[n], nichefit[n], pos,
                                    sp), Tuple(names))
    return nichefitcombine(nichefit)(NamedTuple{Tuple(names)}(results))
end

# A layer's values on the allowed cells, normalised, each checked finite and non-negative.
function _layerweight(layer::AbstractLayer, allowed)
    cellvalues = reshape(parent(layer.matrix), length(allowed))
    w = zeros(Float64, length(allowed))
    for pos in eachindex(allowed)
        allowed[pos] || continue
        x = cellvalues[pos]
        (!ismissing(x) && isfinite(ustrip(x)) && ustrip(x) >= 0) ||
            error("`ByLayer` weights cells by its layer's values, which must be finite and " *
                  "non-negative where it may place; the cell at " *
                  "$(Tuple(CartesianIndices(size(layer.matrix))[pos])) holds $x.")
        w[pos] = Float64(ustrip(x))
    end
    return _normalised(w)
end

# The names of the members of `layer` that `axes` select, each by its axis or its own name; every
# member's where none is given.
function _membernames(layer, axes::Tuple, what::AbstractString)
    isempty(axes) && return collect(Symbol, keys(layer))
    names = Symbol[_membername(layer, a, what) for a in axes]
    allunique(names) ||
        error("the habitat's $what member `$(first(n for n in names if count(==(n), names) > 1))` " *
              "is named twice.")
    return names
end

function _membername(layer, name::Symbol, what::AbstractString)
    haskey(layer, name) && return name
    return _refusemember(layer, repr(name), what)
end

function _membername(layer, axis::Type, what::AbstractString)
    found = [n for (n, m) in pairs(layer) if axisof(m) === axis]
    length(found) == 1 && return only(found)
    isempty(found) && return _refusemember(layer, string(nameof(axis)), what)
    return error("the habitat's $what has $(length(found)) members on `$(nameof(axis))` " *
                 "($(join(("`$n`" for n in found), ", "))): name the one you mean.")
end

function _refusemember(layer, given::AbstractString, what::AbstractString)
    listed = join(("`$n` ($(nameof(axisof(m))))" for (n, m) in pairs(layer)),
                  ", ")
    return error("the habitat's $what has no member $given; it has $listed.")
end

# --- How many -----------------------------------------------------------------

# How many individuals each species places, in species order: the community's `total` split across
# the active species from the community stream, `each` active species' own, or the default total.
function _placecounts(numbers::_Total, ::AbstractString, species, habitat,
                      context, allowed, ::Vector{Problem})
    n = _countof(numbers.value, _areaof(habitat, allowed), "total")
    return _splittotal(n, species, context)
end

function _placecounts(numbers::_Each, ::AbstractString, species, habitat,
                      context, allowed, ::Vector{Problem})
    active = _activespecies(species)
    areaof = _areaof(habitat, allowed)
    counts = zeros(Int, length(species.names))
    for (sp, name) in enumerate(species.names)
        active[sp] || continue
        value = _eachvalue(numbers.value, species.traits[sp, :],
                           context.streams[name], name)
        counts[sp] = _countof(value, areaof, "each")
    end
    return counts
end

function _placecounts(::_DefaultCount, what::AbstractString, species, habitat,
                      context, allowed, problems::Vector{Problem})
    per = _default(Val(:individuals)).value
    k = count(_activespecies(species))
    push!(problems,
          Problem(ProblemNotice(), :default_numbers,
                  "$what was given neither `total` nor `each`, so it places $per individuals per " *
                  "active species, $(per * k) in all, split across them at random."))
    return _splittotal(per * k, species, context)
end

# A community total split across the active species, in species order, by one draw from the
# community stream; the species that are not active get none.
function _splittotal(n::Integer, species, context)
    active = _activespecies(species)
    counts = zeros(Int, length(species.names))
    k = count(active)
    k == 0 && return counts
    counts[active] = rand(context.community, Multinomial(n, k))
    return counts
end

# Which species are active, in species order.
_activespecies(species::SpeciesList) = Vector{Bool}(species.traits.active)

# The area of the allowed cells in a unit asked for, worked out once a unit.
function _areaof(habitat::GridHabitat, allowed)
    areas = Dict{Any, Any}()
    return u -> get!(areas, u) do
        cells = vec(collect(getcellareas(u, habitat)))
        return sum_kbn(cells[allowed])
    end
end

# A number of individuals: a whole number as it is, a density times the area of the allowed cells,
# rounded half to even.
function _countof(x::Real, ::Any, what::AbstractString)
    _checkcount(x, what)
    return Int(x)
end

function _countof(x::Unitful.AbstractQuantity, areaof, what::AbstractString)
    _checkcount(x, what)
    return round(Int, ustrip(NoUnits, x * areaof(unit(inv(x)))))
end

function _countof(x, ::Any, what::AbstractString)
    return _checkcount(x, what)
end

# One species' `each`, from its row of the traits: a number as it is, a column's value (filled where
# missing), a draw from its own stream, or a function's answer.
_eachvalue(x::Number, ::Any, ::AbstractRNG, ::AbstractString) = x

function _eachvalue(d::DrawFrom, ::Any, rng::AbstractRNG, ::AbstractString)
    return _drawn(d, rng)
end

function _eachvalue(f::Function, row, rng::AbstractRNG, name::AbstractString)
    return _eachvalue(f(row), row, rng, name)
end

function _eachvalue(c::TraitColumn, row, rng::AbstractRNG, name::AbstractString)
    hasproperty(row, c.name) ||
        error("`each = $c` reads the species' traits, which have no column `:$(c.name)`; they " *
              "have $(join(("`:$n`" for n in propertynames(row)), ", ")).")
    v = row[c.name]
    ismissing(v) && (v = _eachfill(c, c.fill, row, rng, name))
    return v isa Unitful.AbstractQuantity ? v : _withunit(v, c.unit)
end

function _eachvalue(x, ::Any, ::AbstractRNG, name::AbstractString)
    return error("`each` gave species `$name` a value of type $(typeof(x)); a function given as `each` returns " *
                 "a number of individuals or a density.")
end

# A missing value of a column `each` reads: its fill, drawn from the species' own stream where it is
# a distribution.
function _eachfill(c::TraitColumn, ::Nothing, row, rng, name)
    return error("species `$name` has no value in `:$(c.name)`, and `$c` gives no `fill`.")
end

_eachfill(::TraitColumn, fill::Number, row, rng, name) = fill

function _eachfill(::TraitColumn, fill::UnivariateDistribution, row, rng, name)
    return rand(rng, fill)
end

_eachfill(::TraitColumn, fill::Function, row, rng, name) = fill(row)

# Each species' individuals drawn among the cells, in proportion to `weight(sp)`, from its own
# stream, as rows of species, cell and count. A species with individuals to place and no weight on
# any allowed cell places none, and is named in a warning.
function _drawweighted(weight, counts::AbstractVector{<:Integer}, species,
                       habitat, context, problems::Vector{Problem})
    cellkeys = CartesianIndices(getgridshape(habitat))
    placed = _emptyplaced()
    drawn = zeros(Int, length(cellkeys))
    lost = Pair{String, Int}[]
    for (sp, name) in enumerate(species.names)
        n = counts[sp]
        n == 0 && continue
        probabilities = weight(sp)
        if isnothing(probabilities)
            push!(lost, name => n)
            continue
        end
        rand!(context.streams[name], Multinomial(n, probabilities), drawn)
        for pos in eachindex(drawn)
            iszero(drawn[pos]) && continue
            push!(placed.species, name)
            push!(placed.cell, cellkeys[pos])
            push!(placed.count, drawn[pos])
        end
    end
    isempty(lost) ||
        push!(problems,
              Problem(ProblemWarning(), :nowhere_to_place,
                      "$(_individualsword(sum(last, lost))) were not placed, their species having no " *
                      "weight on any cell they may be placed in: $(_listphrase(lost))."))
    return placed
end

# Species with their individuals, as a phrase, the first ten listed.
function _listphrase(pairs::AbstractVector{<:Pair})
    shown = join(("`$(first(p))` ($(last(p)))" for p in first(pairs, 10)), ", ")
    more = length(pairs) - 10
    return shown * (more > 0 ? " and $more more species" : "")
end

# --- Placing from occurrence records ----------------------------------------------

# Refuse a record's date and its time both, or a `date` tuple that is not a year with a month and a
# day or fewer.
_checkdatecolumns(::Nothing, ::Nothing) = nothing

_checkdatecolumns(::Nothing, ::Symbol) = nothing

_checkdatecolumns(::Symbol, ::Nothing) = nothing

function _checkdatecolumns(date::Tuple, ::Nothing)
    1 <= length(date) <= 3 ||
        error("`date` as a tuple names a record's year, month and day columns, as `(:year, :month)` " *
              "or `(:year, :month, :day)`; got $(repr(date)).")
    return nothing
end

function _checkdatecolumns(::Union{Symbol, Tuple}, ::Symbol)
    return error("`FromOccurrences` dates its records by `date` or by `time`, not both.")
end

# Refuse a density as `total`: a record has no area for one to multiply.
_refusedensity(::Any, ::AbstractString) = nothing

function _refusedensity(x::Unitful.AbstractQuantity, what::AbstractString)
    return error("`FromOccurrences` spreads whole numbers of individuals over its records, which have " *
                 "no area for a density to multiply; give `$what` as a number of individuals, not " *
                 "$x.")
end

# `FromOccurrences`' numbers as given: a community `total`, a number `per_record`, or the default, one.
_occurrencenumbers(::Nothing, ::Nothing) = _DefaultCount()

_occurrencenumbers(total, ::Nothing) = _Total(_checkcount(total, "total"))

function _occurrencenumbers(::Nothing, n)
    (n isa Real && isinteger(n) && n >= 0) ||
        error("`per_record` is a whole number of individuals for each record; got $(repr(n)).")
    return _PerRecord(Int(n))
end

function _occurrencenumbers(total, n)
    return error("give `FromOccurrences` a `total` or a number `per_record`, not both: `total` is the " *
                 "community's, spread over the records, and `per_record` each record's.")
end

# The records of `data` that can be placed, as a `DataFrame` of every column of the data with each
# record's `species` as the species list names it, its `cell` after snapping and, where dated, its
# `time` as elapsed time; and the problems of those that cannot be placed, at the placement's level,
# those whose date cannot be read in a notice.
function _occurrencerows(p::FromOccurrences, data::DataFrame,
                         species::SpeciesList, habitat::GridHabitat,
                         context::PlacementContext)
    _checkoccurrencecolumns(p, data)
    index = _matchindex(species)
    given = data[!, p.species]
    matched = _matchspecies(given, Iterators.repeated(1, length(given)), index,
                            _nameanswers(index, given))
    problems = _matchproblems(matched.dropped, p.level, "record", "not placed",
                              index, :normal)
    named = findall(!isnothing, matched.rows)
    placed = _snapcells(habitat, data[named, p.latitude],
                        data[named, p.longitude], p.snap)
    append!(problems,
            _placingproblems(placed, p.level, p.snap, "record", "not placed",
                             "moved into it"))
    located = findall(!isnothing, placed.cells)
    keep = named[located]
    rows = data[keep, :]
    rows.species = species.names[Vector{Int}(matched.rows[keep])]
    rows.cell = Vector{CartesianIndex{2}}(placed.cells[located])
    times = _recordtimes(p.date, p.time, rows, context)
    isnothing(times) && return (rows = rows, problems = problems)
    dated = .!isnothing.(times)
    undated = count(!, dated)
    iszero(undated) ||
        push!(problems,
              Problem(ProblemNotice(), :undated,
                      "$(_countword(undated, "record")) with no date that could be read " *
                      "$(_were(undated)) not placed."))
    rows = rows[dated, :]
    rows.time = [t for t in times if !isnothing(t)]
    return (rows = rows, problems = problems)
end

# Refuse occurrence records lacking a column the placement reads, naming the columns they have.
function _checkoccurrencecolumns(p::FromOccurrences, data::DataFrame)
    needed = Symbol[p.species, p.latitude, p.longitude]
    isnothing(p.date) || append!(needed, p.date isa Symbol ? (p.date,) : p.date)
    isnothing(p.time) || push!(needed, p.time)
    absent = setdiff(needed, propertynames(data))
    isempty(absent) && return nothing
    return error("the occurrence records have no $(join(("`$c`" for c in absent), ", ")) " *
                 "column; they have $(_columnsphrase(propertynames(data))). Name the columns " *
                 "holding them as `species`, `latitude`, `longitude`, `date` or `time`.")
end

# Each record's time as elapsed time - from its date through the run's epoch and calendar, or as
# given - or `nothing` where it cannot be read; `nothing` for all when the records are undated.
_recordtimes(::Nothing, ::Nothing, ::DataFrame, ::PlacementContext) = nothing

function _recordtimes(date::Symbol, ::Nothing, rows::DataFrame,
                      context::PlacementContext)
    return _datedtimes([_recorddate(d) for d in rows[!, date]], context)
end

function _recordtimes(date::Tuple, ::Nothing, rows::DataFrame,
                      context::PlacementContext)
    parts = [rows[!, c] for c in date]
    return _datedtimes([_partsdate(v...) for v in zip(parts...)], context)
end

function _recordtimes(::Nothing, time::Symbol, rows::DataFrame,
                      context::PlacementContext)
    return [ismissing(t) ? nothing : _placementtime(t, context, r)
            for (r, t) in enumerate(rows[!, time])]
end

# Dates placed through the run's epoch and calendar, each distinct date once, `nothing` staying so.
function _datedtimes(dates::AbstractVector, context::PlacementContext)
    times = Vector{Union{Nothing, typeof(1.0s)}}(nothing, length(dates))
    placed = Dict{Dates.Date, typeof(1.0s)}()
    for (r, date) in enumerate(dates)
        isnothing(date) && continue
        times[r] = get!(placed, date) do
            return _placementtime(date, context, r)
        end
    end
    return times
end

# A record's date, or `nothing` where it has none that can be read: a date as it is, or Darwin
# Core's `eventDate` text, a partial date (`2001`, `2001-03`) or a range (`2001-03/2001-05`) taken at
# its start.
_recorddate(::Missing) = nothing

_recorddate(date::Dates.TimeType) = Dates.Date(date)

function _recorddate(text::AbstractString)
    start = strip(first(split(first(split(text, '/')), 'T')))
    parts = count(==('-'), start)
    parts <= 2 || return nothing
    format = (Dates.dateformat"yyyy", Dates.dateformat"yyyy-mm",
              Dates.dateformat"yyyy-mm-dd")[parts + 1]
    return tryparse(Dates.Date, start, format)
end

function _recorddate(x)
    return error("a `date` column holds dates or dates written as text, and this one holds " *
                 "$(repr(x)); give a year's, a month's and a day's columns as a tuple, " *
                 "`date = (:year, :month, :day)`.")
end

# A record's date from its year, month and day, a month or a day left out taking the start of the
# year or the month; `nothing` where they make no date.
function _partsdate(year, month = missing, day = missing)
    y = _datepart(year)
    m = ismissing(month) ? 1 : _datepart(month)
    d = ismissing(day) ? 1 : _datepart(day)
    (isnothing(y) || isnothing(m) || isnothing(d)) && return nothing
    (1 <= m <= 12 && 1 <= d <= Dates.daysinmonth(y, m)) || return nothing
    return Dates.Date(y, m, d)
end

# One part of a date as a whole number, or `nothing` where it is not one.
_datepart(::Missing) = nothing

_datepart(x::Integer) = Int(x)

_datepart(x::Real) = isinteger(x) ? Int(x) : nothing

_datepart(x::AbstractString) = tryparse(Int, strip(x))

_datepart(::Any) = nothing

# Refuse records dated after a build's start, saying how to place them.
function _checkstarted(rows::DataFrame, context::PlacementContext)
    late = count(>(context.time), rows.time)
    iszero(late) && return nothing
    return error("$(_countword(late, "record")) $(late == 1 ? "is" : "are") dated after the run's " *
                 "start, and a placement at build places the starting state: give `date = " *
                 "nothing` to place every record at the start, or a `do` block that drops or " *
                 "redates them.")
end

# How many records each species has in each cell and, where dated, at each time: the groups' keys,
# sorted, their counts, and each species' range of them.
function _recordgroups(rows::DataFrame, dated::Bool)
    all = dated ? collect(zip(rows.species, rows.cell, rows.time)) :
          collect(zip(rows.species, rows.cell))
    tally = Dict{eltype(all), Int}()
    for k in all
        tally[k] = get(tally, k, 0) + 1
    end
    sorted = sort!(collect(keys(tally)))
    ranges = Dict{String, UnitRange{Int}}()
    for (i, k) in enumerate(sorted)
        r = get(ranges, first(k), i:(i - 1))
        ranges[first(k)] = first(r):i
    end
    return (keys = sorted, counts = [tally[k] for k in sorted],
            ranges = ranges)
end

# The groups' individuals as a placement's table, the empty ones left out.
function _grouptable(groups::NamedTuple, counts::AbstractVector{<:Integer},
                     dated::Bool)
    kept = findall(>(0), counts)
    ks = groups.keys[kept]
    table = (species = String[first(k) for k in ks],
             cell = CartesianIndex{2}[k[2] for k in ks], count = counts[kept])
    return dated ? merge(table, (time = [k[3] for k in ks],)) : table
end

# The individuals the records place: one for each record by default, or `per_record`; or the
# community's `total` spread over every record, from the community stream, so a species with more
# records gets more.
function _occurrencecounts(::_DefaultCount, rows::DataFrame,
                           ::SpeciesList, ::PlacementContext,
                           ::Vector{Problem}, dated::Bool)
    groups = _recordgroups(rows, dated)
    return _grouptable(groups, groups.counts, dated)
end

function _occurrencecounts(numbers::_PerRecord, rows::DataFrame,
                           ::SpeciesList, ::PlacementContext, ::Vector{Problem},
                           dated::Bool)
    groups = _recordgroups(rows, dated)
    return _grouptable(groups, groups.counts .* numbers.value, dated)
end

function _occurrencecounts(numbers::_Total, rows::DataFrame, ::SpeciesList,
                           context::PlacementContext, problems::Vector{Problem},
                           dated::Bool)
    groups = _recordgroups(rows, dated)
    n = Int(numbers.value)
    if isempty(groups.keys)
        n > 0 &&
            push!(problems,
                  Problem(ProblemWarning(), :nowhere_to_place,
                          "$(_individualsword(n)) were not placed: no record could be placed to " *
                          "place them by."))
        return _grouptable(groups, Int[], dated)
    end
    counts = groups.counts
    drawn = rand(context.community, Multinomial(n, counts ./ sum(counts)))
    return _grouptable(groups, drawn, dated)
end

# --- Filling to capacity -----------------------------------------------------------

# A table's shares resolved against the species and the grid: each row's species row, cell and share,
# its time ignored, and the rows that could not be resolved reported as warnings, counted in rows.
function _sharesof(table, species::SpeciesList, habitat::GridHabitat,
                   context::PlacementContext)
    canonical = _canonicaltable(table, habitat, context, times = false)
    shares = [_placedshare(x, row) for (row, x) in enumerate(canonical.count)]
    index = _matchindex(species)
    matched = _matchspecies(canonical.species,
                            Iterators.repeated(1, length(shares)), index,
                            _nameanswers(index, canonical.species))
    problems = _matchproblems(matched.dropped, ProblemWarning(), "row",
                              "not placed", index, :normal)
    allowed = _allowedmask(habitat, context.cells)
    linear = LinearIndices(getgridshape(habitat))
    rows = (species = Int[], cells = CartesianIndex{2}[], shares = Float64[])
    offgrid = outside = 0
    for (sp, cell, share) in zip(matched.rows, canonical.cell, shares)
        isnothing(sp) && continue
        if isnothing(cell)
            offgrid += 1
        elseif !allowed[linear[cell]]
            outside += 1
        else
            push!(rows.species, sp)
            push!(rows.cells, cell)
            push!(rows.shares, share)
        end
    end
    iszero(offgrid) ||
        push!(problems,
              Problem(ProblemWarning(), :off_grid,
                      "$(_countword(offgrid, "row")) naming a place off the grid " *
                      "$(_were(offgrid)) not placed."))
    iszero(outside) ||
        push!(problems,
              Problem(ProblemWarning(), :not_allowed,
                      "$(_countword(outside, "row")) naming a cell the placement may not place " *
                      "in $(_were(outside)) not placed."))
    return (rows = rows, problems = problems)
end

# A row's share as a number, refused where it is not a finite, non-negative one.
function _placedshare(x::Real, row::Int)
    (isfinite(x) && x >= 0) ||
        error("row $row of a placement's table gives a share of $x: a share is a finite, " *
              "non-negative number.")
    return Float64(x)
end

function _placedshare(x, row::Int)
    return error("row $row of a placement's table gives its count as `$(repr(x))`: give a number.")
end

# Each cell filled to `fraction` of what its scarcest resource supports at the share-weighted mean
# demand of the species in it, each species keeping its share, rounded to the nearest whole
# individual: a placement's table, and the cells whose species demand nothing, which no supply can
# bound, reported and left empty.
function _filltocapacity(shares::NamedTuple, species::SpeciesList,
                         habitat::GridHabitat, fraction::Float64)
    rows = shares.rows
    incell = Dict{CartesianIndex{2}, Dict{Int, Float64}}()
    for (sp, cell, share) in zip(rows.species, rows.cells, rows.shares)
        bycell = get!(Dict{Int, Float64}, incell, cell)
        bycell[sp] = get(bycell, sp, 0.0) + share
    end
    supplies = map(_getsupply, values(habitat.supply))
    demands = map(d -> d.resource, values(species.demand))
    table = (species = String[], cell = CartesianIndex{2}[], count = Int[])
    unbounded = 0
    for cell in sort!(collect(keys(incell)))
        held = sort!(collect(incell[cell]))
        total = sum_kbn(last.(held))
        total > 0 || continue
        capacity = _cellcapacity(supplies, demands, held, total, cell)
        if isinf(capacity)
            unbounded += 1
            continue
        end
        for (sp, share) in held
            n = round(Int, share * fraction * capacity / total)
            n > 0 || continue
            push!(table.species, species.names[sp])
            push!(table.cell, cell)
            push!(table.count, n)
        end
    end
    problems = Problem[]
    iszero(unbounded) ||
        push!(problems,
              Problem(ProblemWarning(), :unbounded,
                      "$(_countword(unbounded, "cell")) whose species demand nothing of any " *
                      "resource, which no supply can fill to, $(_were(unbounded)) left empty."))
    return (table = table, problems = problems)
end

# How many individuals of a cell's species, in their shares, its scarcest resource supports: the
# least, over the resources, of the cell's supply over the share-weighted mean demand, a resource the
# species demand nothing of not binding.
function _cellcapacity(supplies::Tuple, demands::Tuple,
                       held::AbstractVector{<:Pair}, total::Float64,
                       cell::CartesianIndex{2})
    capacity = Inf
    for (supply, demand) in zip(supplies, demands)
        mean = sum_kbn([share * demand[sp] for (sp, share) in held]) / total
        iszero(mean) && continue
        capacity = min(capacity, uconvert(NoUnits, supply[cell] / mean))
    end
    return capacity
end

# --- Checking a table ---------------------------------------------------------

# A placement's table in one form, whatever form it came in: `species` as given, `cell` as the
# grid's own key (`nothing` for a place off the grid), `count` (one each where the table has none)
# and, where it is dated and `times` asks for them, `time` as elapsed time.
function _canonicaltable(table, habitat::GridHabitat,
                         context::PlacementContext; times::Bool = true)
    columns = Tables.columns(table)
    names = Tables.columnnames(columns)
    :species in names ||
        error("a placement's table names each row's species in a `species` column; it has " *
              "$(_columnsphrase(names)).")
    species = collect(Tables.getcolumn(columns, :species))
    cells = _gridkeys(habitat, columns, names)
    counts = :count in names ? collect(Tables.getcolumn(columns, :count)) :
             fill(1, length(species))
    (times && :time in names) ||
        return (species = species, cell = cells, count = counts)
    elapsed = [_placementtime(t, context, row)
               for (row, t) in enumerate(Tables.getcolumn(columns, :time))]
    return (species = species, cell = cells, count = counts, time = elapsed)
end

# Each row's cell as the grid's own key, from whichever place the table gives.
function _gridkeys(habitat::GridHabitat, columns, names)
    :cell in names &&
        return _cellkeys(habitat, Tables.getcolumn(columns, :cell))
    :location in names &&
        return getcellat(habitat, collect(Tables.getcolumn(columns, :location)))
    if :latitude in names && :longitude in names
        places = [LatLong(_degrees(lat) * °, _degrees(long) * °)
                  for (lat, long) in zip(Tables.getcolumn(columns, :latitude),
                                         Tables.getcolumn(columns, :longitude))]
        return getcellat(habitat, places)
    end
    return error("a placement's table names each row's place by `latitude` and `longitude`, a " *
                 "`location`, or the grid's own key as a `cell`; it has " *
                 "$(_columnsphrase(names)).")
end

# A `cell` column's keys checked against the grid: a key off it is `nothing`, as a place off it is.
function _cellkeys(habitat::GridHabitat, given)
    shape = CartesianIndices(getgridshape(habitat))
    return Union{Nothing, CartesianIndex{2}}[_cellkey(shape, key, row)
                                             for (row, key) in enumerate(given)]
end

function _cellkey(shape, key::CartesianIndex{2}, ::Int)
    return checkbounds(Bool, shape, key) ? key : nothing
end

_cellkey(::Any, ::Nothing, ::Int) = nothing

function _cellkey(::Any, key, row::Int)
    return error("row $row of a placement's table gives its `cell` as `$(repr(key))`; a `cell` is " *
                 "the grid's own key, as `getcellat` returns it. Give a place as `latitude` and " *
                 "`longitude`, or a `location`.")
end

# A row's time as elapsed time: an elapsed time, or a date placed through the context's epoch and
# calendar.
function _placementtime(t::Unitful.Time, ::PlacementContext, ::Int)
    return float(uconvert(s, t))
end

function _placementtime(date::Dates.TimeType, context::PlacementContext,
                        row::Int)
    isnothing(context.epoch) &&
        error("row $row of a placement's table is dated, but this run has no epoch to place it " *
              "from. Pass `epoch = ` to `build_ecosystem`, or give the rows elapsed times.")
    return float(uconvert(s, _elapsedat(context.epoch, date, context.calendar)))
end

function _placementtime(t, ::PlacementContext, row::Int)
    return error("row $row of a placement's table gives its `time` as `$(repr(t))`: give an " *
                 "elapsed time or a date.")
end

# Several placements' tables, in the one form, as one: dated where any is, the rest at the context's
# time.
function _vcattables(tables::AbstractVector, context::PlacementContext)
    species = reduce(vcat, (t.species for t in tables))
    cells = reduce(vcat, (t.cell for t in tables))
    counts = reduce(vcat, (t.count for t in tables))
    any(t -> haskey(t, :time), tables) ||
        return (species = species, cell = cells, count = counts)
    times = reduce(vcat,
                   (haskey(t, :time) ? t.time :
                    fill(context.time, length(t.species))
                    for t in tables))
    return (species = species, cell = cells, count = counts, time = times)
end

# A placement's table, in the one form, with each row's species matched to its row in the species
# list as a filter matches it, and the rows naming no active species dropped, each kind counted in
# one warning; a count that is not a whole non-negative number stops it. It runs inside the
# placement's shared build, so it looks its names up itself.
function _matchedtable(table::NamedTuple, species::SpeciesList,
                       verbosity::Symbol)
    counts = [_placedcount(n, row) for (row, n) in enumerate(table.count)]
    index = _matchindex(species)
    matched = _matchspecies(table.species, counts, index,
                            _nameanswers(index, table.species))
    keep = .!isnothing.(matched.rows)
    kept = merge(map(column -> column[keep], table),
                 (species = Vector{Int}(matched.rows[keep]),
                  count = counts[keep]))
    return (table = kept,
            problems = _matchproblems(matched.dropped, ProblemWarning(),
                                      "individual", "not placed", index,
                                      verbosity))
end

# A placement's matched table checked against the grid and the time, as what is written: each row's
# species row, storage column and count. A row naming a place off the grid or a cell it may not place
# in is dropped, each kind counted in one warning; a row after the context's time stops it.
function _resolveplaced(placed::PlacedAbundances, habitat::GridHabitat,
                        context::PlacementContext)
    table = placed.table
    allowed = _allowedmask(habitat, context.cells)
    linear = LinearIndices(getgridshape(habitat))
    resolved = (species = Int[], cells = Int[], counts = Int[])
    offgrid = (rows = 0, individuals = 0)
    outside = (rows = 0, individuals = 0)
    dated = haskey(table, :time)
    for row in eachindex(table.species)
        n = table.count[row]
        dated && table.time[row] > context.time &&
            error("row $row of a placement's table is at $(table.time[row]), after " *
                  "$(context.time): a placement at build places the starting state, so every row " *
                  "is at or before the start.")
        sp = table.species[row]
        key = table.cell[row]
        if isnothing(key)
            offgrid = (rows = offgrid.rows + 1,
                       individuals = offgrid.individuals + n)
            continue
        end
        column = linear[key]
        if !allowed[column]
            outside = (rows = outside.rows + 1,
                       individuals = outside.individuals + n)
            continue
        end
        n == 0 && continue
        push!(resolved.species, sp)
        push!(resolved.cells, column)
        push!(resolved.counts, n)
    end
    problems = Problem[]
    iszero(offgrid.rows) ||
        push!(problems,
              Problem(ProblemWarning(), :off_grid,
                      "$(_individualsword(offgrid.individuals)) in $(_rowsword(offgrid.rows)) naming a " *
                      "place off the grid were not placed."))
    iszero(outside.rows) ||
        push!(problems,
              Problem(ProblemWarning(), :not_allowed,
                      "$(_individualsword(outside.individuals)) in $(_rowsword(outside.rows)) naming a " *
                      "cell the placement may not place in - an inactive cell, at build - were " *
                      "not placed."))
    return (species = resolved.species, cells = resolved.cells,
            counts = resolved.counts, problems = problems)
end

# A row's count as a whole number, refused where it is not a non-negative one.
function _placedcount(n::Real, row::Int)
    (isinteger(n) && n >= 0) ||
        error("row $row of a placement's table places $n individuals: a count is a whole, " *
              "non-negative number.")
    return Int(n)
end

function _placedcount(n, row::Int)
    return error("row $row of a placement's table gives its count as `$(repr(n))`: give a whole " *
                 "number.")
end
