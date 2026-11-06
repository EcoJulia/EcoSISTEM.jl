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
[`NoAbundance`](@ref), and several added with `+` ([`SumOfPlacements`](@ref)). Called with
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
        _checktabledata(data)
        return new{typeof(data)}(data, _recordvector(provenance))
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

_reusable(p::SumOfPlacements) = all(_reusable, p.placements)

# Whether a table's data can be read again: a file, a vector of rows, or a table of columns.
_rereadable(::Union{AbstractString, ProvenanceRecord, AbstractVector}) = true

_rereadable(table) = Tables.columnaccess(table)

# A table's data checked as it is given: a file that exists, a record naming one, or a table.
function _checktabledata(path::AbstractString)
    isfile(path) ||
        error("`FromTable` was given the path $(repr(path)), and there is no file there.")
    return nothing
end

function _checktabledata(record::ProvenanceRecord)
    isnothing(record.location) &&
        error("$record names no file: a record `Provenancer.record_provenance` returned, or " *
              "`provenance(path)` read, says where its file is, and `FromTable` takes one of those.")
    return nothing
end

function _checktabledata(table)
    Tables.istable(table) ||
        error("`FromTable` reads a path to a CSV file, a `ProvenanceRecord` of one, or a " *
              "Tables.jl table; got a value of type $(typeof(table)).")
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

# --- Checking a table ---------------------------------------------------------

# A placement's table in one form, whatever form it came in: `species` as given, `cell` as the
# grid's own key (`nothing` for a place off the grid), `count` (one each where the table has none)
# and, where it is dated, `time` as elapsed time.
function _canonicaltable(table, habitat::GridHabitat,
                         context::PlacementContext)
    columns = Tables.columns(table)
    names = Tables.columnnames(columns)
    :species in names ||
        error("a placement's table names each row's species in a `species` column; it has " *
              "$(_columnsphrase(names)).")
    species = collect(Tables.getcolumn(columns, :species))
    cells = _gridkeys(habitat, columns, names)
    counts = :count in names ? collect(Tables.getcolumn(columns, :count)) :
             fill(1, length(species))
    :time in names ||
        return (species = species, cell = cells, count = counts)
    times = [_placementtime(t, context, row)
             for (row, t) in enumerate(Tables.getcolumn(columns, :time))]
    return (species = species, cell = cells, count = counts, time = times)
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
# one warning; a count that is not a whole non-negative number stops it.
function _matchedtable(table::NamedTuple, species::SpeciesList,
                       verbosity::Symbol)
    counts = [_placedcount(n, row) for (row, n) in enumerate(table.count)]
    index = _matchindex(species)
    matched = _matchspecies(table.species, counts, index)
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
