# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a species list is built from, one spec per trait - a tolerance, a demand, a movement - each
# parameter a trait value: one value for every species, a column of the pool, a distribution drawn
# once per species, or a function of the species' row; a spec chosen per species by a key; and the
# filling of the pool's missing values, each drawn value from a stream of its own.

using DataFrames: DataFrames, DataFrame, nrow, eachrow
using Distributions: Distributions, UnivariateDistribution, Normal
using SHA: sha256

"""
    TraitColumn(name::Symbol; unit = nothing, fill = nothing)

A trait value read from a column of a [`SpeciesPool`](@ref), one value per species, for a parameter
of a [`ToleranceSpec`](@ref), a [`DemandSpec`](@ref) or a [`MovementSpec`](@ref), or a demographic
rate of [`build_species`](@ref).

A column of plain numbers, as a CSV file gives, says the unit they are in; a column of quantities
carries its own, and a dimensionless parameter takes plain numbers as they are. When a column is
read with a unit, its mean, spread and range are noted in that unit, so a unit given wrongly - a
temperature in degrees C declared as K - shows; a column needing a unit and given none is refused
with the same summary, to tell which it is in.

# Fields

  - `name`: the column of the pool.
  - `unit`: the unit its plain numbers are in, or `nothing` where its values carry their own or the
    parameter is dimensionless; a unit given for a column of quantities must be theirs.
  - `fill`: what a missing value becomes - `nothing`, the default, refusing one; a value; a
    distribution, drawn from that species' own fill stream (see [`fill_traits`](@ref)); or a
    function of the species' row of the pool returning either. A fill is in the column's unit.
"""
struct TraitColumn{U, F}
    name::Symbol
    unit::U
    fill::F

    function TraitColumn(name::Symbol;
                         unit::Union{Nothing, Unitful.Units} = nothing,
                         fill = nothing)
        _checkfill(fill)
        return new{typeof(unit), typeof(fill)}(name, unit, fill)
    end
end

function Base.show(io::IO, c::TraitColumn)
    print(io, "TraitColumn(", repr(c.name))
    isnothing(c.unit) || print(io, ", unit = ", c.unit)
    isnothing(c.fill) || print(io, ", fill = ", c.fill)
    return print(io, ")")
end

"""
    DrawFrom(dist::UnivariateDistribution; unit = nothing)

A trait value drawn once per species from a distribution, each species' value from a stream of its
own (see [`fill_traits`](@ref)), so the same `seed` gives the same values whatever else is drawn. A
distribution holds bare numbers, so `unit` says what they are: `DrawFrom(Uniform(0.6, 2.4), unit =
km)`. A distribution given as a trait value without one is taken as dimensionless, as `longevity`
or a kernel's `shape` is.

# Fields

  - `dist`: the distribution drawn from.
  - `unit`: the unit of its values, or `nothing` for dimensionless ones.
"""
struct DrawFrom{D <: UnivariateDistribution, U}
    dist::D
    unit::U

    function DrawFrom(dist::UnivariateDistribution;
                      unit::Union{Nothing, Unitful.Units} = nothing)
        return new{typeof(dist), typeof(unit)}(dist, unit)
    end
end

function Base.show(io::IO, d::DrawFrom)
    print(io, "DrawFrom(", d.dist)
    isnothing(d.unit) || print(io, ", unit = ", d.unit)
    return print(io, ")")
end

"""
    AbstractTraitSpec

What [`build_species`](@ref) builds one trait of every species from: a [`ToleranceSpec`](@ref), a
[`DemandSpec`](@ref), a [`MovementSpec`](@ref), or a [`ChosenBy`](@ref) choosing among them per
species. Each parameter is a trait value - one value for every species, a [`TraitColumn`](@ref), a
[`DrawFrom`](@ref), or a function of the species' row of the pool.
"""
abstract type AbstractTraitSpec end

"""
    ToleranceSpec(axis::Type{<:NicheAxis}, D::Type; params...)
    ToleranceSpec(axis::Type{<:NicheAxis}; μ, σ)
    ToleranceSpec(axis::Type{<:NicheAxis}, D::Type, params...)

Each species' tolerance on one niche axis, a response distribution of family `D` - `Normal` by
default - built into a [`NicheTolerance`](@ref) in the axis's own unit. Its parameters are named as
the family names them - `ToleranceSpec(Temperature, Normal, μ = TraitColumn(:tmin, unit = K), σ =
2.0K)`, `ToleranceSpec(Precipitation, Gamma, α = 2.0, θ = TraitColumn(:tp, unit = mm / d))` - which
no positional order can confuse. A family whose fields are not its parameters (`InverseGamma`) takes
them in order instead, as its own constructor does; naming them is refused. Each parameter is
converted by its role: a location as a position, a spread as an interval, a shape as a plain number.

# Arguments

  - `axis`: the niche axis the tolerance is on, which the habitat's regime must match.
  - `D`: the family of the response distribution, `Normal` by default.
  - `params`: its parameters, each a trait value ([`AbstractTraitSpec`](@ref EcoSISTEM.AbstractTraitSpec)),
    by the family's own names, or in its order.
"""
struct ToleranceSpec{A <: NicheAxis, D <: UnivariateDistribution, P <: Tuple} <:
       AbstractTraitSpec
    params::P

    function ToleranceSpec(axis::Type{A}, ::Type{D},
                           params...) where {A <: NicheAxis,
                                             D <: UnivariateDistribution}
        n = length(param_roles_resolved(D, probes = _defaultprobes(D)))
        length(params) == n ||
            error("`$D` takes $n parameters, but $(length(params)) were given.")
        foreach(_checktraitvalue, params)
        return new{A, D, typeof(params)}(params)
    end
end

function ToleranceSpec(axis::Type{<:NicheAxis},
                       ::Type{D}; params...) where {D <: UnivariateDistribution}
    names = _keywordparams(D)
    given = collect(keys(params))
    unknown = setdiff(given, names)
    absent = setdiff(names, given)
    (isempty(unknown) && isempty(absent)) ||
        error("`$D`'s parameters are $(join(("`$n`" for n in names), ", "))" *
              (isempty(unknown) ? "" :
               "; $(join(("`$n`" for n in unknown), ", ")) is not one") *
              (isempty(absent) ? "" :
               "; $(join(("`$n`" for n in absent), ", ")) is missing") * ".")
    return ToleranceSpec(axis, D, (params[n] for n in names)...)
end

function ToleranceSpec(axis::Type{<:NicheAxis}; params...)
    return ToleranceSpec(axis, Normal; params...)
end

# Printed as it was written: by name where the family's parameters can be named, in order where not.
function Base.show(io::IO, s::ToleranceSpec{A, D}) where {A, D}
    print(io, "ToleranceSpec(", nameof(A), ", ", nameof(D))
    names = try
        _keywordparams(D)
    catch
        nothing
    end
    for (i, p) in enumerate(s.params)
        print(io, ", ", isnothing(names) ? "" : "$(names[i]) = ", p)
    end
    return print(io, ")")
end

"""
    DemandSpec(axis::Type{<:NicheAxis}, value)

Each species' demand on one resource, built into a [`Demand`](@ref) in the axis's own unit:
`DemandSpec(SolarRadiation, TraitColumn(:ssr, unit = kJ / d))`.

# Arguments

  - `axis`: the resource axis, which the habitat's supply must match.
  - `value`: the demand, a trait value ([`AbstractTraitSpec`](@ref EcoSISTEM.AbstractTraitSpec)),
    in a unit of the axis's dimension.
"""
struct DemandSpec{A <: NicheAxis, V} <: AbstractTraitSpec
    value::V

    function DemandSpec(axis::Type{A}, value) where {A <: NicheAxis}
        _demandtype(axis)
        _checktraitvalue(value)
        return new{A, typeof(value)}(value)
    end
end

function Base.show(io::IO, s::DemandSpec{A}) where {A}
    return print(io, "DemandSpec(", nameof(A), ", ", s.value, ")")
end

"""
    MovementSpec(mode::AbstractMovement; kernel...)

Each species' movement: a dispersal mode ([`AbstractMovement`](@ref EcoSISTEM.AbstractMovement)),
and the keywords of its kernel, each a trait value -
`MovementSpec(InsectDispersal(), meandist = TraitColumn(:dispersal, unit = km))`,
`MovementSpec(WindDispersal(), meandist = 5.0km, scale = 10.0km)`. [`NoDispersal`](@ref) takes none.
Every species gets the mode and a kernel of its family, as the mode called with those keywords
gives one.

# Arguments

  - `mode`: the dispersal mode, which fixes when individuals move, the kernel's family and what
    becomes of one aimed at a dead cell.
  - `kernel`: the kernel's keywords - a distance (`meandist`, `sd` or `scale`), a `shape`, a
    `thresh` - as the family takes them.
"""
struct MovementSpec{M <: AbstractMovement, P <: NamedTuple} <: AbstractTraitSpec
    mode::M
    params::P

    function MovementSpec(mode::AbstractMovement; kernel...)
        params = NamedTuple(kernel)
        allowed = _kernelkeywords(mode)
        unknown = setdiff(collect(keys(params)), allowed)
        isempty(unknown) ||
            error("$(sprint(show, mode)) takes " *
                  (isempty(allowed) ? "no kernel keywords" :
                   "the kernel keywords $(join(("`$k`" for k in allowed), ", "))") *
                  "; $(join(("`$k`" for k in unknown), ", ")) is not one.")
        foreach(_checktraitvalue, values(params))
        return new{typeof(mode), typeof(params)}(mode, params)
    end
end

function Base.show(io::IO, s::MovementSpec)
    print(io, "MovementSpec(", sprint(show, s.mode))
    isempty(s.params) ||
        print(io, "; ",
              join(("$k = $(sprint(show, v))" for (k, v) in pairs(s.params)),
                   ", "))
    return print(io, ")")
end

"""
    ChosenBy(key, choices::Pair...)

A spec chosen per species: each species' `key` - the value of a column of its row, or what a function
of its row returns - picks the spec of the choice it names, so species of one pool can disperse by
different modes, each with its own kernel and parameters:

```julia
ChosenBy(:syndrome,
         "wind" => MovementSpec(WindDispersal(), meandist = 5.0km, scale = 10.0km),
         "insect" => MovementSpec(InsectDispersal(), meandist = 1.0km))
```

Every species takes exactly one spec; a species whose key no choice names, or is `missing`, is
refused, named. A movement may be chosen so; a tolerance not yet.

# Arguments

  - `key`: a column of the pool, by name, or a function of a species' row of the pool.
  - `choices`: `key => spec` pairs, each spec of the same kind.
"""
struct ChosenBy{K, C <: Tuple} <: AbstractTraitSpec
    key::K
    choices::C

    function ChosenBy(key::Union{Symbol, Function}, choices::Pair...)
        isempty(choices) &&
            error("`ChosenBy` needs at least one `key => spec` choice.")
        specs = map(last, choices)
        all(s -> s isa MovementSpec, specs) ||
            error("`ChosenBy` chooses among `MovementSpec`s; " *
                  "$(join(unique(string.(nameof.(typeof.(specs)))), ", ")) cannot be chosen per " *
                  "species yet.")
        allunique(map(first, choices)) ||
            error("`ChosenBy` names a key twice; each choice needs a key of its own.")
        return new{typeof(key), typeof(choices)}(key, choices)
    end
end

function Base.show(io::IO, c::ChosenBy)
    print(io, "ChosenBy(", repr(c.key))
    for (k, s) in c.choices
        print(io, ", ", repr(k), " => ", s)
    end
    return print(io, ")")
end

# One parameter of one spec as the build reads it: its `name`, which a drawn value's column takes;
# its trait `value`; the `rows` of the species it is for (`Colon()` for all); and the unit it needs,
# `nothing` where nothing says.
struct _TraitSlot
    name::Symbol
    value::Any
    rows::Union{Colon, AbstractVector{Int}}
    expected::Union{Nothing, Unitful.Units}
end

# == Functions ==================================================================================

"""
    fill_traits(pool::SpeciesPool, specs...; seed = nothing, verbosity = :normal)

Fill in the values `specs` will read from `pool`, and return them, building nothing: the first step
[`build_species`](@ref) takes, made public so that the values can be checked first. The result is a
named tuple of two tables, one row per species:

  - `values`: the pool as a table, identity and values, with every missing value a
    [`TraitColumn`](@ref) reads filled from its `fill`, and a column for each value drawn or computed
    for every species - each [`DrawFrom`](@ref), distribution and function - under the name of the
    parameter it is for (`movement_meandist`, `tolerance_σ`);
  - `filled`: per column, `true` where a value was filled in or drawn rather than given.

```julia
niche = ToleranceSpec(Temperature, μ = TraitColumn(:tmin, unit = K, fill = 286.0),
                      σ = TraitColumn(:tmin_sd, unit = K, fill = Uniform(1.0, 3.0)))
filled = fill_traits(pool, niche, seed = 1)
filled.values.tmin_sd     # each species' width, the missing ones drawn
filled.filled.tmin_sd     # which of them were drawn
```

A missing value with no fill is refused, naming the species and the column; removing the species,
or giving the column a `fill`, is the remedy. Each value drawn comes from a stream of its own, seeded
from `seed`, the parameter and the species' key in the pool, so a species draws the same value
however many others are drawn for, and a copy ([`copy_species!`](@ref)) its original's. Left unset,
the draws come from the global generator.

# Arguments

  - `pool`: the species.
  - `specs`: what reads the values - [`ToleranceSpec`](@ref)s, [`DemandSpec`](@ref)s,
    [`MovementSpec`](@ref)s, [`ChosenBy`](@ref)s, tuples or named tuples of them, or `name => value`
    pairs naming a trait value directly.
  - `seed`: an integer fixing every value drawn, or `nothing`.
  - `verbosity`: `:silent` notes nothing, `:normal` how many values of each column were filled.
"""
function fill_traits(pool::SpeciesPool, specs...;
                     seed::Union{Nothing, Integer} = nothing,
                     verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    slots = _TraitSlot[]
    for spec in specs
        _slots!(slots, pool, _slotprefix(spec), spec, Colon())
    end
    filled = _fillpool(pool, slots, seed, verbosity)
    return (values = DataFrames.select(filled.values,
                                       DataFrames.Not([:key, :copykey])),
            filled = filled.filled)
end

# A spec's parameters as slots, named under `prefix`, for the species at `rows`.
function _slots!(slots, pool, prefix::AbstractString,
                 spec::ToleranceSpec{A, D}, rows) where {A, D}
    units = _expectedunits(spec)
    for (n, p, u) in zip(_paramnames(D), spec.params, units)
        push!(slots, _TraitSlot(Symbol(prefix, "_", n), p, rows, u))
    end
    return slots
end

function _slots!(slots, pool, prefix::AbstractString, spec::DemandSpec{A},
                 rows) where {A}
    return push!(slots,
                 _TraitSlot(Symbol(prefix), spec.value, rows,
                            canonicalunit(Resource, A)))
end

function _slots!(slots, pool, prefix::AbstractString, spec::MovementSpec, rows)
    for (k, v) in pairs(spec.params)
        push!(slots,
              _TraitSlot(Symbol(prefix, "_", k), v, rows, _kernelunit(k)))
    end
    return slots
end

function _slots!(slots, pool, prefix::AbstractString, spec::ChosenBy, rows)
    chosen = _chosenspecs(pool, spec)
    for (i, (_, s)) in enumerate(spec.choices)
        mine = [r for r in _rowsof(rows, length(pool)) if chosen[r] == i]
        _slots!(slots, pool, prefix, s, mine)
    end
    return slots
end

function _slots!(slots, pool, prefix::AbstractString, specs::NamedTuple, rows)
    for (k, s) in pairs(specs)
        _slots!(slots, pool, string(prefix, "_", k), s, rows)
    end
    return slots
end

function _slots!(slots, pool, prefix::AbstractString, specs::Tuple, rows)
    for s in specs
        _slots!(slots, pool, string(prefix, "_", _specaxisname(s)), s, rows)
    end
    return slots
end

function _slots!(slots, pool, prefix::AbstractString, pair::Pair{Symbol}, rows)
    return push!(slots, _TraitSlot(first(pair), last(pair), rows, nothing))
end

# A requirement built already reads nothing from the pool. Deliberately untyped, as the fallback for
# whatever `build_species` accepts built.
_slots!(slots, pool, prefix::AbstractString, ::Any, rows) = slots

# The prefix a spec's slots are named under where it is given on its own.
_slotprefix(::ToleranceSpec) = "tolerance"

_slotprefix(::DemandSpec) = "demand"

_slotprefix(::Union{MovementSpec, ChosenBy}) = "movement"

_slotprefix(::Pair{Symbol}) = ""

_slotprefix(specs::Union{Tuple, NamedTuple}) = _slotprefix(first(values(specs)))

# How a spec in a tuple is named among its siblings: by its axis.
function _specaxisname(::Union{ToleranceSpec{A}, DemandSpec{A}}) where {A}
    return string(nameof(A))
end

_specaxisname(::Any) = "spec"

# The rows a slot's `rows` stand for.
_rowsof(::Colon, n::Integer) = 1:n

_rowsof(rows::AbstractVector, ::Integer) = rows

# The names a family's parameters are given by, refusing a family whose fields are not its
# parameters: a probe instance's fields must be its `params`, in order, or naming them would place
# a value in the wrong slot.
function _keywordparams(::Type{D}) where {D}
    d = try
        guess_params(D, probes = _defaultprobes(D))
    catch
        nothing
    end
    names = fieldnames(D)
    ok = !isnothing(d) && length(Distributions.params(d)) == length(names) &&
         all(isequal(getfield(d, f), p)
             for (f, p) in zip(names, Distributions.params(d)))
    ok ||
        error("`$D`'s fields are not its parameters, so they cannot be named; give them in " *
              "order, `ToleranceSpec(axis, $D, params...)`, as `$D`'s constructor takes them.")
    return collect(names)
end

# A family's parameters' names where they can be named, else `param1`, `param2` and so on.
function _paramnames(::Type{D}) where {D}
    n = length(param_roles_resolved(D, probes = _defaultprobes(D)))
    names = try
        _keywordparams(D)
    catch
        nothing
    end
    return isnothing(names) ? [Symbol("param", i) for i in 1:n] : names
end

# The unit each of a tolerance's parameters needs, by its role in the axis's own frame.
function _expectedunits(::ToleranceSpec{A, D}) where {A, D}
    return paramunits(D, _defaultsupport(A), probes = _defaultprobes(D))
end

# The unit a kernel keyword needs: a distance, or none for a shape and a threshold.
_kernelunit(k::Symbol) = k in (:shape, :thresh) ? NoUnits : km

# The kernel keywords a mode's family takes, read off its keyword constructor; none for a mode that
# disperses nothing.
_kernelkeywords(::NoDispersal) = Symbol[]

function _kernelkeywords(::_KernelDispersal{K}) where {K}
    for m in methods(K)
        names = Base.kwarg_decl(m)
        isempty(names) || return names
    end
    return Symbol[]
end

# Refuse what is not a trait value, a vector above all: values per species are never taken by
# position.
function _checktraitvalue(v)
    v isa
    Union{Number, TraitColumn, DrawFrom, UnivariateDistribution, Function} &&
        return nothing
    v isa AbstractVector &&
        error("a trait value is never a vector, whose order nothing ties to the species: add the " *
              "values to the pool with `add_trait!`, keyed by species, and name their column with " *
              "`TraitColumn`.")
    return error("a trait value is one value for every species, a `TraitColumn`, a `DrawFrom`, a " *
                 "distribution or a function of the species' row; got a $(typeof(v)).")
end

# Refuse a fill that is none of the forms a fill takes.
function _checkfill(fill)
    fill isa Union{Nothing, Number, UnivariateDistribution, Function} &&
        return nothing
    return error("a `fill` is a value, a distribution or a function of the species' row; got a " *
                 "$(typeof(fill)).")
end

# The word standing for a fill in a stream's seed - "fill" in ASCII - distinct from the abundance
# split's and a phylogeny's, so no two kinds of draw share a stream.
const _FILL_STREAM = 0x66696c6c

# The generator one species' value of one parameter is drawn from: a stream of its own with a seed,
# keyed by the parameter's name and the species' key; the global one without.
_fillrng(::Nothing, ::Symbol, ::Integer) = Random.default_rng()

function _fillrng(seed::Integer, slot::Symbol, key::Integer)
    return Random.Xoshiro(UInt64[seed % UInt64, _FILL_STREAM, _slotword(slot),
                                 key % UInt64])
end

# A parameter's name as a word of a stream's seed: its SHA-256's first eight bytes, never
# `Base.hash`, whose values change between Julia versions.
_slotword(slot::Symbol) = first(reinterpret(UInt64, sha256(string(slot))[1:8]))

# The pool's values with what `slots` read filled: each column a `TraitColumn` names filled where
# missing, a column per drawn or computed value, and a mask of what was filled or drawn. A species'
# draws come from its key's streams, a copy's from its original's.
function _fillpool(pool::SpeciesPool, slots::AbstractVector{_TraitSlot},
                   seed, verbosity::Symbol)
    table = copy(pool.species)
    original = copy(table)
    mask = DataFrame()
    n = nrow(table)
    streamkeys = Int[coalesce(c, k) for (c, k) in zip(table.copykey, table.key)]
    problems = Problem[]
    unfilled = Pair{Symbol, Vector{String}}[]
    for (name, column) in _columnsread(slots, n)
        name in propertynames(table) ||
            error("the pool has no column `$name` for $(column.spec); its values are " *
                  "$(join(("`$c`" for c in _traitnames(pool)), ", ")).")
        empty = [r for r in column.rows if ismissing(table[r, name])]
        if isempty(empty)
            mask[!, name] = falses(n)
        elseif isnothing(column.spec.fill)
            push!(unfilled, name => table.name[empty])
        else
            filled = Vector{Any}(table[!, name])
            for r in empty
                filled[r] = _fillvalue(column.spec, original[r, :],
                                       table[!, name],
                                       _fillrng(seed, name, streamkeys[r]),
                                       table.name[r])
            end
            table[!, name] = identity.(filled)
            mask[!, name] = [r in empty for r in 1:n]
            push!(problems,
                  Problem(ProblemNotice(), :filled,
                          "`:$name`: $(length(empty)) missing " *
                          "value$(length(empty) == 1 ? "" : "s") filled."))
        end
    end
    isempty(unfilled) ||
        _report(ProblemError(),
                "values are missing, with no `fill` to give them: " *
                join(("`:$c` for $(_namelist(s, verbosity))"
                      for (c, s) in unfilled),
                     "; ") *
                ". Give the column a `fill`, or remove those species from the pool.")
    for (name, drawn) in _drawsread(slots)
        name in propertynames(pool.species) &&
            error("the pool already has a column `$name`, which the values drawn for that " *
                  "parameter would take; rename the column.")
        column = Vector{Any}(missing, n)
        drew = falses(n)
        for (slot, rows) in drawn, r in _rowsof(rows, n)
            column[r], drew[r] = _drawvalue(slot.value, table[r, :],
                                            _fillrng(seed, name, streamkeys[r]))
        end
        table[!, name] = identity.(column)
        mask[!, name] = drew
    end
    _reportpool(problems, verbosity)
    return (values = table, filled = mask)
end

# The columns `slots` read, each once, with the rows read of the pool's `n` and the `TraitColumn`
# read by; two reading one column must agree on its unit and fill.
function _columnsread(slots, n::Integer)
    columns = Dict{Symbol, Any}()
    order = Symbol[]
    for slot in slots
        slot.value isa TraitColumn || continue
        name = slot.value.name
        rows = collect(_rowsof(slot.rows, n))
        if haskey(columns, name)
            seen = columns[name]
            (isequal(seen.spec.unit, slot.value.unit) &&
             isequal(seen.spec.fill, slot.value.fill)) ||
                error("`:$name` is read twice with different units or fills, " *
                      "$(seen.spec) and $(slot.value); give one.")
            columns[name] = (spec = seen.spec,
                             rows = sort!(union(seen.rows, rows)))
        else
            push!(order, name)
            columns[name] = (spec = slot.value, rows = rows)
        end
    end
    return [name => columns[name] for name in order]
end

# The slots whose values are drawn or computed per species, by the column their values take - one
# per parameter, the choices of a `ChosenBy` sharing it, each for its own species.
function _drawsread(slots)
    drawn = Dict{Symbol, Vector{Tuple{_TraitSlot, Any}}}()
    order = Symbol[]
    for slot in slots
        slot.value isa Union{DrawFrom, UnivariateDistribution, Function} ||
            continue
        haskey(drawn, slot.name) || push!(order, slot.name)
        push!(get!(Vector{Tuple{_TraitSlot, Any}}, drawn, slot.name),
              (slot, slot.rows))
    end
    return [name => drawn[name] for name in order]
end

# One missing value of a column filled: its `fill`, drawn from where it is a distribution, called on
# the species' row where it is a function, in the column's unit - a column of quantities taking the
# unit of its values, a column of plain numbers a plain number.
function _fillvalue(column::TraitColumn, row, values, rng,
                    species::AbstractString)
    v = column.fill isa Function ? column.fill(row) : column.fill
    v isa UnivariateDistribution && (v = rand(rng, v))
    v isa Union{Missing, Nothing} &&
        error("the fill of `:$(column.name)` gave `$v` for `$species`; it gives a value.")
    v isa Number ||
        error("the fill of `:$(column.name)` gave a $(typeof(v)) for `$species`; " *
              "it gives a value or a distribution.")
    u = _columnunit(values)
    isnothing(u) && return v isa Unitful.AbstractQuantity ?
           _plainin(column, v) : v
    return v isa Unitful.AbstractQuantity ? v : v * u
end

# A filled quantity as a plain number of its column, which says its unit.
function _plainin(column::TraitColumn, v)
    isnothing(column.unit) &&
        error("the fill of `:$(column.name)` gives a quantity, $v, for a column of plain numbers " *
              "with no `unit`; give the column its unit.")
    return ustrip(column.unit, v)
end

# The unit a column's values carry, or `nothing` for plain numbers.
function _columnunit(values)
    for v in values
        ismissing(v) && continue
        return v isa Unitful.AbstractQuantity ? Unitful.unit(v) : nothing
    end
    return nothing
end

# One species' value of a drawn or computed parameter, with whether it was drawn: a `DrawFrom`'s
# draw in its unit; a distribution's, dimensionless; a function's value, or its draw where it gives a
# distribution or a `DrawFrom`.
_drawvalue(d::DrawFrom, row, rng) = (_drawn(d, rng), true)

_drawvalue(d::UnivariateDistribution, row, rng) = (rand(rng, d), true)

function _drawvalue(f::Function, row, rng)
    v = f(row)
    v isa Union{DrawFrom, UnivariateDistribution} &&
        return (first(_drawvalue(v, row, rng)), true)
    v isa Number ||
        error("a trait value's function gave a $(typeof(v)) for `$(row.name)`; " *
              "it gives a value, a distribution or a `DrawFrom`.")
    return (v, false)
end

_drawn(d::DrawFrom, rng) = _withunit(rand(rng, d.dist), d.unit)

# A drawn number in a `DrawFrom`'s unit, or as it is where it has none.
_withunit(x, ::Nothing) = x

_withunit(x, u::Unitful.Units) = x * u

# Which choice of a `ChosenBy` each species of the pool takes, by position among the choices,
# refusing a species whose key names none or is missing.
function _chosenspecs(pool::SpeciesPool, spec::ChosenBy)
    table = _pooltable(pool)
    keys_ = [k for (k, _) in spec.choices]
    chosen = Vector{Int}(undef, nrow(table))
    unmatched = String[]
    for (r, row) in enumerate(eachrow(table))
        key = spec.key isa Symbol ? _keycolumn(row, spec.key) : spec.key(row)
        i = ismissing(key) ? nothing : findfirst(isequal(key), keys_)
        isnothing(i) ? push!(unmatched, "`$(row.name)` ($(repr(key)))") :
        (chosen[r] = i)
    end
    isempty(unmatched) ||
        error("$(length(unmatched)) species' keys name no choice of $(spec): " *
              join(first(unmatched, 10), ", ") *
              (length(unmatched) > 10 ? " and $(length(unmatched) - 10) more" :
               "") *
              ". Each key needs a choice; the choices are $(join(repr.(keys_), ", ")).")
    return chosen
end

# A species' key read from its row, refusing a column the pool does not have.
function _keycolumn(row, name::Symbol)
    name in propertynames(row) ||
        error("`ChosenBy` reads its key from `:$name`, which the pool has no column of.")
    return row[name]
end

# One parameter's values for the species at `rows`, in a unit of the dimension `expected` asks for -
# plain `Float64`s where it is dimensionless - from the filled values `table`: one value for all, a
# column (its plain numbers given their `unit`, refused without one, and noted in it), or the
# column of values drawn or computed for it. `what` names the parameter in a message.
function _resolvevalue(slot::_TraitSlot, table::DataFrame, rows,
                       what::AbstractString, problems, verbosity::Symbol)
    v = slot.value
    if v isa Number
        return fill(_ofdimension(v, slot.expected, what), length(rows))
    elseif v isa TraitColumn
        return _readcolumn(v, table, rows, slot.expected, what, problems,
                           verbosity)
    end
    return [_ofdimension(x, slot.expected, what)
            for x in table[rows, slot.name]]
end

# A column's values for `rows`, in their unit: plain numbers given the column's `unit`, or taken as
# they are where the parameter is dimensionless; refused, with what the numbers look like, where a
# unit is needed and none is given. A column read with a unit is noted, in it.
function _readcolumn(c::TraitColumn, table::DataFrame, rows, expected, what,
                     problems, verbosity)
    raw = table[rows, c.name]
    carried = _columnunit(raw)
    dimensionless = !isnothing(expected) && expected == NoUnits
    if isnothing(carried)
        if dimensionless
            isnothing(c.unit) ||
                error("$what is dimensionless, so `$(c)` needs no `unit`.")
            return Float64.(raw)
        end
        isnothing(c.unit) &&
            error("$c gives $what plain numbers, so it needs a `unit`" *
                  (isnothing(expected) ? "" :
                   " ($(expected), or another of its dimension)") *
                  ": " * _summaryphrase(raw) * ".")
        quantities = raw .* c.unit
    else
        (isnothing(c.unit) || c.unit == carried) ||
            error("$c says `$(c.unit)`, but its values are in `$carried`.")
        quantities = raw
    end
    checked = [_ofdimension(x, expected, what) for x in quantities]
    push!(problems,
          Problem(ProblemNotice(), :column_read,
                  "`:$(c.name)` read in $(Unitful.unit(first(checked))): " *
                  _summaryphrase(checked) * "."))
    return checked
end

# One value of a parameter, checked against the dimension it needs - and made a plain `Float64` where
# that is none.
function _ofdimension(x, expected, what)
    isnothing(expected) && return x
    if expected == NoUnits
        x isa Unitful.AbstractQuantity && dimension(x) != NoDims &&
            error("$what is dimensionless, but was given $x.")
        return Float64(x isa Unitful.AbstractQuantity ? uconvert(NoUnits, x) :
                       x)
    end
    (x isa Unitful.AbstractQuantity && dimension(x) == dimension(expected)) ||
        error("$what needs a value in $(expected) or another unit of its dimension; got $x.")
    return x
end

# What a column's values look like: mean, standard deviation and range in their unit, and how many
# there are. Computed on the numbers in that unit, since a mean of temperatures in degrees C is not
# a sum Unitful will form.
function _summaryphrase(values)
    present = collect(skipmissing(values))
    isempty(present) && return "no values"
    u = first(present) isa Unitful.AbstractQuantity ?
        Unitful.unit(first(present)) :
        nothing
    numbers = isnothing(u) ? Float64.(present) :
              Float64[ustrip(u, x) for x in present]
    lo, hi = extrema(numbers)
    spread = length(numbers) > 1 ? Distributions.std(numbers) : 0.0
    missingcount = length(values) - length(present)
    shown(x) = _fourfigures(x) * (isnothing(u) ? "" : " $u")
    return "mean $(shown(Distributions.mean(numbers))), sd $(shown(spread)), " *
           "range $(shown(lo)) to $(shown(hi)) ($(length(present)) values" *
           (missingcount > 0 ? ", $missingcount missing" : "") * ")"
end

# A number for a summary, to four significant figures.
_fourfigures(x::Real) = string(round(x, sigdigits = 4))

# --- Building each trait from the filled values -------------------------------
# Each walks a spec as `_slots!` names it, reading every parameter's values for the species at `rows`
# from the filled values `table`; `problems` gathers the notes of the columns read.

# A tolerance: each parameter converted by its role into the axis's own frame, as `NicheTolerance`
# reads bare numbers; a collection member by member; one built already as it is.
function _buildtrait(spec::ToleranceSpec{A, D}, pool, table, prefix, rows,
                     problems, verbosity) where {A, D}
    slots = _slots!(_TraitSlot[], pool, prefix, spec, rows)
    support = _defaultsupport(A)
    roles = param_roles_resolved(D, probes = _defaultprobes(D))
    bare = map(slots, roles) do slot, role
        values = _resolvevalue(slot, table, rows,
                               "`$(slot.name)` of $(spec)", problems, verbosity)
        return [_bareparam(role, support, x) for x in values]
    end
    return NicheTolerance(A, D, bare..., support = support)
end

function _buildtrait(spec::DemandSpec{A}, pool, table, prefix, rows, problems,
                     verbosity) where {A}
    slot = only(_slots!(_TraitSlot[], pool, prefix, spec, rows))
    return _demand(A,
                   _resolvevalue(slot, table, rows, "the demand of $(spec)",
                                 problems, verbosity))
end

function _buildtrait(specs::NamedTuple, pool, table, prefix, rows, problems,
                     verbosity)
    built = map(keys(specs), values(specs)) do k, s
        return _buildtrait(s, pool, table, string(prefix, "_", k), rows,
                           problems, verbosity)
    end
    return SpeciesRequirementCollection(NamedTuple{keys(specs)}(built))
end

function _buildtrait(specs::Tuple, pool, table, prefix, rows, problems,
                     verbosity)
    built = map(specs) do s
        return _buildtrait(s, pool, table,
                           string(prefix, "_", _specaxisname(s)),
                           rows, problems, verbosity)
    end
    return SpeciesRequirementCollection(built)
end

function _buildtrait(built::AbstractTolerance, pool, table, prefix, rows,
                     problems, verbosity)
    _checktolerancespecies(built, length(rows))
    return built
end

function _buildtrait(built::AbstractDemand, pool, table, prefix, rows,
                     problems, verbosity)
    counttypes(built) == length(rows) ||
        error("a demand built already covers $(counttypes(built)) species, but the pool has " *
              "$(length(rows)): it is used as given, so it needs one entry per species.")
    return built
end

# Anything else is not a tolerance or a demand. Deliberately untyped, to say so rather than raise a
# `MethodError`.
function _buildtrait(other, pool, table, prefix, rows, problems, verbosity)
    return error("`$prefix` is a spec - `ToleranceSpec`, `DemandSpec` - a tuple or named tuple of " *
                 "them, or a requirement built already; got a $(typeof(other)).")
end

# One parameter value as a bare number in a tolerance's frame `support`, converted by its role: a
# location as a position, a spread as an interval, a rate as its inverse; a dimensionless one as it
# is.
function _bareparam(role, support, x)
    x isa Unitful.AbstractQuantity || return Float64(x)
    return Float64(_readparam(role, support, Unitful.unit(x), ustrip(x)))
end

# Each species' movement: its mode called with its kernel's values, every species of a `ChosenBy`
# from the spec its key chooses.
function _buildmovement(spec::MovementSpec, pool, table, prefix, rows,
                        problems, verbosity)
    spec.mode isa NoDispersal &&
        return SpeciesDispersal[spec.mode() for _ in rows]
    slots = _slots!(_TraitSlot[], pool, prefix, spec, rows)
    kernelvalues = map(zip(keys(spec.params), slots)) do (k, slot)
        return _resolvevalue(slot, table, rows, "`$k` of $(spec)", problems,
                             verbosity)
    end
    kernelnames = keys(spec.params)
    return SpeciesDispersal[spec.mode(;
                                      NamedTuple{kernelnames}(Tuple(v[i]
                                                                    for v in kernelvalues))...)
                            for i in eachindex(rows)]
end

function _buildmovement(spec::ChosenBy, pool, table, prefix, rows, problems,
                        verbosity)
    chosen = _chosenspecs(pool, spec)
    movement = Vector{SpeciesDispersal}(undef, length(rows))
    for (i, (_, s)) in enumerate(spec.choices)
        mine = findall(r -> chosen[r] == i, rows)
        isempty(mine) && continue
        movement[mine] = _buildmovement(s, pool, table, prefix, rows[mine],
                                        problems, verbosity)
    end
    return movement
end

function _buildmovement(other, pool, table, prefix, rows, problems, verbosity)
    return error("`movement` is a `MovementSpec`, or a `ChosenBy` choosing among them; got a " *
                 "$(typeof(other)).")
end

# The demographic rates, every species' own: births and deaths per unit time, in one unit, and each
# species' longevity; survival one number for every species.
function _buildparams(birth, death, longevity, survival, pool, table, rows,
                      problems, verbosity)
    survival isa Real ||
        error("`survival` is one number for every species; got $(survival).")
    rate = inv(Unitful.yr)
    resolved(name, value, expected) = _resolvevalue(_TraitSlot(name, value,
                                                               Colon(),
                                                               expected),
                                                    table, rows, "`$name`",
                                                    problems, verbosity)
    b = resolved(:birth, birth, rate)
    d = resolved(:death, death, rate)
    u = Unitful.unit(first(b))
    return PopGrowth{typeof(u)}(float.(uconvert.(u, b)),
                                float.(uconvert.(u, d)),
                                resolved(:longevity, longevity, NoUnits),
                                Float64(survival))
end
