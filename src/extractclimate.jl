# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Sampling a climate raster at places and times. One public verb, `extract_values`, and the helpers
# that turn each of its keywords into indices on the right axis:
#
#   place     `_axisindexat` (`cellgeometry.jl`) resolves `lat`/`long` against the spatial lookup
#   time      `_dateindices`, `_monthindices`, `_offsetindices`, `_valueindices`
#   layer     `_codeindices`
#   errors    `_requirecoord`, `_refusekeyword` - what makes a wrong keyword say what the dataset
#             does have, rather than failing inside a lookup
#
# It names no data source, working on any `AbstractClimate` whatever produced it, which is why it sits
# beside the readers rather than inside `ClimatePref`.

using Unitful

using Unitful.DefaultSymbols

using EcoSISTEM.Units

using EcoSISTEM.Units: _monthindex

using IntervalSets: Interval

using DimensionalData

using Dates: Dates

"""
    extract_values(dat::AbstractClimate; lat, long, date = :, month = :, offset = :, code = :)

Extract the value(s) of a climate dataset at the grid cell(s) containing given place(s).

**The dataset is the only positional argument** - everything else is a keyword, so a bare value is
never left to be interpreted by its position. `lat` and `long` are required; every other keyword
defaults to taking its whole dimension.

Each keyword becomes a selector on its **own named dimension**, so a 2-D grid, a monthly stack and a
multi-layer read are one code path rather than three. The result's shape follows the request: a
single place and a single slice give a scalar, and every dimension asked for by a vector is kept.
Precisely: **a dimension is dropped only when every keyword addressing it named exactly one thing**,
so a span, a vector, or a filter such as `month` on a dated axis all keep it, even where just one
slice comes back.

**Vectors cross, they do not zip.** `lat = [a, b], long = [c, d]` samples all four combinations
rather than the two pairs. Pairing a coordinate list with a per-record date belongs to sampling
*occurrence records*, which this deliberately does not do.

# Arguments

  - `dat`: the dataset to sample - any [`AbstractClimate`](@ref), including a [`ClimateRaster`](@ref)
    straight from [`read`](@ref).
  - `lat`, `long`: **required** - where to sample, as a scalar for one place or a vector for several.
    The **containing cell** is returned, so a coordinate anywhere inside a cell gives that cell's
    value, and a coordinate off the grid is an error rather than the nearest cell.
  - `date`: for a dataset on a **calendar** time axis - one `Date` or `DateTime`, a vector of them,
    or a closed interval, `Date(2000, 1, 1) .. Date(2000, 12, 31)`.
  - `month`: a month number (`March`), several (`[June, July]`) or a range (`March:May`). On a
    **monthly climatology** this names the slice; on a **dated** axis it is a *filter* meaning every
    March across all the years present, so it keeps the time dimension even for a single month, since
    how many slices match is not fixed by the request. Refused on a time axis that is not monthly,
    naming `offset` instead.
  - `offset`: for a time axis of plain elapsed-time offsets rather than calendar dates - one value, a
    vector, or an interval, in the axis's own units.
  - `code`: which layer of a multi-layer dataset, by its code (`:bio2`, or the number the catalogue
    names it by). A code is a label, so an unknown one is refused rather than read as a position.

Naming a keyword the dataset has no dimension for is an error that says what the dataset *does* have
and which keyword addresses it.
"""
function extract_values(dat::AbstractClimate; lat = _requirecoord(:lat),
                        long = _requirecoord(:long), date = Colon(),
                        month = Colon(), offset = Colon(), code = Colon())
    A = dat.array
    y = _axisindexat(DimensionalData.dims(A, Y), lat)
    x = _axisindexat(DimensionalData.dims(A, X), long)
    return A[Y(y), X(x), _extraindices(A, date, month, offset, code)...]
end

# Sentinel default for the two required keywords. Its own message rather than `EcoSISTEM._require`'s,
# which offers `DefaultEcosystem()` - a builder that has nothing to do with sampling a raster.
function _requirecoord(field)
    return error("the required keyword `$field` was not passed to `extract_values` - a sample has to " *
                 "be taken *somewhere*, so `lat` and `long` are the only compulsory arguments. Pass " *
                 "a scalar for one place or a vector for several: " *
                 "`extract_values(dat, lat = 55.9°, long = -3.2°)`.")
end

# The `Ti` or `Dim{:layer}` axis this dataset actually has, or `nothing` for a 2-D one. A dataset has
# at most one, which is `_stackaxis`' contract, so this is a lookup rather than a search.
function _extraaxis(A)
    hasdim(A, Ti) && return Ti
    hasdim(A, Dim{:layer}) && return Dim{:layer}
    return nothing
end

# Refuse a keyword that names a dimension this dataset does not have, saying what it **does** have and
# which keyword addresses it. Being able to say that is what the keyword form is for: the dimension
# is known before the value is resolved, so the message can name the alternative.
function _refusekeyword(field::Symbol, A)
    axis = _extraaxis(A)
    isnothing(axis) &&
        return error("`$field` was passed, but this dataset is a plain (lat, long) grid with no " *
                     "third axis to select along - it has only `lat` and `long`.")
    axis === Dim{:layer} &&
        return error("`$field` was passed, but this dataset's third axis is a layer axis, not a " *
                     "time axis. Use `code` to choose a layer " *
                     "(`code = $(repr(first(parent(DimensionalData.lookup(A, Dim{:layer})))))`).")
    return error("`$field` was passed, but this dataset's time axis holds " *
                 "$(eltype(parent(DimensionalData.lookup(A, Ti)))) values, which `$field` cannot " *
                 "address. Use `date` for a calendar axis, `month` for a monthly climatology, and " *
                 "`offset` for an axis of plain elapsed-time offsets.")
end

# `date =` - an exact date, a vector of them, or an interval, resolved against a `TimeType` axis.
# Refused on any other: a `Unitful.Time` axis has no calendar, so a date cannot address it, and
# guessing an epoch is what `MonthOfYearSeries` exists to make a caller state explicitly.
function _dateindices(A, date)
    lookup = DimensionalData.lookup(A, Ti)
    eltype(parent(lookup)) <: Dates.TimeType || _refusekeyword(:date, A)
    return _valueindices(lookup, date)
end

# `offset =` - the mirror of `date`, for an axis of plain elapsed-time offsets rather than calendar
# dates, which is what a model output or a non-calendar series carries. The same three shapes.
function _offsetindices(A, offset)
    lookup = DimensionalData.lookup(A, Ti)
    eltype(parent(lookup)) <: Unitful.Time || _refusekeyword(:offset, A)
    return _valueindices(lookup, offset)
end

# The shared body of `date`/`offset`: one value, several values, or a closed interval.
_valueindices(lookup, v) = DimensionalData.selectindices(lookup, At(v))

function _valueindices(lookup, vs::AbstractVector)
    return [DimensionalData.selectindices(lookup, At(v)) for v in vs]
end

function _valueindices(lookup, iv::Interval)
    return collect(DimensionalData.selectindices(lookup, iv))
end

# `month =` - a month *number* (`March`), several, or a range. Two axes can answer it, and they answer
# it differently, which is the whole reason the keyword exists rather than an axis being sniffed:
#
#   * a **monthly climatology** names its slices by month, so the month is a coordinate; and
#   * a **dated** axis has a month per slice, so "March" means *every* March in it.
#
# The coordinate is read **from the axis**, `lookup[_monthindex(m)]`, and never rebuilt as
# `m * month_mean_duration`, so a partial read such as `month = 2:4` still resolves February to
# February. `_monthindex` is the package's single 1-based conversion, shared with the plot recipes, so
# there is one place an off-by-one could live rather than several.
function _monthindices(A, month)
    lookup = DimensionalData.lookup(A, Ti)
    values = parent(lookup)
    eltype(values) <: Dates.TimeType &&
        return findall(d -> Dates.month(d) in _monthset(month), values)
    eltype(values) <: Unitful.Time || _refusekeyword(:month, A)
    _ismonthly(values) ||
        error("`month` was passed, but this dataset's time axis is not a monthly climatology: its " *
              "coordinates are $(values), which are not whole months of `month_mean_duration`. Use " *
              "`offset` to address an axis of plain elapsed-time offsets by its own coordinates.")
    idx = [DimensionalData.selectindices(lookup, At(values[_monthindex(m)]))
           for m in _monthset(month)]
    return length(idx) == 1 &&
           !(month isa Union{AbstractVector, AbstractRange}) ?
           only(idx) : idx
end

# A month request as a list of 1-based month numbers, whatever shape it arrived in.
_monthset(m::Integer) = (Int(m),)

_monthset(ms) = collect(ms)

# Whether a `Unitful.Time` axis really is a monthly climatology - every coordinate a whole number of
# `month_mean_duration`, within 1 to 12. Checked rather than assumed, because such an axis may equally
# hold plain offsets: a run on a roughly 30-day timestep produces one that *looks* monthly, and that
# is exactly where an unguarded `month = March` would pick a wrong slice in silence.
function _ismonthly(values)
    return all(values) do v
        months = uconvert(Unitful.NoUnits, v / month_mean_duration)
        return isinteger(months) && January <= months <= December
    end
end

# `code =` - a layer code (`:bio2`, or the `Int` the catalogue names it by), one or several. A code is
# a **label**, never a position: reads label the axis by code, so a wrong one is refused rather than
# silently answering with whatever band that number happened to be.
function _codeindices(A, code)
    hasdim(A, Dim{:layer}) || _refusekeyword(:code, A)
    return _valueindices(DimensionalData.lookup(A, Dim{:layer}), code)
end

# Narrow the extra axis by whichever of the four keywords were given, intersecting where more than one
# applies. That the keywords compose is the point: on a dated axis `date = <a span>, month = [June,
# July]` means "those months, within that span".
#
# The result keeps its `Int`-against-`Vector{Int}` character, and that is load-bearing: a lone `Int`
# drops the dimension under `getindex`, so a single slice yields a scalar per place and several yield
# a vector, with no flag anywhere saying which.
function _extraindices(A, date, month, offset, code)
    axis = _extraaxis(A)
    given = filter(!isnothing,
                   (date === Colon() ? nothing : (:date => date),
                    month === Colon() ? nothing : (:month => month),
                    offset === Colon() ? nothing : (:offset => offset),
                    code === Colon() ? nothing : (:code => code)))
    isnothing(axis) && !isempty(given) && _refusekeyword(first(first(given)), A)
    isnothing(axis) && return ()
    isempty(given) && return (axis(Colon()),)
    # The axis kind is checked before the keyword is resolved, so `month` on a layer-axis dataset
    # reaches the tailored refusal rather than a `MethodError` from looking up a `Ti` that is not
    # there.
    resolved = map(given) do (field, value)
        if field === :code
            axis === Dim{:layer} || _refusekeyword(:code, A)
            return _codeindices(A, value)
        end
        axis === Ti || _refusekeyword(field, A)
        field === :date && return _dateindices(A, value)
        field === :month && return _monthindices(A, value)
        return _offsetindices(A, value)
    end
    return (axis(_intersectselections(resolved, first.(given))),)
end

# Combine several keywords' selections on one axis, and decide whether the axis survives.
#
# One rule, the same one the spatial axes follow: **a dimension is dropped only when every keyword
# addressing it named exactly one thing.** So `month = March` alone gives a scalar, while
# `month = March, date = <a whole year>` gives a one-element vector - the caller asked for a span on
# `date`, and a span is a set even where it turns out to hold one slice. Worth stating because it
# looks like an inconsistency and is not: the shape follows the request, never the answer's size.
_intersectselections(sel::Tuple{Any}, _) = only(sel)

function _intersectselections(sels, fields)
    kept = sort(collect(intersect(map(_aslist, sels)...)))
    isempty(kept) &&
        error("$(join(("`$f`" for f in fields), " and ")) select no slice in common - nothing in " *
              "this dataset satisfies all of them at once.")
    return all(s -> s isa Integer, sels) ? only(kept) : kept
end

# A selection as a list, whichever shape it came back in.
_aslist(x::AbstractVector) = x

_aslist(x) = (x,)
