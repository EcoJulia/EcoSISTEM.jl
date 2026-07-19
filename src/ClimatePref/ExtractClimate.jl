# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using EcoSISTEM: LatLong
using AxisArrays
using IndexedTables

# Grid step of a climate array's `dim`-th coordinate axis.
function _axisstep(arr::AxisArray, dim::Int)
    vals = AxisArrays.axes(arr, dim).val
    return vals[2] - vals[1]
end

# Half-cell intervals selecting the grid cell containing (`lat`, `long`), each using its own axis
# step — dim 1 = latitude, dim 2 = longitude.
function _cellselectors(arr::AxisArray, lat::typeof(1.0°), long::typeof(1.0°))
    dlat = _axisstep(arr, 1) / 2
    dlong = _axisstep(arr, 2) / 2
    return (lat - dlat) .. (lat + dlat), (long - dlong) .. (long + dlong)
end

# Positions of `requested` third-axis values within `axisvals`. With `strict`, an absent value is an
# error; otherwise it is skipped (used by the `year` window, which returns whatever exists).
function _findslices(axisvals, requested; strict::Bool)
    idxs = Int[]
    for v in requested
        i = findfirst(isequal(v), axisvals)
        if isnothing(i)
            strict && error("third-axis selector $v not found")
        else
            push!(idxs, i)
        end
    end
    return idxs
end

# Resolve the third-axis selection into a tuple of index/selector args to splat (empty for the 2-D
# `Reference`) plus a flag for whether the result collapses to a single value.
function _sliceselector(dat::AbstractClimate, slice, year)
    arr = dat.array
    if ndims(arr) == 2
        (slice === Colon() && isnothing(year)) ||
            error("$(nameof(typeof(dat))) has no third axis to select from")
        return (), true
    end
    axisvals = AxisArrays.axes(arr, 3).val
    if isnothing(year)
        slice === Colon() && return (Colon(),), false
        slice isa AbstractVector &&
            return (_findslices(axisvals, slice; strict = true),), false
        return (slice,), true      # lone scalar: Int band (positional) or Time (value lookup)
    end
    eltype(axisvals) <: Unitful.Time ||
        error("`year` only applies to a dataset with a time axis")
    # Calendar year `Y` spans months `Y·12 … Y·12+11` in the `month` unit (12month == year),
    # matching how ERA/CERA time axes are built. `:` takes the whole year; a month set filters it.
    # Either way only entries that actually exist are returned (none if the year is out of range).
    if slice === Colon()
        return ((year * 12) * month .. (year * 12 + 11) * month,), false
    end
    months = slice isa AbstractVector ? slice : (slice,)
    requested = [(year * 12 + (m - 1)) * month for m in months]
    return (_findslices(axisvals, requested; strict = false),), false
end

"""
    extractvalues(lat, long, dat::AbstractClimate, slice = :; year = nothing)
    extractvalues(lats, longs, dat::AbstractClimate, slice = :; year = nothing)
    extractvalues(loc::LatLong, dat::AbstractClimate, slice = :; year = nothing)
    extractvalues(locs, dat::AbstractClimate, slice = :; year = nothing)

Extract the value(s) at the grid cell(s) containing (`lat`, `long`) from a climate dataset.

A location is a [`LatLong`](@ref) point; pass a single one, a vector of them, or the equivalent
separate `lat`/`long` scalars or `lats`/`longs` vectors (which build the `LatLong`s and forward).
`slice` chooses along the non-spatial third axis — a month (`Unitful.Time`) or a variable/band index
(`Int`), depending on the dataset — as a single value (one slice), a vector (those slices), or `:`
(all slices, the default); it is ignored for the 2-D [`Reference`](@ref). For a dataset with a
calendar time axis, the `year` keyword selects the entries falling within that calendar year —
however many exist, none if the year lies outside the data — and `slice` is then read as month
number(s) within the year.

A single location with a single slice returns a scalar; a single location with several/all slices
returns a vector; `N` locations return an `N×1` (single slice) or `N×M` matrix.
"""
function extractvalues(loc::LatLong{typeof(1.0°)}, dat::AbstractClimate,
                       slice = Colon(); year = nothing)
    sel3, single = _sliceselector(dat, slice, year)
    latsel, longsel = _cellselectors(dat.array, loc.lat, loc.long)
    cell = dat.array[latsel, longsel, sel3...]
    return single ? first(cell) : cell.data[1, 1, :]
end

function extractvalues(locs::AbstractVector{<:LatLong{typeof(1.0°)}},
                       dat::AbstractClimate, slice = Colon(); year = nothing)
    return transpose(hcat(map(loc -> extractvalues(loc, dat, slice;
                                                   year = year),
                              locs)...))
end

function extractvalues(lat::typeof(1.0°), long::typeof(1.0°),
                       dat::AbstractClimate, slice = Colon(); year = nothing)
    return extractvalues(LatLong(lat, long), dat, slice; year = year)
end

function extractvalues(lats::AbstractVector{typeof(1.0°)},
                       longs::AbstractVector{typeof(1.0°)},
                       dat::AbstractClimate, slice = Colon(); year = nothing)
    length(lats) == length(longs) ||
        error("`lats` and `longs` must have the same length")
    return extractvalues(map(LatLong, lats, longs), dat, slice; year = year)
end
