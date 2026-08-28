# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Reading ERA5 and CERA-20C: `read(ERA, ...)`, `read(CERA, ...)`, and the four helpers only they use
# — `_readera_array` (one variable out of one file), `_readeradir` (a directory of them concatenated
# along time), `_parsecfunit` (a CF `units` attribute to a Unitful unit) and `_wraplong180` (ERA
# longitudes run 0 to 360, and every other reader produces -180 to 180).
#
# Separate from `datasetread.jl` because this is the only code in the package that reads netCDF:
# everything else arrives as GeoTIFF through a catalogued data source.

using Unitful

using Unitful.DefaultSymbols

using EcoSISTEM.Units

using Dates: Dates

using Rasters

using DimensionalData

using DimensionalData.Lookups: Intervals, locus, ForwardOrdered, Irregular,
                               Regular

import Rasters: Projected

import Rasters: X, Y, Ti

# Both of these are imported for a SIDE EFFECT and reference no symbol: ArchGDAL registers the GDAL
# backend Rasters uses for GeoTIFFs, NCDatasets the netCDF one these readers need. Neither may be
# removed in an unused-import pass -- the loss only shows up as a data-dependent read failure.
using NCDatasets

import ArchGDAL

import Base.read

"""
    read(::Type{ERA}, file::String, param::String; cut = nothing)

Read one variable out of a single ERA netCDF file, returning a [`ClimateRaster`](@ref)`{ERA}` —
`ERA` names the *source*, and the raster carries the data.

The file's own time coordinate is used as it stands — a genuine `Dates.DateTime` per layer, decoded
from its CF metadata — and its longitudes are wrapped onto `(-180°, 180°]` to match every other
reader. Use the three-argument method to override the time coordinate instead of trusting it, and the
four-argument one to read and concatenate a whole directory.

# Arguments

  - `file`: the netCDF file to read.
  - `param`: the name of the variable to take out of it.
  - `cut`: an optional region to crop to.
"""
function Base.read(::Type{ERA}, file::AbstractString, param::AbstractString;
                   cut = nothing)
    return _timeseriesraster(ERA, _applycut(_readera_array(file, param), cut))
end

"""
    read(::Type{ERA}, file::String, param::String, dim::Vector{<:Unitful.Time}; cut = nothing)

As the two-argument method, but with the file's own time coordinate overridden. Returns a
[`ClimateRaster`](@ref)`{ERA}`.

# Arguments

  - `file`, `param`: as above.
  - `dim`: one time per layer, replacing what the file says. For files whose CF time metadata is
    missing or wrong, or which are to be labelled in elapsed time rather than calendar dates.
  - `cut`: an optional region to crop to.
"""
function Base.read(::Type{ERA}, file::AbstractString, param::AbstractString,
                   dim::Vector{<:Unitful.Time}; cut = nothing)
    aa = _readera_array(file, param)
    world = DimArray(aa.data, (dims(aa, Y), dims(aa, X), Ti(collect(dim))))
    return _timeseriesraster(ERA, _applycut(world, cut))
end

"""
    read(::Type{ERA}, dir::String, file::String, param::String,
         dim::Vector{<:AbstractVector{<:Unitful.Time}}; cut = nothing)

Read one variable out of every matching ERA netCDF file in a directory, concatenating them along
time, and return a [`ClimateRaster`](@ref)`{ERA}`. The multi-file wrapper over the three-argument
method above.

# Arguments

  - `dir`: the directory to search.
  - `file`: matched against each filename; every file containing it is read, in `readdir` order.
  - `param`: the name of the variable to take out of each.
  - `dim`: the time coordinate **per file** — a vector of time vectors, one per matched file.
  - `cut`: an optional region to crop to.
"""
function Base.read(::Type{ERA}, dir::AbstractString, file::AbstractString,
                   param::AbstractString,
                   dim::Vector{<:AbstractVector{<:Unitful.Time}}; cut = nothing)
    return _timeseriesraster(ERA, _readeradir(dir, file, param, dim, cut = cut))
end

"""
    read(::Type{CERA}, dir::String, file::String, param::String; cut = nothing)

Read one variable out of the CERA-20C netCDF files in a directory — the archive is one file per
decade — concatenate them along time, and return a [`ClimateRaster`](@ref)`{CERA}`.

Unlike the `ERA` readers this needs no `dim`: the monthly time coordinate for each decade is
generated here, because the archive's layout is known. Longitudes are wrapped onto `(-180°, 180°]` as
for `ERA`.

# Arguments

  - `dir`: the directory to search.
  - `file`: matched against each filename.
  - `param`: the name of the variable to take out of each.
  - `cut`: an optional region to crop to.
"""
function Base.read(::Type{CERA}, dir::AbstractString, file::AbstractString,
                   param::AbstractString; cut = nothing)
    # one decade-length monthly time vector per file (CERA-20C is archived a decade per file)
    times = [collect((1901year + 1month_mean_duration):(1month_mean_duration):(1910year))]
    for i in 2:12
        push!(times,
              1900year .+ ifelse(i == 12,
                     collect(((i - 1) * 120month_mean_duration + 1month_mean_duration):(1month_mean_duration):((i - 1) * 120month_mean_duration + 1year)),
                     collect(((i - 1) * 120month_mean_duration + 1month_mean_duration):(1month_mean_duration):(i * 10year))))
    end
    return _timeseriesraster(CERA,
                             _readeradir(dir, file, param, times, cut = cut))
end

# Parse a CF `units` attribute, as read off a netCDF file ("J m**-2", "m**3 m**-3"), into a Unitful
# unit. CF spells exponentiation `**` where Unitful spells it `^`, and a bare space is an implicit
# product; a missing or blank attribute means no known unit. The same `unit_context` as
# `LayerCatalogue.jl` uses for the shipped CSV tables, so the two agree on what a unit string means.
function _parsecfunit(s::AbstractString)
    isempty(s) && return NoUnits
    return uparse(replace(s, "**" => "^", " " => "*"),
                  unit_context = [Unitful, Units])
end

# Roll a `(Y, X, Ti)` ERA array from the 0 to 360 degree longitude convention onto (-180, 180],
# reordering the data columns so that longitude stays ascending. Data already in range is unchanged.
#
# `order` and `span` are stated rather than left `Auto`, which is what DimensionalData defaults to for
# a plain `Vector`-backed lookup because it cannot cheaply infer either from an arbitrary vector. An
# `Auto` span here is not resolved lazily later — it is unindexable: `_slicespan` has no method for
# `AutoSpan` at all, so any later indexing of the result, `_applycut`'s crop included, throws a
# `MethodError`. `order` is `ForwardOrdered` because `sortperm` has just made it so, and `span` is
# `Irregular` with explicit edge bounds rather than `Regular` so that nothing assumes perfectly
# uniform spacing survived the wrap.
function _wraplong180(aa::DimArray)
    xdim = dims(aa, X)
    lonvals = parent(DimensionalData.lookup(aa, X))
    wrapped = mod.(ustrip.(lonvals) .+ 180, 360) .- 180
    perm = sortperm(wrapped)
    sortedlon = wrapped[perm] .* °
    lo = length(sortedlon) > 1 ? sortedlon[2] - sortedlon[1] :
         zero(eltype(sortedlon))
    hi = length(sortedlon) > 1 ? sortedlon[end] - sortedlon[end - 1] :
         zero(eltype(sortedlon))
    newx = X(Projected(sortedlon, sampling = Intervals(_locus(xdim)),
                       crs = Rasters.crs(xdim), order = ForwardOrdered(),
                       span = Irregular((first(sortedlon) - lo / 2,
                                         last(sortedlon) + hi / 2))))
    return rebuild(aa, data = aa.data[:, perm, :],
                   dims = (dims(aa, Y), newx, dims(aa, Ti)))
end

# Read one ERA netCDF file into a longitude-wrapped array, shared by both `read(ERA, ...)` methods
# above. NCDatasets, through Rasters, handles the CF coordinate variables, the `units` attribute and
# `scale_factor`/`add_offset`/`_FillValue` on its own, so the time axis arrives as real calendar
# dates and `_rastertodimarray` preserves it rather than rebuilding a synthetic one. The physical
# unit comes from the variable's `units` attribute, through `_parsecfunit`.
#
# `source = NCDsource()` is forced because a file downloaded from the CDS carries no `.nc` extension,
# so Rasters' extension-based backend guess falls back to GDAL — which reads the data but drops both
# the CF coordinates, leaving integer indices, and the `units` attribute.
function _readera_array(file::String, param::String)
    ras = Raster(file, name = Symbol(param), source = Rasters.NCDsource())
    u = _parsecfunit(string(get(Rasters.metadata(ras), "units", "")))
    # ERA5 longitudes run 0–360°; roll them onto (-180, 180] to match the other (tif) readers
    return _wraplong180(_rastertodimarray(ras, unit = u))
end

# Read the files in `dir` matching `file` via the single-file `read(ERA, ...)`, one per time-vector
# in `times`, and concatenate them along time (dim 3). Shared by the directory `ERA`/`CERA` readers.
# Not `dims`, which would shadow `DimensionalData.dims` in a body that also passes `dims` as a
# keyword to `cat`.
function _readeradir(dir::String, file::String, param::String, times;
                     cut = nothing)
    filenames = _searchdir(dir, file)
    arrays = [read(ERA, joinpath(dir, filenames[i]), param, times[i],
                   cut = cut).array
              for i in eachindex(times)]
    return cat(arrays..., dims = 3)
end
