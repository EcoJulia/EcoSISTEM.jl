# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations — `ClimatePref` submodule
#
# (One live deprecation cannot be moved here: the `xmin`/`xmax`/`ymin`/`ymax`
# keyword branch inside `readfile(file; …)` in `src/datasetread.jl` is part of that
# live method's body.)
# ===========================================================================

# Declared here, where it is used. A `using` is module-scoped, so an import stated in a neighbouring
# file would satisfy this one invisibly and break the moment that file changed.
using Statistics
using DimensionalData: rebuild, lookup

# **Three groups of shims moved to `EcoSISTEMRasterDataSourcesExt`** — the five per-source wrapper
# constructors (`Worldclim_bioclim` and friends, which forward to a `ClimateRaster{WorldClim{…}}`),
# `readworldclim`, and `readCHELSA_monthly`. Each names a `RasterDataSources` type in its signature
# or its target, so none can be defined without that package; the names stay exported by the
# submodule regardless.
# Everything below names no dataset and stays: the ERA/CERA readers, the four-argument `readfile`,
# `readCRUTS` (which forwards to a moved method but names nothing itself), and the two
# naming-standardisation redirects.

# ---------------------------------------------------------------------------
# Names whose methods live in `EcoSISTEMRasterDataSourcesExt`
#
# Declared here, method-less, for three reasons: the submodule already `export`s them, so the
# binding must exist; `@deprecate` in the extension resolves `ClimatePref.<name>` at macro-expansion
# time and needs something to attach to; and a docstring on a name whose only method is in an
# extension is invisible to `@autodocs` unless it stays in the parent.
# Calling one without `RasterDataSources` loaded is a `MethodError` naming the function — the same
# behaviour as `retrieve_era5` without `PyCall`.
# ---------------------------------------------------------------------------
function Worldclim_bioclim end
function CHELSA_bioclim end
function Landcover end
function Worldclim_monthly end
function CHELSA_monthly end

"""
    readworldclim(T::Type{WorldClim{Climate}}, files; cut = nothing)

Deprecated — use `read(WorldClim{Climate}, layers; …)` instead, which downloads via `getraster`
and reads through the same machinery. Retained to read an already-resolved set of monthly climate
raster file paths.
"""
function readworldclim end

"""
    readCHELSA_monthly(dir::String, var_name::String; scale = 1, fn = mean, cut = nothing)

Deprecated — use `read(CHELSA{Climate}, dir, var_name; scale, fn, cut)` instead.
"""
function readCHELSA_monthly end

# ---------------------------------------------------------------------------
# Readers → the unified `read`/`cut` API
# ---------------------------------------------------------------------------
# Deprecated positional-extent form `readfile(file, xmin, xmax, ymin, ymax)`; convert to `cut` and forward
# to the keyword `readfile(file; cut)` method in `src/datasetread.jl`.
function readfile(file::String, xmin, xmax, ymin, ymax)
    Base.depwarn("`readfile(file, xmin, xmax, ymin, ymax)` is deprecated; pass " *
                 "`cut = Extent(Y = (ymin, ymax), X = (xmin, xmax))` (e.g. from `boundingbox`) " *
                 "instead.", :readfile)
    return readfile(file,
                    cut = Extents.Extent(Y = (ymin, ymax), X = (xmin, xmax)))
end

# ---------------------------------------------------------------------------
# Source-named readers → `read(T, ...)` dispatch on the result type
#
# `readCRUTS`/`readCHELSA_monthly`/`readERA`/`readCERA` named themselves after their source instead
# of using multiple dispatch — the same anti-pattern `_thirdaxis` had (fixed in `datasetread.jl`).
# `CRUTS`/`ERA`/`CERA` are already real result types (`src/Climate.jl`), reused here as the `read`
# dispatch tag (`read(CRUTS, dir, var_name)`, `read(ERA, file, param)`) — Julia doesn't distinguish
# "source type" from "result type" for dispatch, the same `parse(Float64, "3.14")`-shaped idiom.
# `readCHELSA_monthly` dispatches on the existing `CHELSA{Climate}` source type instead, matching
# `Base.read(T::Type{<:RDS.RasterDataSource}, ...)` in the same file. `readfile` is *not* included —
# it has no result-wrapper type to hang a dispatch tag on, and shadowing `Base.read(::AbstractString)`
# (raw-bytes read) with a type-punned return value would be confusing, not clarifying.
# ---------------------------------------------------------------------------
"""
    readCRUTS(dir::String, var_name::String; cut = nothing)

Deprecated — use `read(CRUTS, dir, var_name; cut)` instead.
"""
function readCRUTS(dir::String, var_name::String; cut = nothing)
    Base.depwarn("`readCRUTS` is deprecated; use `read(CRUTS, dir, var_name; cut)`.",
                 :readCRUTS)
    return read(CRUTS, dir, var_name, cut = cut)
end

"""
    readERA(file::String, param::String; cut = nothing)

Deprecated — use `read(ERA, file, param; cut)` instead.
"""
function readERA(file::String, param::String; cut = nothing)
    Base.depwarn("`readERA(file, param; cut)` is deprecated; use `read(ERA, file, param; cut)`.",
                 :readERA)
    return read(ERA, file, param, cut = cut)
end

"""
    readERA(file::String, param::String, dim::Vector{<:Unitful.Time}; cut = nothing)

Deprecated — use `read(ERA, file, param, dim; cut)` instead.
"""
function readERA(file::String, param::String, dim::Vector{<:Unitful.Time};
                 cut = nothing)
    Base.depwarn("`readERA(file, param, dim; cut)` is deprecated; use " *
                 "`read(ERA, file, param, dim; cut)`.", :readERA)
    return read(ERA, file, param, dim, cut = cut)
end

"""
    readERA(dir::String, file::String, param::String,
            dim::Vector{<:AbstractVector{<:Unitful.Time}}; cut = nothing)

Deprecated — use `read(ERA, dir, file, param, dim; cut)` instead.
"""
function readERA(dir::String, file::String, param::String,
                 dim::Vector{<:AbstractVector{<:Unitful.Time}}; cut = nothing)
    Base.depwarn("`readERA(dir, file, param, dim; cut)` is deprecated; use " *
                 "`read(ERA, dir, file, param, dim; cut)`.", :readERA)
    return read(ERA, dir, file, param, dim, cut = cut)
end

"""
    readCERA(dir::String, file::String, param::String; cut = nothing)

Deprecated — use `read(CERA, dir, file, param; cut)` instead.
"""
function readCERA(dir::String, file::String, param::String; cut = nothing)
    Base.depwarn("`readCERA` is deprecated; use `read(CERA, dir, file, param; cut)`.",
                 :readCERA)
    return read(CERA, dir, file, param, cut = cut)
end

# ---------------------------------------------------------------------------
# `compressLC` → `compress_landcover`: the `LC` abbreviation is expanded, and the name is `snake_case`
# like every other function. The released name is `compressLC`; the intermediate
# `compressLandCover` never shipped, so this points straight at the final name and owes no second shim.
# ---------------------------------------------------------------------------
@deprecate compressLC compress_landcover

# ---------------------------------------------------------------------------
# Naming standardisation (v0.5.0): `camelCase` phylogenetic-fit functions → `snake_case`.
# ---------------------------------------------------------------------------
@deprecate fitBrownian fitbrownian

# ---------------------------------------------------------------------------
# `fitLambda`/`fitlambda` DELETED (v0.5.0) — not deprecated, and no shim
# ---------------------------------------------------------------------------
# Pagel's-lambda fitting is gone with no redirect, because it never worked on Julia 1.x and so no
# user code can depend on a result it never returned. Two independent faults: `varcovar` threw
# `UndefVarError` on every call from Julia 0.7 (it used `indmax`, removed in 0.7), and `fitLambda`
# itself threw `ArgumentError: matrix contains Infs or NaNs` — it optimised **with bounds** on
# lambda, so `Optim` used `Fminbox`, which needs a gradient and takes one by finite differences,
# while the objective's `log(abs(det(x[1] * V)))` underflows for a scaled n x n covariance
# (`det(0.1 * V)` carries a factor `0.1^n`).
#
# `fitbrownian` shares that expression and survives **only** because its optimisation is unbounded:
# Nelder-Mead needs no gradient, so it never evaluates the objective where it blows up.
# A working lambda model needs a numerically stable log-determinant (`logdet`, or a Cholesky), which
# is a modelling decision rather than a repair — see `NEWS.md`.
# ---------------------------------------------------------------------------
# `extractvalues` → `extract_values` (v0.5.0): `snake_case` like every other public function, and a
# keyword signature in place of four positional methods.
#
# **A hand-written shim rather than `@deprecate`, because the *signature* changed too.** A bare
# `@deprecate extractvalues extract_values` forwards the arguments, so the old positional call would
# warn and then die on a `MethodError` that says nothing about how to fix it. It cannot forward
# honestly either: `slice` meant a month, a band index or a layer code depending on the dataset, and
# `year` needs an epoch to become a date range — the ambiguity the keywords exist to remove is
# exactly what makes the translation impossible to do for the caller. So it names the replacement
# for each old argument and stops.
# ---------------------------------------------------------------------------
"""
    extractvalues(args...; kwargs...)

Deprecated — use [`extract_values`](@ref), which takes the dataset positionally and everything else
by keyword.
"""
function extractvalues(args...; kwargs...)
    Base.depwarn("`extractvalues` is deprecated; use `extract_values`.",
                 :extractvalues)
    return error("`extractvalues(lat, long, dat, slice; year)` has been replaced by " *
                 "`extract_values(dat; lat, long, date, month, offset, code)`, and the arguments " *
                 "do not map across automatically. `lat`/`long` keep their meaning; `slice` " *
                 "becomes `month` (a monthly climatology), `offset` (an elapsed-time axis) or " *
                 "`code` (a layer), depending on which axis the dataset has; and `year = 2000` " *
                 "becomes `date = Date(2000, 1, 1) .. Date(2000, 12, 31)`.")
end
export extractvalues

# ---------------------------------------------------------------------------
# **`upresolution` / `downresolution` / `downresolution!` — DEPRECATED, and moved here whole.**
#
# **Nothing in the package calls these.** The only caller is the testset that exists to test them:
# read-time coarsening is done by block aggregation in `datasetread.jl`, and putting two layers onto
# one grid is what a `StudyArea` is for.
#
# **Rasters.jl has no equivalent**, measured on a 6 by 6 grid of `10i + j`:
# - `downresolution(x, 2)` and `Rasters.aggregate(mean, x, 2)` give the **same size** `(3, 3)` and
# **different values** — `[11.0, 13.0, 15.0]` against `[16.5, 18.5, 20.5]`. This takes an
# *overlapping* window centred on every `rescale`th cell, collapsing to a single cell at the
# edges; `aggregate` takes *disjoint tiled blocks*. Same shape, different meaning — the
# dangerous kind of near-miss, and why no forwarding shim is possible.
# - `upresolution(x, 2)` gives **(11, 11)** against `disaggregate`'s **(12, 12)**: this is
# **bilinear interpolation between cell centres** (`(n - 1) * rescale + 1` points), not the
# block replication `disaggregate` performs (`n * rescale`).
# Note the docstring below claims the values "are repeated". That is **wrong** — they are
# interpolated. Left uncorrected deliberately: it is the released text of a deprecated function.
#
# The implementations stay here, working, until the removal decision is taken. Each warns in its
# innermost method, so a call through an `ERA`/`ClimateRaster` wrapper warns once, not twice.
# ---------------------------------------------------------------------------

"""
    upresolution(data::Union{ERA, ClimateRaster}, rescale::Int64; fn)

Function to increase the resolution of a climate dataset, by a factor, `rescale`.

# Arguments

  - `data`: the dataset to refine — an [`ERA`](@ref) or a [`ClimateRaster`](@ref).
  - `rescale`: how many output cells each input cell becomes along **each** axis, so `2` gives four
    times as many cells.
  - `fn`: how to fill them. Splitting a cell invents no information: the values are repeated,
    so this is a change of *grid*, not of detail.
"""
function upresolution end

function upresolution(era::ClimateRaster{ERA}, rescale::Int64)
    array = upresolution(era.array, rescale)
    return ClimateRaster(ERA, array)
end

function upresolution(bioclim::ClimateRaster{T}, rescale::Int64) where {T}
    array = upresolution(bioclim.array, rescale)
    return ClimateRaster(T, array)
end

function upresolution(aa::DimensionalData.AbstractDimArray{T, 3} where {T},
                      rescale::Int64)
    Base.depwarn("`upresolution` is deprecated and unused by the package. " *
                 "`Rasters.disaggregate` is the nearest replacement but is NOT " *
                 "equivalent: this interpolates bilinearly between cell centres, " *
                 "giving `(n - 1) * rescale + 1` cells, where `disaggregate` " *
                 "replicates each cell into a block, giving `n * rescale`.",
                 :upresolution)
    grid = size(aa)
    grid = (grid[1] .* rescale .- (rescale - 1),
            grid[2] .* rescale .- (rescale - 1), grid[3])
    array = Array{typeof(aa[1]), 3}(undef, grid)
    aa_ax = Base.axes(aa)
    map(aa_ax[3]) do t
        for x in 1:(size(aa, 1) - 1)
            for y in 1:(size(aa, 2) - 1)
                for dx in 0:rescale
                    fx = dx / rescale
                    ix = rescale * x - (rescale - 1)
                    for dy in 0:rescale
                        fy = dy / rescale
                        iy = rescale * y - (rescale - 1)
                        array[ix + dx, iy + dy, t] = aa[x, y, t] * (1 - fx) *
                                                     (1 - fy) +
                                                     aa[x + 1, y, t] * fx *
                                                     (1 - fy) +
                                                     aa[x, y + 1, t] *
                                                     (1 - fx) * fy +
                                                     aa[x + 1, y + 1, t] * fx *
                                                     fy
                    end
                end
            end
        end
    end

    # Carry the input's axis names through unchanged (rather than hardcoding a lat/long/time
    # convention), so the resampled result keeps whatever axes it came in with.
    d1, d2 = dims(aa, 1), dims(aa, 2)
    v1, v2 = parent(lookup(d1)), parent(lookup(d2))
    return DimArray(array,
                    (rebuild(d1, range(v1[1], v1[end], grid[1])),
                     rebuild(d2, range(v2[1], v2[end], grid[2])),
                     dims(aa, 3)))
end

function upresolution(aa::DimensionalData.AbstractDimArray{T, 2} where {T},
                      rescale::Int64)
    Base.depwarn("`upresolution` is deprecated and unused by the package. " *
                 "`Rasters.disaggregate` is the nearest replacement but is NOT " *
                 "equivalent: this interpolates bilinearly between cell centres, " *
                 "giving `(n - 1) * rescale + 1` cells, where `disaggregate` " *
                 "replicates each cell into a block, giving `n * rescale`.",
                 :upresolution)
    grid = size(aa) .* rescale .- (rescale - 1)
    array = Matrix{typeof(aa[1])}(undef, grid)
    for x in 1:(size(aa, 1) - 1)
        for y in 1:(size(aa, 2) - 1)
            for dx in 0:rescale
                fx = dx / rescale
                ix = rescale * x - (rescale - 1)
                for dy in 0:rescale
                    fy = dy / rescale
                    iy = rescale * y - (rescale - 1)
                    array[ix + dx, iy + dy] = aa[x, y] * (1 - fx) * (1 - fy) +
                                              aa[x + 1, y] * fx * (1 - fy) +
                                              aa[x, y + 1] * (1 - fx) * fy +
                                              aa[x + 1, y + 1] * fx * fy
                end
            end
        end
    end

    # Carry the input's axis names through unchanged (see the 3-D method).
    d1, d2 = dims(aa, 1), dims(aa, 2)
    v1, v2 = parent(lookup(d1)), parent(lookup(d2))
    return DimArray(array,
                    (rebuild(d1, range(v1[1], v1[end], grid[1])),
                     rebuild(d2, range(v2[1], v2[end], grid[2]))))
end

"""
    downresolution(data::Union{ERA, ClimateRaster{WorldClim{BioClim}, <: DimensionalData.AbstractDimArray}}, rescale::Int64; fn)

Function to decrease the resolution of a climate dataset, by a factor, `rescale`, and aggregation function, `fn`. The aggregation function has a default setting of taking the mean value.

# Arguments

  - `data`: the dataset to coarsen — an [`ERA`](@ref) or a [`ClimateRaster`](@ref).
  - `rescale`: how many input cells are combined along **each** axis, so `2` combines blocks of
    four.
  - `fn`: how to combine them, `mean` by default. A mean is meaningless for **class codes** —
    use a nearest-class rule there, which is what the resampling path does for a layer on a
    a `TypologyAxis`.
"""
function downresolution end

function downresolution(era::ClimateRaster{ERA}, rescale::Int64;
                        fn::Function = mean)
    array = downresolution(era.array, rescale, fn = fn)
    return ClimateRaster(ERA, array)
end

function downresolution(bioclim::ClimateRaster{T}, rescale::Int64;
                        fn::Function = mean) where {T}
    array = downresolution(bioclim.array, rescale, fn = fn)
    return ClimateRaster(T, array)
end

function downresolution(aa::DimensionalData.AbstractDimArray{T, 3} where {T},
                        rescale::Int64;
                        fn::Function = mean)
    Base.depwarn("`downresolution` is deprecated and unused by the package. " *
                 "`Rasters.aggregate` is the nearest replacement but is NOT " *
                 "equivalent: this aggregates an *overlapping* window centred on " *
                 "every `rescale`th cell (collapsing at the edges), where " *
                 "`aggregate` takes disjoint tiled blocks.", :downresolution)
    grid = size(aa)
    grid = ceil.(Int64, (grid[1] / rescale, grid[2] / rescale, grid[3]))
    array = Array{typeof(aa[1]), 3}(undef, grid)
    map(1:grid[3]) do t
        for x in Base.axes(array, 1)
            for y in Base.axes(array, 2)
                xcoords = filter(k -> 1 ≤ rescale * (x - 1) + 1 + k ≤
                                      size(aa, 1),
                                 round(Int, -rescale / 2):round(Int,
                                                                rescale / 2))
                ycoords = filter(k -> 1 ≤ rescale * (y - 1) + 1 + k ≤
                                      size(aa, 2),
                                 round(Int, -rescale / 2):round(Int,
                                                                rescale / 2))
                xrange = min(-minimum(xcoords), maximum(xcoords))
                yrange = min(-minimum(ycoords), maximum(ycoords))
                xcoords = (rescale * (x - 1) + 1) .+ ((-xrange):xrange)
                ycoords = (rescale * (y - 1) + 1) .+ ((-yrange):yrange)
                array[x, y, t] = fn(filter(!isnan, aa[xcoords, ycoords, t]))
            end
        end
    end
    # Carry the input's axis names through unchanged (rather than hardcoding a lat/long/time
    # convention), so the resampled result keeps whatever axes it came in with.
    d1, d2 = dims(aa, 1), dims(aa, 2)
    v1, v2 = parent(lookup(d1)), parent(lookup(d2))
    return DimArray(array,
                    (rebuild(d1, range(v1[1], v1[end], grid[1])),
                     rebuild(d2, range(v2[1], v2[end], grid[2])),
                     dims(aa, 3)))
end

function downresolution(aa::DimensionalData.AbstractDimArray{T, 2} where {T},
                        rescale::Int64;
                        fn::Function = mean)
    Base.depwarn("`downresolution` is deprecated and unused by the package. " *
                 "`Rasters.aggregate` is the nearest replacement but is NOT " *
                 "equivalent: this aggregates an *overlapping* window centred on " *
                 "every `rescale`th cell (collapsing at the edges), where " *
                 "`aggregate` takes disjoint tiled blocks.", :downresolution)
    grid = size(aa)
    grid = ceil.(Int64, (grid[1] / rescale, grid[2] / rescale))
    array = Matrix{typeof(aa[1])}(undef, grid)
    for x in Base.axes(array, 1)
        for y in Base.axes(array, 2)
            xcoords = filter(k -> 1 ≤ rescale * (x - 1) + 1 + k ≤ size(aa, 1),
                             round(Int, -rescale / 2):round(Int, rescale / 2))
            ycoords = filter(k -> 1 ≤ rescale * (y - 1) + 1 + k ≤ size(aa, 2),
                             round(Int, -rescale / 2):round(Int, rescale / 2))
            xrange = min(-minimum(xcoords), maximum(xcoords))
            yrange = min(-minimum(ycoords), maximum(ycoords))
            xcoords = (rescale * (x - 1) + 1) .+ ((-xrange):xrange)
            ycoords = (rescale * (y - 1) + 1) .+ ((-yrange):yrange)
            array[x, y] = fn(filter(!isnan, aa[xcoords, ycoords]))
        end
    end

    # Carry the input's axis names through unchanged (see the 3-D method).
    d1, d2 = dims(aa, 1), dims(aa, 2)
    v1, v2 = parent(lookup(d1)), parent(lookup(d2))
    return DimArray(array,
                    (rebuild(d1, range(v1[1], v1[end], grid[1])),
                     rebuild(d2, range(v2[1], v2[end], grid[2]))))
end

"""
    downresolution!(resized_array::Matrix{T}, array::Matrix{T}, rescale::Int64, fn)
    downresolution!(resized_array::Array{T, 3}, array::Matrix{T}, dim::Int64, rescale::Int64, fn)

Function to decrease the resolution of a climate dataset in place, by a factor, `rescale`, and aggregation function, `fn`. The aggregation function has a default setting of taking the mean value.

# Arguments

  - `resized_array`: the destination, written in place — its size is what fixes the output grid.
  - `array`: the source values, one 2-D slice.
  - `dim`: which slice of a 3-D destination to write, for the second form only.
  - `rescale`: how many input cells are combined along each axis.
  - `fn`: how to combine them, applied to the non-`NaN` values of each block.
"""
function downresolution! end

function downresolution!(resized_array::Matrix{T}, array::Matrix{T},
                         dim::Int64, rescale::Int64;
                         fn::Function = mean) where {T}
    Base.depwarn("`downresolution!` is deprecated and unused by the package; " *
                 "see `downresolution`.", :downresolution!)
    dim == 1 || error("Accessing invalid 3rd dimension of 2d array")
    return downresolution!(resized_array, array, rescale, fn = fn)
end

function downresolution!(resized_array::Matrix{T}, array::Matrix{T},
                         rescale::Int64; fn::Function = mean) where {T}
    Base.depwarn("`downresolution!` is deprecated and unused by the package; " *
                 "see `downresolution`.", :downresolution!)
    Threads.@threads for i in eachindex(resized_array)
        (y, x) = EcoSISTEM.convert_coords(i, size(resized_array, 1))
        ycoords = filter(k -> 1 ≤ rescale * (y - 1) + 1 + k ≤ size(array, 1),
                         round(Int, -rescale / 2):round(Int, rescale / 2))
        xcoords = filter(k -> 1 ≤ rescale * (x - 1) + 1 + k ≤ size(array, 2),
                         round(Int, -rescale / 2):round(Int, rescale / 2))
        yrange = min(-minimum(ycoords), maximum(ycoords))
        xrange = min(-minimum(xcoords), maximum(xcoords))
        ycoords = (rescale * (y - 1) + 1) .+ ((-yrange):yrange)
        xcoords = (rescale * (x - 1) + 1) .+ ((-xrange):xrange)
        resized_array[y, x] = fn(filter(!isnan, array[ycoords, xcoords]))
    end
end

function downresolution!(resized_array::Array{T, 3}, array::Matrix{T},
                         dim::Int64, rescale::Int64;
                         fn::Function = mean) where {T}
    Base.depwarn("`downresolution!` is deprecated and unused by the package; " *
                 "see `downresolution`.", :downresolution!)
    new_dims = size(resized_array, 1) * size(resized_array, 2)
    Threads.@threads for i in 1:new_dims
        (y, x) = EcoSISTEM.convert_coords(i, size(resized_array, 1))
        ycoords = filter(k -> k .<= size(array, 1),
                         (rescale * y - (rescale - 1)):(rescale * y))
        xcoords = filter(k -> k .<= size(array, 2),
                         (rescale * x - (rescale - 1)):(rescale * x))
        resized_array[y, x, dim] = fn(filter(!isnan, array[ycoords, xcoords]))
    end
end
export upresolution, downresolution, downresolution!
