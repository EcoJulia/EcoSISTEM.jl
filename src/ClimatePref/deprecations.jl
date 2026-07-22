# SPDX-License-Identifier: LGPL-3.0-or-later

# ---------------------------------------------------------------------------
# Readers → the unified `read`/`cut` API
# ---------------------------------------------------------------------------
"""
    readworldclim(T::Type{WorldClim{Climate}}, files; cut = nothing)

Deprecated — use `read(WorldClim{Climate}, layers; …)` instead, which downloads via `getraster`
and reads through the same machinery. Retained to read an already-resolved set of monthly climate
raster file paths.
"""
function readworldclim(T::Type{WorldClim{Climate}}, files; cut = nothing)
    Base.depwarn("`readworldclim` is deprecated; use `read(WorldClim{Climate}, layers; …)`.",
                 :readworldclim)
    return _readsource(T, _filelist(files); cut = cut)
end

# Deprecated positional-extent form `readfile(file, xmin, xmax, ymin, ymax)`; convert to `cut` and forward
# to the keyword `readfile(file; cut)` method in `ReadData.jl`.
function readfile(file::String, xmin, xmax, ymin, ymax)
    Base.depwarn("`readfile(file, xmin, xmax, ymin, ymax)` is deprecated; pass " *
                 "`cut = LatLong(ymin .. ymax, xmin .. xmax)` (e.g. from `boundingbox`) " *
                 "instead.", :readfile)
    return readfile(file; cut = LatLong(ymin .. ymax, xmin .. xmax))
end
