# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A file downloaded once and cached outside the repository, under `EcoSISTEM.assetdir`.

"""
    CachedAsset(owner::Type, url::AbstractString)

An immutable descriptor for a file downloaded from `url` and cached under
`EcoSISTEM.assetdir(owner = owner)`. Nothing is downloaded at construction — call
[`assetpath`](@ref) to get the local path, downloading it into the cache first if it is not there.

# Fields

  - `owner`: the type the download belongs to, which names its own subdirectory of the cache so that
    one type's downloads cannot collide with another's.
  - `url`: where to fetch it from.
"""
struct CachedAsset
    owner::Type
    url::String
    CachedAsset(owner::Type, url::AbstractString) = new(owner, String(url))
end

"""
    EcoSISTEM.assetdir(mod::Module = EcoSISTEM; owner::Union{Type, Nothing} = nothing)

Path to a subdirectory of EcoSISTEM's Scratch.jl space, for storing downloaded data outside the
repository. Creating the directory is part of answering.

One EcoSISTEM-owned space with a subdirectory per package, rather than a space each, keeps the whole
cache under one lifecycle: created on first use, and reclaimed by `Pkg.gc()` when EcoSISTEM is
removed. `RasterDataSources` is put under it the same way, by the extension's `__init__` setting
`RASTERDATASOURCES_PATH`.

# Arguments

  - `mod`: whose subdirectory to return, defaulting to EcoSISTEM's own.
  - `owner`: a type to nest a further subdirectory under, so each owner's downloads are separate —
    see [`CachedAsset`](@ref).
"""
function assetdir(mod::Module = EcoSISTEM;
                  owner::Union{Type, Nothing} = nothing)
    dir = joinpath(get_scratch!(EcoSISTEM, "assets"), string(nameof(mod)))
    isnothing(owner) || (dir = joinpath(dir, string(nameof(owner))))
    return mkpath(dir)
end

"""
    assetpath(asset::CachedAsset)
    assetpath(relative_path::AbstractString)

Return a local path in EcoSISTEM's asset cache, downloading the file first if it is not there yet.

# Arguments

  - `asset`: a [`CachedAsset`](@ref), fetched from its `url` into its owner's subdirectory if the
    cache does not already hold it.
  - `relative_path`: a path resolved under [`assetdir`](@ref) directly, for a cached file whose
    name is already known and which needs no descriptor.
"""
function assetpath(asset::CachedAsset)
    path = joinpath(assetdir(owner = asset.owner), basename(asset.url))
    isfile(path) && return path
    # Downloaded to a sibling temporary file and renamed in, never written straight to `path`.
    # `isfile(path)` is the only cache check, so a half-written file there would be served as valid
    # forever, failing with an opaque reader error that re-running never clears. `Downloads.download`
    # already unwinds cleanly on a raised error, Ctrl-C included; the rename closes the one window it
    # cannot, a process killed outright part-way through. A rename within a directory is atomic, and
    # concurrent downloads get their own `tempname`, so `force` means last-writer-wins on identical
    # content.
    temp = tempname(dirname(path))
    try
        Downloads.download(asset.url, temp)
        mv(temp, path, force = true)
    finally
        rm(temp, force = true)          # a no-op once the rename has succeeded
    end
    return path
end
assetpath(relative_path::AbstractString) = joinpath(assetdir(), relative_path)

# Search a directory `path` for the files whose names contain `key` -- a plain `occursin` test, not a
# glob or a regular expression. It knows nothing about climate data: both callers, `clearcache!` and
# `gettimes`, are listing a cache folder, which is why it sits beside the other filesystem helpers
# and is private.
_searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))
