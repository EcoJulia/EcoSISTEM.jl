# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations - main `EcoSISTEM` module
#
# Every deprecated public API is collected here, in sections, and included late in
# `EcoSISTEM.jl`, after all the types it shims. Each shim warns - so downstream code gets a
# migration message rather than a silent `MethodError` - and forwards to the current API. Mirrored
# by `test/test_deprecations.jl`. The shims whose signatures name a `RasterDataSources` type are in
# `ext/EcoSISTEMRasterDataSourcesExt/deprecations.jl`, since they cannot be defined without it;
# it holds none at present.
#
# ## Every section says which release deprecated it
#
# Each section header ends with a `Deprecated in vX.Y.Z` line saying which release the shim belongs
# to, and **the sections are ordered newest release first**, so the oldest deprecations sit at the
# foot of the file and dropping a release's worth is deleting from the bottom up to the first
# section carrying a newer label. `clean_Deprecations.jl` asserts every section has a label naming
# a release it knows and that the labels never get newer going down the file.
#
# Two things a whole-section delete does not reach, and both fail quietly:
#
#   - an `export` or `public` declaration for a deprecated name that sits anywhere but beside its
#     shim. A stale `export` of a deleted binding is a precompile *warning*, not an error.
#   - `test/test_deprecations.jl`, which mirrors this file and must lose the same sections.
# ===========================================================================

# ---------------------------------------------------------------------------
# The kernels take their distances by keyword, each named for what it is, since the positional
# distance is a Gaussian kernel's mean but neither the mean nor the scale of a long-tailed one. Each
# shim must build the kernel its positional form denotes to the bit, or every lookup and canonical
# reference built from one moves: the long-tailed distance is the scale `(2 * dist) / sqrt(pi)`.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
@deprecate(GaussianKernel(dist::Unitful.Length{Float64}, thresh::Float64),
           GaussianKernel(meandist = dist, thresh = thresh), false)
@deprecate(LongTailKernel(dist::Unitful.Length{Float64}, shape::Float64,
                          thresh::Float64),
           LongTailKernel(scale = (2 * uconvert(km, dist)) / sqrt(pi),
                          shape = shape, thresh = thresh), false)

# ---------------------------------------------------------------------------
# A species list's movement is one `SpeciesDispersal` per species, a dispersal mode with the species'
# kernel, where it was one object for the whole list. Each list-wide movement type is a function
# returning that vector, in the mode it named, each kernel's type its family and each species' flag
# its dead-cell fate; given to `build_species` as `movement`, each name stands for its mode.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
function BirthOnlyMovement(kernels::AbstractVector{<:AbstractKernel},
                           disperse_safely::AbstractVector{Bool} = fill(true,
                                                                        length(kernels)))
    Base.depwarn("`BirthOnlyMovement(kernels, disperse_safely)` is deprecated: give each " *
                 "species its movement, as `SeedDispersal(family; disperse_safely).(kernels)`.",
                 :BirthOnlyMovement)
    return _listmovement(SeedDispersal, kernels, disperse_safely)
end

function AlwaysMovement(kernels::AbstractVector{<:AbstractKernel},
                        disperse_safely::AbstractVector{Bool} = fill(true,
                                                                     length(kernels)))
    Base.depwarn("`AlwaysMovement(kernels, disperse_safely)` is deprecated: give each " *
                 "species its movement, as `AlwaysDispersal(family; disperse_safely).(kernels)`.",
                 :AlwaysMovement)
    return _listmovement(AlwaysDispersal, kernels, disperse_safely)
end

function NoMovement(kernels::AbstractVector{<:AbstractKernel},
                    ::AbstractVector{Bool} = fill(true, length(kernels)))
    Base.depwarn("`NoMovement(kernels)` is deprecated: give each species its movement, " *
                 "`NoDispersal()()`.", :NoMovement)
    return fill(NoDispersal()(), length(kernels))
end

# A list-wide movement's species, each in `mode` with its own kernel and dead-cell fate.
function _listmovement(mode, kernels, disperse_safely)
    length(disperse_safely) == length(kernels) ||
        error("`disperse_safely` has $(length(disperse_safely)) entries but there are " *
              "$(length(kernels)) dispersal kernels: it is a per-species setting, so it needs one " *
              "entry per species.")
    return [SpeciesDispersal(mode(typeof(k), disperse_safely = s), k)
            for (k, s) in zip(kernels, disperse_safely)]
end

function _speciesmovement(::typeof(BirthOnlyMovement), args...)
    Base.depwarn("`movement = BirthOnlyMovement` is deprecated: give `movement = " *
                 "InsectDispersal()`, which disperses the same.",
                 :build_species)
    return _speciesmovement(InsectDispersal(), args...)
end

function _speciesmovement(::typeof(AlwaysMovement), args...)
    Base.depwarn("`movement = AlwaysMovement` is deprecated: give `movement = AlwaysDispersal()`.",
                 :build_species)
    return _speciesmovement(AlwaysDispersal(), args...)
end

function _speciesmovement(::typeof(NoMovement), args...)
    Base.depwarn("`movement = NoMovement` is deprecated: give `movement = NoDispersal()`.",
                 :build_species)
    return _speciesmovement(NoDispersal(), args...)
end

export BirthOnlyMovement, AlwaysMovement, NoMovement

# ---------------------------------------------------------------------------
# `InputRecord` is `ProvenanceRecord`: a record describes a file however it enters a run, an
# output of the run that made it included, so its name no longer calls it an input. The binding
# stays `public`, as the old name was.
#
# Deprecated in v0.10.0.
# ---------------------------------------------------------------------------
Base.@deprecate_binding InputRecord ProvenanceRecord false
public InputRecord
