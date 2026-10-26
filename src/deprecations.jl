# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations - main `EcoSISTEM` module
#
# Every deprecated public API is collected here, in sections, and included late in
# `EcoSISTEM.jl`, after all the types it shims. Each shim warns - so downstream code gets a
# migration message rather than a silent `MethodError` - and forwards to the current API. Mirrored
# by `test/test_deprecations.jl`. The shims whose signatures name a `RasterDataSources` type are in
# `ext/EcoSISTEMRasterDataSourcesExt/deprecations.jl`, since they cannot be defined without it.
# There are none at present.
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
