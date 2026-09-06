# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Land cover ------------------------------------------------------------------------------------
#
# The two operations keyed on `EarthEnv{LandCover}`. Both are documented in the parent, on
# method-less stubs in `src/extensions.jl`.
#
# `landcoverclass` is the only **public, non-deprecated** name in the package whose implementation
# is here, so it is where a user meets this dependency without having asked for climate data:
# `SetLandCover(:open_water)` reaches it, while `SetLandCover(7)` takes the untyped path and does not.

# `dominant_class` on the twelve bands: the code is the band's position, since EarthEnv's bands are
# labelled `1:12` in class order, and a cell with no data in any band comes out absent rather than
# as the first class. The result carries `DerivedData{T}` and no code, as any dominant class does.
function EcoSISTEM.compress_landcover(landcover::ClimateRaster{T}) where
    {T <: EarthEnv{<:LandCover}}
    return EcoSISTEM.dominant_class(landcover)
end

function EcoSISTEM.landcoverclass(name::Symbol)
    aliases = layerinfo(EarthEnv{LandCover}, name).aliases
    return parse(Int, only(filter(a -> all(isdigit, a), aliases)))
end
