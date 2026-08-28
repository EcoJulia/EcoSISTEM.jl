# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Land cover ------------------------------------------------------------------------------------
#
# The two operations keyed on `EarthEnv{LandCover}`. Both are documented in the parent, on
# method-less stubs in `src/extensions.jl`.
#
# `landcoverclass` is the only **public, non-deprecated** name in the package whose implementation
# is here, so it is where a user meets this dependency without having asked for climate data:
# `SetLandCover(:open_water)` reaches it, while `SetLandCover(7)` takes the untyped path and does not.

function EcoSISTEM.compress_landcover(landcover::ClimateRaster{T}) where
    {T <: EarthEnv{<:LandCover}}
    array = landcover.array
    yx = dims(array, (Y, X))
    codes = Array{Int64}(undef, size(array, 1), size(array, 2))
    Threads.@threads for i in Base.axes(array, 1)
        for j in Base.axes(array, 2)
            codes[i, j] = argmax(view(array, i, j, :))
        end
    end

    # Deliberately **no** `code`: the result is not any one of the bands it came from, so there is
    # nothing honest to put there.
    # And **`DerivedData{T}`, not `T`**, for the same reason: a dominant-class layer is computed
    # *from* land cover and is not land cover, so claiming the dataset would be a provenance claim
    # the values do not support. The lineage is carried by the parameter instead.
    #
    # **Nothing here declares that the result is categorical, and nothing should.** That is what
    # the spec's `axis = LandCoverTypology` says — a `TypologyAxis`, so class codes — and it says it
    # where the caller can see it. Stamping `valuetype = :categorical` on the raster here would be
    # working around the shipped table having no row for a derived layer; the axis has no such gap.
    return ClimateRaster(_derivedfrom(T), DimArray(codes, yx))
end

function EcoSISTEM.landcoverclass(name::Symbol)
    aliases = layerinfo(EarthEnv{LandCover}, name).aliases
    return parse(Int, only(filter(a -> all(isdigit, a), aliases)))
end
