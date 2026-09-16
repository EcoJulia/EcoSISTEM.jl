# SPDX-License-Identifier: LGPL-3.0-or-later

# Field-by-field equality of built layers, shared by `test_layerpayload.jl` and `SmallMPItest.jl`.
# `include`d into each test module; not named `test_*.jl`, so `runtests.jl` neither runs it directly
# nor expects a matching `src/` file. The including module must have `DimensionalData` and `Rasters`
# loaded.
#
# Nothing defines `==` for a layer, and DimensionalData's `==` ignores reference dims, names,
# metadata, the CRS and the lookup kind, so equality is written out. `typeof` catches the lookup
# kinds and every type parameter; the CRS is compared on its own because two CRSs of one type
# differ only in value.
function samedimarray(a, b)
    return typeof(a) === typeof(b) && isequal(parent(a), parent(b)) &&
           dims(a) == dims(b) &&
           map(Rasters.crs, dims(a)) == map(Rasters.crs, dims(b)) &&
           refdims(a) == refdims(b) &&
           DimensionalData.name(a) == DimensionalData.name(b) &&
           DimensionalData.metadata(a) == DimensionalData.metadata(b) &&
           isequal(Rasters.missingval(a), Rasters.missingval(b))
end

samechange(a::EcoSISTEM.NoLayerChange, b) = b isa EcoSISTEM.NoLayerChange

function samechange(a::EcoSISTEM.SeriesLayerChange, b)
    return typeof(a) === typeof(b) && isequal(a.slices, b.slices) &&
           a.times == b.times && a.origin == b.origin && a.atend == b.atend &&
           a.calendar == b.calendar && samedimarray(a.baseline, b.baseline)
end

function samelayer(a::EcoSISTEM.LayerCollection, b)
    return typeof(a) === typeof(b) && keys(a) == keys(b) &&
           all(samelayer(x, y) for (x, y) in zip(values(a), values(b)))
end

function samelayer(a, b)
    return typeof(a) === typeof(b) && samedimarray(a.matrix, b.matrix) &&
           a.size == b.size && samechange(a.change, b.change)
end
