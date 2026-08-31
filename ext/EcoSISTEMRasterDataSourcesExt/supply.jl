# SPDX-License-Identifier: LGPL-3.0-or-later

# --- The one supply constructor keyed on a dataset ---------------------------------------------------
#
# Every other `Supply{A}` constructor stays in `src/Layer.jl`: they take a `Matrix` or a
# `DimArray` and name no source. This one reads a WorldClim BioClim raster, so it is here.

# The raw raster is a per-area rate at its own native reporting period (e.g. `L*m^-2*yr^-1`
# for `bio12`) - `cancel` multiplies by the grid's (area-preserving, see `_cellsize`) cell area and
# states the result in the axis's canonical unit, whatever the native time unit; no separate
# pre-`uconvert` step needed. `bioclim.array` is kept as its real `(Y, X)` `DimArray`, not stripped,
# so the supply shares the regime's real grid.
function EcoSISTEM.Supply{A}(bioclim::ClimateRaster{WorldClim{BioClim}}) where {A <:
                                                                                NicheAxis}
    cellareas = _cellareas(bioclim.array)
    return Supply{A}(cancel.(bioclim.array, cellareas, A))
end
