# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Conformance to EcoBase, not our own API: every method this package defines on an `EcoBase` generic,
# ordered by type. Its siblings are `BaseInterface.jl` and `DiversityInterface.jl`.
#
# Two types answer, and only two. `StudyGrid` is the package's only `EcoBase.AbstractGrid`, and a
# `GridHabitat` answers by handing over the grid it was built on. **A layer is not a grid**, and must
# not be given these methods: it knows a cell size but not where its cells are, so `xmin` would have
# to be invented and `xcellsize` would answer in whatever unit the layer happened to hold.
#
# **The one trap in this file: EcoBase's `indices` and `coordinates` are `(x, y)` COLUMNS, the
# opposite order to this package**, which is `(y, x)` throughout. `EcoBase.convert_to_image` reads `indices(grd, 1)` as the
# matrix *column*. Getting it the wrong way round transposes every grid plotted through EcoBase, and
# a square test grid cannot show it.
#
# Every definition is written qualified, so the file needs no `import` - only the `using` below,
# which binds the module name.

using EcoBase

# == StudyGrid - the AbstractGrid interface =========================================================
# **`xmin`/`ymin` are the smallest coordinate LABEL, not the grid's outer edge.** EcoBase derives
# `xmax = xmin + xcellsize * (xcells - 1)` - the label of the last cell - so the two must be in the
# same space. These grids label cells by their lower corner (`Intervals(Start)`), and these methods
# stay in that space rather than re-anchoring to a centre or an outer bound.
EcoBase.xmin(grid::StudyGrid) = minimum(DimensionalData.lookup(grid.x))
EcoBase.ymin(grid::StudyGrid) = minimum(DimensionalData.lookup(grid.y))

# **Unitful, and in the grid's OWN unit.** A geographic grid's cells are an angle across, so
# converting to a length here would give `° km^-1`: not a length, not an angle, and wrong without
# complaining. Whoever wants metres asks `getcellsizes` for them.
EcoBase.xcellsize(grid::StudyGrid) = _axisstep(grid.x)
EcoBase.ycellsize(grid::StudyGrid) = _axisstep(grid.y)

EcoBase.xcells(grid::StudyGrid) = length(grid.x)
EcoBase.ycells(grid::StudyGrid) = length(grid.y)
# Column 1 is x and column 2 is y, per EcoBase - see the header. Rows are in the package's own cell
# order, column-major over `(Y, X)` with y fastest, so row `i` describes the same cell as column `i`
# of a `GridLandscape`'s abundance matrix.
function EcoBase.indices(grid::StudyGrid)
    ny, nx = length(grid.y), length(grid.x)
    out = Matrix{Int}(undef, ny * nx, 2)
    i = 0
    for c in 1:nx, r in 1:ny
        i += 1
        out[i, 1] = _ascendingrank(grid.x, c)
        out[i, 2] = _ascendingrank(grid.y, r)
    end
    return out
end

EcoBase.indices(grid::StudyGrid, idx::Integer) = EcoBase.indices(grid)[:, idx]

# The real coordinates of each cell, in the same column and cell order as `indices` above. These are
# the cells' own labels - their lower corners - rather than midpoints reconstructed here, so they
# line up exactly with `xrange`/`yrange`.
function EcoBase.coordinates(grid::StudyGrid)
    ylk, xlk = DimensionalData.lookup(grid.y), DimensionalData.lookup(grid.x)
    ny, nx = length(ylk), length(xlk)
    out = Matrix{eltype(ylk)}(undef, ny * nx, 2)
    i = 0
    for c in 1:nx, r in 1:ny
        i += 1
        out[i, 1] = xlk[c]
        out[i, 2] = ylk[r]
    end
    return out
end

# == GridHabitat - the AbstractPlaces interface, answered from its StudyGrid ========================
# The habitat's location data - its `StudyGrid` - which is what makes a `GridHabitat` an
# `EcoBase.AbstractPlaces` with real coordinates rather than one that fakes its own.
EcoBase.getcoords(habitat::GridHabitat) = habitat.area.builtgrid
function EcoBase.coordinates(habitat::GridHabitat)
    return EcoBase.coordinates(getcoords(habitat))
end
EcoBase.nplaces(habitat::GridHabitat) = countsubcommunities(habitat.regime)
EcoBase.placenames(habitat::GridHabitat) = _getsubcommunitynames(habitat)
