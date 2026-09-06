# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Conformance to EcoBase, not our own API: every method this package defines on an `EcoBase` generic,
# ordered by type. Its siblings are `BaseInterface.jl` and `DiversityInterface.jl`.
#
# Two types answer, and only two. `StudyGrid` is the package's only `EcoBase.AbstractRegularGrid`,
# and a `GridHabitat` answers by handing over the grid it was built on. **A layer is not a grid**,
# and must not be given these methods: it knows a cell size but not where its cells are, so `xmin`
# would have to be invented and `xcellsize` would answer in whatever unit the layer happened to
# hold.
#
# **`indices` and `coordinates` report `(y, x)` columns, the order this package uses throughout,
# and the grid says so**: `coordinateorder` declares `YThenX()`, and EcoBase reorders for anyone who
# asks `XThenY()`, its own `convert_to_image` included. EcoBase's default is `XThenY()`, so the
# declaration is what keeps a grid plotted through it the right way up; a square test grid cannot
# show it either way round.
#
# Every definition is written qualified, so the file needs no `import` - only the `using` below,
# which binds the module name.

using EcoBase

# == StudyGrid - the AbstractRegularGrid interface ==================================================
# **`xmin`/`ymin` are the smallest coordinate LABEL, not the grid's outer edge.** EcoBase builds the
# whole coordinate surface from this label, the cell size and the cell count - the edges, `xmax`,
# `xrange` - so all three must be in the same space. These methods stay in the space the grid labels
# its cells in rather than re-anchoring to a centre or an outer bound.
EcoBase.xmin(grid::StudyGrid) = minimum(DimensionalData.lookup(grid.x))
EcoBase.ymin(grid::StudyGrid) = minimum(DimensionalData.lookup(grid.y))

# **And that space is the cell's lower corner** (`Intervals(Start)`), which is what this declares.
# EcoBase assumes a centre unless told, so without it every edge it derives - and every heatmap
# drawn from those edges - sits half a cell low, silently and on every grid at once.
EcoBase.cellanchor(::StudyGrid) = EcoBase.CellCorner()

# **And its columns come `y` first**, as everything in this package does. EcoBase assumes `x` first
# unless told; a habitat and an ecosystem answer the same through the grid they hold.
EcoBase.coordinateorder(::StudyGrid) = EcoBase.YThenX()

# **Unitful, and in the grid's OWN unit.** A geographic grid's cells are an angle across, so
# converting to a length here would give `° km^-1`: not a length, not an angle, and wrong without
# complaining. Whoever wants metres asks `getcellsizes` for them.
EcoBase.xcellsize(grid::StudyGrid) = _axisstep(grid.x)
EcoBase.ycellsize(grid::StudyGrid) = _axisstep(grid.y)

EcoBase.xcells(grid::StudyGrid) = length(grid.x)
EcoBase.ycells(grid::StudyGrid) = length(grid.y)
# Column 1 is y and column 2 is x, the order the grid declares. Rows are in the package's own cell
# order, column-major over `(Y, X)` with y fastest, so row `i` describes the same cell as column `i`
# of a `GridLandscape`'s abundance matrix.
function EcoBase.indices(grid::StudyGrid)
    ny, nx = length(grid.y), length(grid.x)
    out = Matrix{Int}(undef, ny * nx, 2)
    i = 0
    for c in 1:nx, r in 1:ny
        i += 1
        out[i, 1] = _ascendingrank(grid.y, r)
        out[i, 2] = _ascendingrank(grid.x, c)
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
        out[i, 1] = ylk[r]
        out[i, 2] = xlk[c]
    end
    return out
end

# == GridHabitat - the AbstractPlaces interface, answered from its StudyGrid ========================
# The habitat's location data - its `StudyGrid` - which is what makes a `GridHabitat` an
# `EcoBase.AbstractPlaces` with real coordinates rather than one that fakes its own.
#
# **Four methods, and four is all a habitat owes.** EcoBase answers every gridded question - `xmin`,
# `xrange`, `xedges`, `indices`, `cellanchor` and the rest - for anything holding gridded location
# data, and for an assemblage of such places, so a habitat and the `Ecosystem` over it both get the
# whole surface from the grid below without this file restating a line of it.
EcoBase.getcoords(habitat::GridHabitat) = habitat.area.builtgrid
function EcoBase.coordinates(habitat::GridHabitat)
    return EcoBase.coordinates(getcoords(habitat))
end
EcoBase.nplaces(habitat::GridHabitat) = countsubcommunities(habitat.regime)
EcoBase.placenames(habitat::GridHabitat) = _getsubcommunitynames(habitat)
