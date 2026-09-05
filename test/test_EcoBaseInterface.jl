# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Needs `RasterDataSources` (a weak dependency) for the geographic fixtures, so it must go through
# `Pkg.test`:
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["test_EcoBaseInterface.jl"])'

module TestEcoBaseInterface

using EcoSISTEM
using EcoSISTEM: StudyGrid, CellNames, _axisstep
using EcoSISTEM.Units
using EcoBase
using DimensionalData: DimensionalData, DimArray, X, Y
using Rasters
using RasterDataSources
using Unitful, Unitful.DefaultSymbols
using Test

include("rasterfixtures.jl")

# **Non-square on purpose, and not evenly divisible either** - 11 rows of y against 7 columns of
# x. Every `(y, x)` mix-up this package has had was invisible on a square grid.
const YEXTENT, XEXTENT, CELL = 11.0km, 7.0km, 1.0km
function _projected()
    area = StudyArea(extent = (YEXTENT, XEXTENT), cellsize = CELL,
                     verbosity = :silent)
    return GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                       supply = UniformSpec(10.0kJ / m^2 / day,
                                            axis = SolarRadiation),
                       area = area)
end

# A real geographic grid - 0.5° cells over 7 × 11 cells - which is the case the old implementation
# could not express at all.
function _geographic(; lat = (50.0:0.5:53.0) .* °, long = (0.0:0.5:5.0) .* °)
    raster = _testraster(WorldClim{BioClim},
                         fill(291.0K, length(lat),
                              length(long)),
                         lat = lat, long = long)
    area = StudyArea(regime = _reg(raster), verbosity = :silent)
    return GridHabitat(regime = ConstructedRasterSpec(() -> raster,
                                                      axis = Temperature),
                       supply = UniformSpec(10.0kJ / m^2 / day,
                                            axis = SolarRadiation),
                       area = area)
end

@testset "a habitat's grid is its own location data" begin
    hab = _projected()
    grid = EcoBase.getcoords(hab)
    @test grid isa StudyGrid
    # **The type parameter reaches the `Ecosystem`**, which is the whole reason `StudyArea` is
    # parameterised: an `Any` here would put a dynamic read on every coordinate question.
    @test isconcretetype(typeof(grid))
    @test hab.area isa StudyArea{typeof(grid)}
    @test hab isa EcoSISTEM.AbstractHabitat{<:Any, <:Any, typeof(grid)}
    # ...and it really is EcoBase's location data, not merely shaped like it.
    @test grid isa EcoBase.AbstractGrid
    @test grid isa EcoBase.AbstractLocationData

    # **The same dimension objects `active` is indexed by, not a copy** - the rule that stops a
    # grid and the array it describes drifting apart.
    yd, xd = DimensionalData.dims(hab.active, (Y, X))
    @test grid.y === yd
    @test grid.x === xd
end

@testset "the seven grid questions, answered from the coordinates" begin
    hab = _projected()
    grid = EcoBase.getcoords(hab)
    ny, nx = Base.size(hab.active)
    @test (ny, nx) == (11, 7)

    @test EcoBase.xcells(grid) == nx
    @test EcoBase.ycells(grid) == ny
    # **Unitful, in the grid's own unit** - the old implementation returned a bare `Float64`.
    @test EcoBase.xcellsize(grid) == CELL
    @test EcoBase.ycellsize(grid) == CELL
    @test EcoBase.cellsize(grid) == (CELL, CELL)
    # `xmin`/`ymin` are read, not fabricated. This grid starts at zero, so the assertion that
    # carries the information is the *geographic* one below - here the point is that the units and
    # EcoBase's own derivations survive.
    @test EcoBase.xmin(grid) == 0.0km
    @test EcoBase.xmax(grid) == (nx - 1) * CELL
    @test EcoBase.ymax(grid) == (ny - 1) * CELL
    @test collect(EcoBase.xrange(grid)) == collect((0:(nx - 1)) .* CELL)
    @test length(EcoBase.yrange(grid)) == ny
end

# **The measurement that made B3 a defect fix rather than a tidy-up.** The old
# `Float64(regime.size / km)` did not throw on a geographic grid - it returned `1.0 ° km^-1`, a unit
# that means nothing, *silently*, and the `DimensionError` only surfaced two functions later inside
# EcoBase's own `xmax`, where it read as the caller's fault.
@testset "a geographic grid answers in degrees, and its derivations work" begin
    hab = _geographic()
    grid = EcoBase.getcoords(hab)
    @test EcoBase.xcellsize(grid) == 0.5°
    @test EcoBase.ycellsize(grid) == 0.5°
    @test Unitful.dimension(EcoBase.xcellsize(grid)) === Unitful.NoDims
    # The origin is the grid's own, not `0` - the second half of the same defect.
    @test EcoBase.ymin(grid) == 50.0°
    @test EcoBase.xmin(grid) == 0.0°
    # And EcoBase's derived questions answer, rather than throwing: they are plain arithmetic over
    # what the grid reports, so they work on angles unchanged.
    @test EcoBase.xmax(grid) == 0.0° + (EcoBase.xcells(grid) - 1) * 0.5°
    @test length(collect(EcoBase.xrange(grid))) == EcoBase.xcells(grid)
    @test !isnothing(grid.crs)
end

@testset "indices and coordinates use EcoBase's (x, y) column order" begin
    hab = _projected()
    grid = EcoBase.getcoords(hab)
    ny, nx = Base.size(hab.active)
    idx = EcoBase.indices(grid)
    @test Base.size(idx) == (ny * nx, 2)
    # **Column 1 is x, column 2 is y** - the opposite order to the rest of this package, and
    # EcoBase's convention rather than a slip: `convert_to_image` uses `indices(grd, 1)` as the
    # matrix *column*. The old implementation handed out `(y, x)` rows, so anything plotting through
    # EcoBase transposed the grid - invisible on a square one.
    @test maximum(idx[:, 1]) == nx
    @test maximum(idx[:, 2]) == ny
    @test EcoBase.indices(grid, 1) == idx[:, 1]
    @test EcoBase.indices(grid, 2) == idx[:, 2]
    # Rows are in the package's own cell order - column-major over `(Y, X)`, y fastest - so cell 2
    # is the next one *down* the first column, not along the first row.
    @test idx[1, :] == [1, 1]
    @test idx[2, :] == [1, 2]
    @test idx[ny + 1, :] == [2, 1]

    coords = EcoBase.coordinates(grid)
    @test Base.size(coords) == (ny * nx, 2)
    @test coords[1, :] == [0.0km, 0.0km]
    @test coords[2, :] == [0.0km, CELL]
    @test coords[ny + 1, :] == [CELL, 0.0km]
end

# A raster's `Y` commonly runs north to south, so its array rows descend while EcoBase's `yrange`
# ascends by construction. Ranking rather than passing the row straight out is what keeps such a grid
# the right way up - and it is a branch that fires on no grid this package builds, so it needs a
# fixture that is deliberately reversed.
@testset "a descending Y axis is ranked, not passed through" begin
    hab = _geographic(lat = (53.0:-0.5:50.0) .* °)
    grid = EcoBase.getcoords(hab)
    ylk = DimensionalData.lookup(grid.y)
    @test !issorted(ylk)
    # `ymin` is the smallest coordinate whichever way the array runs...
    @test EcoBase.ymin(grid) == 50.0°
    @test EcoBase.ycellsize(grid) == 0.5°
    # ...and array row 1, which holds the *largest* latitude, ranks last.
    idx = EcoBase.indices(grid)
    @test idx[1, 2] == EcoBase.ycells(grid)
    @test idx[EcoBase.ycells(grid), 2] == 1

    # **The invariant that ties `indices` and `coordinates` together**, and the one that would
    # catch a rank applied to one but not the other: a cell's coordinate is exactly the entry its own
    # index picks out of `yrange`. It holds whichever way the array runs, which is the whole point.
    coords = EcoBase.coordinates(grid)
    yr, xr = collect(EcoBase.yrange(grid)), collect(EcoBase.xrange(grid))
    @test all(coords[i, 2] ≈ yr[idx[i, 2]] for i in Base.axes(idx, 1))
    @test all(coords[i, 1] ≈ xr[idx[i, 1]] for i in Base.axes(idx, 1))
end

@testset "cells are named by where they are" begin
    hab = _projected()
    names = EcoBase.placenames(hab)
    ny, nx = Base.size(hab.active)

    @test names isa CellNames
    @test names isa AbstractVector{String}
    @test length(names) == ny * nx
    @test length(names) == EcoSISTEM.countsubcommunities(hab.regime)
    # **Y first, then X**, and a half-open interval, so the name says exactly which ground it is.
    @test names[1] == "[0.0, 1.0) × [0.0, 1.0) km"
    @test names[2] == "[1.0, 2.0) × [0.0, 1.0) km"
    @test names[ny + 1] == "[0.0, 1.0) × [1.0, 2.0) km"
    # **Distinct, which the whole point depends on**: a name that repeated would be no better
    # than the indices it replaced.
    @test length(Set(names)) == length(names)

    # **Lazy, and that is a requirement rather than an optimisation.**  The eager equivalent of a
    # 1.2 million-cell grid's names is ~33 MB, on a habitat whose `active` mask was deliberately
    # squeezed to 0.14 MB. Nothing in the simulation reads a cell's name.
    @test Base.summarysize(names) < Base.summarysize(collect(names))
    # ...and `reshape`, which is what Diversity does to them, stays a view.
    @test parent(reshape(names, 1, length(names))) === names

    # Diversity reaches the same names through its own accessor.
    @test EcoSISTEM._getsubcommunitynames(hab) == names
end

@testset "a geographic cell names its axes" begin
    names = EcoBase.placenames(_geographic())
    # `N`/`E` say which axis is which and which way is positive - a negative longitude is still
    # `°E`, just west of the meridian.
    @test names[1] == "[50.0, 50.5)°N × [0.0, 0.5)°E"
    @test occursin("°N", names[end])
    @test occursin("°E", names[end])
    @test length(Set(names)) == length(names)
end

# Enough decimal places to tell neighbouring cells apart, taken from the step rather than fixed:
# a 30-arcsecond grid needs more than a 1 km one, and rounding to a common number would make whole
# rows of cells share a name.
@testset "the precision follows the cell size" begin
    fine = _geographic(lat = (50.0:0.01:50.1) .* °, long = (0.0:0.01:0.1) .* °)
    names = EcoBase.placenames(fine)
    @test names[1] == "[50.0, 50.01)°N × [0.0, 0.01)°E"
    @test length(Set(names)) == length(names)
end

@testset "a grid must be able to say where its cells are" begin
    hab = _projected()
    grid = EcoBase.getcoords(hab)
    # Bare cell indices cannot state a position or a cell size, which is this type's whole job.
    bare = DimArray(zeros(3, 4),
                    (Y(DimensionalData.NoLookup()),
                     X(DimensionalData.NoLookup())))
    @test_throws ErrorException StudyGrid(nothing, bare)
    # ...and neither can numbers dressed as coordinates.
    unitless = DimArray(zeros(3, 4), (Y(1.0:3.0), X(1.0:4.0)))
    @test_throws ErrorException StudyGrid(nothing, unitless)
    # The habitat's own array is of course accepted, and rebuilds the same grid.
    rebuilt = StudyGrid(grid.crs, hab.active)
    @test rebuilt.y === grid.y
    @test rebuilt.x === grid.x
end

# **`stage` is NOT redundant with the area's type parameter**, and this is the case that shows
# it. `StudyArea(habitat)` copies an `AsBuilt` report - it describes a habitat that exists - into a
# *fresh* area that nothing has yet been built on. So the report's stage and the area's `builtgrid`
# answer two different questions, and neither can be derived from the other.
@testset "a copied area keeps AsBuilt but carries no grid" begin
    hab = _projected()
    @test hab.area.report.stage isa EcoSISTEM.AsBuilt
    @test !isnothing(hab.area.builtgrid)

    copied = StudyArea(hab, verbosity = :silent)
    @test copied.report.stage isa EcoSISTEM.AsBuilt
    @test isnothing(copied.builtgrid)
    @test copied isa StudyArea{Nothing}
    # An area that has merely been investigated likewise has a decided grid but nothing built on it.
    fresh = StudyArea(extent = (YEXTENT, XEXTENT), cellsize = CELL,
                      verbosity = :silent)
    @test fresh isa StudyArea{Nothing}
    @test fresh.report.stage isa EcoSISTEM.AsInvestigated
    # ...and building on either gives back a habitat that does carry one.
    rebuilt = GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                          supply = UniformSpec(10.0kJ / m^2 / day,
                                               axis = SolarRadiation),
                          area = copied)
    @test rebuilt.area.builtgrid isa StudyGrid
    # And the caller's area is unchanged - building yields a new one rather than mutating theirs.
    @test isnothing(copied.builtgrid)
end

# **A layer is values ON a grid, and must not be a grid itself.** Putting the seven methods on
# `AbstractRegime` is what lets them fabricate an origin and a cell size it does not have.
@testset "a layer is not an EcoBase grid" begin
    hab = _projected()
    @test !(hab.regime isa EcoBase.AbstractGrid)
    @test !(EcoSISTEM.AbstractLayer <: EcoBase.AbstractGrid)
    @test_throws Exception EcoBase.xcellsize(hab.regime)
end

end
