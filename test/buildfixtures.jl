# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Fixtures shared by the tests of the *building* path — `GridHabitat`, `build_species`,
# `build_ecosystem`, and the raster and study-area machinery underneath them.
#
# Not a test file. It is `include`d, like `rasterfixtures.jl`, and defines no testsets: naming it
# in `Pkg.test(test_args = …)` would run nothing. The fixtures live here rather than in any one test
# file because the building path is tested from several, and they are needed by every piece.
#
# `_area`'s leading `;` is load-bearing: without it `kw...` splats positionally and every caller
# gets a `MethodError` about `Pair`s. See `CLAUDE.md`.

# A synthetic single-polygon shapefile (a rectangle over `(x1,y1)-(x2,y2)` in `sr`, or with no CRS
# metadata at all if `sr = nothing`), written to a temp dir for `ShapeSpec` tests — no checked-in
# binary fixture needed.
function _testshapefile(x1, y1, x2, y2; sr = nothing)
    path = joinpath(mktempdir(), "test.shp")
    poly = ArchGDAL.createpolygon([[
                                      (x1, y1),
                                      (x2, y1),
                                      (x2, y2),
                                      (x1, y2),
                                      (x1, y1)
                                  ]])
    ArchGDAL.create(path, driver = ArchGDAL.getdriver("ESRI Shapefile")
                    ) do dataset
        makelayer = layer -> begin
            ArchGDAL.addfeature(layer) do feature
                return ArchGDAL.setgeom!(feature, poly)
            end
            ArchGDAL.copy(layer, dataset = dataset)
        end
        if isnothing(sr)
            ArchGDAL.createlayer(makelayer, name = "test",
                                 geom = ArchGDAL.wkbPolygon)
        else
            ArchGDAL.createlayer(makelayer, name = "test",
                                 geom = ArchGDAL.wkbPolygon,
                                 spatialref = sr)
        end
    end
    return path
end

# Shared supply spec for the synthetic/data builders — a flat solar rate. Every `GridHabitat`
# call needs an explicit `supply` spec: a bare quantity is not accepted.
const SUP = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
# `_reg` (wrap a raster as a layer spec) is shared, in `rasterfixtures.jl`.
# A study area, silenced: these tests assert on the grid they get, not on how it was announced, so
# `:silent` keeps their log expectations about the *build* rather than the area decision. Every
# constraint (`regime`/`supply`/`within`/`cellsize`/`crs`/`extent`/`align`) passes straight through.
_area(; kw...) = StudyArea(; verbosity = :silent, kw...)
# The common shape by far: an environment on an area decided from the regime alone. Where the supply
# is also meant to shape the grid, or the area is itself what is being tested, it is written out.
function _env(regime, supply; kw...)
    return GridHabitat(regime = regime, supply = supply,
                       area = _area(; regime = regime, kw...))
end
# The former `build_species` defaults, now that `tolerance`/`demand` are required — reused by the
# many calls that only exercise the other keywords (movement, rates, abundance, …).
const TOL = (298.0K, 2.0K)
const DEM = 10.0 * EcoSISTEM.canonicalunit(EcoSISTEM.Resource, SolarRadiation)
