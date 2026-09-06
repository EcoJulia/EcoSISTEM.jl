# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `src/materialise.jl` - what inspection shows against what the builder builds.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_materialise.jl"])'

module TestMaterialise

using Test
using EcoSISTEM
using EcoSISTEM: hasdata, landcoverclass
using EcoSISTEM: materialise
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using Rasters
using RasterDataSources
using ArchGDAL
using Distributions: Normal
using Extents: Extent
include("rasterfixtures.jl")
include("buildfixtures.jl")

# **`[ONE-PATH]`: `GridHabitat` puts in the habitat exactly what `materialise` shows.**
# This is a *structural* guard, not a numerical coincidence. A builder running its own near-copy of
# the inspection chain (`_materialiseon` -> `_assemble` -> `_resolve_regime`/`_resolve_supply`) drifts
# from it: three times over, in un-united dims, a hardcoded `Intervals(Start())`, and a `NicheSpec`
# that `materialise` built and the builder could not. Each
# was found by accident. The builder now calls `materialise`, so the assertions below can only fail
# if the two are pulled apart again.
#
# Both roles, both kinds of spec (data-backed and synthetic), on both kinds of positioned area -
# because the two paths differed *per kind*, so a single case would prove almost nothing.
# **`NicheSpec` is deliberately absent**: it is stochastic and unseeded (A19), so two
# materialisations of one spec disagree with each other, never mind with the builder.
@testset "what `materialise` shows is what `GridHabitat` builds" begin
    data = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)),
                axis = Temperature)
    watersrc = _reg(_bngraster(WorldClim{BioClim}, fill(4.0mm / day, 9, 9)),
                    axis = Precipitation)
    sun = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
    warm = UniformSpec(291.0K, axis = Temperature)

    for (label, area) in ("projected" => _area(regime = data),
                          "synthetic" => _area(extent = (40.0km, 40.0km),
                                cellsize = 10.0km))
        # A synthetic area cannot take a data-backed layer at all, so it tests the synthetic pair.
        regime = label == "projected" ? data : warm
        supply = label == "projected" ? watersrc : sun
        env = GridHabitat(regime = regime, supply = supply, area = area)
        seenreg = materialise(regime, area, role = EcoSISTEM.Condition)
        seensup = materialise(supply, area, role = EcoSISTEM.Resource)

        @test env.regime.matrix == seenreg.matrix
        @test env.supply.matrix == seensup.matrix
        # The **dims**, not just the values: A17 was two identical value arrays whose coordinates
        # were described differently, which is exactly what a value-only check misses.
        @test dims(env.regime.matrix, (Y, X)) == dims(seenreg.matrix, (Y, X))
        @test dims(env.supply.matrix, (Y, X)) == dims(seensup.matrix, (Y, X))
        @test env.regime.size == seenreg.size
        @test typeof(env.supply) === typeof(seensup)
    end

    # ...including the *mixed* multi-layer regime, where one member is generated at the grid's shape
    # and the other sampled onto it - the arity and the names both survive the round trip.
    area = _area(regime = data)
    env = GridHabitat(regime = (temperature = data, extra = warm),
                      supply = sun, area = area)
    seen = materialise((temperature = data, extra = warm), area,
                       role = EcoSISTEM.Condition)
    @test keys(EcoSISTEM.NamedTuple(env.regime)) ==
          keys(EcoSISTEM.NamedTuple(seen))
    @test env.regime.temperature.matrix == seen.temperature.matrix
    @test env.regime.extra.matrix == seen.extra.matrix
end

# A `RasterFileSpec` is the lazy route for a raster file that belongs to no dataset. What is pinned:
# it gives the layer `in_memory_raster` gives for the same file, it reads through the cache, the
# read is windowed, and it serves as a combine child. The fixture is a 5 x 7 WGS84 GeoTIFF of one
# degree cells written here, non-square so a transposition would show.
@testset "RasterFileSpec: a file read lazily, through the cache" begin
    dir = mktempdir()
    path = joinpath(dir, "field.tif")
    ArchGDAL.create(path, driver = ArchGDAL.getdriver("GTiff"), width = 7,
                    height = 5, nbands = 1, dtype = Float32) do ds
        # ArchGDAL writes `(x, y)`; the value encodes both so a swapped axis is visible.
        ArchGDAL.write!(ds, Float32[280 + i + 10j for i in 1:7, j in 1:5], 1)
        ArchGDAL.setgeotransform!(ds, [10.0, 1.0, 0.0, 55.0, 0.0, -1.0])
        return ArchGDAL.setproj!(ds,
                                 ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
    end
    spec = RasterFileSpec(path, axis = Temperature, unit = K)
    @test spec isa RasterFileSpec{Temperature}
    @test spec.source === EcoSISTEM.SyntheticData
    @test occursin("RasterFileSpec(", repr(spec))
    @test occursin("unit = K", repr(spec))
    @test occursin("axis = Temperature", repr(spec))
    # A URL is deferred to a cached download, as on `ShapeSpec`; nothing is fetched here.
    @test RasterFileSpec("https://example.org/a.tif",
                         axis = EcoSISTEM.NicheAxis).path isa
          EcoSISTEM.CachedAsset

    # The file shapes the grid, and the layer is the one the eager route gives.
    area = StudyArea(regime = spec, verbosity = :silent)
    @test size(area.report.active) == (5, 7)
    @test area.report.cellsize == 1.0°
    lazy = materialise(spec, area)
    eager = materialise(EcoSISTEM.in_memory_raster(readfile(path, unit = K),
                                                   axis = Temperature), area)
    @test lazy.matrix == eager.matrix
    @test eltype(lazy.matrix) <: Unitful.Temperature
    # Read once, through the cache, keyed on the path and the unit.
    keys_ = collect(keys(area.report.cache.reads))
    @test length(keys_) == 1
    @test only(keys_).code == path
    @test only(keys_).readkw.unit == K

    # A `within` box windows the read before the pixels are fetched, and narrows the grid.
    win = StudyArea(regime = spec,
                    within = Extent(Y = (51.2°, 53.8°), X = (11.5°, 14.5°)),
                    verbosity = :silent)
    @test size(win.report.active) == (3, 4)
    @test !isnothing(only(keys(win.report.cache.reads)).readkw.cut)

    # A combine child, and a mask through one.
    doubled = ConstructedRasterSpec(r -> r .* 2, spec, axis = Temperature)
    @test materialise(doubled, area).matrix == 2 .* lazy.matrix
    mask = ConstructedRasterSpec(r -> r .> 300K, spec,
                                 axis = EcoSISTEM.NicheAxis)
    masked = StudyArea(regime = spec, within = mask, verbosity = :silent)
    @test all(masked.report.active)
    @test size(masked.report.active) == (4, 7)

    # `scale` coarsens on read: 2 x 2 blocks of the 5 x 7 file become one cell each, and each holds
    # the block mean. Blocks are anchored at the south-west corner as the study grid's cells are, so
    # the partial row dropped is the **northern** one (file row 1) and the partial column the eastern
    # (file column 7). Values are `280 + col + 10 row` with rows counted from the top of the file:
    # the southern block (rows 4-5) comes first once latitude ascends. Memoised on disk, as a
    # dataset's is.
    coarse = RasterFileSpec(path, axis = Temperature, unit = K, scale = 2)
    @test occursin("scale = 2", repr(coarse))
    @test_throws ErrorException RasterFileSpec(path, axis = Temperature,
                                               scale = 0)
    carea = StudyArea(regime = coarse, verbosity = :silent)
    @test size(carea.report.active) == (2, 3)
    @test carea.report.cellsize == 2.0°
    @test materialise(coarse, carea).matrix ≈
          [326.5 328.5 330.5; 306.5 308.5 310.5]K
    @test only(keys(carea.report.cache.reads)).readkw.scale == 2
    @test isfile(EcoSISTEM._aggcachepath(path, 2,
                                         EcoSISTEM._reducer(nothing,
                                                            Temperature),
                                         K))

    # A file of class codes must not be averaged: on a `TypologyAxis` the reducer is the most
    # frequent class, ties to the smallest code, and `fn` overrides it. Rows are given top to bottom
    # as the file stores them; the flip to ascending latitude puts the lower block row first.
    cpath = joinpath(dir, "classes.tif")
    ArchGDAL.create(cpath, driver = ArchGDAL.getdriver("GTiff"), width = 4,
                    height = 4, nbands = 1, dtype = Float32) do ds
        codes = Float32[1 1 2 3; 1 2 2 3; 5 5 7 7; 5 7 7 9]
        ArchGDAL.write!(ds, permutedims(codes), 1)      # ArchGDAL writes `(x, y)`
        ArchGDAL.setgeotransform!(ds, [10.0, 1.0, 0.0, 54.0, 0.0, -1.0])
        return ArchGDAL.setproj!(ds,
                                 ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
    end
    classes = RasterFileSpec(cpath, axis = LandCoverTypology, scale = 2)
    cls = StudyArea(regime = classes, verbosity = :silent)
    @test materialise(classes, cls).matrix == [5 7; 1 2]     # {2,3,2,3} ties to 2
    @test isfile(EcoSISTEM._aggcachepath(cpath, 2, EcoSISTEM._majorityclass,
                                         NoUnits))
    biggest = RasterFileSpec(cpath, axis = LandCoverTypology, scale = 2,
                             fn = maximum)
    @test occursin("fn = maximum", repr(biggest))
    @test materialise(biggest,
                      StudyArea(regime = biggest, verbosity = :silent)).matrix ==
          [7 9; 2 3]
    # The same file on a continuous axis is averaged, as before.
    asvalues = RasterFileSpec(cpath, axis = EcoSISTEM.NicheAxis, scale = 2)
    @test materialise(asvalues,
                      StudyArea(regime = asvalues, verbosity = :silent)).matrix ≈
          [5.5 7.5; 1.25 2.5]

    # **Route R1 of the regridding pipeline.** A study grid that is a whole multiple of the file's
    # cells, aligned to them, is an exact block aggregation of the file - what the report has always
    # called `LayerAggregated(f)` - and not a resample. Three things pin it: the values are the block
    # means to the last bit; the categorical file aggregates by majority the same way; and coarsening
    # on read (`scale = f`) then materialising on the file's own lattice gives the **bit-identical**
    # layer, because both routes are one `_blockaggregate`. That last one is the duplicate-against-
    # original test: two implementations of one coarsening can only be trusted while they are compared.
    M = ustrip.(K, lazy.matrix)                          # the 5 x 7 file, ascending latitude
    for f in (2,)
        area_f = StudyArea(regime = spec, cellsize = float(f) * °,
                           verbosity = :silent)
        @test only(area_f.report.layers).kind == EcoSISTEM.LayerAggregated(f)
        got = materialise(spec, area_f).matrix
        want = [sum(M[((r - 1) * f + 1):(r * f), ((c - 1) * f + 1):(c * f)]) /
                f^2 for r in 1:(5 ÷ f), c in 1:(7 ÷ f)]
        @test ustrip.(K, got) == want
        pre = RasterFileSpec(path, axis = Temperature, unit = K, scale = f)
        area_pre = StudyArea(regime = pre, verbosity = :silent)
        @test only(area_pre.report.layers).kind == EcoSISTEM.LayerKeptExactly()
        @test materialise(pre, area_pre).matrix == got
    end
    cgrid = StudyArea(regime = RasterFileSpec(cpath, axis = LandCoverTypology),
                      cellsize = 2.0°, verbosity = :silent)
    @test materialise(RasterFileSpec(cpath, axis = LandCoverTypology), cgrid).matrix ==
          [5 7; 1 2]

    # And a habitat builds on it - geographic, so it can be inspected but not simulated.
    h = GridHabitat(regime = spec,
                    supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                         axis = SolarRadiation),
                    area = area, topology = Torus())
    @test h.regime.matrix == lazy.matrix
end

end
