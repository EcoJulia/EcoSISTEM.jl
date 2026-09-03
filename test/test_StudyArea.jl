# SPDX-License-Identifier: LGPL-3.0-or-later

module TestStudyArea

using EcoSISTEM
# `[C7-VIS]` C: these are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: SteadyLayerChange
using EcoSISTEM: LayerCache, ReadKey, _asraster, Problem, _originaligned,
                 _resamplecost, _choosealign, _gridfootprint, _footprintmessage,
                 _analyse, _snapbounds, _crsname, materialise, _materialiseon,
                 _sampledata
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using DimensionalData: DimensionalData, DimArray, X, Y
using Rasters
using RasterDataSources
using Extents: Extent
using Test

include("rasterfixtures.jl")
using EcoSISTEM: hasdata, landcoverclass
using EcoSISTEM: materialise
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using ArchGDAL
using Distributions: Normal
include("buildfixtures.jl")

@testset "_analyse" begin
    @testset "a single layer is adopted exactly" begin
        # Nothing specified: the layer's own CRS, its own cell size, and itself as the alignment
        # layer - so it is copied rather than resampled. This is the case that must never degrade.
        bng = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
        p = _analyse((regime = bng,))
        @test p.crssource isa EcoSISTEM.AdoptedFromLayers
        @test p.cellsize == 2500.0m
        @test p.cellsizesource isa EcoSISTEM.TakenFromAlignedLayer
        @test p.align === :regime
        @test size(p.active) == (9, 9)
        @test only(p.layers).kind isa EcoSISTEM.LayerKeptExactly
        @test isempty(p.problems)
        @test p.footprint.cells == 81
    end

    @testset "a whole multiple aggregates exactly rather than resampling" begin
        bng = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
        p = _analyse((regime = bng,), cellsize = 5000.0m)
        plan = only(p.layers)
        @test plan.kind isa EcoSISTEM.LayerAggregated
        @test plan.kind.factor == 2
        @test p.cellsizesource isa EcoSISTEM.GivenByUser
    end

    @testset "a finer target upsamples, and says so" begin
        bng = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
        p = _analyse((regime = bng,), cellsize = 1250.0m)
        @test only(p.layers).kind isa EcoSISTEM.LayerResampled
        @test :upsampling in [pr.code for pr in p.problems]
    end

    @testset "an explicit `align` must name a layer that was given" begin
        bng = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
        @test_throws ErrorException _analyse((regime = bng,), align = :supply)
    end

    @testset "a geographic target warns that it cannot be simulated" begin
        wgs = _reg(_testraster(WorldClim{BioClim}, fill(291.0K, 5, 5)))
        p = _analyse((regime = wgs,))
        @test :geographic in [pr.code for pr in p.problems]
    end

    @testset "a synthetic area needs both extent and cellsize" begin
        p = _analyse(NamedTuple(), extent = (4km, 12km), cellsize = 1km)
        @test isnothing(p.crs)
        @test p.crssource isa EcoSISTEM.NoRealWorldPosition
        @test size(p.active) == (4, 12)     # (y, x), as everywhere else
        @test isempty(p.layers)
        @test_throws ErrorException _analyse(NamedTuple(), extent = (4km, 12km))
    end
end

# The whole point of windowing, and the property that makes it safe to have: a study area decided
# from a *windowed* read must equal one decided from the whole layer. The window is computed before
# any pixel is read - the target CRS comes from the file header via `sourcecrs` - and padded by
# `_WINDOW_PAD`, because snapping the grid outwards onto cell boundaries, `_atleast2` and resampling
# all need ground beyond the mask. A window cut exactly to the mask is *not* enough: it comes up one
# cell short at each snapped edge. Guarded: it needs the real global layers.
if !Sys.iswindows()
    @testset "windowing the reads does not change the answer" begin
        scot = EcoSISTEM.boundingbox("Scotland", coverage = AllTerritories())
        for (src, code) in ((WorldClim{BioClim}, :bio1),   # scale 1
            (EarthEnv{LandCover}, 7))      # scale 10 - aggregation blocks
            # An already-read raster cannot be windowed, so this is the unwindowed reference.
            whole = EcoSISTEM._read(SourceSpec(src, code))
            ref = investigate_study_area(regime = ConstructedSpec(() -> whole,
                                                                  axis = EcoSISTEM.NicheAxis),
                                         within = scot)
            win = investigate_study_area(regime = SourceSpec(src, code),
                                         within = scot)
            @test size(win.active) == size(ref.active)
            @test win.active == ref.active
            @test win.cellsize == ref.cellsize
            @test [l.kind for l in win.layers] == [l.kind for l in ref.layers]
        end
    end
end

# A multi-layer read is assembled from **per-layer** cache entries, so a whole-dataset request and a
# single-layer one share what they have in common: without that, a study area with an
# `EarthEnv{LandCover}` regime and a `SourceSpec(EarthEnv{LandCover}, 7)` supply reads file 7 twice.
# Narrowing the *key* instead would be wrong - it must describe what was asked for, or one request is
# served another's data. Decomposing is possible only because a spec's code list is explicit.
if !Sys.iswindows()
    @testset "a multi-layer read caches its layers individually" begin
        cache = LayerCache()
        whole = _asraster(SourceSpec(EarthEnv{LandCover}), cache)
        @test length(cache) == 12                 # one entry per layer, not one for the request
        # ...so asking for a single layer afterwards is a hit, adding nothing.
        one = _asraster(SourceSpec(EarthEnv{LandCover}, 7), cache)
        @test length(cache) == 12
        @test size(one.array) == size(whole.array)[1:2]

        # Assembling the stack here rather than in `_readmultilayer` must not change it: same data,
        # and the same canonical *names* on the layer axis (`EarthEnv` code 7 has always shown as
        # `:cultivated_and_managed`, never as `7`).
        direct = EcoSISTEM._read(SourceSpec(EarthEnv{LandCover}))
        lax(a) = parent(DimensionalData.lookup(a.array,
                                               DimensionalData.Dim{:layer}))
        @test lax(whole) == lax(direct)
        @test first(lax(whole)) === :needleleaf_trees
        @test isequal(parent(whole.array), parent(direct.array))
    end
end

@testset "StudyArea and investigate_study_area" begin
    bng() = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))

    @testset "the two faces share one analysis" begin
        # The investigator must describe the grid that would actually be built, not a reconstruction.
        r = investigate_study_area(regime = bng())
        a = StudyArea(regime = bng(), verbosity = :silent)
        @test r.crs == a.report.crs
        @test r.cellsize == a.report.cellsize
        @test r.align == a.report.align
        @test size(r.active) == size(a.report.active)
        @test [l.kind for l in r.layers] == [l.kind for l in a.report.layers]
    end

    @testset "the report reads clearly and is programmable" begin
        r = investigate_study_area(regime = bng())
        text = sprint(show, MIME"text/plain"(), r)
        @test occursin("EPSG:27700", text)          # not GeoFormatTypes.EPSG{1}((27700,))
        @test occursin("9 × 9 cells", text)
        @test occursin("kept exactly", text)
        @test occursin("no problems found", text)
        # ...and the same facts are reachable without parsing the text.
        @test isempty(r.problems)
        @test only(r.layers).kind isa EcoSISTEM.LayerKeptExactly
        @test r.footprint.cells == 81
    end

    @testset "refinement inherits everything not overridden, including the cache" begin
        a = StudyArea(regime = bng(), verbosity = :silent)
        b = StudyArea(a, cellsize = 5000.0m, verbosity = :silent)
        @test b.report.cache === a.report.cache      # no re-reading between attempts
        @test b.report.specs.regime === a.report.specs.regime
        @test b.report.cellsize == 5000.0m
        @test only(b.report.layers).kind isa EcoSISTEM.LayerAggregated

        # `nothing` clears an inherited constraint, where `missing` would have kept it.
        @test isnothing(StudyArea(b, cellsize = nothing,
                                  verbosity = :silent).report.constraints.cellsize)
        @test StudyArea(b, verbosity = :silent).report.constraints.cellsize ==
              5000.0m
    end

    @testset "a report can be committed to, reusing its cache" begin
        # The investigate -> commit loop: the report is a valid `base`, so accepting what it proposes
        # needs no re-reading and no restating of the arguments.
        r = investigate_study_area(regime = bng())
        a = StudyArea(r, verbosity = :silent)
        @test a.report.cache === r.cache
        @test a.report.crs == r.crs
        @test a.report.cellsize == r.cellsize
        @test size(a.report.active) == size(r.active)
        # Constraints are inherited *as given*, so a derived CRS is re-derived rather than being
        # relabelled as though the user had supplied it.
        @test a.report.crssource isa EcoSISTEM.AdoptedFromLayers

        # A report also refines, and can be investigated further before committing.
        r2 = investigate_study_area(r, cellsize = 5000.0m)
        @test r2.cache === r.cache
        @test only(r2.layers).kind isa EcoSISTEM.LayerAggregated
        @test StudyArea(r2, verbosity = :silent).report.cellsize == 5000.0m
    end

    @testset "both entry points treat missing/nothing identically" begin
        # The sentinels are resolved in one shared place, so investigating a refinement and
        # committing to it must agree: `missing` inherits, `nothing` clears.
        base = investigate_study_area(regime = bng(), cellsize = 5000.0m,
                                      crs = Rasters.EPSG(27700))
        @test base.constraints.cellsize == 5000.0m
        @test base.crssource isa EcoSISTEM.GivenByUser

        # `missing` (the default) carries both constraints through.
        kept = investigate_study_area(base)
        @test kept.constraints.cellsize == 5000.0m
        @test kept.constraints.crs == Rasters.EPSG(27700)

        # `nothing` clears them, so each is derived afresh - the cell size falls back to the
        # aligned layer's own, and the CRS goes back to being adopted rather than given.
        cleared = investigate_study_area(base, cellsize = nothing,
                                         crs = nothing)
        @test isnothing(cleared.constraints.cellsize)
        @test cleared.cellsize == 2500.0m
        @test cleared.cellsizesource isa EcoSISTEM.TakenFromAlignedLayer
        @test cleared.crssource isa EcoSISTEM.AdoptedFromLayers

        # ...and `StudyArea` does exactly the same with the same arguments.
        a = StudyArea(base, cellsize = nothing, crs = nothing,
                      verbosity = :silent)
        @test a.report.cellsize == cleared.cellsize
        @test a.report.cellsizesource == cleared.cellsizesource
        @test a.report.crssource == cleared.crssource
    end

    @testset "no compulsory arguments, but not any combination" begin
        @test StudyArea(extent = (4km, 12km), cellsize = 1km,
                        verbosity = :silent) isa StudyArea
        @test_throws ErrorException StudyArea()
        @test_throws ErrorException StudyArea(extent = (4km, 12km))
    end

    @testset "verbosity is presentation only" begin
        # `:silent` says nothing, `:normal` announces what it guessed, `:verbose` shows the report -
        # and all three build the identical grid.
        @test_logs StudyArea(regime = bng(), verbosity = :silent)
        @test_logs (:info,) match_mode=:any StudyArea(regime = bng(),
                                                      verbosity = :normal)
        quiet = StudyArea(regime = bng(), verbosity = :silent)
        loud = StudyArea(regime = bng(), verbosity = :normal)
        @test size(quiet.report.active) == size(loud.report.active)
        @test_throws ErrorException StudyArea(regime = bng(),
                                              verbosity = :chatty)
    end

    @testset "a StudyArea shows its essentials on one line" begin
        text = sprint(show, StudyArea(regime = bng(), verbosity = :silent))
        @test occursin("9 × 9 cells of 2500.0 m", text)
        @test occursin("EPSG:27700", text)
        @test occursin("81 active", text)
    end
end

@testset "materialise onto a decided grid" begin
    bng() = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))

    @testset "a data layer lands on the area's own grid" begin
        area = StudyArea(regime = bng(), verbosity = :silent)
        layer = materialise(bng(), area)
        # A **layer**, not a bare array: the array type is an implementation detail and does not
        # cross the public boundary. The values are reached through `matrix`, as on any layer.
        @test layer isa EcoSISTEM.AbstractLayer
        @test size(layer.matrix) == size(area.report.active)
        # Real coordinates, so it plots on true axes and can be compared with any other layer
        # cell-wise.
        # **The layer's coordinates are the area's, unit and all**, and the equality is asserted
        # directly. Nothing here may strip a unit off one side: the grid carries its own, so an
        # assertion that had to `ustrip` one side would be hiding a disagreement rather than
        # testing for one.
        northings = parent(DimensionalData.lookup(layer.matrix, Y))
        @test unit(eltype(northings)) == m
        @test northings == parent(DimensionalData.lookup(area.report.active, Y))
    end

    @testset "a synthetic layer is generated at the grid's shape" begin
        # The case that is otherwise invisible: a gradient has no coordinates of its own, so the only
        # way to check its direction and range is to put it on the grid and look.
        area = StudyArea(regime = bng(), verbosity = :silent)
        grad = materialise(GradientSpec(0.0K, 100.0K, axis = Temperature),
                           area)
        @test size(grad.matrix) == size(area.report.active)
        @test extrema(ustrip.(grad.matrix)) == (0.0, 100.0)
    end

    @testset "`role` applies the final conversion, and nothing else does" begin
        area = StudyArea(regime = bng(), verbosity = :silent)
        spec = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
        # Role-free: the per-area rate as written, held in a regime layer.
        plain = materialise(spec, area)
        @test plain isa EcoSISTEM.AbstractLayer{EcoSISTEM.Condition}
        @test unit(first(plain.matrix)) == unit(1.0kJ / (km^2 * day))
        # As a Resource: a supply, converted to an absolute per-cell rate against the area's real
        # cell size. The role decides the layer's *role parameter*, not merely its numbers.
        resource = materialise(spec, area, role = EcoSISTEM.Resource)
        @test resource isa Supply{SolarRadiation}
        @test dimension(first(resource.matrix)) == dimension(1.0kJ / day)
    end

    # **The claim `materialise` makes is that inspection and building agree**, so assert it rather
    # than trusting that two call paths stay in step: the layer a user inspects must equal the one
    # the builder puts in the environment, cell for cell and type for type.
    @testset "what you see is what the builder gets" begin
        area = StudyArea(regime = bng(), verbosity = :silent)
        seen = materialise(bng(), area, role = EcoSISTEM.Condition)
        built = GridHabitat(regime = bng(),
                            supply = UniformSpec(1.0kJ / (km^2 * day),
                                                 axis = SolarRadiation),
                            area = area).regime
        @test typeof(seen) == typeof(built)
        @test seen.matrix == built.matrix
    end

    # **A synthetic area has no latitude, so no latitude correction.** A raw-value accessor that
    # dresses a CRS-less grid's cell *indices* as degrees sends a supply down the geographic branch:
    # a 1 kJ/(km^2*day) rate over 10 km cells then inspects as 100.365 kJ/day where the builder -
    # which multiplies by a bare `cellsize^2` - says 100.0. Asserted against the **builder** rather
    # than against 100.0, so the two cannot drift apart in either direction.
    @testset "a synthetic area gets no latitude correction" begin
        area = StudyArea(extent = (90.0km, 90.0km), cellsize = 10.0km,
                         verbosity = :silent)
        spec = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
        seen = materialise(spec, area, role = EcoSISTEM.Resource)
        built = GridHabitat(regime = UniformSpec(298.0K,
                                                 axis = Temperature),
                            supply = spec, area = area).supply
        @test first(seen.matrix) == 100.0kJ / day
        @test seen.matrix == built.matrix

        # ...and a geographic area still gets one, so the guard has not simply turned it off.
        geo = StudyArea(regime = _reg(_testraster(WorldClim{BioClim},
                                                  fill(291.0K, 9, 9),
                                                  lat = (50.0:1.0:58.0) .* °,
                                                  long = (0.0:1.0:8.0) .* °)),
                        verbosity = :silent)
        @test !all(isone,
                   EcoSISTEM._areafactor(EcoSISTEM._cellintervals(geo.report.active,
                                                                  Y)))
    end
end

@testset "when a ConstructedSpec's combine runs" begin
    # A 2.5 km fixture and a 5 km study grid over the same ground, so every target cell is built
    # from a block of source cells and the two orderings have something to disagree about.
    east = (245000.0:2500.0:252500.0) .* m
    north = (640000.0:2500.0:647500.0) .* m
    fine(values) = _testraster(WorldClim{BioClim}, values, lat = north,
                               long = east, crs = Rasters.EPSG(27700))
    # Growing-season water and growing-season length, the pair the whole option exists for: the
    # water is uniform, the season alternates, so their ratio varies where neither layer is
    # remarkable on its own.
    water = fine(fill(100.0, 4, 4))
    days = fine(repeat([50.0 100.0 50.0 100.0], 4, 1))
    ratio(a, b) = ClimateRaster(WorldClim{BioClim}, a.array ./ b.array)

    # The cache is pre-loaded with the fixtures rather than reading anything: `_asraster` hits
    # it before `_read`, so the two `SourceSpec`s below name real layers but never touch the disk.
    # That is what lets an ordering test use grids chosen for the arithmetic rather than whatever
    # WorldClim happens to ship.
    cache = LayerCache()
    wspec = SourceSpec(WorldClim{BioClim}, :bio12)
    dspec = SourceSpec(WorldClim{BioClim}, :bio1)
    cache.reads[ReadKey(wspec)] = water
    cache.reads[ReadKey(dspec)] = days
    # `_reg` only where the raster is used as a *layer*; the two cache entries above stay raw
    # rasters, because a cache holds reads and a read is a raster.
    # The area itself is kept, not just its grid: `_materialisefield` takes a `StudyArea`, and
    # rebuilding one inside a testset is both wasteful and a different object.
    area = StudyArea(regime = _reg(water), cellsize = 5000.0m,
                     verbosity = :silent)
    target = area.report.active

    @testset "the stage decides which grid the combine sees" begin
        # The contract itself, independent of any resampler: how many cells the combine is handed.
        seen = Ref(0)
        function counting(stage)
            spec = ConstructedSpec(wspec, axis = EcoSISTEM.NicheAxis,
                                   combinestage = stage) do layer
                seen[] = length(layer.array)
                return layer
            end
            _materialiseon(spec, target, cache)
            return seen[]
        end
        @test counting(CombineOnSourceGrid()) == length(water.array)
        @test counting(CombineOnTargetGrid()) == length(target)
        @test length(target) < length(water.array)   # or the two would not differ
    end

    @testset "a nonlinear combine gives different answers either way" begin
        late = _materialiseon(ConstructedSpec(ratio, wspec, dspec,
                                              axis = EcoSISTEM.NicheAxis,
                                              combinestage = CombineOnTargetGrid()),
                              target, cache)
        early = _materialiseon(ConstructedSpec(ratio, wspec, dspec,
                                               axis = EcoSISTEM.NicheAxis,
                                               combinestage = CombineOnSourceGrid()),
                               target, cache)
        # Division is cell-wise, so the old "cell-wise => safe to combine late" reading would call
        # both of these correct. They are not the same number: early is the area-weighted mean of
        # the two cells' daily rates, late the total water over the mean season length.
        @test parent(early.array) != parent(late.array)
        # Each is exactly its own ordering, spelled out - early combines then sags onto the grid...
        @test parent(early.array) ==
              parent(_sampledata(ratio(water, days), target, name = "layer",
                                 categorical = false))
        # ...late samples each layer and combines there.
        @test parent(late.array) ==
              parent(ratio(_materialiseon(wspec, target, cache),
                           _materialiseon(dspec, target, cache)).array)
    end

    # **`materialise` must give what the builder gives.** If it read the layers on their own grid
    # and combined there - `CombineOnSourceGrid`'s ordering - **whatever the spec said**, then on the
    # default stage inspecting a layer would show a different field from the one the model uses. The
    # testset above is what makes that observable: the two orderings genuinely differ here.
    @testset "materialise honours combinestage" begin
        # `_materialisefield` reads through the *area's own* cache, not the one the fixtures
        # preloaded - without this it would try a real read and divide `bio12` by `bio1`'s affine °C.
        area.report.cache.reads[ReadKey(wspec)] = water
        area.report.cache.reads[ReadKey(dspec)] = days
        for stage in (CombineOnTargetGrid(), CombineOnSourceGrid())
            spec = ConstructedSpec(ratio, wspec, dspec,
                                   axis = EcoSISTEM.NicheAxis,
                                   combinestage = stage)
            viaspec = _materialiseon(spec, target, cache)
            viaarea = EcoSISTEM._materialisefield(spec, area)
            @test isequal(parent(viaspec.array), parent(viaarea.values))
        end
        # ...and the two stages still differ, so the equality above is not vacuous.
        onto = _materialiseon(ConstructedSpec(ratio, wspec, dspec,
                                              axis = EcoSISTEM.NicheAxis,
                                              combinestage = CombineOnTargetGrid()),
                              target, cache)
        src = _materialiseon(ConstructedSpec(ratio, wspec, dspec,
                                             axis = EcoSISTEM.NicheAxis,
                                             combinestage = CombineOnSourceGrid()),
                             target, cache)
        @test !isequal(parent(onto.array), parent(src.array))
    end

    # **A synthetic layer has no grid of its own, so on the early path it adopts the one the data
    # layers agree on.** That is the whole rule: a generated layer needs shape and orientation, never
    # coordinates - it can take any grid it is given, it just cannot supply one.
    @testset "an early combine generates a synthetic layer at the data's grid" begin
        mixed = ConstructedSpec(wspec,
                                UniformSpec(2.0, axis = EcoSISTEM.NicheAxis),
                                axis = EcoSISTEM.NicheAxis,
                                combinestage = CombineOnSourceGrid()) do w, u
            return w .* u
        end
        out = _materialiseon(mixed, target, cache)
        # Against the same combine done by hand: the synthetic layer at the *data's* grid, doubled.
        @test parent(out.array) ==
              parent(_sampledata(ClimateRaster(WorldClim{BioClim},
                                               water.array .* 2.0), target,
                                 name = "l", categorical = false))
        # ...and the synthetic layer really was generated at the source grid, not the target: had it
        # been made at the target the combine would have been 2×2, and this is 4×4 before sampling.
        @test length(water.array) > length(target)

        # **All synthetic is the one case that cannot work** - nothing supplies a grid, and none
        # can be invented.
        allsynth = ConstructedSpec(UniformSpec(1.0, axis = EcoSISTEM.NicheAxis),
                                   UniformSpec(2.0, axis = EcoSISTEM.NicheAxis),
                                   axis = EcoSISTEM.NicheAxis,
                                   combinestage = CombineOnSourceGrid()) do a, b
            return a
        end
        @test_throws ErrorException _materialiseon(allsynth, target, cache)
        # ...while a spec with **no** layers is a thunk supplying its own raster, and must still work:
        # `all` of an empty collection is `true`, which is exactly how that broke once.
        @test_nowarn _materialiseon(_reg(water), target, cache)
    end

    @testset "the default is unchanged" begin
        # Everything already written keeps combining on the target grid, bit for bit.
        default = _materialiseon(ConstructedSpec(ratio, wspec, dspec,
                                                 axis = EcoSISTEM.NicheAxis),
                                 target,
                                 cache)
        late = _materialiseon(ConstructedSpec(ratio, wspec, dspec,
                                              axis = EcoSISTEM.NicheAxis,
                                              combinestage = CombineOnTargetGrid()),
                              target, cache)
        @test parent(default.array) == parent(late.array)
    end

    # Sampling the combine's *result* is what finally gives the declared **axis** something to do:
    # it is read before the resample rather than after it, so a derived class-code layer is taken by
    # nearest class instead of interpolated. Codes 1 and 3 in alternating columns make the
    # difference visible - interpolating them invents a 2, which is a class that occurs nowhere.
    # There is deliberately no separate `valuetype` keyword: the axis says it, and so cannot
    # contradict it.
    @testset "a declared axis reaches the resampler" begin
        codes = fine(repeat([1.0 3.0 1.0 3.0], 4, 1))
        function collapsed(axis)
            spec = ConstructedSpec(() -> codes, axis = axis,
                                   combinestage = CombineOnSourceGrid())
            out = _materialiseon(spec, target, cache)
            return Set(filter(!isnan, vec(parent(out.array))))
        end
        @test issubset(collapsed(LandCoverTypology), Set([1.0, 3.0]))
        # On a non-typology axis the same layer is interpolated and lands between the two classes.
        @test !issubset(collapsed(EcoSISTEM.NicheAxis), Set([1.0, 3.0]))
    end

    @testset "a combine's layers must share one grid" begin
        # The guard is on the read path, so it covers `_analyse` as well as the early collapse -
        # a `ConstructedSpec` states its own extent by combining its layers on their own grid.
        other = LayerCache()
        other.reads[ReadKey(wspec)] = water
        other.reads[ReadKey(dspec)] = _testraster(WorldClim{BioClim},
                                                  fill(1.0, 3, 3),
                                                  lat = (640000.0:5000.0:650000.0) .*
                                                        m,
                                                  long = (245000.0:5000.0:255000.0) .*
                                                         m,
                                                  crs = Rasters.EPSG(27700))
        spec = ConstructedSpec(ratio, wspec, dspec, axis = EcoSISTEM.NicheAxis)
        err = try
            _asraster(spec, other)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("must be on one grid", err)
        @test occursin("bio1", err)          # and it names the layer that does not fit
        # A pair that does share a grid is unaffected.
        @test _asraster(spec, cache) isa ClimateRaster
    end
end

@testset "the decision core" begin
    @testset "cell-boundary alignment" begin
        # A target starting a whole number of source cells along is aligned; half a cell along is not.
        @test _originaligned(0.0m, 0.0m, 100.0m)
        @test _originaligned(0.0m, 500.0m, 100.0m)
        @test _originaligned(245000.0m, 245500.0m, 250.0m)
        @test !_originaligned(0.0m, 50.0m, 100.0m)
        @test !_originaligned(245000.0m, 245125.0m, 250.0m)

        # On an angular lattice the same question is integer divisibility, so it needs no tolerance.
        arcsec(n) = uconvert(°, float(n) * EcoSISTEM.Units.arcsecond)
        @test _originaligned(-180.0°, arcsec(-647910), arcsec(30))   # 3 whole cells along
        @test !_originaligned(-180.0°, arcsec(-647895), arcsec(30))  # half a cell along
    end

    # The case that matters is a bound **already on** a cell boundary: it must snap to itself, or the
    # study area silently grows by a spurious cell - and *which* cell you get then depends on how much
    # of the layer was read, which is the read-extent variance the exact-lattice work exists to remove.
    # In float degrees `(bound - origin) / step` lands either side of the whole number it should be, so
    # `floor` and `ceil` disagree; over a global 30 arcsec lattice that happened for **7664 of 21601**
    # bounds. Exact origins and steps do not fix it - `range` elements in degrees cannot be exact - but
    # integer arcseconds do. Swept in aggregate: one invariant, not 21601 assertions.
    @testset "a bound already on a cell boundary snaps to itself" begin
        arcsec(n) = uconvert(°, float(n) * EcoSISTEM.Units.arcsecond)
        o, s = -180.0°, arcsec(30)
        onlattice = [arcsec(-648000 + 30k) for k in 0:21600]
        @test all(EcoSISTEM._snap(b, o, s, floor) ==
                  EcoSISTEM._snap(b, o, s, ceil)
                  for b in onlattice)
        # ...and it snaps to *itself*, not merely to something self-consistent.
        @test all(EcoSISTEM._snap(b, o, s, floor) == b for b in onlattice)

        # A bound genuinely off the lattice must still be grown outwards, never shrunk inwards.
        offlattice = [arcsec(-648000 + 30k + d)
                      for k in 0:2000, d in (7, 15, 23)]
        @test all(EcoSISTEM._snap(b, o, s, floor) <= b <=
                  EcoSISTEM._snap(b, o, s, ceil)
                  for b in offlattice)
        # Compared in arcseconds rather than by subtracting the two snapped degree values: both sit
        # near -180°, so differencing them cancels most of their significance and the result cannot
        # equal the step exactly. "Exactly one cell apart" is only an exact statement on the lattice.
        @test all(EcoSISTEM._arcsecs(EcoSISTEM._snap(b, o, s, ceil)) -
                  EcoSISTEM._arcsecs(EcoSISTEM._snap(b, o, s, floor)) == 30
                  for b in offlattice)
        # A projected lattice has no arcsecond form and keeps the float path.
        @test EcoSISTEM._snap(245125.0m, 245000.0m, 250.0m, floor) == 245000.0m
        @test EcoSISTEM._snap(245125.0m, 245000.0m, 250.0m, ceil) == 245250.0m
    end

    @testset "what a grid choice costs each layer" begin
        # Same CRS, same size, aligned: the layer is copied untouched.
        @test _resamplecost(true, 100.0m, 100.0m, true) isa
              EcoSISTEM.LayerKeptExactly
        # Same size but offset boundaries: every target cell straddles two source cells.
        @test _resamplecost(true, 100.0m, 100.0m, false) isa
              EcoSISTEM.LayerResampled
        # A whole multiple, aligned: exact aggregation, and the factor is reported.
        agg = _resamplecost(true, 100.0m, 500.0m, true)
        @test agg isa EcoSISTEM.LayerAggregated
        @test agg.factor == 5
        # A whole multiple but offset: aggregation would cross source cells, so it is a resample.
        @test _resamplecost(true, 100.0m, 500.0m, false) isa
              EcoSISTEM.LayerResampled
        # Not a whole multiple: resampled even when aligned.
        @test _resamplecost(true, 100.0m, 250.0m, true) isa
              EcoSISTEM.LayerResampled
        # Finer target than source - upsampling, which invents detail; the reason says so.
        up = _resamplecost(true, 500.0m, 100.0m, true)
        @test up isa EcoSISTEM.LayerResampled
        @test occursin("finer", up.reason)
        # A different CRS always costs a resample, whatever the sizes.
        diff = _resamplecost(false, 100.0m, 100.0m, true)
        @test diff isa EcoSISTEM.LayerResampled
        @test occursin("different CRS", diff.reason)
    end

    # Angular grids must be compared as whole arcseconds, not as degrees. Steps reach here already
    # snapped to the arcsecond lattice by `_snaparcsec`, but exact *equality* on those degree floats is
    # not enough: the question `_resamplecost` asks is whether one step is a whole multiple of another,
    # and the float ratio gets that wrong even for exactly-snapped inputs - 21 arcsec into 7 comes out
    # as 3.0000000000000004, and a sweep of every pair up to 200 arcsec finds 107 such cases. Each one
    # would be reported as a lossy `:resampled` when it is in fact an exact `:aggregated`. Integer
    # arcseconds answer it outright, so the cases below are the ones a float implementation fails.
    @testset "whole-multiple angular steps aggregate exactly" begin
        arcsec(n) = (n / 3600)°
        for (src, tgt, factor) in ((21, 7, nothing), (7, 21, 3), (3, 27, 9), (5, 25, 5),
            (13, 39, 3), (30, 90, 3), (30, 600, 20))
            cost = _resamplecost(true, arcsec(src), arcsec(tgt), true)
            if isnothing(factor)          # target finer than source: still a resample
                @test cost isa EcoSISTEM.LayerResampled
            else
                @test cost isa EcoSISTEM.LayerAggregated
                @test cost.factor == factor
            end
        end
        # Equality and the not-a-multiple verdict are unaffected by the change, so pin them too.
        @test _resamplecost(true, arcsec(30), arcsec(30), true) isa
              EcoSISTEM.LayerKeptExactly
        @test _resamplecost(true, arcsec(20), arcsec(30), true) isa
              EcoSISTEM.LayerResampled
        # A step off the arcsecond lattice has no integer form, so it falls back to the float ratio
        # rather than being silently rounded onto a neighbouring lattice point.
        @test isnothing(EcoSISTEM._arcsecs((0.4 / 3600)°))
        @test isnothing(EcoSISTEM._arcsecs(100.0m))     # projected: arcseconds are meaningless
        @test EcoSISTEM._arcsecs(arcsec(90)) === 90
    end

    # A cell size may be *written* in whatever angular unit reads best - the whole reason `arcminute`
    # and `arcsecond` exist - and must still take the exact integer path rather than silently falling
    # back to floating point. Unitful cannot answer "is this an angle?" from the dimension (angles
    # are `NoDims`, like a bare number) nor from convertibility (`uconvert(arcsecond, 0.5)` succeeds,
    # reading the number as radians); an angle is what is dimensionless *and* carries a unit.
    @testset "any angular unit is understood, not just degrees" begin
        @test EcoSISTEM._arcsecs(30arcsecond) == 30
        @test EcoSISTEM._arcsecs(10arcminute) == 600
        @test EcoSISTEM._arcsecs((1 / 6)°) == 600         # ...and all agree with each other
        @test EcoSISTEM._arcsecs(30.0arcsecond) ==
              EcoSISTEM._arcsecs(30arcsecond)
        # `30arcsecond` is an *`Int`* quantity, and `eps` has no `Integer` method - this is the
        # literal a caller actually writes, and it threw a bare `MethodError` before `float` was added.
        @test EcoSISTEM._arcsecs(1°) == 3600
        # Not everything dimensionless is an angle, and a length never is.
        @test isnothing(EcoSISTEM._arcsecs(0.5))
        @test isnothing(EcoSISTEM._arcsecs(2500.0u"km"))
        # A radian is an angle but not a whole number of arcseconds, so it is recognised and rejected.
        @test EcoSISTEM._isangle(1.0u"rad")
        @test isnothing(EcoSISTEM._arcsecs(1.0u"rad"))
        # An entirely arcsecond-denominated lattice snaps exactly, and keeps the origin's unit.
        lo = EcoSISTEM._snap(47arcsecond, 0arcminute, 30arcsecond, floor)
        hi = EcoSISTEM._snap(47arcsecond, 0arcminute, 30arcsecond, ceil)
        @test EcoSISTEM._arcsecs(lo) == 30 && EcoSISTEM._arcsecs(hi) == 60
        @test unit(lo) === arcminute
    end

    # Cell sizes are shown the way a person would write them. A geographic grid is laid out in whole
    # arcseconds, so `0.16666666666666666°` is really 10 arcmin; a derived projected size is floored
    # to one significant figure so it is a whole number of km or m rather than 14218.48446569977 m.
    @testset "cell sizes read as whole units" begin
        cst = EcoSISTEM._cellsizetext
        @test cst((1 / 6)°) == "10 ′"              # WorldClim
        @test cst((1 / 120)°) == "30 ″"            # EarthEnv / CHELSA
        @test cst(1.0°) == "1°"
        @test cst(2500.0u"m") == "2500 m"          # no pointless `.0`
        @test cst(10000.0u"m") == "10 km"          # ...and whole km read as km
        @test cst(5km) == "5 km"
        # Anything off the arcsecond lattice keeps its exact value rather than being dressed up.
        @test cst(0.1234°) == string(0.1234°)

        # Floored, never rounded to nearest: a grid must never come out coarser than the size
        # measured off the projection, or resolution the source carried is silently discarded.
        f1 = EcoSISTEM._floor1sf
        @test f1(14218.48446569977u"m") == 10.0km
        @test f1(927.66u"m") == 900.0u"m"
        @test f1(19999.0u"m") == 10.0km            # the worst case: ~4× the cells
        @test f1(8674.82u"m") == 8.0km
        for x in (14218.48446569977, 927.66, 19999.0, 8674.82, 1437.0, 3.7)
            @test f1(x * u"m") ≤ x * u"m"           # the invariant that matters
        end
    end

    @testset "choosing the layer to align to" begin
        bng, wgs = Rasters.EPSG(27700), Rasters.EPSG(4326)
        # Only layers already in the target CRS are candidates.
        facts = [(name = :regime, crs = wgs, step = 1.0m),
            (name = :supply, crs = bng, step = 5.0m)]
        @test _choosealign(facts, bng).name == :supply

        # Among candidates, the finest wins...
        facts = [(name = :regime, crs = bng, step = 100.0m),
            (name = :supply, crs = bng, step = 25.0m)]
        @test _choosealign(facts, bng).name == :supply

        # ...and an exact tie falls back to the order given, so the choice is deterministic.
        facts = [(name = :regime, crs = bng, step = 25.0m),
            (name = :supply, crs = bng, step = 25.0m)]
        @test _choosealign(facts, bng).name == :regime

        # No layer in the target CRS: nothing can be preserved exactly.
        facts = [(name = :regime, crs = wgs, step = 1.0m)]
        @test isnothing(_choosealign(facts, bng))
    end

    @testset "snapping bounds onto the alignment grid" begin
        # Bounds grow outwards to the alignment layer's own cell boundaries, so every target cell is
        # a whole number of source cells; shrinking instead could drop data the mask asked for.
        ragged = Extent(Y = (1010.0m, 1990.0m), X = (2010.0m, 2990.0m))
        snapped = Extent(Y = (1000.0m, 2000.0m), X = (2000.0m, 3000.0m))
        @test _snapbounds(ragged, (0.0m, 0.0m), 100.0m) == snapped
        # Already-aligned bounds are left alone.
        @test _snapbounds(snapped, (0.0m, 0.0m), 100.0m) == snapped
        # With no alignment layer there is nothing to snap to.
        @test _snapbounds(ragged, nothing, 100.0m) == ragged
    end

    @testset "grid footprint" begin
        fp = _gridfootprint(137, 81)
        @test fp.cells == 137 * 81
        # One environment layer is Float64 per cell; abundances are Int64 per species per cell.
        @test fp.layer == 137 * 81 * 8
        @test fp.perspecies == 137 * 81 * 8
        msg = _footprintmessage(fp)
        @test occursin("11097 cells", msg)
        @test occursin("per environment layer", msg)
        @test occursin("per species", msg)
        # A grid big enough to matter reports in sensible units rather than raw bytes.
        @test occursin("MiB", _footprintmessage(_gridfootprint(2000, 2000)))
    end

    @testset "a CRS used far outside its area of use is reported" begin
        bng = Rasters.EPSG(27700)
        # PROJ's database knows what ground each CRS was defined for.
        use = EcoSISTEM._areaofuse(bng)
        @test !isnothing(use)
        @test use.west≈-9.01 atol=0.5          # the UK, not the world
        @test use.north≈61.01 atol=0.5
        @test occursin("United Kingdom", use.name)

        # On its home ground, silence - a study area routinely overhangs its CRS's declared box a
        # little (BNG's stops at the coastline), and warning on that would be noise.
        quiet = EcoSISTEM.Problem[]
        EcoSISTEM._crsareaofuse!(quiet, bng,
                                 Extent(Y = (7.0e5m, 1.0e6m),
                                        X = (2.0e5m, 4.0e5m)))
        @test isempty(quiet)

        # Stretched across continental Europe it still *computes*, silently distorting distance and
        # area - exactly the failure that should not pass unremarked.
        loud = EcoSISTEM.Problem[]
        EcoSISTEM._crsareaofuse!(loud, bng,
                                 Extent(Y = (0.0m, 4.0e6m),
                                        X = (-1.0e6m, 3.0e6m)))
        @test only(loud).code === :crs_area_of_use
        @test occursin("outside the area", only(loud).message)
        @test occursin("EPSG:27700", only(loud).message)
        # ...and it names a better CRS rather than only complaining.
        @test occursin("crs = Rasters.EPSG(", only(loud).message)
    end

    # The grid convention, and what happens to the remainder - one subject in two halves, so they
    # share a fixture. **Non-square throughout**: the cell count is computed per axis, so a square
    # fixture cannot see a y/x swap in the `ceil` - the same class of bug as `[NICHE-YX]`.
    @testset "the target grid starts on the data's own edge" begin
        # 7 × 11 cells of 2.5 km, so the data spans 17.5 km north-south by 27.5 km east-west.
        raw = _bngraster(WorldClim{BioClim}, fill(291.0K, 7, 11),
                         north = (640000.0:2500.0:655000.0) .* m,
                         east = (245000.0:2500.0:270000.0) .* m)
        ns = _reg(raw)
        # A synthesised grid labels its cells by their lower corner and lays down exactly as many as
        # cover the span - `ceil(span / cellsize)`, with no spare row. Both counts are checked, and
        # they differ, so a transposition could not pass.
        p = _analyse((regime = ns,), cellsize = 5000.0m)
        @test size(p.active) == (3, 5)          # 17.5/5 -> 3 whole cells, 27.5/5 -> 5
        lk = DimensionalData.lookup(p.active, Y)
        @test DimensionalData.Lookups.locus(DimensionalData.Lookups.sampling(lk)) isa
              DimensionalData.Lookups.Start
        # **The property, not the size**: every surviving cell is *wholly* inside the data. It is
        # asserted through `_cellintervals` rather than by rebuilding edges from the coordinates,
        # which would re-make the assumption this change removed.
        e = EcoSISTEM._dimsextent(raw.array, nothing)
        yx = EcoSISTEM._cellintervals(p.active)
        @test all(iv -> e.Y[1] <= iv[1] && iv[2] <= e.Y[2], yx.lat)
        @test all(iv -> e.X[1] <= iv[1] && iv[2] <= e.X[2], yx.long)
        # ...and the first cell begins exactly on the data's own edge, rather than straddling it by
        # half a cell as the `Center` labelling did.
        @test first(yx.lat)[1] == e.Y[1]
        @test first(yx.long)[1] == e.X[1]

        # The coverage test itself must be blind to the labelling, or the rule would travel with
        # the convention it was chosen alongside. The *same cells*, described both ways.
        L = DimensionalData.Lookups
        mk(loc, D, v) = D(Rasters.Projected(collect(v),
                                            crs = Rasters.EPSG(27700),
                                            order = L.ForwardOrdered(),
                                            span = L.Regular(5000.0m),
                                            sampling = L.Intervals(loc)))
        norths = range(e.Y[1], step = 5000.0m, length = 4)
        easts = range(e.X[1], step = 5000.0m, length = 6)
        gs = Rasters.Raster(zeros(4, 6),
                            (mk(L.Start(), Y, norths),
                             mk(L.Start(), X, easts)))
        gc = Rasters.Raster(zeros(4, 6),
                            (mk(L.Center(), Y, norths .+ 2500.0m),
                             mk(L.Center(), X, easts .+ 2500.0m)))
        @test EcoSISTEM._fullycovered(raw, gs) ==
              EcoSISTEM._fullycovered(raw, gc)
        # 3 rows of 5 km fit in 17.5 km and 5 columns in 27.5 km; the fourth and sixth overhang.
        @test count(EcoSISTEM._fullycovered(raw, gs)) == 3 * 5
    end

    @testset "simulate_safely decides what happens to the remainder" begin
        # 9 × 9 of 2.5 km - 22.5 km square, which divides into 4 km and 6 km cells with a remainder
        # of 62.5% and 75% of a cell. **Those two are the cases where the rules disagree**: both
        # are more than half covered, so the centre test keeps them and the coverage test does not.
        sq = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
        for (cell, safe, loose) in ((4000.0m, 5, 6), (6000.0m, 3, 4))
            a = _analyse((regime = sq,), cellsize = cell)
            @test size(a.active) == (safe, safe)
            @test a.simulate_safely
            @test isempty([p.code
                           for p in a.problems if
                           p.code === :partly_covered])
            b = _analyse((regime = sq,), cellsize = cell,
                         simulate_safely = false)
            @test size(b.active) == (loose, loose)
            @test !b.simulate_safely
            # The looser rule announces itself, costed in cells - it is opting in to simulating
            # ground the data does not describe, so it must not be silent.
            note = only([p for p in b.problems if p.code === :partly_covered])
            @test note.severity isa EcoSISTEM.ProblemNotice
            @test occursin("partly backed by data", note.message)
            # Every cell the safe rule keeps is one the loose rule keeps too: it is strictly
            # tighter, never a different selection.
            @test count(a.active) < count(b.active)
        end

        # Where the data divides exactly there is no remainder and the two rules agree - the case a
        # tolerance has to survive, since the last cell's far edge is reached by different
        # arithmetic on each side and can differ by an ULP.
        for cell in (2500.0m, 5000.0m)
            @test size(_analyse((regime = sq,), cellsize = cell).active) ==
                  size(_analyse((regime = sq,), cellsize = cell,
                                simulate_safely = false).active)
        end

        # Cells larger than the data itself leave nothing wholly covered, and the error names the
        # flag rather than blaming the mask and the CRS as the general one does.
        err = try
            _analyse((regime = sq,), cellsize = 40.0km)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("simulate_safely = false", err.msg)

        # Tri-state like every other constraint: `missing` inherits from `base`, `nothing` clears
        # it back to the default, and a `Bool` is used as given.
        base = _analyse((regime = sq,), cellsize = 4000.0m,
                        simulate_safely = false)
        @test size(base.active) == (6, 6)
        inherited = StudyArea(base, verbosity = :silent)
        @test !inherited.report.simulate_safely
        @test size(inherited.report.active) == (6, 6)
        @test size(StudyArea(base, simulate_safely = true,
                             verbosity = :silent).report.active) == (5, 5)
        # `nothing` is not "false" - it discards the inherited value, so the default applies.
        cleared = StudyArea(base, simulate_safely = nothing,
                            verbosity = :silent)
        @test cleared.report.simulate_safely
        @test size(cleared.report.active) == (5, 5)
    end

    # The rule reaches a layer that was never named on the area, which is what makes it a
    # property of the *simulation* rather than of the grid decision.
    @testset "simulate_safely applies to a layer introduced at build time" begin
        big = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
        # **7 cells of 2 km, not 5 of 2.5 km.** The obvious fixture does *not* discriminate: its
        # footprint ends exactly on a target cell's centre, which the half-open centre test excludes,
        # so both rules answered 4 and the test passed for the wrong reason. Measured: this one
        # ends at 652 750 m, four-fifths of the way through the third 5 km cell - covered by the
        # centre test, not wholly covered.
        small = _reg(_bngraster(WorldClim{BioClim}, fill(295.0K, 7, 7),
                                north = (639750.0:2000.0:651750.0) .* m,
                                east = (244750.0:2000.0:256750.0) .* m))
        area = StudyArea(regime = big, cellsize = 5000.0m, verbosity = :silent)
        @test size(area.report.active) == (4, 4)
        # Only the 2 × 2 block the small layer wholly covers survives.
        @test count(EcoSISTEM._coverage(materialise(small, area))) == 4
        # Without the flag the partly-covered ring comes back, because the centre test keeps it.
        loose = StudyArea(regime = big, cellsize = 5000.0m,
                          simulate_safely = false, verbosity = :silent)
        @test count(EcoSISTEM._coverage(materialise(small, loose))) == 9
    end

    # The regression that `extras_examples` caught: `_fullycovered` compared footprints in the
    # **target's** CRS, so a near-global WGS84 layer re-expressed in British National Grid came back
    # as a 125 km strip of eastings and most of Scotland read as outside the data.
    @testset "a layer in another CRS is judged in its own coordinates" begin
        # A BNG study area, and a global WGS84 layer that obviously contains it.
        bng = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 11),
                              east = (245000.0:2500.0:270000.0) .* m))
        area = StudyArea(regime = bng, cellsize = 5000.0m, verbosity = :silent)
        # **Both bounds matter, and I got the longitude wrong once.** Latitude stops at 90° rather
        # than reaching past it, because a cell edge beyond the pole makes PROJ refuse the (wrong-way)
        # transform outright instead of returning the nonsense this test is about. And longitude must
        # end exactly at **180°**, not past it: measured, labels `-175:10:175` put the last cell's
        # far edge at **185°**, and that wrap alone changed the bogus envelope from
        # `256 685...400 000 m` (which excludes most of the grid) to `-101 403...400 000 m` (which
        # contains all of it) - so the fixture silently stopped reproducing the bug at all.
        global_ = _testraster(WorldClim{BioClim}, fill(290.0K, 14, 36),
                              lat = (-50.0:10.0:80.0) .* °,
                              long = (-180.0:10.0:170.0) .* °)
        # Every cell is covered - the layer spans the world.
        @test all(EcoSISTEM._fullycovered(global_, area.report.active))
        # **The regression itself, executable** - not a threshold in metres, which would only be a
        # number I chose. Asked the wrong way round the footprint is computed in the *target's* CRS,
        # and the very same per-cell test then drops cells the layer plainly covers. This is what
        # deleted 1859 of `ScottishCultivatedLand.jl`'s 3168 active cells, and it is what would come
        # back if `_fullycovered` were ever repointed at `_dimsextent(A, tcrs)`.
        wrongway = EcoSISTEM._dimsextent(global_.array,
                                         Rasters.crs(area.report.active))
        @test !all(EcoSISTEM._cellsinbox(EcoSISTEM._cellintervals(area.report.active,
                                                                  Y),
                                         EcoSISTEM._cellintervals(area.report.active,
                                                                  X), wrongway))

        # ...and a layer that genuinely covers only part of the area still says so, so the fix is not
        # simply "always true".
        corner = _testraster(WorldClim{BioClim}, fill(290.0K, 3, 3),
                             lat = (55.7:0.05:55.8) .* °,
                             long = (-4.5:0.05:-4.4) .* °)
        part = EcoSISTEM._fullycovered(corner, area.report.active)
        @test any(part) && !all(part)
    end

    @testset "Problem's severity is a type, not a name" begin
        @test Problem(EcoSISTEM.ProblemWarning(), :upsampling, "...").severity isa
              EcoSISTEM.ProblemWarning
        # An unrecognised severity is refused by the signature, where it is written, rather than by
        # a check inside the constructor listing the names it accepts.
        @test_throws MethodError Problem(:oops, :upsampling, "...")
    end
end

@testset "LayerCache" begin
    @testset "read identity, not spec identity, is the key" begin
        # Two separately-constructed specs for the same layer must share one cache entry - keying on
        # the spec objects would read twice.
        a = SourceSpec(WorldClim{BioClim}, :bio1)
        b = SourceSpec(WorldClim{BioClim}, :bio1)
        @test ReadKey(a) == ReadKey(b)
        @test hash(ReadKey(a)) == hash(ReadKey(b))

        # Different layer, different dataset and different read keywords are all distinct reads.
        @test ReadKey(a) != ReadKey(SourceSpec(WorldClim{BioClim}, :bio12))
        @test ReadKey(a) !=
              ReadKey(SourceSpec(EarthEnv{LandCover}, :cultivated_and_managed))
        @test ReadKey(SourceSpec(WorldClim{Climate}, :wind, month = 1:12)) !=
              ReadKey(SourceSpec(WorldClim{Climate}, :wind, month = 1))

        # A heap-allocated read keyword (a range) must still hash by value - the `===`/`objectid`
        # fallback would miss every hit here, which is why `==`/`hash` are defined explicitly.
        k1 = ReadKey(SourceSpec(WorldClim{Climate}, :wind, month = 1:12))
        k2 = ReadKey(SourceSpec(WorldClim{Climate}, :wind, month = 1:12))
        @test k1 == k2
        @test hash(k1) == hash(k2)
        @test length(Dict(k1 => 1, k2 => 2)) == 1
    end

    @testset "a cached read happens once" begin
        cache = LayerCache()
        @test length(cache) == 0

        spec = SourceSpec(WorldClim{BioClim}, :bio1)
        first = _asraster(spec, cache)
        @test length(cache) == 1
        @test haskey(cache, spec)

        # The same read, and an equal-but-distinct spec, both hit rather than re-reading.
        @test _asraster(spec, cache) === first
        @test _asraster(SourceSpec(WorldClim{BioClim}, :bio1), cache) === first
        @test length(cache) == 1

        # A different layer is a different read.
        _asraster(SourceSpec(WorldClim{BioClim}, :bio12), cache)
        @test length(cache) == 2
    end

    @testset "a ConstructedSpec caches its children, not itself" begin
        # Its `combine` is an opaque closure, so a top-level entry could never be shared; the reads
        # underneath it are where the cost is, and those are cached.
        cache = LayerCache()
        spec = ConstructedSpec(SourceSpec(WorldClim{BioClim}, :bio1),
                               axis = EcoSISTEM.NicheAxis) do layer
            return layer
        end
        _asraster(spec, cache)
        @test length(cache) == 1
        @test haskey(cache, SourceSpec(WorldClim{BioClim}, :bio1))
    end

    # **A raster is not a spec.** Letting an in-memory raster pass straight through is the one way
    # a layer could reach a simulation without declaring what it means: a raster's axis is resolved
    # from its layer *code* through the shipped catalogue, so a derived raster (which has none)
    # silently becomes `EcoSISTEM.NicheAxis`, with no
    # keyword anywhere to correct it.
    @testset "a raster is refused as a layer, and the spec form is not" begin
        cache = LayerCache()
        raster = _asraster(SourceSpec(WorldClim{BioClim}, :bio1), LayerCache())
        @test_throws ErrorException _asraster(raster, cache)
        @test_throws ErrorException _asraster(raster)
        # Both entry points refuse it - deciding a grid and building a layer are separate paths,
        # and a message on only one of them would leave the other failing as a `MethodError`.
        # Building a layer is `materialise`: the builder has no entry point of its own, so this is
        # the one a raster passed as a `regime`/`supply` reaches.
        @test_throws ErrorException materialise(raster,
                                                StudyArea(regime = _reg(raster),
                                                          verbosity = :silent))
        # ...and the error says what to do rather than naming an internal function.
        msg = try
            _asraster(raster, cache)
        catch e
            sprint(showerror, e)
        end
        @test occursin("ConstructedSpec", msg) && occursin("axis", msg)

        # The wrapped form works, and still needs no cache entry: a nullary combine reads nothing.
        @test _asraster(_reg(raster), cache) === raster
        @test length(cache) == 0
    end

    # **The contract is in the signature**, so `methods(GridHabitat)` and the rendered docs
    # both show what a builder takes. The price, taken deliberately: Julia does not dispatch on
    # keyword types, so a wrong `regime`/`supply` is refused by a `TypeError` naming the keyword and
    # printing the union, not by a message offering a remedy - and no fallback method could improve
    # on that, since a second method with the same positional signature replaces the first.
    @testset "regime and supply are typed in the signature" begin
        raster = _asraster(SourceSpec(WorldClim{BioClim}, :bio1), LayerCache())
        sun = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
        area = StudyArea(regime = _reg(raster), verbosity = :silent)

        # Each of the three keyword-taking faces refuses at the keyword.
        @test_throws TypeError StudyArea(regime = raster, verbosity = :silent)
        @test_throws TypeError investigate_study_area(regime = raster)
        @test_throws TypeError GridHabitat(regime = _reg(raster),
                                           supply = raster, area = area)
        # ...including a bare value, the commonest form inherited from v0.4.0.
        @test_throws TypeError GridHabitat(regime = _reg(raster),
                                           supply = 3.0m / s, area = area)

        # `materialise`'s spec is **positional**, so it *could* refuse by `MethodError` - and it
        # must not. A positional argument can carry a fallback method where an annotated keyword
        # cannot, so the two cases below get the tailored message and the remedy instead of a method
        # list. Both messages are checked, not just the fact of an error: a bare value and a raster
        # fail for *different* reasons and name different fixes.
        @test_throws "has no niche axis" materialise(3.0K, area)
        @test_throws "ConstructedSpec(() -> raster; axis = SomeAxis)" materialise(raster,
                                                                                  area)

        # A signature cannot see inside a container, and a multi-layer regime legitimately is a
        # tuple - so a bad *element* is still caught later, with the tailored message.
        err = try
            GridHabitat(regime = (raster, _reg(raster)), supply = sun,
                        area = area)
        catch e
            sprint(showerror, e)
        end
        @test occursin("ConstructedSpec", err)

        # `Varying` is an `AbstractSpec`, which is what keeps `LayerInput` this short.
        @test Varying <: EcoSISTEM.AbstractSpec
        @test StudyArea(regime = Varying(_reg(raster),
                                         IncrementBy(0.5K / year)),
                        verbosity = :silent) isa StudyArea
        # ...and a pre-built supply is still accepted, because a `Supply{A}` carries its own axis.
        @test EcoSISTEM.AbstractSupply <:
              Union{EcoSISTEM.LayerInput, EcoSISTEM.AbstractSupply}
    end
end

@testset "Varying at the build boundary" begin
    bng() = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
    warming = IncrementBy(0.5K / year)
    sun = UniformSpec(1.0kJ / km^2 / day, axis = SolarRadiation)

    @testset "a wrapper is invisible to the study area" begin
        # The regression this guards is silent, not loud. `_probecrs`'s `::Any` fallback returns
        # `nothing` for an unrecognised value, which makes `_probetargetcrs` decline for the *whole*
        # area, which disables read windowing - every layer then read at full global extent, with no
        # error at all. A change says nothing about the grid, so the two areas must agree exactly.
        naked = StudyArea(regime = bng(), verbosity = :silent)
        wrapped = StudyArea(regime = Varying(bng(), warming),
                            verbosity = :silent)
        @test wrapped.report.crs == naked.report.crs
        @test !isnothing(wrapped.report.crs)
        @test wrapped.report.cellsize == naked.report.cellsize
        @test wrapped.report.align == naked.report.align
        @test size(wrapped.report.active) == size(naked.report.active)
        @test [l.kind for l in wrapped.report.layers] ==
              [l.kind for l in naked.report.layers]

        # ...and the same holds for a wrapped *supply*, which passes through the identical
        # analyse-time funnel (`_expandspecs` -> `_shapesgrid` -> `_probecrs`).
        supplied = StudyArea(regime = bng(),
                             supply = Varying(bng(), warming),
                             verbosity = :silent)
        @test supplied.report.crs == naked.report.crs

        # `materialise` ignores a declared change too: it shows what the spec becomes on the grid,
        # and the change is hung on the layer by `GridHabitat` afterwards.
        @test isequal(materialise(Varying(bng(), warming), naked).matrix,
                      materialise(bng(), naked).matrix)
    end

    @testset "the declared change lands on the built layer" begin
        area = StudyArea(regime = bng(), verbosity = :silent)
        env = GridHabitat(regime = Varying(bng(), warming), supply = sun,
                          area = area)
        @test env.regime.change isa SteadyLayerChange
        EcoSISTEM._layerupdate!(env.regime, 1.0year, 1.0year)
        @test all(≈(291.5K), env.regime.matrix)

        # A naked spec is untouched - the wrapper is the only way to declare a change.
        plain = GridHabitat(regime = bng(), supply = sun, area = area)
        @test plain.regime.change isa NoLayerChange
    end

    @testset "each layer of a collection keeps its own change" begin
        # **Named**: both layers come from the same `bng()` fixture and so share an axis.
        area = StudyArea(regime = (a = bng(), b = bng()), verbosity = :silent)
        env = GridHabitat(regime = (a = Varying(bng(), warming), b = bng()),
                          supply = sun, area = area)
        @test length(values(env.regime)) == 2
        @test env.regime.a.change isa SteadyLayerChange
        @test env.regime.b.change isa NoLayerChange
        EcoSISTEM._layerupdate!(env.regime, 1.0year, 1.0year)
        @test all(≈(291.5K), env.regime.a.matrix)
        @test all(≈(291.0K), env.regime.b.matrix)
    end
end

# A layer whose accumulation period is *another layer* - `gsp` over `gsl` - is an amount as a
# regime and a rate as a supply. Only the supply reading needs anything: the spec is rewritten into
# the division, and that is what makes `gsp` usable as a water supply at all.
@testset "a per-cell accumulation period desugars in supply position" begin
    E = EcoSISTEM
    gsp = SourceSpec(CHELSA{BioClimPlus}, :gsp)

    @testset "the rewrite" begin
        out = E._desugarsupply(gsp)
        @test out isa ConstructedSpec
        @test [l.code for l in out.layers] == [:gsp, :gsl]
        # `Precipitation`, not `GrowingSeasonPrecipitation`: a Resource-role axis says *which
        # resource*, and water is water. It also matches what `_wrapsupply` builds regardless.
        @test out.axis == Precipitation
        # Not a preference - division is cell-wise but nonlinear, so it must precede regridding.
        @test out.combinestage isa CombineOnSourceGrid

        # Driven by the catalogue's `percell=`, not by the code name: a layer without one is
        # returned *identically*, so the rewrite cannot fire on anything that does not declare it.
        plain = SourceSpec(WorldClim{Climate}, :prec)
        @test E._desugarsupply(plain) === plain
        # Idempotent, which is what lets `StudyArea` and `GridHabitat` both desugar safely.
        @test E._desugarsupply(out) === out
        # A tuple desugars member-wise and keeps its names.
        pair = E._desugarsupply((water = gsp, sun = plain))
        @test pair.water isa ConstructedSpec && pair.sun === plain
    end

    # The arithmetic, on fixtures rather than downloads: the cache is pre-loaded under each
    # spec's own `ReadKey`, so these name real layers but never touch the disk - the same trick the
    # combine-stage tests use, and what lets the grids be chosen for checkable numbers.
    @testset "the division, and the empty-season cell" begin
        east = (245000.0:2500.0:250000.0) .* m
        north = (640000.0:2500.0:645000.0) .* m
        fix(v) = _testraster(CHELSA{BioClimPlus}, v, lat = north, long = east,
                             crs = Rasters.EPSG(27700))
        # 100 mm of season water; a season that is 50 days in one column and 100 in another, so the
        # ratio differs where neither layer is remarkable alone. One cell has *no* season at all.
        water = fix(fill(100.0, 3, 3))
        days = fix([50.0 100.0 0.0; 50.0 100.0 0.0; 50.0 100.0 0.0])

        cache = LayerCache()
        cache.reads[ReadKey(gsp)] = water
        cache.reads[ReadKey(SourceSpec(CHELSA{BioClimPlus}, :gsl))] = days

        out = E._desugarsupply(gsp)
        # `_reg` for the *layer* use only - the cache entries above stay raw rasters.
        target = StudyArea(regime = _reg(water), cellsize = 2500.0m,
                           verbosity = :silent).report.active
        got = _materialiseon(out, target, cache)
        vals = parent(got.array)

        @test all(≈(2.0), vals[:, 1])       # 100 / 50
        @test all(≈(1.0), vals[:, 2])       # 100 / 100
        # A cell with no growing season has no growing-season water, and must read `NaN` so the
        # coverage check marks it inactive.
        @test all(isnan, vals[:, 3])
        # **This is why `_perperiod` exists and must not be "simplified" back to `./`.** Plain
        # division gives `Inf`, not `NaN` - only `0/0` is NaN - and `Inf` is far worse, because
        # nothing masks it: the cell would report an *infinite* water supply and be perfectly valid.
        @test isinf(100.0 / 0.0) && !isnan(100.0 / 0.0)
        @test isnan(EcoSISTEM._perperiod(100.0, 0.0))
        # ...and it is a guard on the divisor, not a blanket NaN: ordinary cells are untouched.
        @test EcoSISTEM._perperiod(100.0, 50.0) == 2.0
        # It carries units, so the guard works on the real unitful path too, not just the fixtures.
        @test isnan(ustrip(EcoSISTEM._perperiod(100.0u"L/m^2", 0.0u"d")))
        @test unit(EcoSISTEM._perperiod(100.0u"L/m^2", 2.0u"d")) == u"L/m^2/d"
    end
end

# **The reason `build_habitat` can reseed from a habitat at all: the habitat's own as-built grid
# survives the round trip.** A supply the study area never saw can cost cells, and that narrowing is
# recorded *nowhere else* - so re-deriving a `StudyArea` from the report's `specs` answers the
# original question and hands back the wider mask. Measured before the fix at 48 cells against 46.
# `report.stage` is what makes the two distinguishable, and this is its first real consumer.
@testset "an as-built area round-trips, narrowing and all" begin
    reg = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)),
               axis = Temperature)
    # A supply with holes the area cannot know about, because it is named only to the habitat.
    holed = fill(1.0kJ / (km^2 * day), 9, 9)
    holed[1, 1] = NaN * unit(eltype(holed))
    holed[2, 3] = NaN * unit(eltype(holed))
    area = StudyArea(regime = reg, verbosity = :silent)
    h = @test_logs (:warn, r"can only remove cells") match_mode=:any GridHabitat(regime = reg,
                                                                                 supply = _reg(_bngraster(WorldClim{BioClim},
                                                                                                          holed),
                                                                                               axis = SolarRadiation),
                                                                                 area = area)
    @test h.area.report.stage isa EcoSISTEM.AsBuilt
    @test count(h.area.report.active) == count(area.report.active) - 2

    # No keyword: the report is taken **verbatim**, narrowing, problems and stage included.
    back = StudyArea(h, verbosity = :silent)
    @test count(back.report.active) == count(h.area.report.active)
    @test back.report.stage isa EcoSISTEM.AsBuilt
    @test length(back.report.problems) == length(h.area.report.problems)
    # ...and it survives being copied again, which resetting the stage would have broken.
    @test count(StudyArea(back, verbosity = :silent).report.active) ==
          count(h.area.report.active)

    # **Any keyword drops through to the ordinary path**, deliberately and without a warning:
    # overriding is an explicit act. The re-derived area answers the *area's* question again, so the
    # narrowing is gone and the stage falls back to `AsInvestigated`.
    redone = StudyArea(h, simulate_safely = false, verbosity = :silent)
    @test redone.report.stage isa EcoSISTEM.AsInvestigated
    @test count(redone.report.active) == count(area.report.active)

    # An `AsInvestigated` base is untouched by any of this - nothing was narrowed, so re-deriving
    # reaches the same answer and there is nothing a copy could preserve.
    @test count(StudyArea(area, verbosity = :silent).report.active) ==
          count(area.report.active)

    # ...and the whole point: `build_habitat` inherits the narrowed grid, not the wider one.
    rebuilt = build_habitat(h, verbosity = :silent)
    @test count(rebuilt.active) == count(h.active)
end
# **An as-built area must still be able to READ DATA, and for a while it could not.** Three
# decisions composed into a hole: a `GridHabitat` discards its report's reads (`cache === nothing`,
# deliberately - they are consumed inputs and keeping them pins every raster the build touched); the
# copy constructor above took the report *verbatim*, cache included; and `_materialisefield` hands
# `area.report.cache` straight to `_asraster`/`_combineon`, which have methods for a `LayerCache`
# and none for `nothing`.
# **What made it invisible**: rebuilding with a *synthetic* spec worked perfectly, and so did
# inheriting everything - a built layer is passed through, never materialised. Only a **data-backed**
# spec on a habitat-derived area reached the cache at all, and that was the one combination no test
# covered. All nine gates were green over it.
@testset "an as-built area can still read data" begin
    reg = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)),
               axis = Temperature)
    sun = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
    area = StudyArea(regime = reg, verbosity = :silent)
    h = GridHabitat(regime = reg, supply = sun, area = area)

    # The habitat's own report keeps its reads discarded - that is 6e's decision and it stands.
    @test isnothing(h.area.report.cache)
    # But a `StudyArea` is a thing you *build on*, so the copy gets a fresh, empty cache: the
    # reads really are gone, and "nothing cached yet" is both true and usable.
    @test StudyArea(h, verbosity = :silent).report.cache isa
          EcoSISTEM.LayerCache

    # The regression itself - this raised `MethodError: no method matching _combineon(..., ::Nothing)`.
    values = fill(2.0kJ / (km^2 * day), 9, 9)
    datasupply = _reg(_bngraster(WorldClim{BioClim}, values),
                      axis = SolarRadiation)
    reread = build_habitat(h, supply = datasupply, verbosity = :silent)
    @test reread.supply isa Supply{SolarRadiation}
    @test all(≈(2.0kJ / (km^2 * day) *
                first(EcoSISTEM.getcellareas(reread))), reread.supply.matrix)
    # ...on the same grid it inherited, which is the point of reseeding from the habitat.
    @test count(reread.active) == count(h.active)

    # ...and the synthetic case still works, which is what passed while the above did not.
    @test build_habitat(h, supply = sun, verbosity = :silent) isa GridHabitat
end

end
