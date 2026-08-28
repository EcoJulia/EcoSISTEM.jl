# SPDX-License-Identifier: LGPL-3.0-or-later

module TestGridHabitat

using EcoSISTEM
# `public` but not exported, so `using EcoSISTEM` alone does not bring it into scope.
using EcoSISTEM: getcellareas
using Unitful
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using RasterDataSources
using DimensionalData: DimArray, X, Y, Ti
using EcoSISTEM: hasdata, landcoverclass
using EcoSISTEM: materialise
using Unitful, Unitful.DefaultSymbols
using DimensionalData: DimensionalData, DimArray, X, Y, dims
using Rasters
using ArchGDAL
using Distributions: Normal
using Extents: Extent
include("rasterfixtures.jl")
include("buildfixtures.jl")

# **Six testsets removed 2026-08-13** — `simplehabitat`, `simplenichehabitat`, `erahabitat`,
# `worldclimhabitat`, `bioclimhabitat` and `landcoverhabitat`. They tested **deprecated builders**,
# whose home is `test_deprecations.jl`; that file already covers the `*AE → *habitat` chain and both
# gradient families, so no current behaviour lost coverage. They also read ERA, twelve months of
# WorldClim climate, BioClim and EarthEnv land cover between them, which is what made this file
# expensive.
#
# `getraster(WorldClim{Climate}, :wind, month = 1:12)` and `getraster(EarthEnv{LandCover})` went
# with them: only those testsets needed them, and `test_datasetread.jl` fetches both regardless — so the
# CI cache is unchanged and only the processing is saved.
if !Sys.iswindows()
    # `bio1` is what the one surviving testset here reads.
    getraster(WorldClim{BioClim}, :bio1)

    @testset "GridHabitat from source type" begin
        supply = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
        # a SourceSpec reads one layer and attaches its unit from the shipped layer table
        # (defaulted here), so bio1 comes back as a continuous temperature (K) regime with no
        # manual unit juggling
        spec = SourceSpec(WorldClim{BioClim}, 1)
        @test spec isa SourceSpec
        @test spec.unit ==
              EcoSISTEM.layerunit(WorldClim{BioClim}, 1)
        # The grid is decided first, by the study area; the builder only samples onto it.
        area = StudyArea(regime = spec, verbosity = :silent)
        senv = GridHabitat(regime = spec, supply = supply, area = area)
        @test senv isa GridHabitat
        @test EcoSISTEM.iscontinuous(senv.regime)
        @test eltype(senv.regime.matrix) <: Unitful.Temperature

        # **The `(source, layer)` pair form is refused**, at every entry point and with the same
        # message: a tuple `regime` means a *multi-layer* environment, so a pair inside one could be
        # told from a layer by nothing but nesting depth. The remedy is `SourceSpec`, asserted just
        # below to work.
        for bad in (((WorldClim{BioClim}, 1), (WorldClim{BioClim}, 1)),
                    (SourceSpec(WorldClim{BioClim}, 1),
                     (WorldClim{BioClim}, 1)))
            msg = try
                GridHabitat(regime = bad, supply = supply, area = area)
                ""
            catch e
                sprint(showerror, e)
            end
            @test occursin("SourceSpec", msg) && occursin("multi-layer", msg)
        end
        # …and deciding a grid refuses it too, rather than accepting what the builder will not.
        @test_throws ErrorException StudyArea(regime = ((WorldClim{BioClim}, 1),
                                                        (WorldClim{BioClim}, 1)),
                                              verbosity = :silent)

        # a tuple of SourceSpecs is accepted too (same widened builder)
        # **Named**, because both are `bio1` and so on one axis — derived names cannot tell two
        # `Temperature` layers apart, and the collection says so rather than guessing.
        smenv = GridHabitat(regime = (first = SourceSpec(WorldClim{BioClim}, 1),
                                      second = SourceSpec(WorldClim{BioClim},
                                                          1)),
                            supply = supply, area = area)
        @test length(values(smenv.regime)) == 2
    end
end

# R3/R6, decided with the user 2026-08-05: masking is the `active` array's job and nothing else's.
# A habitat keeps every cell's supply, including the inactive ones, so a cell that is reactivated
# later still has resource — and `totalsupply` applies the mask itself rather than relying on
# the values having been destroyed.
#
# The bug this prevents is not "a reactivated cell is barren". `death_resource` is `E / K`, so a
# zeroed supply in a cell holding individuals gives `Inf` and a death probability of exactly 1: the
# cell would kill everything that dispersed into it, every step, silently.
# **Converted 2026-08-13 off a hand-built `active` matrix.** What this needs is *an inactive
# cell*, not an arbitrary mask, and a `CircleMaskSpec` gives one honestly: cells outside the disc are
# genuinely inactive, decided by the study area rather than asserted afterwards. That keeps the
# testset working when `GridHabitat` becomes the sole way to make a `GridHabitat`.
@testset "an inactive cell keeps its supply; `available` applies the mask" begin
    # A 4 × 4 km grid of 1 km cells, with a disc about the centre leaving the corners outside it.
    area = StudyArea(extent = (4.0km, 4.0km), cellsize = 1.0km,
                     within = CircleMaskSpec(radius = 1.9km),
                     verbosity = :silent)
    habitat = GridHabitat(regime = UniformSpec(298.0K,
                                               axis = Temperature),
                          supply = UniformSpec(2.0kJ / (km^2 * day),
                                               axis = SolarRadiation),
                          area = area)

    inactive = count(!, habitat.active)
    usable = count(habitat.active)
    # The mask must actually exclude something, or everything below passes vacuously.
    @test inactive > 0

    # The values survive in **every** cell, including the inactive ones — this is what makes
    # reactivation possible at all, and what stops `E / K` becoming `Inf` in a cell holding
    # individuals.
    persupply = 2.0kJ / (km^2 * day) * first(getcellareas(habitat))
    @test all(habitat.supply.matrix .≈ persupply)
    @test count(iszero, habitat.supply.matrix) == 0

    # …while "available" counts only the usable cells, not all sixteen.
    @test EcoSISTEM.totalsupply(habitat).SolarRadiation ≈ usable * persupply
    @test EcoSISTEM.totalsupply(habitat).SolarRadiation <
          sum(habitat.supply.matrix)

    # A single supply is reported as a one-entry `NamedTuple`, keyed by its **axis** —
    # `:SolarRadiation`, never a role label such as `:supply`. The key must say *what the resource
    # is*, because that is what makes it match the name the species' demand derives for the same
    # axis; a role label would name the same thing on both sides and pair nothing.
    @test EcoSISTEM.totalsupply(habitat) isa NamedTuple
    @test keys(EcoSISTEM.totalsupply(habitat)) == (:SolarRadiation,)

    # A data gap is still zeroed, and that is a different job: `NaN` is "no data", which no mask
    # can rescue, and it would otherwise reach the resource arithmetic. Asserted on the built
    # habitat's own supply, which is the only route by which one arrives.
    # `_zerogaps` **returns** a cleaned layer rather than mutating in place — it rebuilds the
    # layer so the change's stored slices are cleaned with it — so the result is what to assert on.
    habitat.supply.matrix[2, 2] = NaN * unit(persupply)
    cleaned = EcoSISTEM._zerogaps(habitat.supply)
    @test !any(isnan, cleaned.matrix)
    @test iszero(cleaned.matrix[2, 2])
end

# The multi-resource shape had **no** coverage through this accessor before 2026-08-05 — every
# assertion elsewhere in this file uses a single-supply habitat — so the path that returns one total
# per resource was exercised by nothing.
@testset "available supply is one named total per resource" begin
    # Converted off a hand-built mask 2026-08-13, as above — the disc supplies the inactive cells.
    area = StudyArea(extent = (4.0km, 4.0km), cellsize = 1.0km,
                     within = CircleMaskSpec(radius = 1.9km),
                     verbosity = :silent)
    solar = 2.0kJ / (km^2 * day)
    rain = 3.0Unitful.L / (km^2 * day)
    named = GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                        supply = (sunlight = UniformSpec(solar,
                                                         axis = SolarRadiation),
                                  water = UniformSpec(rain,
                                                      axis = Precipitation)),
                        area = area)
    usable = count(named.active)
    cell = first(getcellareas(named))
    avail = EcoSISTEM.totalsupply(named)

    # Keyed by the caller's own names, so a resource is addressed by what it is rather than by
    # where it happened to sit in the tuple.
    @test keys(avail) == (:sunlight, :water)
    @test avail.sunlight ≈ usable * solar * cell
    @test avail.water ≈ usable * rain * cell

    # …and it is *concretely* typed. A `Vector` could only have held these at
    # `Quantity{Float64}` — solar `kJ/day` and water `L/day` have no common concrete type — so every
    # element would have been boxed and the element type would have said nothing.
    @test avail isa NamedTuple
    @test isconcretetype(typeof(avail))
    @test unit(avail.sunlight) != unit(avail.water)

    # A collection built from a plain `Tuple` is named by its members' **axes** where those are
    # distinguishable — so these two supplies are `(:SolarRadiation, :Precipitation)`. Positional
    # names are never invented: two members on one axis is an error asking for names, not a silent
    # fallback. `totalsupply`'s keys follow `supply_names`, which is what the second assertion pins.
    posn = GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                       supply = (UniformSpec(solar,
                                             axis = SolarRadiation),
                                 UniformSpec(rain,
                                             axis = Precipitation)),
                       area = area)
    @test keys(EcoSISTEM.totalsupply(posn)) ==
          (:SolarRadiation, :Precipitation)
    @test keys(EcoSISTEM.totalsupply(posn)) == keys(posn.supply)
end

@testset "GridHabitat" begin
    extent = (5km, 5km)     # with cellsize 1km -> a 5×5 grid
    cellsize = 1km
    flat = UniformSpec(298.0K, axis = Temperature)
    # A synthetic study area needs no layers at all: its geometry is entirely `extent` + `cellsize`.
    square = _area(extent = extent, cellsize = cellsize)

    @testset "cell-count derivation" begin
        u = GridHabitat(regime = flat, supply = SUP, area = square)
        @test size(u.regime.matrix) == (5, 5)
        # extent not a whole number of cells warns — while the *area* is being decided, which is
        # where the geometry now lives
        @test_logs (:warn,) match_mode=:any _area(extent = (5km, 5km),
                                                  cellsize = 2km)

        # A square extent can't distinguish x from y (both dimensions have the
        # same count), which is exactly why a historical X/Y mixup between this
        # synthetic path and the data-driven GridHabitat path went undetected.
        # `extent` is (y, x), so a non-square extent must produce a matrix of
        # the same shape -- dimension 1 = y -- matching the convention used
        # throughout the package, not (nx, ny).
        nonsquare = GridHabitat(regime = flat, supply = SUP,
                                area = _area(extent = (4km, 12km),
                                             cellsize = cellsize))
        @test size(nonsquare.regime.matrix) == (4, 12)
        @test size(nonsquare.active) == (4, 12)
    end

    @testset "regime spec type selects the kind of environment" begin
        # a uniform temperature spec -> flat continuous temperature regime
        u = GridHabitat(regime = flat, supply = SUP, area = square)
        @test u isa GridHabitat
        @test EcoSISTEM.iscontinuous(u.regime)
        @test u.supply isa Supply{SolarRadiation}

        # a temperature gradient/peaked spec -> continuous regime
        g = GridHabitat(regime = GradientSpec(274.0K, 303.0K,
                                              axis = Temperature),
                        supply = SUP, area = square)
        @test EcoSISTEM.iscontinuous(g.regime)
        gp = GridHabitat(regime = PeakedSpec(274.0K, 303.0K,
                                             axis = Temperature),
                         supply = SUP, area = square)
        @test EcoSISTEM.iscontinuous(gp.regime)

        # a NicheSpec -> categorical niches. On a **non-square** grid deliberately: the niche
        # field's cluster finder indexed its neighbours the other way round, which only a grid whose
        # dimensions differ can see (see `test_Layer.jl`'s sweep).
        niche = GridHabitat(regime = NicheSpec(3,
                                               axis = EcoSISTEM.TypologyAxis),
                            supply = SUP,
                            area = _area(extent = (4km, 12km),
                                         cellsize = cellsize))
        @test !EcoSISTEM.iscontinuous(niche.regime)
        @test size(niche.regime.matrix) == (4, 12)

        # a precipitation supply spec -> a `Supply{Precipitation}` (no rainfall→water magic; the
        # water supply happens only because the supply spec declared the `Precipitation` axis)
        rain = GridHabitat(regime = GradientSpec(0.0mm / day,
                                                 100.0mm / day,
                                                 axis = Precipitation),
                           supply = UniformSpec(1.0mm / day,
                                                axis = Precipitation),
                           area = square)
        @test rain.supply isa Supply{Precipitation}
    end

    @testset "supply unit picks type; `within` honoured" begin
        water = GridHabitat(regime = flat,
                            supply = UniformSpec(1.0mm / day,
                                                 axis = Precipitation),
                            area = square)
        @test water.supply isa Supply{Precipitation}

        # Carbon is a third resource family, not a flag: a mass-per-area-per-time supply spec
        # becomes a `Supply{CarbonFlux}` because its spec *declares* `CarbonFlux`, by the same
        # axis lookup that picks solar and water, with no carbon-specific code on the build path.
        carbon = GridHabitat(regime = flat,
                             supply = UniformSpec(1.0g / (km^2 * day),
                                                  axis = CarbonFlux),
                             area = square)
        @test carbon.supply isa Supply{CarbonFlux}
        # 1 g km⁻² day⁻¹ over a 1 km² cell is 1 g/day, absolutely
        @test all(EcoSISTEM._getsupply(carbon.supply) .≈ 1.0g / day)

        active = fill(true, (5, 5))
        active[1, 1] = false
        env = GridHabitat(regime = flat, supply = SUP,
                          area = _area(extent = extent,
                                       cellsize = cellsize,
                                       within = active))
        @test env.active == active
    end

    @testset "a multi-variable synthetic environment" begin
        # Every combination of one/several regimes and one/several supplies on a synthetic area.
        # The 2-regime cases were refused outright until 2026-08-05 — `_shapesgrid`'s `::Any`
        # fallback reads "data-backed", so handing it a *tuple* of synthetic specs classified the
        # whole regime as data-backed and produced an error blaming the caller for passing data they
        # had not passed. Both halves are pinned here: the count, and that the names survive.
        rain = UniformSpec(2.0mm / day, axis = Precipitation)
        water = UniformSpec(2.0mm / day, axis = Precipitation)

        one = GridHabitat(regime = flat, supply = SUP, area = square)
        @test one.regime isa ContinuousLayer
        @test one.supply isa Supply{SolarRadiation}

        multisupply = GridHabitat(regime = flat, supply = (SUP, water),
                                  area = square)
        @test length(values(multisupply.supply)) == 2

        multiregime = GridHabitat(regime = (flat, rain), supply = SUP,
                                  area = square)
        @test multiregime.regime isa LayerCollection
        @test length(values(multiregime.regime)) == 2
        @test size(first(values(multiregime.regime)).matrix) == (5, 5)

        both = GridHabitat(regime = (flat, rain), supply = (SUP, water),
                           area = square)
        @test length(values(both.regime)) == 2
        @test length(values(both.supply)) == 2

        # A `NamedTuple` keeps the caller's names on both sides, which is what makes the members
        # addressable rather than positional.
        named = GridHabitat(regime = (temperature = flat, rain = rain),
                            supply = (sunlight = SUP, water = water),
                            area = square)
        @test keys(named.regime) == (:temperature, :rain)
        @test keys(named.supply) == (:sunlight, :water)

        # Neither guard was weakened to let the above through: a 1-element tuple is still refused
        # (it says "several layers" and does not describe several), and a data-backed spec still
        # cannot be built on an area with no CRS to place it by.
        @test_throws ErrorException GridHabitat(regime = (flat,),
                                                supply = SUP,
                                                area = square)
        @test_throws ErrorException GridHabitat(regime = SourceSpec(WorldClim{BioClim},
                                                                    :bio1),
                                                supply = SUP,
                                                area = square)
    end

    @testset "omitting a required input errors" begin
        # regime is required
        @test_throws ErrorException GridHabitat(supply = SUP,
                                                area = square)
        # supply is required
        @test_throws ErrorException GridHabitat(regime = flat,
                                                area = square)
        # so is the area: the grid is decided before anything is built on it, never inside the build
        @test_throws ErrorException GridHabitat(regime = flat,
                                                supply = SUP)
        # a study area with no layers needs both extent and cellsize to have any geometry at all
        @test_throws ErrorException _area(extent = extent)
        @test_throws ErrorException _area(cellsize = cellsize)
    end
end
# **The test that ties step 16 (per-cell cell area) to step 8 (the space resource), and the only
# supply for which the two are the *same number*.** A `SurfaceSpec()` is areal `1.0` — a whole cell
# per cell — so its per-cell supply **is** the cell's area. On a geographic grid that must therefore
# fall towards the pole, and by exactly the latitude factor.
#
# **Northern hemisphere: decreasing with latitude.** Cells shrink as you go north, so the supply
# does too — worth stating, because "varies with latitude" is easy to assert in the wrong direction
# and this is the assertion that would catch a sign error in `_areafactor`.
@testset "a space supply falls with latitude on a geographic grid" begin
    # Scotland, 54.6°N–60.9°N, south row first.
    lats = collect(range(54.6, 60.9, length = 8)) .* °
    longs = collect(range(-6.0, -2.0, length = 5)) .* °
    temp = _testraster(WorldClim{BioClim}, fill(290.0K, 8, 5), lat = lats,
                       long = longs)
    env = _env(_reg(temp, axis = Temperature), SurfaceSpec())

    supply = Array(env.supply.matrix)
    profile = supply[:, 1]                     # one column: the latitude profile
    @test length(profile) == 8

    # Strictly decreasing northwards — not merely "different", which a constant would also satisfy
    # if the correction were dropped, and not merely "unequal at the ends", which a sign error passes.
    @test issorted(profile, rev = true)
    @test all(>(zero(eltype(profile))), profile)
    # …and constant across longitude, since the factor depends on latitude alone.
    for r in 1:8
        @test all(≈(supply[r, 1]), supply[r, :])
    end

    # The magnitude is the latitude factor itself, not merely "some variation": over this extent
    # `f` runs 1.087 → 0.910, so the southern cell is ~19% larger than the northern one.
    # **These two numbers moved half a cell north on 2026-08-13**, when `_testraster` was corrected
    # to the `Intervals(Start)` its readers actually produce: the first cell now spans 54.6°–55.5°
    # rather than 54.15°–55.05°, and the last 60.9°–61.8°. Both recomputed from the fixture's own
    # `_cellintervals`, not read off the failure.
    @test profile[1] / profile[end] ≈ 1.0871 / 0.9099 rtol=1.0e-3

    # And a **projected** grid of the same shape is flat — the correction is geography, not an
    # artefact of the space axis.
    bng = _bngraster(WorldClim{BioClim}, fill(290.0K, 8, 5),
                     north = (640000.0:2500.0:657500.0) .* m,
                     east = (245000.0:2500.0:255000.0) .* m)
    flat = Array(_env(_reg(bng, axis = Temperature), SurfaceSpec()).supply.matrix)
    @test all(≈(flat[1, 1]), flat)
end
@testset "GridHabitat from data" begin
    temp = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))

    @testset "continuous regime + supply from spec/raster" begin
        env = _env(_reg(temp), SUP)
        @test env isa GridHabitat
        @test EcoSISTEM.iscontinuous(env.regime)
        @test size(env.regime.matrix) == (5, 5)
        @test env.supply isa Supply{SolarRadiation}

        # A hand-built raster carries no layer code, so it has **no axis**, and a supply's meaning
        # must never be inferred from `mm/day`'s dimension — two different quantities can share a
        # unit. Declaring it is required, and is what `ConstructedSpec`'s `axis =` keyword is for.
        water = _testraster(WorldClim{BioClim}, fill(50.0mm / day, 5, 5))
        @test _env(_reg(temp),
                   ConstructedSpec(() -> water, axis = Precipitation)).supply isa
              Supply{Precipitation}
        # …and offering it without one is refused rather than guessed. A `TypeError` rather than
        # an `ErrorException`: `supply` is **typed in the signature**, so the refusal comes from the
        # keyword itself. That is the accepted price of the signature stating what it takes — Julia
        # does not dispatch on keyword types, so a friendlier message cannot also be offered.
        @test_throws TypeError _env(_reg(temp), water)
    end

    @testset "land cover regime is continuous (a class's own % cover, not a categorical code)" begin
        landcover = _testraster(EarthEnv{LandCover},
                                Float64.(repeat(1:5, 1, 5)))
        env = _env(_reg(landcover), SUP)
        @test EcoSISTEM.iscontinuous(env.regime)
    end

    @testset "grid: native by default, an explicit cell size re-grids and is classified" begin
        @test size(_env(_reg(temp), SUP).regime.matrix) == (5, 5)  # native

        # A change of resolution is not a blanket build-time warning: the study area classifies what
        # the target grid costs each layer, so a whole multiple of the native step is reported as an
        # *exact* aggregation, which costs nothing and deserves no warning …
        bng = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)))
        agg = investigate_study_area(regime = bng, cellsize = 5.0km)
        @test only(agg.layers).kind isa EcoSISTEM.LayerAggregated
        @test only(agg.layers).kind.factor == 2
        @test isempty(agg.problems)
        # … while a genuinely lossy one says so, and names why.
        up = investigate_study_area(regime = bng, cellsize = 1.0km)
        @test only(up.layers).kind isa EcoSISTEM.LayerResampled
        @test :upsampling in [p.code for p in up.problems]
        # At `:normal` the area announces both, so the information the old warning carried is still
        # emitted — just per layer, and only when it is actually true.
        @test_logs (:info,) match_mode=:any StudyArea(regime = bng,
                                                      cellsize = 5.0km)
        # 22.5 km of ground at 10 km cells is two whole cells and a 2.5 km remainder that fills no
        # cell. Was 3 × 3, under a labelling where the grid straddled the data's edge by half a
        # cell at the origin and the outer cells were only 39% covered — and were nonetheless given
        # a whole cell's worth of supply. `simulate_safely = false` restores that count.
        @test size(GridHabitat(regime = bng, supply = SUP,
                               area = _area(regime = bng,
                                            cellsize = 10.0km)).regime.matrix) ==
              (2, 2)
    end

    @testset "mask rules: hasdata / landcoverclass" begin
        # `hasdata` — the public combine rule masking a source's non-missing (non-NaN) cells.
        withnan = fill(290.0K, 5, 5)
        withnan[1, 1] = NaN * K
        dm = hasdata(_testraster(WorldClim{BioClim}, withnan))
        # A **raster**, like every other combine result — a mask is one whose element type is
        # `Bool`, not a different container. No public function hands out a bare array, which is
        # what `test_EcoSISTEM.jl` asserts with no carve-outs.
        @test dm isa ClimateRaster
        @test eltype(dm) == Bool
        @test !dm.array[1, 1] && dm.array[2, 2]

        # one cell of every land-cover class, addressed throughout by name (never a raw numeric
        # code) so this test can't repeat the exact bug it's guarding against — a hardcoded
        # `water = [4]` silently meant "other_trees", not water. `landcoverclass` looks each up.
        classes = (:needleleaf_trees, :evergreen_broadleaf_trees,
                   :deciduous_broadleaf_trees, :other_trees, :shrubs,
                   :herbaceous,
                   :cultivated_and_managed, :regularly_flooded, :urban_builtup,
                   :snow_ice, :barren, :open_water)
        codes = collect(landcoverclass.(classes))
        landcover = _testraster(EarthEnv{LandCover},
                                Float64.(reshape(codes, 3, 4)),
                                lat = (0.0:1.0:2.0) .* °,
                                long = (0.0:1.0:3.0) .* °)
        classgrid = Array(landcover.array)
        cellof(name) = findfirst(==(landcoverclass(name)), classgrid)

        # The land-mask predicate (exclude open water) written as the docs recommend it, by name.
        lm = classgrid .!= landcoverclass(:open_water)
        @test count(lm) == 11                # only open_water excluded
        @test !lm[cellof(:open_water)]
        @test _env(_reg(landcover), SUP, within = lm) isa GridHabitat

        # The nature-mask predicate: exclude built/bare/frozen/water surfaces too.
        natureexcluded = landcoverclass.((:open_water, :urban_builtup, :barren,
                                          :snow_ice))
        nm = [c ∉ natureexcluded for c in classgrid]
        @test count(nm) == 8
        @test !nm[cellof(:open_water)] && !nm[cellof(:urban_builtup)] &&
              !nm[cellof(:barren)] && !nm[cellof(:snow_ice)]

        # farmland = false additionally excludes cultivated_and_managed.
        nofarmexcluded = (natureexcluded...,
                          landcoverclass(:cultivated_and_managed))
        nm_nofarm = [c ∉ nofarmexcluded for c in classgrid]
        @test count(nm_nofarm) == 7
        @test !nm_nofarm[cellof(:cultivated_and_managed)]
    end

    @testset "circle mask" begin
        env = _env(_reg(temp), SUP, within = CircleMaskSpec(radius = 100.0km))
        @test env.active isa DimensionalData.AbstractDimArray{Bool}

        # explicit `centre` (a LatLong point, degrees — not the same unit as `radius`) selects a
        # different cell than the default grid-centre
        off_centre = _env(_reg(temp), SUP,
                          within = CircleMaskSpec(radius = 100.0km,
                                                  centre = LatLong(0.0°, 0.0°)))
        @test off_centre.active != env.active
        @test off_centre.active[1, 1]        # the corner cell, now inside the circle
    end

    @testset "lazy land-cover mask via ConstructedSpec" begin
        # A land-cover mask is a `ConstructedSpec` combine over the real EarthEnv{LandCover} data:
        # construction is pure — no read/download happens until `GridHabitat` materialises it.
        # Neither combine names an array type: a raster broadcasts and yields a raster, so the
        # mask is `Bool`-valued but still a raster, exactly like a layer combine's result.
        # Every EarthEnv read here is cut to Europe, and the reason is a hard ceiling rather than
        # tidiness. A `within` mask cannot be windowed by the study area: it is what DECIDES the
        # grid, so with no bound of its own the whole world is fetched, compressed and compared --
        # twice here -- and then a habitat is built on every land cell on earth.
        #
        # Measured peak resident for this file, against a GitHub runner's 16 GB: 2.1 GB without this
        # testset at all, and 14.6 GB with it unbounded. So this testset was the whole of the file's
        # cost, and it is what shut the Linux runner down and what left macOS with so little free
        # memory that Rasters refused a later 60 MB read in this same file.
        #
        # Europe is chosen for what it contains rather than for its size: the assertions below need
        # open water, and they need the four classes the nature mask excludes on top of water to be
        # present somewhere, which coasts, cities, Alpine barren and Scandinavian snow supply. A
        # region without them would make the second test pass vacuously.
        #
        # `cut` is the right lever and two others are not, both measured rather than assumed.
        # `cellsize` on the study area does nothing here (12.5 GB, and an error): the mask is
        # materialised at native resolution to decide the extent, so a coarser target grid arrives
        # too late. `extent` is refused outright -- it is a size rather than a bounding box, and
        # naming one beside data layers that carry their own extent is an error by construction.
        europe = EcoSISTEM.boundingbox("Europe", islands = true)
        landmask = ConstructedSpec(
                                   SourceSpec(EarthEnv{LandCover},
                                              cut = europe),
                                   axis = EcoSISTEM.NicheAxis) do lc
            compress_landcover(lc) .!= landcoverclass(:open_water)
        end
        naturemask = ConstructedSpec(
                                     SourceSpec(EarthEnv{LandCover},
                                                cut = europe),
                                     axis = EcoSISTEM.NicheAxis) do lc
            excluded = landcoverclass.((:open_water, :urban_builtup, :barren,
                                        :snow_ice, :cultivated_and_managed))
            compress_landcover(lc) .∉ Ref(excluded)
        end
        @test landmask isa EcoSISTEM.ConstructedSpec
        @test naturemask isa EcoSISTEM.ConstructedSpec
        @test hasdata isa Function

        # Materialisation — real (small single-class) EarthEnv{LandCover} network read, the
        # regression oracle for the whole chain: spec -> `read(EarthEnv{LandCover})` ->
        # `compress_landcover` -> combine rule -> `_samplemask`.
        cultivated = SourceSpec(EarthEnv{LandCover}, :cultivated_and_managed,
                                cut = europe)
        lenv = _env(cultivated, SUP, within = landmask)
        @test lenv isa GridHabitat
        @test 0 < count(lenv.active) < length(lenv.active)

        nenv = _env(cultivated, SUP, within = naturemask)
        @test nenv isa GridHabitat
        # the nature mask (also excluding cultivated) excludes strictly more than the land mask
        @test count(nenv.active) < count(lenv.active)
    end

    @testset "shapemask" begin
        # No CRS metadata (no .prj) — treated as already WGS84 lat/long. Polygon covers lat/long
        # [-0.5,2.5] x [-0.5,2.5]. **A cell is tested at its representative *midpoint***, and under
        # the `Intervals(Start)` locus a reader produces, `temp`'s cells labelled 0…4 have midpoints
        # 0.5, 1.5, 2.5, 3.5, 4.5. So cells 0 and 1 are inside and cell 2's midpoint falls **exactly
        # on the polygon edge**, which GDAL's `contains` excludes — the same rule the projected
        # shapefile case below already documents.
        path = _testshapefile(-0.5, -0.5, 2.5, 2.5)
        env = _env(_reg(temp), SUP, within = ShapeSpec(path))
        @test env.active isa DimensionalData.AbstractDimArray{Bool}
        # The mask sets the extent, so the grid *is* those cells rather than the data's full 5 × 5
        # with an active corner — and because a native-resolution mask crops the data's own grid
        # instead of re-gridding, the survivors keep the source's own labels exactly.
        @test size(env.active) == (2, 2)
        @test all(env.active)
        @test parent(DimensionalData.lookup(env.active, Y)) == [0.0, 1.0] .* °
        @test parent(DimensionalData.lookup(env.active, X)) == [0.0, 1.0] .* °

        # British National Grid (EPSG:27700) — reprojects to ~lat 55.72-55.81, long -4.39 to -4.23.
        bng = ArchGDAL.importEPSG(27700)
        bngpath = _testshapefile(250000.0, 650000.0, 260000.0, 660000.0,
                                 sr = bng)
        # **Relabelled, not repositioned** (2026-08-13). This fixture's point is a *spatial
        # coincidence* — a ~10 km polygon landing in exactly one of two 83 km cells — so when the
        # locus moved, the labels moved with it to keep the cells over the **same ground**: what
        # `Center` called 55.75°/56.5° and −4.6°/−4.3° is what `Start` calls 55.375°/56.125° and
        # −4.75°/−4.45°. The cells, and so the assertion below, are unchanged. Elsewhere in this
        # file the opposite call was right — a fixture describing ground generically simply covers
        # different ground now — but here preserving the coincidence *is* the fixture's meaning.
        scotland = _testraster(WorldClim{BioClim}, fill(290.0K, 2, 2),
                               lat = [55.375, 56.125] .* °,
                               long = [-4.75, -4.45] .* °)
        scotenv = _env(_reg(scotland), SUP, within = ShapeSpec(bngpath))
        @test scotenv.active == [false true; false false]

        # `ShapeSpec` also accepts a URL, downloaded (via `EcoSISTEM.CachedAsset`)
        # only when materialised — a `file://` URL needs no network access, so this isn't flaky.
        # A single-file format (GeoJSON, unlike a bare `.shp`'s sidecar files) so one `CachedAsset`
        # download is genuinely everything `_shapegeoms` needs to read it.
        geojsonpath = joinpath(mktempdir(), "test.geojson")
        poly = ArchGDAL.createpolygon([[
                                          (-0.5, -0.5),
                                          (2.5, -0.5),
                                          (2.5, 2.5),
                                          (-0.5, 2.5),
                                          (-0.5, -0.5)
                                      ]])
        ArchGDAL.create(geojsonpath,
                        driver = ArchGDAL.getdriver("GeoJSON")) do dataset
            makelayer = layer -> begin
                ArchGDAL.addfeature(layer) do feature
                    return ArchGDAL.setgeom!(feature, poly)
                end
                ArchGDAL.copy(layer, dataset = dataset)
            end
            ArchGDAL.createlayer(makelayer, name = "test",
                                 geom = ArchGDAL.wkbPolygon)
        end
        urlspec = ShapeSpec("file://" * geojsonpath)
        @test urlspec.path isa EcoSISTEM.CachedAsset
        urlenv = _env(_reg(temp), SUP, within = urlspec)
        # same polygon, same mask-set extent as the plain-path case above
        @test size(urlenv.active) == (2, 2)
        @test all(urlenv.active)
        # tidy up the shared cache: this test's fixture has no business persisting there
        rm(EcoSISTEM.assetpath(urlspec.path))
    end

    @testset "geodesy: area-preserving cell size + top/bottom report" begin
        spanning = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5),
                               lat = (0.0:10.0:40.0) .* °,
                               long = (0.0:10.0:40.0) .* °)
        @test_logs (:info,) match_mode=:any _env(_reg(spanning), SUP)
    end

    @testset "projected grids: target CRS staged rule + crs= keyword" begin
        bng = _bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9))

        # A projected reference keeps its own metric grid, and `cellsize` there means a real physical
        # cell side — exactly, not an approximation.
        penv = _env(_reg(bng, axis = Temperature), SUP, cellsize = 5.0km)
        @test EcoSISTEM._uniformcellside(EcoSISTEM.getregime(build_ecosystem(build_species(2,
                                                                                           tolerance = (291.0K,
                                                                                                        2.0K),
                                                                                           toleranceaxis = Temperature,
                                                                                           demand = DEM,
                                                                                           demandaxis = SolarRadiation,
                                                                                           seed = 1),
                                                                             penv,
                                                                             seed = 1))) ≈
              5.0km
        # 22.5 km of ground (9 cells of 2.5 km, edge to edge) at 5 km cells → 4 whole cells each
        # way, starting on the data's own lower edge, with 2.5 km left over that fills no cell.
        # Was 5 × 5 while the grid was labelled by cell *centres* and laid one extra cell: the
        # first cell then straddled the data's edge and the last was half outside it.
        @test size(penv.regime.matrix) == (4, 4)

        # A single projected CRS among otherwise-geographic inputs is adopted for the whole grid, so
        # combining BNG land cover with WGS84 climate keeps the metric grid rather than degrees —
        # even though the geographic layer is the *first* (reference) one. The WGS84 layer must
        # actually cover the same ground (central Scotland) for the layers to be combinable at all.
        scotgeo = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5),
                              lat = (55.6:0.1:56.0) .* °,
                              long = (-4.5:0.1:-4.1) .* °)
        pair = (regime1 = _reg(scotgeo), regime2 = _reg(bng))
        mixed = GridHabitat(regime = pair, supply = SUP,
                            area = _area(regime = pair))
        @test map(EcoSISTEM.iscontinuous, values(mixed.regime)) == (true, true)
        # The adopted grid is the projected one: its cell side is a real length, not a degree step.
        @test EcoSISTEM._crsunit(EcoSISTEM._targetcrs((scotgeo, bng), nothing,
                                                      nothing)) == u"m"

        # A physical `cellsize` over a purely geographic grid fails closed rather than approximating
        # degrees as kilometres — and the error names a concrete UTM `crs` to paste in.
        err = try
            _area(regime = _reg(temp), cellsize = 5.0km)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("crs = Rasters.EPSG(", sprint(showerror, err))

        # An explicit `crs` reprojects the target: a WGS84 source built onto a UTM grid is metric, so
        # `cellsize` is honoured there too.
        uenv = _env(_reg(temp, axis = Temperature), SUP, cellsize = 200.0km,
                    crs = Rasters.EPSG(32631))
        @test uenv isa GridHabitat
        @test EcoSISTEM._uniformcellside(EcoSISTEM.getregime(build_ecosystem(build_species(2,
                                                                                           tolerance = (290.0K,
                                                                                                        2.0K),
                                                                                           toleranceaxis = Temperature,
                                                                                           demand = DEM,
                                                                                           demandaxis = SolarRadiation,
                                                                                           seed = 1),
                                                                             uenv,
                                                                             seed = 1))) ≈
              200.0km
    end

    @testset "CRS suggestions come from PROJ's database, not arithmetic" begin
        # The trap this locks in: "Mercator" appears in the name of three *different* projections.
        # Only the cylindrical variants distort area unboundedly; Transverse and Oblique Mercator are
        # near-area-true in their own zone, and excluding them would discard British National Grid
        # and every UTM zone — i.e. almost every good local answer.
        @test EcoSISTEM._areadistorting("Mercator (variant A)")
        @test EcoSISTEM._areadistorting("Popular Visualisation Pseudo Mercator")
        @test !EcoSISTEM._areadistorting("Transverse Mercator")
        @test !EcoSISTEM._areadistorting("Hotine Oblique Mercator (variant A)")
        @test !EcoSISTEM._areadistorting("Lambert Azimuthal Equal Area")

        # The database is real and large enough to be worth consulting.
        @test length(EcoSISTEM._crscandidates()) > 1000

        # A national extent gets its own national grid rather than a generic UTM zone…
        scotland = Extent(Y = (54.63°, 60.86°), X = (-8.65°, -0.72°))
        @test EcoSISTEM._suggestcrs(scotland).code == 27700
        # …and every suggestion genuinely contains what was asked for, and is not area-distorting.
        for b in (scotland,
            Extent(Y = (-34.8°, 37.4°), X = (-17.5°, 51.4°)),   # Africa: many UTM zones
            Extent(Y = (-34.0°, -33.5°), X = (150.9°, 151.3°))) # Sydney
            s = EcoSISTEM._suggestcrs(b)
            @test !isnothing(s)
            @test s.west <= ustrip(°, b.X[1]) && s.east >= ustrip(°, b.X[2])
            @test s.south <= ustrip(°, b.Y[1]) && s.north >= ustrip(°, b.Y[2])
            @test !EcoSISTEM._areadistorting(s.method)
        end

        # The advice is paste-able and names the CRS, so it can be judged rather than trusted.
        advice = EcoSISTEM._crsadvice(scotland)
        @test occursin("crs = Rasters.EPSG(27700)", advice)
        @test occursin("British National Grid", advice)
    end

    @testset "resolution inference: aligned layer's own, else unanimous, else require cellsize" begin
        # With `cellsize` omitted the target adopts the alignment layer's own resolution — and where
        # no layer can be aligned to, the layers in the target CRS must at least agree on one.
        bng = _bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9))
        @test _env(_reg(bng), SUP) isa GridHabitat

        # Two same-CRS layers at different resolutions: the finer is aligned to, so the coarser is
        # the one that gets resampled — and the area says which, and why.
        fine = _bngraster(WorldClim{BioClim}, fill(292.0K, 17, 17),
                          east = (245000.0:1250.0:265000.0) .* m,
                          north = (640000.0:1250.0:660000.0) .* m)
        pair = (regime1 = _reg(bng), regime2 = _reg(fine))
        report = investigate_study_area(regime = pair)
        @test report.cellsize == 1250.0m         # the finer layer's own step
        # A multi-variable regime is expanded into one numbered layer per element, so each is
        # planned (and can be aligned to) in its own right.
        @test [l.name for l in report.layers] == [:regime1, :regime2]
        @test report.align === :regime2
        @test GridHabitat(regime = pair, supply = SUP,
                          area = _area(regime = pair)) isa GridHabitat
        # …and an explicit `cellsize` overrides that choice.
        @test _area(regime = pair, cellsize = 5.0km).report.cellsize == 5.0km

        # A `crs` no input layer uses leaves no layer to align to, so nothing states a resolution in
        # the target's units. Rather than refuse, the cell is *measured across the projection*:
        # putting WGS84 climate on a metric grid is an entirely reasonable request.
        wgs = _reg(_testraster(WorldClim{BioClim}, fill(290.0K, 9, 9),
                               lat = (55.0:0.1:55.8) .* °,
                               long = (-4.5:0.1:-3.7) .* °))
        across = _area(regime = wgs, crs = Rasters.EPSG(27700))
        @test isnothing(across.report.align)          # nothing can be kept exactly
        # 0.1° of latitude is ~11.1 km, and of longitude at 55°N ~6.4 km, so the area-preserving
        # side is ~sqrt(11.1 × 6.4) ≈ 8.4 km — measured against the real projection, not a constant.
        # A measured size is then floored to one significant figure, because nothing aligns to a
        # synthesised grid anyway and 8438.6 m is not a cell size anyone would choose to read.
        @test across.report.cellsizesource isa EcoSISTEM.RoundedFromMeasurement
        @test across.report.cellsize == 8.0km
        @test unit(across.report.cellsize) === u"km"
        # Floored, never nearest: the grid must not end up coarser than the projection measured,
        # or resolution the source carried is thrown away.
        @test across.report.cellsize ≤ 8.5km
        # Both sizes are reported, so the rounding is visible rather than silent.
        rounding = only(filter(p -> p.code === :cellsize_rounded,
                               across.report.problems))
        @test rounding.severity isa EcoSISTEM.ProblemNotice
        @test occursin(r"was \d+\.\d+ m", rounding.message)   # the measured size…
        @test occursin("8 km", rounding.message)              # …and the one adopted

        # …and the derived grid genuinely builds.
        @test GridHabitat(regime = wgs, supply = SUP, area = across) isa
              GridHabitat

        # Only when nothing can state a resolution at all does it still refuse.
        @test isnothing(EcoSISTEM._stepacross(_testraster(WorldClim{BioClim},
                                                          fill(290.0K, 1, 1),
                                                          lat = [55.0] .* °,
                                                          long = [-4.0] .* °),
                                              Rasters.EPSG(27700)))
    end

    @testset "a NicheSpec builds on a POSITIONED area, and both paths agree" begin
        # Inspection and building must generate a `NicheSpec`'s field the same way: they are two
        # routes to one answer and have drifted apart before. Where `materialise` reaches
        # `_specfield` (which knows `NicheSpec`) and the builder reaches `_syntheticsupplyfield`
        # (which does not), a spec `materialise` shows happily cannot be built, and the failure is a
        # bare `MethodError: AbstractFloat(::NicheSpec)`. The other NicheSpec test
        # above uses a *synthetic* area, which is why nothing caught it — a positioned area is the
        # only thing that reaches the builder's branch.
        #
        # **Every area here is deliberately NON-SQUARE**, and they were all square until a second
        # bug was found underneath this one: the niche field's cluster finder called
        # `get_neighbours(M, y, x)` where the first argument is the *row*, so it read the transposed
        # neighbourhood on a square grid and threw *"Coordinates outside grid"* on any other. A
        # square fixture cannot see either half of that. (`test_Layer.jl` sweeps the shapes; these
        # are the same property reached through the two public paths.)
        geo = _reg(_testraster(WorldClim{BioClim}, fill(290.0K, 7, 5),
                               lat = (0.0:1.0:6.0) .* °,
                               long = (0.0:1.0:4.0) .* °))
        bng = _reg(_bngraster(WorldClim{BioClim}, fill(290.0K, 9, 7),
                              north = (640000.0:2500.0:660000.0) .* m,
                              east = (245000.0:2500.0:260000.0) .* m))
        for area in (_area(regime = geo),                            # geographic
            _area(regime = bng),                            # projected
            _area(extent = (40km, 70km), cellsize = 10km))  # synthetic
            spec = NicheSpec(3, axis = EcoSISTEM.TypologyAxis)
            # …and the fixture really is non-square, so this cannot rot back into a square one.
            @test size(area.report.active, 1) != size(area.report.active, 2)
            # Both paths give a *categorical* layer — the thing most at risk here, because the
            # builder's raster route asks `iscategorical`, which answers `false` for a code-less
            # synthetic raster and would silently have built a continuous regime.
            shown = materialise(spec, area)
            @test shown isa CategoricalLayer
            built = GridHabitat(regime = spec, supply = SUP, area = area)
            @test built.regime isa CategoricalLayer
            @test !EcoSISTEM.iscontinuous(built.regime)
            # Same grid and the same cell size on both —  an **angle** on the geographic area,
            # which is the second bug this fixes: the cell size was threaded into `_randomniches`
            # only to be discarded, and its `::Unitful.Length` annotation refused a degree grid.
            @test size(shown.matrix) == size(built.regime.matrix)
            @test shown.size == built.regime.size
            @test EcoSISTEM._yxcompatible(EcoSISTEM._yx(built.regime),
                                          EcoSISTEM._yx(shown))
            # Class codes, in range, and never canonicalised as if they were measurements.
            @test eltype(built.regime.matrix) <: Integer
            @test issubset(unique(built.regime.matrix), 1:3)
        end
        # **Not** asserted: that the two show the *same* niches. `_randomniches` is stochastic
        # and unseeded, so each call draws its own pattern — a known limit of inspecting a
        # `NicheSpec`, recorded in the plan rather than fixed here.
    end

    @testset "geographic grids build (with a warning) but cannot be simulated" begin
        # A degree grid is fine to *build* — useful for inspecting/plotting data in its own
        # coordinates — and the study area warns, naming a projected CRS to use instead.
        report = investigate_study_area(regime = _reg(temp))
        geo = only(filter(p -> p.code === :geographic, report.problems))
        @test occursin("crs = Rasters.EPSG(", geo.message)
        @test occursin("cannot be simulated", geo.message)
        # At `:normal` that problem really is emitted, rather than only collected.
        @test_logs (:warn, r"cannot be simulated") match_mode=:any StudyArea(regime = _reg(temp))

        genv = _env(_reg(temp), SUP)
        @test genv isa GridHabitat

        # …but assembling an ecosystem on it is refused, because dispersal assumes uniform cells
        sp = build_species(2, tolerance = (290.0K, 2.0K),
                           toleranceaxis = Temperature, demand = DEM,
                           demandaxis = SolarRadiation,
                           seed = 1)
        err = try
            build_ecosystem(sp, genv, seed = 1)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("cannot be simulated", sprint(showerror, err))

        # a projected environment assembles fine
        penv = _env(_reg(_bngraster(WorldClim{BioClim}, fill(290.0K, 9, 9)),
                         axis = Temperature),
                    SUP)
        @test build_ecosystem(sp, penv, seed = 1) isa Ecosystem
        # and so does a *synthetic* one, which has no CRS at all but a genuine metric cell size
        senv = GridHabitat(regime = UniformSpec(290.0K,
                                                axis = Temperature),
                           supply = SUP,
                           area = _area(extent = (5km, 5km),
                                        cellsize = 1km))
        @test build_ecosystem(sp, senv, seed = 1) isa Ecosystem
        # the guard is on the Ecosystem constructor too, not only on `build_ecosystem`
        @test_throws ErrorException Ecosystem(sp, genv,
                                              NicheSuitability{Temperature,
                                                               typeof(1.0K)}())
    end

    @testset "masks on a projected target" begin
        # `CircleMaskSpec`/`ShapeSpec` both assumed degree coordinates before A2; on a projected
        # target they must work in the target's own metric plane instead.
        bng = _bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9))

        # Default centre = the grid's own midpoint (CRS-agnostic): a 6 km disc over a 20 × 20 km grid
        # of 2.5 km cells selects a real subset, not everything and not nothing.
        cenv = _env(_reg(bng), SUP, within = CircleMaskSpec(radius = 6.0km))
        @test 0 < count(cenv.active) < length(cenv.active)
        # Symmetric about the grid centre, as a disc on a uniform metric grid must be. Compared as
        # plain arrays: `reverse` on a `DimArray` reverses its coordinate lookup too, so comparing the
        # `DimArray`s directly would compare (reversed) dims rather than just the mask values.
        @test Array(cenv.active) == reverse(Array(cenv.active), dims = 1)
        @test Array(cenv.active) == reverse(Array(cenv.active), dims = 2)

        # An explicit `centre` is a WGS84 `LatLong` whatever the target CRS, so it is transformed
        # into the target's coordinates — giving a *different* (off-centre) mask.
        offc = _env(_reg(bng), SUP,
                    within = CircleMaskSpec(radius = 6.0km,
                                            centre = LatLong(55.80°, -4.30°)))
        @test offc.active != cenv.active

        # A BNG-CRS shapefile against a BNG target: `_shapegeoms` reprojects features into the
        # target CRS (rather than always WGS84), so the point-in-polygon test happens in metres.
        # The polygon covers the lower-left quadrant. Under the `Intervals(Start)` locus the cells
        # labelled 640000…660000 have **midpoints** 641250…661250, so the four inside
        # [640000, 650000] are indices 1:4 — and none of them lands on the polygon edge, unlike the
        # WGS84 case above.
        bngsr = ArchGDAL.importEPSG(27700)
        path = _testshapefile(245000.0, 640000.0, 255000.0, 650000.0,
                              sr = bngsr)
        senv = _env(_reg(bng), SUP, within = ShapeSpec(path))
        # As in "shapemask", the mask now sets the extent: the grid is exactly those cells, keeping
        # the source grid's own 2.5 km registration — and the labels are the cells' lower corners.
        @test size(senv.active) == (4, 4)
        @test all(senv.active)
        @test parent(DimensionalData.lookup(senv.active, Y)) ==
              [640000.0, 642500.0, 645000.0, 647500.0] .* m
        @test parent(DimensionalData.lookup(senv.active, X)) ==
              [245000.0, 247500.0, 250000.0, 252500.0] .* m
    end

    @testset "a supply's own gaps restrict the grid, and are zeroed only in the habitat" begin
        # A data-driven supply's `NaN`s must reach the coverage check rather than being zeroed at
        # construction: zeroing makes a missing supply cell indistinguishable from a genuinely zero
        # one, and a supply-shaped hole then restricts nothing. A supply's gaps restrict the grid
        # exactly as a regime's do.
        reg = _bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9))
        supvals = fill(600.0mm / day, 9, 9)
        supvals[1, :] .= NaN * mm / day      # a missing edge row …
        supvals[:, 1] .= NaN * mm / day      # … and a missing edge column
        holed = ConstructedSpec(() -> _bngraster(WorldClim{BioClim}, supvals),
                                axis = Precipitation)

        # Naming the supply lets it shape the grid: its gaps are then part of the coverage the area
        # recuts to, exactly as a regime's are.
        area = _area(regime = _reg(reg), supply = holed)
        env = GridHabitat(regime = _reg(reg), supply = holed, area = area)
        # The 17 missing cells (a row + a column, sharing a corner) are dropped, and the recut then
        # tightens the 9 × 9 grid to the 8 × 8 block that is left.
        @test size(env.active) == (8, 8)
        @test all(env.active)
        # …but the habitat's supply is still `NaN`-free: the gaps survive only long enough to be seen.
        @test !any(isnan, ustrip.(env.supply.matrix))

        # An *unnamed* supply cannot move the frame — it can only mark cells inactive, and says so.
        unshaped = @test_logs (:warn, r"can only remove cells") match_mode=:any GridHabitat(regime = _reg(reg),
                                                                                            supply = holed,
                                                                                            area = _area(regime = _reg(reg)))
        @test size(unshaped.active) == (9, 9)          # frame unchanged
        @test count(unshaped.active) == 81 - 17        # but the gaps are inactive

        # A supply with no gaps must not restrict anything — the 9 × 9 grid stays whole.
        whole = ConstructedSpec(() -> _bngraster(WorldClim{BioClim},
                                                 fill(600.0mm / day, 9, 9)),
                                axis = Precipitation)
        @test size(_env(_reg(reg), whole).active) == (9, 9)
    end

    @testset "supply constructors preserve NaN and leave the caller's array alone" begin
        # Both halves of the old behaviour were wrong: the zeroing destroyed the distinction between
        # "no data" and "no resource", and it did so by mutating the array it was handed.
        rate = 1.0 *
               EcoSISTEM.canonicalunit(EcoSISTEM.Resource, SolarRadiation)
        arr = fill(rate, 3, 3)
        arr[1, 1] = NaN * rate
        # Real coordinates, because a layer must be able to state its own grid — `NoLookup` dims
        # are refused at construction now (`[CELL-DO]` 4b). Incidental to what this testset checks,
        # which is NaN handling, so the fixture simply says what its cells are.
        supply = Supply{SolarRadiation}(EcoSISTEM._sizedyx(arr, 1.0km))
        @test isnan(ustrip(supply.matrix[1, 1]))
        @test isnan(ustrip(arr[1, 1]))     # the caller's own array is untouched
    end

    @testset "a mask reaching past the data is reported, not silently shrunk" begin
        bng = _bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9))
        bngsr = ArchGDAL.importEPSG(27700)

        # Twice the layer's own 20 × 20 km footprint: the grid cannot exceed the data, so the clamp is
        # warned about, and the result is the whole (fully active) data grid.
        big = _testshapefile(235000.0, 630000.0, 275000.0, 670000.0, sr = bngsr)
        report = investigate_study_area(regime = _reg(bng),
                                        within = ShapeSpec(big))
        @test :mask_clamped in [p.code for p in report.problems]
        wide = GridHabitat(regime = _reg(bng), supply = SUP,
                           area = StudyArea(report, verbosity = :silent))
        @test size(wide.active) == (9, 9)
        @test all(wide.active)
        # …and at `:normal` the area really says so rather than only recording it.
        @test_logs (:warn, r"clamped to the data") match_mode=:any StudyArea(regime = _reg(bng),
                                                                             within = ShapeSpec(big))

        # A mask covering the data exactly must *not* warn: the envelope is the layer's own
        # footprint, which is coverage, not a cut.
        # **These four numbers ARE the data's extent and moved with the locus** — under
        # `Intervals(Start)` the 9 cells of 2.5 km labelled 640000…660000 span [640000, 662500), so
        # the exact envelope is that, not the `centre ± half` box (638750…661250) a cell-centre
        # reading would give.
        exact = _testshapefile(245000.0, 640000.0, 267500.0, 662500.0,
                               sr = bngsr)
        exactreport = investigate_study_area(regime = _reg(bng),
                                             within = ShapeSpec(exact))
        @test :mask_clamped ∉ [p.code for p in exactreport.problems]
        @test size(_env(_reg(bng), SUP, within = ShapeSpec(exact)).active) ==
              (9, 9)

        # A mask that misses the data altogether is a clear error, not an empty grid.
        away = _testshapefile(500000.0, 100000.0, 510000.0, 110000.0,
                              sr = bngsr)
        @test_throws Exception _area(regime = _reg(bng),
                                     within = ShapeSpec(away))
    end

    @testset "GridHabitat on a decided StudyArea" begin
        bng = _bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9))

        # The area decides the grid; the builder only samples onto it.
        area = StudyArea(regime = _reg(bng), verbosity = :silent)
        env = GridHabitat(regime = _reg(bng), supply = SUP, area = area)
        @test env isa GridHabitat
        @test size(env.active) == size(area.report.active)
        @test size(env.regime.matrix) == size(area.report.active)

        # A synthetic area needs no layers and no CRS, and keeps the (y, x) convention.
        synth = StudyArea(extent = (4km, 12km), cellsize = 1km,
                          verbosity = :silent)
        senv = GridHabitat(regime = UniformSpec(298.0K,
                                                axis = Temperature),
                           supply = SUP, area = synth)
        @test size(senv.active) == (4, 12)

        # **`extent` with data layers must be refused, never silently ignored.** It says how *big*
        # an area is, never *where*, so once layers position the grid it cannot be honoured — and
        # dropping it quietly lets a global dataset produce a grid nobody asked for. The message
        # names `within`, which is what actually positions an area.
        @test_throws "cannot be combined with data layers" StudyArea(regime = _reg(bng),
                                                                     extent = (10km,
                                                                               10km),
                                                                     cellsize = 2.5km,
                                                                     verbosity = :silent)

        # A data layer needs a positioned area to be sampled onto — that mismatch is named rather
        # than left as a `MethodError`.
        @test_throws ErrorException GridHabitat(regime = _reg(bng),
                                                supply = SUP,
                                                area = synth)

        # **The mirror case is NOT an error.** A *synthetic* regime on a *positioned* area is
        # supported, because a generated layer needs shape and orientation rather than coordinates —
        # which is exactly why the **supply** side allows it too. Refusing it would be inconsistent
        # in a way the package itself could not survive: `build_habitat()` defaults its regime to a
        # synthetic `UniformSpec`, so the package's own "sensible defaults" entry point could not be
        # used with a real study area at all.
        syntheticreg = GridHabitat(regime = UniformSpec(298.0K,
                                                        axis = Temperature),
                                   supply = SUP, area = area)
        @test syntheticreg isa GridHabitat
        @test all(==(298.0K), syntheticreg.regime.matrix)
        # And it does not narrow coverage: a generated layer has no gaps, so `active` is still
        # whatever the area and the data-backed supply decided.
        @test count(syntheticreg.active) == count(area.report.active)
        # …the same route `DefaultEcosystem()` takes.
        @test build_habitat(supply = SUP, area = area) isa GridHabitat

        # The frame is fixed: a layer the area never saw can remove cells, but not resize the grid.
        holed = fill(291.0K, 9, 9)
        holed[1, :] .= NaN * K
        env2 = @test_logs (:warn, r"can only remove cells") match_mode=:any GridHabitat(regime = _reg(_bngraster(WorldClim{BioClim},
                                                                                                                 holed)),
                                                                                        supply = SUP,
                                                                                        area = area)
        @test size(env2.active) == size(area.report.active)   # unchanged frame
        @test count(env2.active) < count(area.report.active)  # fewer active cells
    end

    @testset "cell size dispatches on coordinate units" begin
        deg = (0.0:1.0:2.0) .* °
        km_axis = (0.0:1.0:2.0) .* km
        # projected (length) coordinates -> planar, exactly sqrt(Δlat × Δlong)
        @test EcoSISTEM._cellsize(km_axis, km_axis) ≈ 1.0km
        # geographic (degree) coordinates -> spherical, shrinking toward the poles
        @test EcoSISTEM._cellsize((58.0:1.0:60.0) .* °, deg) <
              EcoSISTEM._cellsize(deg, deg)
    end

    # `_cellsize` is one scalar *side*, but a geographic cell's **area** falls with latitude, so
    # every areal supply is multiplied by `_cellareas` instead — a column indexed by latitude.
    # **Deliberately an 11 × 7 grid.** A square grid cannot see a Y/X mix-up and a projected one
    # cannot see a missing/spurious `Ref`, so the only shape that catches either is geographic and
    # non-square.
    @testset "cell areas vary with latitude on a geographic grid" begin
        lats = collect(range(50.25°, 59.75°, length = 11))
        longs = collect(range(-5.75°, 1.75°, length = 7))
        # `_cellareas` reads an array's own dims now, asking each axis for what it means —
        # midpoints for the nominal size, intervals for the latitude correction (`[LOCUS-BLIND]`) —
        # so the fixture is a raster rather than two coordinate vectors.
        grid(la, lo) = _testraster(WorldClim{BioClim},
                                   fill(291.0K, length(la), length(lo)),
                                   lat = la, long = lo).array
        areas = EcoSISTEM._cellareas(grid(lats, longs))

        # A column, not a full matrix and not a scalar: it must broadcast across X (and across a
        # monthly stack's `Ti`) without ever materialising one value per cell.
        @test size(areas) == (11, 1)
        # Cells shrink towards the pole. `issorted(rev = true)` rather than a spot check, because a
        # reversed or mis-indexed column would still pass on its endpoints alone.
        @test issorted(vec(areas), rev = true)
        # …and a south-up grid gives the mirror image, so nothing has to be reversed by hand.
        @test vec(EcoSISTEM._cellareas(grid(reverse(lats), longs))) ≈
              reverse(vec(areas))

        # **The exact property, and the one the scalar could not have.** The per-cell areas sum to
        # the true spherical area of the band the grid covers. *Not* "the centre cell equals the
        # old scalar": at the centre latitude the factor is `sin(h)/h`, short of 1 by ~1.1e-5 for
        # these cells, because the nominal `Δφ × 111.32 km/°` overstates a meridional arc.
        # **The band comes from the lookup's own cell edges, not from `centre ± half`**
        # (`[LOCUS-BLIND]`). It was `sin(max(lats) + h) − sin(min(lats) − h)`, which silently assumed
        # the coordinates were cell *centres* — so when `_testraster` was corrected to the
        # `Intervals(Start)` its readers really produce, the test's own expectation described
        # different ground from the cells it was checking. Asking `_cellintervals` is exact and
        # cannot encode a locus again.
        R = uconvert(km, EcoSISTEM.LONGITUDE_DEGREE_LENGTH)  # km/° *is* a radius: ° is dimensionless
        ivs = EcoSISTEM._cellintervals(grid(lats, longs), Y)
        band = R^2 * ustrip(u"rad", abs(longs[2] - longs[1]) * length(longs)) *
               (sin(maximum(last, ivs)) - sin(minimum(first, ivs)))
        @test sum(areas) * length(longs) ≈ band rtol=1.0e-12

        # **Why this was invisible for so long**: the old scalar is right in aggregate and badly
        # wrong per cell, so anything checking a total — or a uniform grid — would have passed.
        # Measured on this grid: **1.3% out in total, 1.2%–15.5% out per cell**, an order of
        # magnitude apart. Asserted on the **worst** cell rather than on `areas[1]`, which under
        # the corrected locus happens to be only 9.4% out — a threshold that picked one arbitrary
        # cell was always going to be fragile against a change in where the cells sit.
        old = EcoSISTEM._cellsize(lats, longs)^2 * length(lats) * length(longs)
        percell = old / length(lats) / length(longs)
        @test isapprox(old, band, rtol = 2.0e-2)          # total: nearly right
        @test maximum(abs(percell - a) / a for a in vec(areas)) > 0.15  # per cell: not

        # A projected grid is exact already, so the factor is a plain `1.0` and the area stays a
        # scalar — no array is introduced where the old code had none.
        pl = collect(range(0.0km, 50.0km, length = 11))
        pg = collect(range(0.0km, 30.0km, length = 7))
        pgrid = _bngraster(WorldClim{BioClim}, fill(291.0K, 11, 7),
                           north = pl, east = pg).array
        @test EcoSISTEM._areafactor(EcoSISTEM._cellintervals(pgrid, Y)) === 1.0
        @test EcoSISTEM._cellareas(pgrid) isa Number
        @test EcoSISTEM._cellareas(pgrid) ≈ EcoSISTEM._cellsize(pl, pg)^2

        # **Locus-blindness, the property that stops a labelling change propagating**
        # (`[LOCUS-BLIND]`): the *same cells* described `Intervals(Start)` and `Intervals(Center)`
        # must give identical answers from every accessor. A raw `parent(lookup(...))` would not.
        L = DimensionalData.Lookups
        mk(loc, v) = Y(Rasters.Projected(collect(v), crs = Rasters.EPSG(27700),
                                         order = L.ForwardOrdered(),
                                         span = L.Regular(10.0km),
                                         sampling = L.Intervals(loc)))
        st = DimArray(zeros(3),
                      (mk(L.Start(), range(0.0km, step = 10.0km,
                                           length = 3)),))
        ct = DimArray(zeros(3),
                      (mk(L.Center(), range(5.0km, step = 10.0km,
                                            length = 3)),))
        @test EcoSISTEM._cellintervals(st, Y) == EcoSISTEM._cellintervals(ct, Y)
        @test EcoSISTEM._cellcentres(st, Y) == EcoSISTEM._cellcentres(ct, Y)
    end
end
@testset "multi-raster environments" begin
    temp = _testraster(WorldClim{BioClim}, fill(290.0K, 5, 5))
    rain = _testraster(WorldClim{BioClim}, fill(50.0mm / day, 5, 5))
    landcover = _testraster(EarthEnv{LandCover}, Float64.(repeat(1:5, 1, 5)))

    @testset "two regimes + two supplies" begin
        env = _env((_reg(temp, axis = Temperature),
                    _reg(rain, axis = Precipitation)),
                   (UniformSpec(1.0kJ / (km^2 * day),
                                axis = SolarRadiation),
                    UniformSpec(1.0mm / day, axis = Precipitation)))
        @test length(values(env.regime)) == 2
        @test length(values(env.supply)) == 2
        # Named by axis — `_reg(temp)`/`_reg(rain)` resolve to distinguishable ones.
        @test Base.size(first(values(env.regime)).matrix) == (5, 5)
    end

    @testset "temperature + land cover, and a 3-tuple (both continuous)" begin
        env = _env((warmth = _reg(temp), cover = _reg(landcover)), SUP)
        @test map(EcoSISTEM.iscontinuous, values(env.regime)) == (true, true)
        @test length(values(_env((warmth = _reg(temp), wet = _reg(rain),
                                  cover = _reg(landcover)),
                                 SUP).regime)) == 3
    end

    @testset "arity is limited only by what is written" begin
        four = _env((a = _reg(temp), b = _reg(rain), c = _reg(landcover),
                     d = _reg(temp)), SUP)
        @test length(values(four.regime)) == 4
        # a 1-tuple is not a multi-layer environment: pass the spec on its own
        @test_throws ErrorException _env((_reg(temp),), SUP)
    end

    @testset "every supply is built, however many there are" begin
        env = _env((a = _reg(temp), b = _reg(rain), c = _reg(landcover)),
                   (p = UniformSpec(1.0kJ / (km^2 * day),
                                    axis = SolarRadiation),
                    q = UniformSpec(2.0kJ / (km^2 * day),
                                    axis = SolarRadiation),
                    r = UniformSpec(3.0kJ / (km^2 * day),
                                    axis = SolarRadiation)))
        @test length(values(env.supply)) == 3
    end

    @testset "named layers keep their names" begin
        env = _env((temperature = _reg(temp), rainfall = _reg(rain)), SUP)
        @test keys(env.regime) == (:temperature, :rainfall)
        @test Base.size(env.regime.temperature.matrix) == (5, 5)
    end

    @testset "default active = AND of non-missing cells" begin
        withnan = fill(290.0K, 5, 5)
        withnan[1, 1] = NaN * K
        env = _env((warmth = _reg(_testraster(WorldClim{BioClim}, withnan)),
                    wet = _reg(rain)), SUP)
        @test env.active isa DimensionalData.AbstractDimArray{Bool}
        @test !env.active[1, 1] && env.active[2, 2]
    end
end
@testset "SolarRadiation: differently-scaled values reconcile to one canonical unit" begin
    # WorldClim's srad ships in kJ·m⁻²·day⁻¹, CHELSA's rsds* in MJ·m⁻²·day⁻¹ — this is the bug
    # `canonicalunit(::Type{<:SolarRadiationAxis})` fixes: a tolerance built (by default) in the canonical unit
    # and a regime built from an honestly-tagged-but-differently-scaled raster value must still
    # reconcile, not silently disagree by a factor of 1000.
    sp = build_species(2,
                       tolerance = (500.0kJ / (m^2 * day),
                                    100.0kJ / (m^2 * day)),
                       demand = DEM, demandaxis = SolarRadiation,
                       toleranceaxis = SolarRadiation, seed = 1)
    @test sp.tolerance isa NicheTolerance{SolarRadiation}
    # projected (BNG) fixture, since this test simulates — see the multi-regime testset above
    mjraster = _bngraster(WorldClim{Climate},
                          fill(0.5MJ / (m^2 * day), 9, 9))
    env = _env(_reg(mjraster, axis = SolarRadiation), SUP)
    # the regime is canonicalised to kJ·m⁻²·day⁻¹ regardless of the MJ input
    canon = typeof(1.0 * EcoSISTEM.canonicalunit(SolarRadiation))
    @test eltype(env.regime) == canon
    eco = build_ecosystem(sp, env, seed = 1)
    @test eco.nichefit isa NicheSuitability{SolarRadiation, canon}
    @test_nowarn simulate!(eco, 2month_mean_duration, 1month_mean_duration)
end
# **A declared axis must decide a supply's type, never the value's dimension.** `m/s` and `mm/day`
# are both `𝐋𝐓⁻¹`, so a dimension-driven supply path turns `UniformSpec(3.0m/s, axis = WindSpeed)`
# into a water supply — a wrong answer rather than an error. Asking the axis makes that
# unrepresentable, and this test is what holds the supply path to it.
# **An already-built layer is accepted in either role**, by `materialise`'s `AbstractLayer`
# method, which checks it is on this area's grid and passes it through. That one method serves both
# roles: a built layer carries its axis in its own type, so it has said what it means — the property
# a bare raster lacks and the whole refusal turns on.
@testset "an already-built regime is accepted, on its own grid" begin
    # The axis is named on the wrapper, so the assertion below is that a built layer *keeps*
    # it — the property that makes accepting one safe.
    reg = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)),
               axis = Temperature)
    sun = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
    area = StudyArea(regime = reg, verbosity = :silent)
    built = EcoSISTEM.materialise(reg, area, role = EcoSISTEM.Condition)
    @test built isa EcoSISTEM.AbstractRegime

    # **A17: two layers built on ONE area must agree they are on one grid.** The study area's
    # template is deliberately unitless — `Rasters.resample` needs a bare `to =` target — and the data
    # path re-attaches the unit afterwards in `_targetyx`. A **synthetic** layer never reprojects, so
    # without `_unitedyx` in `_materialisefield` it keeps the bare numbers, and `_yxcompatible` then
    # reads two descriptions of the *same* grid as two different grids: coordinate values identical
    # to the last digit, rejected with "`supply`'s grid does not match the regime's grid" on the unit
    # alone.
    builtsun = EcoSISTEM.materialise(sun, area, role = EcoSISTEM.Resource)
    @test EcoSISTEM._yxcompatible(dims(built.matrix, (Y, X)),
                                  dims(builtsun.matrix, (Y, X)))
    @test unit(eltype(DimensionalData.lookup(builtsun.matrix, Y))) ==
          unit(eltype(DimensionalData.lookup(built.matrix, Y)))
    # …and the pair actually builds, which is the failure A17 was reported as.
    @test_nowarn GridHabitat(regime = built, supply = builtsun,
                             area = area)

    env = GridHabitat(regime = built, supply = sun, area = area)
    @test env.regime isa EcoSISTEM.AbstractRegime
    # Used **as it stands**: the same values, not resampled and not canonicalised a second time.
    @test env.regime.matrix == built.matrix
    # …and it keeps the axis it was built with, which is where its meaning lives.
    @test EcoSISTEM._layeraxis(env.regime) === Temperature

    # **A layer built on a different grid is refused, not silently adopted.** It matters more here
    # than for a supply: `GridHabitat` checks the supply and mask *against the regime*, so a wrongly
    # gridded regime would become the reference everything else is validated against.
    other = StudyArea(extent = (40.0km, 40.0km), cellsize = 10.0km,
                      verbosity = :silent)
    elsewhere = EcoSISTEM.materialise(UniformSpec(291.0K, axis = Temperature),
                                      other, role = EcoSISTEM.Condition)
    err = try
        GridHabitat(regime = elsewhere, supply = sun, area = area)
    catch e
        sprint(showerror, e)
    end
    @test occursin("on its own grid", err)

    # …and a synthetic area refuses it too, but on the **grid**, not on the area's kind: `built`
    # is 9 × 9 on the BNG grid and `other` is 4 × 4, so it is the size check that turns it away.
    err2 = try
        GridHabitat(regime = built, supply = sun, area = other)
    catch e
        sprint(showerror, e)
    end
    @test occursin("on its own grid", err2)

    # **A layer built on the synthetic area itself IS accepted** — the correction 6f made, and a
    # correction rather than a loosening: the blanket refusal was written when a synthetic grid had
    # `NoLookup` dims and so could not say where anything was, and `[S4]` gave it real coordinates.
    # **The supply must be built on `other` itself, or this tests the wrong thing.** A built layer
    # whose grid does *not* match is caught by the shape check above, so this is the only case the
    # relaxed rule reaches — and it is exactly the one `build_habitat` needs, reseeding a habitat
    # from its own as-built area.
    fitting = EcoSISTEM.materialise(sun, other, role = EcoSISTEM.Resource)
    ownsun = GridHabitat(regime = UniformSpec(291.0K, axis = Temperature),
                         supply = fitting, area = other)
    @test ownsun.supply.matrix == fitting.matrix
    # …and still copied in rather than aliased, which the built-layer path owes either way.
    @test parent(ownsun.supply.matrix) !== parent(fitting.matrix)

    # **The general claim is untouched: arbitrary data still cannot be sampled onto a synthetic
    # grid**, which has no CRS to place it by. That is what the relaxation must not have reached.
    ownreg = EcoSISTEM.materialise(UniformSpec(291.0K, axis = Temperature),
                                   other, role = EcoSISTEM.Condition)
    @test_nowarn GridHabitat(regime = ownreg, supply = sun, area = other)
    err3 = try
        GridHabitat(regime = reg, supply = sun, area = other)
    catch e
        sprint(showerror, e)
    end
    @test occursin("synthetic study area", err3)
    @test occursin("no CRS", err3)
end
# **A habitat OWNS its layers: a pre-built one is COPIED in, never aliased.**
#
# Why this is a correctness test and not tidiness: `_applychange!` writes `layer.matrix .= …` on
# **every timestep**, so a layer stored by reference would make the caller's own object live
# simulation state — and two habitats handed the same layer would share it. `_zerogaps` names this
# hazard but rebuilds only where there are *gaps* to clean, so it cannot be relied on to copy: a
# gap-free layer passes straight through it.
#
# The spec path was never affected — `materialise` builds fresh arrays every call — which is why
# this is asserted on the **built-layer** path specifically.
@testset "a pre-built layer is copied into the habitat, not aliased" begin
    reg = _reg(_bngraster(WorldClim{BioClim}, fill(291.0K, 9, 9)),
               axis = Temperature)
    area = StudyArea(regime = reg, verbosity = :silent)
    built = EcoSISTEM.materialise(reg, area, role = EcoSISTEM.Condition)
    sun = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)
    builtsun = EcoSISTEM.materialise(sun, area, role = EcoSISTEM.Resource)

    env = GridHabitat(regime = built, supply = builtsun, area = area)

    # The values are still what was handed in — copying must not change what is built.
    @test env.regime.matrix == built.matrix
    @test env.supply.matrix == builtsun.matrix

    # …but the habitat's layers are its own objects, sharing no storage with the caller's.
    @test env.regime !== built
    @test env.supply !== builtsun
    @test parent(env.regime.matrix) !== parent(built.matrix)
    @test parent(env.supply.matrix) !== parent(builtsun.matrix)

    # The property that matters, stated as the mutation it prevents: a timestep writing the
    # habitat's regime must not reach back into the layer the caller still holds.
    before = copy(parent(built.matrix))
    parent(env.regime.matrix) .= 999.0K
    @test parent(built.matrix) == before

    beforesun = copy(parent(builtsun.matrix))
    parent(env.supply.matrix) .*= 2
    @test parent(builtsun.matrix) == beforesun

    # …and two habitats built from one layer are independent of each other.
    a = GridHabitat(regime = built, supply = builtsun, area = area)
    b = GridHabitat(regime = built, supply = builtsun, area = area)
    @test parent(a.regime.matrix) !== parent(b.regime.matrix)
    beforeb = copy(parent(b.regime.matrix))
    parent(a.regime.matrix) .= 123.0K
    @test parent(b.regime.matrix) == beforeb

    # The **change** is deliberately shared, not copied: nothing mutates one in place
    # (`_repointseries!` rebuilds, `_cleanstored` copies before cleaning, `setchange!` replaces the
    # field), so copying it would clone a `SeriesLayerChange`'s whole slice stack for nothing.
    @test env.regime.change === built.change
end
@testset "a non-resource axis cannot be built as a supply" begin
    square = _area(extent = (5km, 5km), cellsize = 1km)
    flat = UniformSpec(298.0K, axis = Temperature)

    # The control: a real resource axis still builds, so the refusal below is about resource-hood
    # and not about synthetic supplies having broken.
    solar = GridHabitat(regime = flat,
                        supply = UniformSpec(10.0 *
                                             EcoSISTEM.canonicalunit(EcoSISTEM.Resource,
                                                                     SolarRadiation) /
                                             m^2,
                                             axis = SolarRadiation),
                        area = square)
    @test EcoSISTEM.axisof(solar.supply) === SolarRadiation

    # **Asserted on the dimensional collision itself, not against a unit-keyed lookup table.**
    # There is no such table on either side, and there should not be: stating the *reason* a unit
    # cannot decide a quantity's meaning outlives any lookup that happens to encode it.
    #
    # The collision is not in `m/s` itself — that is `𝐋𝐓⁻¹`, while `VolumeFlow` is `𝐋³𝐓⁻¹`. It
    # appears one step later, because the supply path multiplies by the **cell area** first: a
    # 3 m/s wind over a 1 km² cell becomes 2.592e14 L/day, which is **dimensionally identical** to
    # rainfall over the same cell. Any lookup keyed on dimension must confuse the two.
    windflow = EcoSISTEM.cancel(3.0m / s, 1.0km^2)
    rainflow = EcoSISTEM.cancel(2.0mm / day, 1.0km^2)
    @test dimension(windflow) == dimension(rainflow)
    @test dimension(windflow) == dimension(1.0Unitful.L / day)
    # …and it really is a wind speed that got there, not a quantity that was water all along.
    @test dimension(3.0m / s) != dimension(2.0mm / day * 1.0km^2)

    # …and the supply side must not consult the unit at all. The error must name the axis, since
    # "this is not a resource" is the whole content of the failure.
    err = try
        GridHabitat(regime = flat,
                    supply = UniformSpec(3.0m / s, axis = WindSpeed),
                    area = square)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test occursin("WindSpeed", err)
    @test occursin("not a consumable resource", err)

    # And a bare value — the form that could only work by guessing from the unit — is refused by
    # the **signature**, `supply::Union{LayerInput, AbstractSupply}`, so the error names the keyword
    # and prints what it does accept. The trade is deliberate: a signature that documents itself,
    # at the cost of a message tailored to this particular mistake.
    err2 = try
        GridHabitat(regime = flat, supply = 3.0m / s, area = square)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err2)
    @test occursin("in keyword argument supply", err2)
    @test occursin("AbstractSpec", err2)
    # The tailored message survives where a signature cannot reach — inside a container, since a
    # multi-layer supply legitimately is a tuple.
    err3 = try
        GridHabitat(regime = flat, supply = (3.0m / s, SUP),
                    area = square)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test occursin("no niche axis", err3)
end

end
