# SPDX-License-Identifier: LGPL-3.0-or-later

module TestCollections

using EcoSISTEM
using Diversity
# `nichefitcombine` is `public`, not exported.
using EcoSISTEM: nichefitcombine
using EcoSISTEM.Units
using DimensionalData
using Distributions: Normal
using RasterDataSources
using Rasters
using Unitful, Unitful.DefaultSymbols
using Test

include("rasterfixtures.jl")

# Wrap an already-in-memory `ClimateRaster` as a lazy layer spec, and build on the grid the layers
# themselves decide — the same escape hatch `test_GridHabitat.jl` uses, since a multi-layer regime has
# to be data-backed (a synthetic area has no CRS, so there is nothing to place the data by).
_area(; kw...) = StudyArea(; verbosity = :silent, kw...)
function _env(regime, supply; kw...)
    return GridHabitat(regime = regime, supply = supply,
                       area = _area(; regime = regime, kw...))
end
const SUP = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation)

# The hot-path allocation measurements, behind a function barrier: `eco` arrives as a typed
# argument, so the test's own (necessarily type-unstable) local cannot put an allocation on the
# measurement, and the first call warms the same specialisation the second is measured in.
function _suitalloc(eco)
    EcoSISTEM.suitability(eco, 1, 1)
    return @allocated EcoSISTEM.suitability(eco, 1, 1)
end
function _adjustalloc(eco)
    EcoSISTEM._resourceadjustment(eco, eco.habitat.supply, 1, 1)
    return @allocated EcoSISTEM._resourceadjustment(eco, eco.habitat.supply, 1,
                                                    1)
end

@testset "Collection backing" begin
    @testset "folds and zipmaps are unrolled and allocation-free" begin
        t2 = (1.0, 2.0)
        t10 = ntuple(i -> Float64(i), 10)
        fold2() = EcoSISTEM._fold(*, t2)
        fold10() = EcoSISTEM._fold(*, t10)
        zip3() = EcoSISTEM._zipmap(*, t2, t2, t2)
        fold2(), fold10(), zip3()
        @test fold2() == 2.0
        @test fold10() == prod(t10)
        @test zip3() == (1.0, 8.0)
        @test (@allocated fold10()) == 0
        @test (@allocated zip3()) == 0
        @test @inferred(fold10()) isa Float64
        # collections that must line up member for member are rejected when they do not
        @test_throws ErrorException EcoSISTEM._zipmap(+, (1.0, 2.0),
                                                      (1.0, 2.0, 3.0))
    end

    @testset "the fold and the zipmap stay free at every arity" begin
        # `map` stops being unrolled at around seven elements and starts allocating; the recursion
        # must not, at any arity a caller can reach. Pinned at the arities the design was measured
        # at rather than only at the two the hot paths happen to use today.
        for n in (2, 3, 6, 10)
            a = ntuple(i -> Float64(i), n)
            function zipped()
                return EcoSISTEM._zipmap(a, a, a) do x, y, z
                    return x * y * z
                end
            end
            folded() = EcoSISTEM._fold(+, zipped())
            zipped()
            folded()
            @test (@allocated zipped()) == 0
            @test (@allocated folded()) == 0
            @test @inferred(folded()) isa Float64
        end
    end
end

@testset "Layer collections" begin
    grid = (5, 5)
    cellsize = 1.0km
    temp = EcoSISTEM.simpleregime(298.0K, cellsize, grid, Temperature)
    rain = EcoSISTEM.simpleregime(50.0mm / day, cellsize, grid,
                                  Precipitation)
    solar = Supply{SolarRadiation}(fill(10.0kJ / day, grid...))
    water = Supply{Precipitation}(fill(10.0Unitful.L / day, grid...))

    # **A collection built from a plain `Tuple` is named by its members' AXES**, where those are
    # distinguishable — so this pair is `(:Temperature, :Precipitation)`, not `(:one, :two)`. Two
    # structures that pair correctly therefore derive the *same* names with the caller writing
    # nothing, which is what `_checknaming` compares.
    # `.one`/`.two` survive only where they *are* the names — see the fallback testset below.
    # Nothing is owed for the change: the `:one`/`:two` scheme is absent from v0.4.0.
    @testset "a positional collection is named by its members' axes" begin
        lc = LayerCollection((temp, rain))
        @test lc isa EcoSISTEM.LayerCollection{EcoSISTEM.Condition}
        @test lc.Temperature === temp
        @test lc.Precipitation === rain
        @test values(lc) === (temp, rain)
        @test keys(lc) == (:Temperature, :Precipitation)
        @test NamedTuple(lc) == (Temperature = temp, Precipitation = rain)
        @test propertynames(lc) == (:Temperature, :Precipitation)
        # A `FieldError`, because member access is plain `getfield` on the backing, whose own
        # message already names the alternatives.
        @test_throws FieldError lc.three
    end

    # **The fallback, and it is what keeps every existing positional call site working.** Two
    # Members on the same axis cannot be told apart by name, so the collection **errors and asks
    # for names** rather than falling back to `:one`/`:two`. The model itself is never refused — a
    # repeated axis is legitimate (`CombiningFit`'s own docstring pairs a summer and a winter
    # temperature) — only the *guessing* is. This fell back silently until 2026-08-18; the
    # fallback made sense when most members had no axis at all, and step A removed that.
    @testset "a repeated axis errors, asking for names" begin
        msg = try
            LayerCollection((temp, temp))
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("Temperature", msg) &&
              occursin("cannot be told apart", msg)
        # …and naming them is the remedy the message gives, so it must work.
        named = LayerCollection((summer = temp, winter = temp))
        @test keys(named) == (:summer, :winter)
        @test named.summer === temp
        @test named.winter === temp
    end

    @testset "named access uses the caller's own names" begin
        lc = LayerCollection((temperature = temp, rainfall = rain))
        @test lc.temperature === temp
        @test lc.rainfall === rain
        @test keys(lc) == (:temperature, :rainfall)
        @test NamedTuple(lc) == (temperature = temp, rainfall = rain)
        # a named collection is a different type from a positional one, so the two never mix
        @test lc isa LayerCollection{EcoSISTEM.Condition, <:Any, <:NamedTuple}
    end

    @testset "a bare layer answers as a one-member collection" begin
        @test values(temp) === (temp,)
        @test values(solar) === (solar,)
        # **Named by its own AXIS**, not by a role label. It read `(:regime,)`/`(:supply,)` until
        # 2026-08-18, and those could never match the species side's `(:tolerance,)`/`(:demand,)` —
        # so single-member pairings had to be *skipped* by the naming check. Both sides now derive
        # the same name from the same axis, and the skip is gone with the labels.
        @test keys(temp) == (:Temperature,)
        @test keys(solar) == (:SolarRadiation,)
    end

    @testset "the role is dispatched on, and must agree" begin
        # **The role guarantee lives on the side builders**, which is the one place a mislabelled
        # side would actually do damage. It cannot live in the container interface, which is
        # role-blind by design — `values(x)` works on either side.
        @test_throws MethodError EcoSISTEM._regimeside(solar)
        @test_throws MethodError EcoSISTEM._supplyside(temp)
        @test_throws ErrorException LayerCollection((temp, solar))
    end

    @testset "arity is limited only by what is written" begin
        four = LayerCollection((a = temp, b = rain, c = temp, d = rain))
        @test length(values(four)) == 4
        @test four.d === rain
        @test map(EcoSISTEM.iscontinuous, values(four)) == ntuple(_ -> true, 4)
        @test length(values(LayerCollection((a = solar, b = water,
                                             c = solar)))) == 3
    end

    @testset "the standard container interface" begin
        lc = LayerCollection((temp, rain))
        # Indexing by position **and** by name, iteration, and the rest — forwarded to the
        # backing `NamedTuple` (`src/collections.jl`).
        @test lc[1] === temp
        @test lc[:Precipitation] === rain
        @test length(lc) == 2
        @test keys(lc) == (:Temperature, :Precipitation)
        @test values(lc) === (temp, rain)
        @test collect(x for x in lc) == [temp, rain]
        @test haskey(lc, :Temperature)
        @test !haskey(lc, :Nope)
        @test get(lc, :Nope, nothing) === nothing
        @test NamedTuple(lc) === (Temperature = temp, Precipitation = rain)
        @test propertynames(lc) == keys(lc)
        # **A leaf answers identically as a one-member container**, which is the whole point:
        # nothing downstream has to know whether it holds one layer or several.
        @test keys(temp) == (:Temperature,)
        @test values(temp) === (temp,)
        @test length(temp) == 1
        @test temp[1] === temp
        @test temp[:Temperature] === temp
        @test NamedTuple(temp) === (Temperature = temp,)
        # **Iterating a mixed collection must not allocate** — that claim is what the old
        # "no container interface" rule rested on, and it was measured backwards: the fold it
        # protected allocated where the loop does not. Guard it so it cannot regress.
        count_members(c) = (n = 0; for _ in c
                                n += 1
                            end; n)
        count_members(lc)
        @test @allocated(count_members(lc)) == 0
    end
end

@testset "Species-side collections" begin
    n = 5
    temptol = NicheTolerance(Temperature, Normal, fill(298.0K, n),
                             fill(2.0K, n))
    raintol = NicheTolerance(Precipitation, Normal, fill(50.0mm / day, n),
                             fill(5.0mm / day, n))

    @testset "tolerances" begin
        tc = SpeciesRequirementCollection((temptol, raintol))
        # Axis-derived, exactly as the layer side is — which is *why* a tolerance collection and
        # the regime collection it pairs with check out against one another unaided.
        @test tc.Temperature === temptol
        @test values(tc) === (temptol, raintol)
        @test keys(tc) == (:Temperature, :Precipitation)
        named = SpeciesRequirementCollection((temperature = temptol,
                                              rainfall = raintol))
        @test named.rainfall === raintol
        @test keys(named) == (:temperature, :rainfall)
        @test values(temptol) === (temptol,)
    end

    @testset "demands, and one totalE column per resource" begin
        dc = SpeciesRequirementCollection((Demand{SolarRadiation}(fill(10.0kJ /
                                                                       day, n)),
                                           Demand{Precipitation}(fill(1.0Unitful.L /
                                                                      day,
                                                                      n))))
        @test length(values(dc)) == 2
        @test EcoSISTEM.numdemands(typeof(dc)) == 2
        # `length` is the **member** count, like any container. The *species* count is Diversity's
        # `counttypes`, and the two must not be confused on the Resource side.
        @test length(dc) == 2
        @test counttypes(dc) == n
        three = SpeciesRequirementCollection((a = Demand{SolarRadiation}(fill(10.0kJ /
                                                                              day,
                                                                              n)),
                                              b = Demand{Precipitation}(fill(1.0Unitful.L /
                                                                             day,
                                                                             n)),
                                              c = Demand{SolarRadiation}(fill(5.0kJ /
                                                                              day,
                                                                              n))))
        @test EcoSISTEM.numdemands(typeof(three)) == 3
    end

    @testset "nichefits combine by a whole-tuple function" begin
        f = MultiplicativeFit((NicheSuitability(temptol),
                               NicheSuitability(raintol)))
        @test f isa EcoSISTEM.CombiningFit
        @test length(values(f)) == 2
        # Named by AXIS, not by position — a nichefit carries its own axis since D3(a), so a
        # `CombiningFit` derives the same names its paired tolerance and regime do. It read
        # `(:one, :two)` until then, because every fit answered the root and so clashed.
        @test keys(f) == (:Temperature, :Precipitation)
        @test nichefitcombine(f)((Temperature = 2.0, Precipitation = 3.0)) ==
              6.0
        @test nichefitcombine(AdditiveFit((NicheSuitability(temptol),
                                           NicheSuitability(raintol))))((Temperature = 2.0,
                                                                         Precipitation = 3.0)) ==
              5.0
        # a lone fit combines by taking the only result
        @test nichefitcombine(NicheSuitability(temptol))((only = 4.0,)) == 4.0
        # an arbitrary rule reads by name, and so cannot be reordered by mistake
        custom = CombiningFit((temperature = NicheSuitability(temptol),
                               rainfall = NicheSuitability(raintol))) do s
            return s.temperature + 2 * s.rainfall
        end
        @test nichefitcombine(custom)((temperature = 1.0, rainfall = 3.0)) ==
              7.0
        @test custom.rainfall isa NicheSuitability
    end
end

@testset "Named layers reach the built ecosystem" begin
    n = 5
    # A projected (British National Grid) fixture: a degree grid cannot be simulated on, since
    # dispersal needs one uniform cell size.
    temp = _bngraster(WorldClim{BioClim}, fill(298.0K, 9, 9))
    rain = _bngraster(WorldClim{BioClim}, fill(50.0mm / day, 9, 9))
    env = _env((temperature = _reg(temp, axis = Temperature),
                rainfall = _reg(rain, axis = Precipitation)), SUP)
    @test keys(env.regime) == (:temperature, :rainfall)
    @test env.regime.temperature isa EcoSISTEM.ContinuousRegime

    species = build_species(n,
                            tolerance = (temperature = (298.0K, 2.0K),
                                         rainfall = (50.0mm / day,
                                                     5.0mm / day)),
                            toleranceaxis = (Temperature, Precipitation),
                            demand = 10.0kJ / day, demandaxis = SolarRadiation)
    @test keys(species.tolerance) == (:temperature, :rainfall)

    eco = build_ecosystem(species, env)
    @test eco.habitat.regime.temperature isa EcoSISTEM.ContinuousRegime
    @test keys(eco.nichefit) == (:temperature, :rainfall)
    @test EcoSISTEM.suitability(eco, 1, 1) > 0
    @test_nowarn simulate!(eco, 2month_mean_duration, 1month_mean_duration)

    @testset "paired structures must line up, and the error names both sides" begin
        positional = build_species(n,
                                   tolerance = ((298.0K, 2.0K),
                                                (50.0mm / day, 5.0mm / day)),
                                   toleranceaxis = (Temperature, Precipitation),
                                   demand = 10.0kJ / day,
                                   demandaxis = SolarRadiation)
        # named on one side, positional on the other is its own error
        @test_throws ErrorException build_ecosystem(positional, env)
        # so is the same arity under different names, in a different order
        swapped = build_species(n,
                                tolerance = (rainfall = (50.0mm / day,
                                                         5.0mm / day),
                                             temperature = (298.0K, 2.0K)),
                                toleranceaxis = (Precipitation, Temperature),
                                demand = 10.0kJ / day,
                                demandaxis = SolarRadiation)
        @test_throws ErrorException build_ecosystem(swapped, env)
        # and so is a plain arity mismatch
        onetol = build_species(n, tolerance = (298.0K, 2.0K),
                               toleranceaxis = Temperature,
                               demand = 10.0kJ / day,
                               demandaxis = SolarRadiation)
        @test_throws ErrorException build_ecosystem(onetol, env)
    end

    @testset "a four-layer regime builds and simulates" for named in (false, true)
        specs = (a = _reg(temp, axis = Temperature),
                 b = _reg(rain, axis = Precipitation),
                 c = _reg(temp, axis = Temperature),
                 d = _reg(rain, axis = Precipitation))
        tols = ((298.0K, 2.0K), (50.0mm / day, 5.0mm / day), (298.0K, 2.0K),
                (50.0mm / day, 5.0mm / day))
        axes = (Temperature, Precipitation, Temperature, Precipitation)
        keys4 = (:summer, :wet, :winter, :dry)
        # **Both branches are named now, and they must be**: two of the four layers are on
        # `Temperature` and two on `Precipitation`, so derived names cannot tell them apart and the
        # collection says so. What the `named` sweep still varies is *which* names — the caller's
        # own against the fixture's `(:a, :b, :c, :d)` — which is what the pairing check consumes.
        four = _env(named ? NamedTuple{keys4}(values(specs)) : specs, SUP)
        @test length(values(four.regime)) == 4
        @test keys(four.regime) == (named ? keys4 : keys(specs))
        spp4 = build_species(n,
                             tolerance = named ? NamedTuple{keys4}(tols) :
                                         NamedTuple{keys(specs)}(tols),
                             toleranceaxis = axes, demand = 10.0kJ / day,
                             demandaxis = SolarRadiation)
        eco4 = build_ecosystem(spp4, four)
        named && @test eco4.habitat.regime.winter isa EcoSISTEM.ContinuousRegime
        @test EcoSISTEM.suitability(eco4, 1, 1) > 0
        @test_nowarn simulate!(eco4, 2month_mean_duration, 1month_mean_duration)
    end

    @testset "every supply is built, however many there are" begin
        env3 = _env((a = _reg(temp, axis = Temperature),
                     b = _reg(rain, axis = Precipitation),
                     c = _reg(temp, axis = Temperature)),
                    (p = UniformSpec(1.0kJ / (km^2 * day),
                                     axis = SolarRadiation),
                     q = UniformSpec(2.0kJ / (km^2 * day),
                                     axis = SolarRadiation),
                     r = UniformSpec(3.0kJ / (km^2 * day),
                                     axis = SolarRadiation)))
        @test length(values(env3.supply)) == 3
        # The collection is parameterised on its **role** first, so the arity is read off the
        # second parameter: `numdemands` sizes the `totalE` cache from the type alone, without an
        # instance, and it must keep doing that now that the two collections have merged.
        @test EcoSISTEM.numdemands(EcoSISTEM.SpeciesRequirementCollection{EcoSISTEM.Resource,
                                                                          Tuple{SolarRadiation,
                                                                                Precipitation,
                                                                                SolarRadiation},
                                                                          @NamedTuple{a::Int,
                                                                                      b::Int,
                                                                                      c::Int}}) ==
              3
    end
end

@testset "Named layers are what StudyArea reports and aligns by" begin
    # Step 5's visible payoff: a named multi-layer `regime` is reported and aligned by the caller's
    # own names, where a positional one falls back to `:regime1`/`:regime2`.
    coarse = _bngraster(WorldClim{BioClim}, fill(298.0K, 9, 9))
    fine = _bngraster(WorldClim{BioClim}, fill(50.0mm / day, 17, 17),
                      east = (245000.0:1250.0:265000.0) .* m,
                      north = (640000.0:1250.0:660000.0) .* m)
    named = (temperature = _reg(coarse, axis = Temperature),
             rainfall = _reg(fine, axis = Precipitation))

    # aligning to each named layer picks that layer's own step
    @test _area(regime = named, align = :temperature).report.cellsize ==
          2500.0m
    @test _area(regime = named, align = :rainfall).report.cellsize == 1250.0m
    # the per-layer report is keyed by those names too
    @test Set(l.name for l in _area(regime = named).report.layers) ==
          Set((:temperature, :rainfall))
    # and the numbered fallback is gone once names are given — with both listed in the error
    err = try
        _area(regime = named, align = :regime1)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin(":temperature", err.msg) && occursin(":rainfall", err.msg)

    # a positional tuple still numbers, as before
    positional = (_reg(coarse, axis = Temperature),
                  _reg(fine, axis = Precipitation))
    @test Set(l.name for l in _area(regime = positional).report.layers) ==
          Set((:regime1, :regime2))
    @test _area(regime = positional, align = :regime2).report.cellsize ==
          1250.0m
end

@testset "Hot paths stay allocation-free" begin
    # `_suitability` and `_resourceadjustment` run per cell per species, so the folds have to cost
    # nothing at either arity — this is the gate that keeps `map`/`reduce` out of them.
    temp = _bngraster(WorldClim{BioClim}, fill(298.0K, 9, 9))
    rain = _bngraster(WorldClim{BioClim}, fill(50.0mm / day, 9, 9))
    # Named, because the third layer repeats `Temperature` and derived names cannot tell the
    # two apart — the arity is what this fixture is varying, not the naming.
    # A named tuple does not slice, so the arity sweep takes its first `n` members by name.
    _firstn(nt::NamedTuple, n) = NamedTuple{keys(nt)[1:n]}(values(nt)[1:n])
    specs = (a = _reg(temp, axis = Temperature),
             b = _reg(rain, axis = Precipitation),
             c = _reg(temp, axis = Temperature))
    tols = (a = (298.0K, 2.0K), b = (50.0mm / day, 5.0mm / day),
            c = (298.0K, 2.0K))
    axes = (Temperature, Precipitation, Temperature)

    @testset "suitability over $nlayers layers" for nlayers in (2, 3)
        env = _env(_firstn(specs, nlayers), SUP)
        species = build_species(5, tolerance = _firstn(tols, nlayers),
                                toleranceaxis = axes[1:nlayers],
                                demand = 10.0kJ / day,
                                demandaxis = SolarRadiation)
        eco = build_ecosystem(species, env)
        @test _suitalloc(eco) == 0
        @test @inferred(EcoSISTEM.suitability(eco, 1, 1)) isa Float64
    end

    @testset "resource adjustment over $nres supplies" for nres in (2, 3)
        sup = (a = UniformSpec(1.0kJ / (km^2 * day), axis = SolarRadiation),
               b = UniformSpec(1.0mm / day, axis = Precipitation),
               c = UniformSpec(2.0kJ / (km^2 * day), axis = SolarRadiation))
        dem = (a = 10.0kJ / day, b = 1.0Unitful.L / day, c = 20.0kJ / day)
        env = _env(_reg(temp, axis = Temperature),
                   _firstn(sup, nres))
        species = build_species(5, tolerance = (298.0K, 2.0K),
                                toleranceaxis = Temperature,
                                demand = _firstn(dem, nres),
                                demandaxis = (SolarRadiation, Precipitation,
                                              SolarRadiation)[1:nres])
        eco = build_ecosystem(species, env)
        EcoSISTEM.update_resource_usage!(eco)
        @test size(eco.cache.totalE, 2) == nres
        @test _adjustalloc(eco) == 0
        @test @inferred(EcoSISTEM._resourceadjustment(eco, eco.habitat.supply,
                                                      1, 1)) isa
              Tuple{Float64, Float64}
    end
end

end
