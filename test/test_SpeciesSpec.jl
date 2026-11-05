# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/SpeciesSpec.jl` and `build_species(pool)` in `src/actions.jl`: the trait values and
# specs a species list is built from, the filling of a pool's missing values from streams of each
# species' own, the units a column is read in, and species built from a pool - the same, value for
# value, as the released form builds from the same values, and kept whole by `makeunique`. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_SpeciesSpec.jl"])'

module TestSpeciesSpec

using DataFrames
using Distributions
using Diversity: UniqueTypes
using EcoSISTEM
using EcoSISTEM: AbstractTraitSpec, makeunique
using EcoSISTEM.Units
using Test
using Unitful
using Unitful.DefaultSymbols

include("wfofixtures.jl")

# The two specs every build here needs and no test is about.
const TEMPERATURE = ToleranceSpec(Temperature, μ = 285.0K, σ = 2.0K)
const SOLAR = DemandSpec(SolarRadiation, 1.0e9kJ / day)

# A build that says nothing, the defaults given.
function quietly(pool; kw...)
    return build_species(pool; tolerance = TEMPERATURE, demand = SOLAR,
                         movement = MovementSpec(InsectDispersal(),
                                                 meandist = 1.0km),
                         birth = 0.6 / year, death = 0.6 / year,
                         longevity = 1.0,
                         survival = 0.2, verbosity = :silent, kw...)
end

@testset "trait values and specs say what they read" begin
    @test string(TraitColumn(:tmin, unit = K, fill = 285.0)) ==
          "TraitColumn(:tmin, unit = K, fill = 285.0)"
    drawn = string(DrawFrom(Uniform(0.6, 2.4), unit = km))
    @test startswith(drawn, "DrawFrom(Uniform") &&
          endswith(drawn, ", unit = km)")
    @test_throws "a `fill` is a value" TraitColumn(:x, fill = "zero")
    # A family's parameters by name, placed in its order whatever order they are written in.
    gamma = ToleranceSpec(Precipitation, Gamma, θ = 3.0mm / day, α = 2.0)
    @test gamma.params == (2.0, 3.0mm / day)
    @test gamma isa AbstractTraitSpec
    @test ToleranceSpec(Precipitation, Gamma, 2.0, 3.0mm / day).params ==
          gamma.params
    @test string(gamma) ==
          "ToleranceSpec(Precipitation, Gamma, α = 2.0, θ = 3.0 mm d⁻¹)"
    @test ToleranceSpec(Temperature, μ = 285.0K, σ = 2.0K) isa
          ToleranceSpec{Temperature, Normal}
    @test_throws "`ν` is not one" ToleranceSpec(Temperature, μ = 285.0K,
                                                ν = 2.0K)
    @test_throws "`σ` is missing" ToleranceSpec(Temperature, μ = 285.0K)
    @test_throws "takes 2 parameters" ToleranceSpec(Temperature, Normal, 285.0K)
    # A family whose fields are not its parameters is named in order only.
    @test_throws "give them in order" ToleranceSpec(Precipitation, InverseGamma,
                                                    α = 2.0, θ = 3.0mm / day)
    @test ToleranceSpec(Precipitation, InverseGamma, 2.0, 3.0mm / day) isa
          ToleranceSpec{Precipitation, InverseGamma}
    @test string(ToleranceSpec(Precipitation, InverseGamma, 2.0, 3.0mm / day)) ==
          "ToleranceSpec(Precipitation, InverseGamma, 2.0, 3.0 mm d⁻¹)"
    # Nothing takes values by position.
    @test_throws "never a vector" ToleranceSpec(Temperature,
                                                μ = [285.0K, 286.0K],
                                                σ = 2.0K)
    @test_throws "not a consumable resource" DemandSpec(Temperature, 1.0K)
    @test_throws "`speed` is not one" MovementSpec(InsectDispersal(),
                                                   speed = 1.0km)
    @test_throws "no kernel keywords" MovementSpec(NoDispersal(),
                                                   meandist = 1.0km)
    @test string(MovementSpec(WindDispersal(), meandist = 5.0km,
                              scale = 10.0km)) ==
          "MovementSpec(WindDispersal(WaldKernel, disperse_safely = false); " *
          "meandist = 5.0 km, scale = 10.0 km)"
    # A movement may be chosen per species; a tolerance not yet; each key once.
    @test_throws "cannot be chosen per species yet" ChosenBy(:x,
                                                             "a" => TEMPERATURE)
    @test_throws "names a key twice" ChosenBy(:x,
                                              "a" =>
                                                  MovementSpec(NoDispersal()),
                                              "a" =>
                                                  MovementSpec(NoDispersal()))
end

@testset "missing values are filled, each species from a stream of its own" begin
    table = DataFrame(species = ["w", "x", "y", "z"],
                      sd = [2.0, missing, 3.0, missing],
                      lo = [1.0, missing, 3.0, 4.0])
    pool = SpeciesPool(table)
    specs = (:spread => TraitColumn(:sd, unit = K, fill = Uniform(1.0, 3.0)),
             :lower => TraitColumn(:lo, fill = row -> 10.0),
             :draw => DrawFrom(Uniform(0.0, 1.0), unit = km))
    filled = fill_traits(pool, specs..., seed = 3, verbosity = :silent)
    @test filled.values.sd[[1, 3]] == [2.0, 3.0]
    @test all(1.0 .<= filled.values.sd[[2, 4]] .<= 3.0)
    @test filled.values.lo == [1.0, 10.0, 3.0, 4.0]
    @test eltype(filled.values.sd) == Float64
    @test filled.filled.sd == [false, true, false, true]
    @test filled.filled.lo == [false, true, false, false]
    @test all(filled.filled.draw) && unit(eltype(filled.values.draw)) == km
    # The same seed gives the same values, and a species' value is its own: removing another
    # species, which leaves every key as it was, leaves it unchanged.
    @test isequal(fill_traits(pool, specs..., seed = 3, verbosity = :silent),
                  filled)
    remove_species!(pool, "w", verbosity = :silent)
    fewer = fill_traits(pool, specs..., seed = 3, verbosity = :silent)
    @test fewer.values.sd == filled.values.sd[2:4]
    @test fewer.values.draw == filled.values.draw[2:4]
    @test fill_traits(pool, specs..., seed = 4, verbosity = :silent).values.draw !=
          fewer.values.draw
    # A copy draws its original's values.
    copy_species!(pool, "x", as = "x copy")
    copied = fill_traits(pool, specs..., seed = 3, verbosity = :silent)
    @test copied.values.draw[end] == copied.values.draw[1]
    @test copied.values.sd[end] == copied.values.sd[1]
    # A function may give a distribution, drawn from the species' stream.
    spread = fill_traits(pool,
                         :g => TraitColumn(:sd, fill = row -> Normal(2.0, 0.1)),
                         seed = 3, verbosity = :silent)
    @test all(1.5 .< spread.values.sd[spread.filled.sd] .< 2.5)
    # A missing value with no fill is refused, naming the species and the column.
    @test_throws "`:lo` for `x`" fill_traits(pool, :lo => TraitColumn(:lo),
                                             verbosity = :silent)
    @test_throws "no column `nope`" fill_traits(pool, :n => TraitColumn(:nope),
                                                verbosity = :silent)
    @test_logs (:info, r"`:sd`: 3 missing values filled") fill_traits(pool,
                                                                      :s =>
                                                                          TraitColumn(:sd,
                                                                                      fill = 2.0))
end

@testset "a column is read in its unit, and refused without one where one is needed" begin
    pool = SpeciesPool(DataFrame(species = ["a", "b", "c"],
                                 tmin = [12.0, 15.5, 9.25],
                                 shape = [2.0, 3.0, 4.0]))
    plain = ToleranceSpec(Temperature, μ = TraitColumn(:tmin), σ = 2.0K)
    @test_throws r"needs a `unit`.*mean 12.25, sd 3.132, range 9.25 to 15.5 \(3 values\)" quietly(pool,
                                                                                                  tolerance = plain)
    # Read in its unit, a column is noted in it, which shows a wrong unit for what it is.
    kelvin = ToleranceSpec(Temperature, μ = TraitColumn(:tmin, unit = K),
                           σ = 2.0K)
    @test_logs (:info, r"`:tmin` read in K: mean 12.25 K") match_mode = :any quietly(pool,
                                                                                     tolerance = kelvin,
                                                                                     verbosity = :normal)
    # Degrees C convert as positions; a spread converts as an interval.
    celsius = quietly(pool,
                      tolerance = ToleranceSpec(Temperature,
                                                μ = TraitColumn(:tmin,
                                                                unit = °C),
                                                σ = 2.0K))
    @test mean.(celsius.tolerance.dists) ≈ [285.15, 288.65, 282.4]
    @test std.(celsius.tolerance.dists) ≈ fill(2.0, 3)
    spreadc = quietly(pool,
                      tolerance = ToleranceSpec(Temperature,
                                                μ = TraitColumn(:tmin,
                                                                unit = °C),
                                                σ = 2.0°C))
    @test std.(spreadc.tolerance.dists) ≈ fill(2.0, 3)
    @test_throws "needs a value in K" quietly(pool,
                                              tolerance = ToleranceSpec(Temperature,
                                                                        μ = TraitColumn(:tmin,
                                                                                        unit = km),
                                                                        σ = 2.0K))
    # A dimensionless parameter takes plain numbers, and refuses a unit.
    longtail = quietly(pool,
                       movement = MovementSpec(VertebrateDispersal(),
                                               shape = TraitColumn(:shape),
                                               scale = 1.0km))
    @test [m.kernel.shape for m in longtail.movement] == [2.0, 3.0, 4.0]
    @test_throws "needs no `unit`" quietly(pool,
                                           movement = MovementSpec(VertebrateDispersal(),
                                                                   shape = TraitColumn(:shape,
                                                                                       unit = km),
                                                                   scale = 1.0km))
    @test_throws "one number for every species" quietly(pool, survival = [0.2])
end

@testset "species built from a pool are the released form's, value for value" begin
    old = build_species(5, tolerance = (285.0K, 3.0K),
                        toleranceaxis = Temperature,
                        demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                        seed = 3)
    new = @test_logs((:info, r"`movement` defaulted"),
                     (:info, r"`birth` defaulted"),
                     (:info, r"`death` defaulted"),
                     (:info, r"`longevity` defaulted"),
                     (:info, r"`survival` defaulted"),
                     build_species(SpeciesPool(5),
                                   tolerance = ToleranceSpec(Temperature,
                                                             μ = 285.0K,
                                                             σ = 3.0K),
                                   demand = SOLAR, seed = 3))
    @test new.names == old.names && new.abun == old.abun
    @test new.tolerance.dists == old.tolerance.dists
    @test new.demand.resource == old.demand.resource
    @test new.movement == old.movement &&
          typeof(new.movement) == typeof(old.movement)
    @test typeof(new.params) == typeof(old.params)
    @test (new.params.birth, new.params.death, new.params.longevity,
           new.params.survival) ==
          (old.params.birth, old.params.death, old.params.longevity,
           old.params.survival)
    @test_throws "`tolerance` was not passed" build_species(SpeciesPool(2),
                                                            demand = SOLAR)
    @test_throws "`demand` was not passed" build_species(SpeciesPool(2),
                                                         tolerance = TEMPERATURE)
    # A tolerance built already is taken as it is.
    @test quietly(SpeciesPool(5), tolerance = old.tolerance).tolerance ===
          old.tolerance
end

@testset "a pool's species built with every kind of value" begin
    table = DataFrame(species = ["a", "b", "c", "d"],
                      tmin = [280.0, 285.0, missing, 290.0],
                      tp = [1.0, 2.0, 3.0, 4.0],
                      ssr = [1.0e9, 2.0e9, 3.0e9, 4.0e9],
                      syndrome = ["wind", "insect", "insect", "vertebrate"])
    pool = SpeciesPool(table)
    deactivate_species!(pool, "d", verbosity = :silent)
    built() = quietly(pool,
                      tolerance = (temperature = ToleranceSpec(Temperature,
                                                               μ = TraitColumn(:tmin,
                                                                               unit = K,
                                                                               fill = 285.0),
                                                               σ = DrawFrom(Uniform(1.0,
                                                                                    3.0),
                                                                            unit = K)),
                                   rainfall = ToleranceSpec(Precipitation,
                                                            Gamma, α = 2.0,
                                                            θ = TraitColumn(:tp,
                                                                            unit = mm /
                                                                                   day))),
                      demand = (light = DemandSpec(SolarRadiation,
                                                   TraitColumn(:ssr,
                                                               unit = kJ / day)),
                                water = DemandSpec(Precipitation,
                                                   1.0e3Unitful.L / day)),
                      movement = ChosenBy(:syndrome,
                                          "wind" =>
                                              MovementSpec(WindDispersal(),
                                                           meandist = 5.0km,
                                                           scale = 10.0km),
                                          "insect" =>
                                              MovementSpec(InsectDispersal(),
                                                           meandist = DrawFrom(Uniform(0.6,
                                                                                       2.4),
                                                                               unit = km)),
                                          "vertebrate" =>
                                              MovementSpec(VertebrateDispersal(),
                                                           shape = 2.0,
                                                           meandist = 2.0km)),
                      longevity = row -> row.tp / 2, seed = 7)
    species = built()
    @test isequal(species.traits, built().traits)
    @test keys(species.tolerance) == (:temperature, :rainfall)
    @test mean(species.tolerance.temperature.dists[3]) == 285.0
    @test all(1.0 .<= std.(species.tolerance.temperature.dists) .<= 3.0)
    @test species.params.longevity == [0.5, 1.0, 1.5, 2.0]
    # Each species disperses by the mode its key chooses, with its own kernel.
    @test [nameof(typeof(m.mode)) for m in species.movement] ==
          [
        :WindDispersal,
        :InsectDispersal,
        :InsectDispersal,
        :VertebrateDispersal
    ]
    @test ismissing(species.traits.movement_meandist[1]) &&
          all(0.6km .<= species.traits.movement_meandist[2:3] .<= 2.4km)
    @test species.filled.tmin == [false, false, true, false]
    # A deactivated species starts with none; the others share 1000 each.
    @test species.abun[4] == 0 && sum(species.abun) == 3000
    # A key no choice names is refused, named.
    odd = SpeciesPool(DataFrame(species = ["a", "b"],
                                syndrome = ["wind", "ballistic"]))
    @test_throws r"`b` \(\"ballistic\"\)" quietly(odd,
                                                  movement = ChosenBy(:syndrome,
                                                                      "wind" =>
                                                                          MovementSpec(WindDispersal(),
                                                                                       meandist = 5.0km,
                                                                                       scale = 10.0km)))
end

@testset "a species list keeps its pool's identity, names and records" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        pool = SpeciesPool(DataFrame(species = [
                                         "Acacia alpha",
                                         "Mimosa alpha",
                                         "Acacia zeta"
                                     ],
                                     v = 1:3),
                           authority = wfo, merge = MergeWith(first),
                           verbosity = :silent)
        copy_species!(pool, "Acacia zeta", as = "zeta copy")
        species = quietly(pool)
        @test species.names == ["Acacia alpha", "Acacia zeta", "zeta copy"]
        @test isequal(species.traits.id,
                      ["wfo-0000000001", "wfo-0000000010", missing])
        @test species.traits.copykey[3] == species.traits.key[2]
        @test species.authority === wfo
        # A name the pool's data used finds its species: a synonym merged into another included.
        @test species.namemap["Mimosa alpha"] ==
              species.namemap["Acacia alpha"] == 1
        @test species.namemap["zeta copy"] == 3
        @test any(r -> r.doi == "10.5281/zenodo.20782718", species.inputs)
        # `makeunique` replaces the similarity and nothing else: every other field of the list, and
        # of the ecosystem, is the original's. Compared field by field over the type, so a field
        # added to either is checked without being named here; a drawn longevity gives `filled`
        # something to lose.
        drawn = quietly(pool, longevity = DrawFrom(Uniform(0.5, 1.5)),
                        seed = 1)
        @test drawn.authority === wfo && any(any, eachcol(drawn.filled))
        habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
        eco = build_ecosystem(drawn, habitat, seed = 1)
        stripped = makeunique(eco)
        @test stripped.spplist.types isa UniqueTypes
        # The fields of `a`, but `except`, holding something other than `b`'s.
        function changed(a, b, except)
            return filter(fieldnames(typeof(a))) do name
                return name != except &&
                       !isequal(getfield(a, name), getfield(b, name))
            end
        end
        @test isempty(changed(stripped.spplist, drawn, :types))
        @test isempty(changed(stripped, eco, :spplist))
    end
    # Species built otherwise have an identity of their own.
    plain = build_species(3, tolerance = (285.0K, 3.0K),
                          toleranceaxis = Temperature,
                          demand = 1.0e9kJ / day, demandaxis = SolarRadiation)
    @test plain.traits.key == [1, 2, 3] && all(plain.traits.active)
    @test plain.namemap == Dict("1" => 1, "2" => 2, "3" => 3)
end

@testset "species built from a pool run as the released form's do, and take an arrival" begin
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    # A demand the default habitat can carry, so the run has something to compare.
    tiny = 1.0e-4kJ / day
    niche = ToleranceSpec(Temperature, μ = 298.0K, σ = 2.0K)
    old = build_species(3, tolerance = (298.0K, 2.0K),
                        toleranceaxis = Temperature,
                        demand = tiny, demandaxis = SolarRadiation,
                        dispersal = 1.0km,
                        seed = 1)
    pool = SpeciesPool(DataFrame(species = ["1", "2", "3"],
                                 life = [0.5, 1.0, 2.0],
                                 need = [0.5, 1.0, 1.5] .* 1.0e-4))
    new = quietly(pool, tolerance = niche,
                  demand = DemandSpec(SolarRadiation, tiny),
                  seed = 1)
    runs = map((old, new)) do species
        eco = build_ecosystem(species, habitat, seed = 1)
        simulate!(eco, 6month_mean_duration, 1month_mean_duration)
        return eco.abundances.matrix
    end
    @test sum(last(runs)) > 0
    @test first(runs) == last(runs)
    # Each species' own longevity runs too, and changes the run where demands differ - with every
    # species' demand alike, ϵ̄ = 1 and no exponent moves it.
    needs = DemandSpec(SolarRadiation, TraitColumn(:need, unit = kJ / day))
    lives = map((1.0, TraitColumn(:life))) do longevity
        species = quietly(pool, tolerance = niche, demand = needs,
                          longevity = longevity, seed = 1)
        eco = build_ecosystem(species, habitat, seed = 1)
        simulate!(eco, 6month_mean_duration, 1month_mean_duration)
        return eco
    end
    eco = last(lives)
    @test sum(eco.abundances.matrix) > 0
    @test eco.abundances.matrix != first(lives).abundances.matrix
    simulate!(eco, 1month_mean_duration, 1month_mean_duration,
              intervention = Intervention(AtTime(7.0month_mean_duration),
                                          AllCells(),
                                          AddSpecies(abundance = 10,
                                                     name = "arrival")))
    @test eco.spplist.names[end] == "arrival"
    @test eco.spplist.traits.name[end] == "arrival"
    @test eco.spplist.traits.key[end] == 4 &&
          eco.spplist.namemap["arrival"] == 4
    @test eco.spplist.params.longevity == [0.5, 1.0, 2.0, 2.0]
end

end
