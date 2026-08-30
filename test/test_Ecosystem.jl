# SPDX-License-Identifier: LGPL-3.0-or-later

module TestEcosystem

using EcoSISTEM
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: getnichefit, update!
using EcoSISTEM: getregime
using Test
using Unitful, Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Diversity
using DimensionalData

include("TestCases.jl")

# A stand-in for any non-`UniqueTypes` similarity structure (a phylogeny, in practice), so the
# refusal can be checked without building a tree.
struct NotUniqueTypes <: Diversity.AbstractTypes end

@testset "Ecosystem" begin
    # NB not `@test_nowarn`: the fixture is still built through the deprecated `simplehabitat`,
    # so under `--depwarn=yes` (what `Pkg.test` runs with) it warns exactly once per session.
    eco = Test1Ecosystem()
    @test eco isa Ecosystem
    @test sum(eco.abundances.matrix, dims = 2)[:, 1] == eco.spplist.abun
    # tolerances, regime and nichefit line up — arity, names, unit and kind, member by member
    @test isnothing(EcoSISTEM._checkaligned("t",
                                            EcoSISTEM._toleranceside(eco.spplist.tolerance),
                                            EcoSISTEM._regimeside(eco.habitat.regime),
                                            EcoSISTEM._nichefitside(eco.nichefit)))
    @test isnothing(EcoSISTEM._checkaligned("d",
                                            EcoSISTEM._demandside(eco.spplist.demand),
                                            EcoSISTEM._supplyside(eco.habitat.supply)))

    sppl = SpeciesList{typeof(eco.spplist.tolerance),
                       typeof(eco.spplist.demand),
                       typeof(eco.spplist.movement), UniqueTypes,
                       typeof(eco.spplist.params)}(eco.spplist.names,
                                                   eco.spplist.tolerance,
                                                   eco.spplist.abun,
                                                   eco.spplist.demand,
                                                   UniqueTypes(length(eco.spplist.names)),
                                                   eco.spplist.movement,
                                                   eco.spplist.params,
                                                   eco.spplist.native)
    eco = Ecosystem(sppl, eco.habitat, eco.nichefit)
    # Adding a species mid-run is the `AddSpecies` operation, applied either by declaring an
    # intervention or — as here — through `applyinterventions!`, the supported imperative route.
    # `step` must be one the run has actually reached: the selection RNG is keyed on it.
    before = length(eco.spplist.names)
    @test_nowarn EcoSISTEM.applyinterventions!(eco,
                                               Intervention(EveryStep(),
                                                            AllCells(),
                                                            AddSpecies(abundance = 10)),
                                               0.0s, 1month_mean_duration, 1)
    @test length(eco.spplist.names) == before + 1

    @testset "get functions" begin
        # Test Simulation get functions
        @test_nowarn getnichefit(eco)
        @test getnichefit(eco) == eco.nichefit
        @test_nowarn getregime(eco)
        @test getregime(eco) == eco.habitat.regime
        @test_nowarn EcoSISTEM.getgridarea(eco)
        @test EcoSISTEM.getgridarea(eco) ≈
              size(eco.abundances.matrix, 2) .* eco.habitat.regime.size^2
        # `getgridshape` reports the grid's shape in cells. It replaces `getdimension`; the old
        # `getgridsize` gave one cell's side length and is deprecated onto that meaning, so the
        # released name keeps answering correctly rather than silently changing type.
        @test_nowarn EcoSISTEM.getgridshape(eco)
        @test EcoSISTEM.getgridshape(eco) == size(eco.habitat.regime.matrix)
        # **No `@test_nowarn` and no log assertion**, deliberately: these are deprecated, so they
        # warn — but only when `--depwarn=yes`, which `Pkg.test` sets and a direct `include` does
        # not. Asserting either way makes the test pass in one environment and fail in the other.
        # What matters is that the shims still *work*, so the values are what is pinned.
        @test getdispersaldist(eco, 1) == eco.spplist.movement.kernels[1].dist
        @test getdispersaldist(eco, "1") == eco.spplist.movement.kernels[1].dist
        @test getdispersalvar(eco, 1) ==
              (eco.spplist.movement.kernels[1].dist)^2 * pi / 4
        @test getdispersalvar(eco, "1") == getdispersalvar(eco, 1)

        # …and what supersedes them: the kernel itself, from either an ecosystem or the species
        # list `build_species` returned, by index or by name.
        @test_nowarn EcoSISTEM.speciesdispersal(eco, 1)
        @test EcoSISTEM.speciesdispersal(eco, 1) ===
              eco.spplist.movement.kernels[1]
        @test EcoSISTEM.speciesdispersal(eco, "1") ===
              EcoSISTEM.speciesdispersal(eco, 1)
        @test EcoSISTEM.speciesdispersal(eco.spplist, 1) ===
              EcoSISTEM.speciesdispersal(eco, 1)
        # The round trip the return type buys: it is exactly `AddSpecies`' `dispersal` keyword.
        @test AddSpecies(dispersal = EcoSISTEM.speciesdispersal(eco, 1),
                         abundance = 10).dispersal ===
              eco.spplist.movement.kernels[1]
        @test EcoSISTEM.getlookup(eco, "1") == eco.lookup[1]
        @test EcoSISTEM.getlookup(eco, 1) == eco.lookup[1]
        # `setchange!` installs a steady increment, checked against the regime's own unit per
        # second — here a dimensionless (land-cover) regime, so `s⁻¹`. A rate of the wrong
        # dimension for *this* layer is rejected, where the old `resetrate!` methods accepted any
        # `𝐓⁻¹` or `𝚯 𝐓⁻¹` rate on any regime whatsoever.
        @test_nowarn EcoSISTEM.setchange!(eco.habitat.regime,
                                          IncrementBy(0.1 / s))
        @test eco.habitat.regime.change ==
              EcoSISTEM.SteadyLayerChange{typeof(0.1 / s)}(0.1 / s)
        @test_throws ErrorException EcoSISTEM.setchange!(eco.habitat.regime,
                                                         IncrementBy(0.1K / s))
        # `resettime!` resets the simulation **clock**, which is what puts a stored series back to
        # its first slice. It must not reach for a cursor on the regime: a series is indexed by
        # elapsed time, so there is no cursor to reset and a regime holding none would have no
        # method at all.
        EcoSISTEM._advanceclock!(eco, 1month_mean_duration)
        @test EcoSISTEM.simulationtime(eco) > 0.0s
        resettime!(eco)
        @test EcoSISTEM.simulationtime(eco) == 0.0s
    end
    @testset "movement types" begin
        # Test other movement types
        eco = Test1Ecosystem()
        # No topology here since step 19 — it belongs to the habitat, not the movement.
        mov = AlwaysMovement(fill(LongTailKernel(10.0km, 10.0, 1e-10),
                                  length(eco.spplist.names)))
        sppl = SpeciesList{typeof(eco.spplist.tolerance),
                           typeof(eco.spplist.demand), typeof(mov),
                           typeof(eco.spplist.types),
                           typeof(eco.spplist.params)}(eco.spplist.names,
                                                       eco.spplist.tolerance,
                                                       eco.spplist.abun,
                                                       eco.spplist.demand,
                                                       eco.spplist.types, mov,
                                                       eco.spplist.params,
                                                       eco.spplist.native)
        @test_nowarn Ecosystem(sppl, eco.habitat, eco.nichefit)
    end
    @testset "diversity" begin
        # Test Diversity get functions
        @test_nowarn getmetaabundance(eco)
        @test_nowarn getpartition(eco)
        @test getpartition(eco) == eco.habitat
        @test_nowarn gettypes(eco)
        @test gettypes(eco) == eco.spplist
        @test_nowarn getordinariness!(eco)
        # **Recomputed every call, never cached.** Asserting that a *cache* was populated would be
        # asserting the hazard: nothing that mutates the ecosystem can invalidate one, so after
        # `_addspecies!` a cached matrix is the wrong shape. Measured: caching saved
        # 16,720 B on a 217 MB `norm_sub_alpha` call.
        @test getordinariness!(eco) == getordinariness!(eco)
    end
    @testset "collection-trait incompatibility message" begin
        # `iscontinuous` of a collection is a `Vector{Bool}`; the constructor's incompatibility errors must
        # format that rather than crash on a `Vector{Bool}` in a `?:` (the old bug threw a `TypeError`).
        @test EcoSISTEM._kindlabel(true) == "continuous"
        @test EcoSISTEM._kindlabel(false) == "categorical"
        @test EcoSISTEM._kindlabel([true, false]) == "[continuous, categorical]"

        # Build a valid two-variable (temperature + rainfall) collection environment + species, then pass a
        # mismatched single nichefit so the alignment check fails on the collection path.
        numSpecies = 4
        grid = (5, 5)
        area = 10000.0km^2
        # Two regimes *and* two supplies — `GridHabitat` takes a tuple on either side, so the
        # eltypes still differ ([K, mm/d]) without hand-assembling the collections.
        side = sqrt(area)
        habitat = GridHabitat(regime = (UniformSpec(10.0K,
                                                    axis = Temperature),
                                        UniformSpec(10.0mm / day,
                                                    axis = Precipitation)),
                              supply = (UniformSpec(1000.0kJ / km^2 / day,
                                                    axis = SolarRadiation),
                                        UniformSpec(100.0mm / day,
                                                    axis = Precipitation)),
                              area = StudyArea(extent = (side, side),
                                               cellsize = side / grid[1],
                                               verbosity = :silent))
        tolerance = SpeciesRequirementCollection((NicheTolerance(Temperature,
                                                                 Normal,
                                                                 fill(10.0K,
                                                                      numSpecies),
                                                                 fill(0.1K,
                                                                      numSpecies)),
                                                  NicheTolerance(Precipitation,
                                                                 Uniform,
                                                                 fill(1.0mm /
                                                                      day,
                                                                      numSpecies),
                                                                 fill(5.0mm /
                                                                      day,
                                                                      numSpecies))))
        abun = rand(Multinomial(1000, numSpecies))
        movement = BirthOnlyMovement(GaussianKernel.(fill(1.0km, numSpecies),
                                                     10e-4))
        native = fill(true, numSpecies)
        resource = SpeciesRequirementCollection((Demand{SolarRadiation}(fill(2.0kJ /
                                                                             day,
                                                                             numSpecies)),
                                                 Demand{Precipitation}(fill(2.0Unitful.L /
                                                                            day,
                                                                            numSpecies))))
        param = EqualPop(0.6 / month_mean_duration, 0.6 / month_mean_duration,
                         1.0, 0.0, 1000.0)
        sppl = SpeciesList(numSpecies, tolerance, abun, resource, movement,
                           param,
                           native)
        # the matching nichefit is `MultiplicativeFit((...))`; a single `NicheSuitability` mismatches
        badrel = NicheSuitability{Temperature, typeof(1.0K)}()
        err = try
            Ecosystem(sppl, habitat, badrel)
            nothing
        catch e
            e
        end
        @test err isa ErrorException                       # not a TypeError from the message itself
        # the merged check names both sides and the counts, rather than only "incompatible"
        @test occursin("species tolerance", err.msg)
        @test occursin("trait nichefit", err.msg)
        @test !occursin("non-boolean", err.msg)

        # a per-member disagreement names the offending *layer*, which comparing two `Vector`s
        # with `==` (the old `tematch`/`nfmatch`) never could
        swapped = MultiplicativeFit((NicheSuitability{Precipitation,
                                                      typeof(1.0mm / day)}(),
                                     NicheSuitability{Temperature,
                                                      typeof(1.0K)}()))
        swaperr = try
            Ecosystem(sppl, habitat, swapped)
            nothing
        catch e
            e
        end
        @test swaperr isa ErrorException
        # **The offending layers are named by their AXES**, not by their positions — that is D25
        # earning its keep on the error path: the message says *which measurements* disagree rather
        # than which slots, and the reader does not have to count members to find out.
        # **D3(a) moved which check catches this.** A nichefit now carries its own axis, so the
        # swapped fit derives `(:Precipitation, :Temperature)` against the tolerance's
        # `(:Temperature, :Precipitation)` and the **naming** check fires first — where before, the
        # fit had no axis, both sides fell back to positional names, and it reached the per-member
        # check. Caught earlier and by a better question; the axes are still what is named.
        @test occursin("Temperature", swaperr.msg) &&
              occursin("Precipitation", swaperr.msg)
    end

    @testset "GridLandscape takes its grid from a StudyGrid" begin
        # A landscape is positioned by the `StudyGrid` it is built on, which is what lets its flat
        # `location` dimension name each cell by extent. The `(Y, X)` order this used to guard by
        # `MethodError` is now enforced one level up, by `StudyGrid` itself, so it cannot be got
        # wrong here at all.
        area = StudyArea(extent = (2.0km, 3.0km), cellsize = 1.0km,
                         verbosity = :silent)
        habitat = GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                              supply = UniformSpec(1.0e11kJ / km^2 / day,
                                                   axis = SolarRadiation),
                              area = area, topology = Torus())
        grid = EcoSISTEM.getcoords(habitat)
        abun = reshape(collect(1:24), 4, 6)
        names = string.(1:4)
        gl = GridLandscape(abun, names, grid)
        @test gl isa GridLandscape
        # the flat view now says where each cell is, not merely which index it has
        @test occursin("km",
                       string(DimensionalData.lookup(gl.dimmatrix, :location)[1]))
        # and a grid of the wrong size is still refused
        @test_throws DimensionMismatch GridLandscape(reshape(collect(1:20), 4,
                                                             5),
                                                     names, grid)
    end
end

@testset "_addtypes! keeps the names, and refuses what it cannot extend" begin
    # Real names are essential: `build_species` yields `"1"`, `"2"`, … , so the bug this guards —
    # `UniqueTypes(count)` renaming **every existing species** to its index on each `_addspecies!` —
    # cannot show itself on a synthetic list.
    named = ["oak", "ash", "elm"]
    types = EcoSISTEM._uniquetypes(named)
    grown = EcoSISTEM._addtypes!(types, "yew")
    @test gettypenames(grown, true) == [named; "yew"]
    @test counttypes(grown, true) == 4

    # And a species list whose types are not `UniqueTypes` — a phylogeny, in practice — is refused
    # with the reason rather than falling off the end of dispatch as a `MethodError`.
    @test_throws ErrorException EcoSISTEM._addtypes!(NotUniqueTypes(), "yew")
end

# --- `CachedEcosystem`, which is declared in `Ecosystem.jl` ---

# Build a small, fully deterministic ecosystem. With a fixed `seed` the initial
# abundances (a fixed vector) and the per-species RNG streams are identical every
# time, so runs are byte-for-byte reproducible.
function cache_test_eco(seed)
    numSpecies = 8
    grid = (4, 4)
    area = 100.0km^2
    individuals = 1_000
    totalK = 10000.0kJ / km^2 / day

    resource = Demand{SolarRadiation}(fill(10.0kJ / day, numSpecies))
    param = EqualPop(0.2 / year, 0.2 / year, 1.0, 0.0, 1.0)
    kernel = fill(GaussianKernel(2.0km, 1.0e-3), numSpecies)
    movement = BirthOnlyMovement(kernel)
    tolerance = NicheTolerance(Temperature, Normal,
                               fill(274.0K, numSpecies),
                               fill(0.5K, numSpecies))
    native = fill(true, numSpecies)
    abun = fill(div(individuals, numSpecies), numSpecies)

    sppl = SpeciesList(numSpecies, tolerance, abun, resource, movement, param,
                       native)
    # The grid is decided first, then built on: 100 km² over 4 × 4 is 2.5 km cells — the same
    # `sqrt(area / prod(dimension))` the deprecated builder computed internally.
    studyarea = StudyArea(extent = (sqrt(area), sqrt(area)),
                          cellsize = sqrt(area) / grid[1], verbosity = :silent)
    habitat = GridHabitat(regime = UniformSpec(274.0K,
                                               axis = Temperature),
                          supply = UniformSpec(totalK,
                                               axis = SolarRadiation),
                          area = studyarea)
    nichefit = NicheSuitability{Temperature, typeof(1.0K)}()
    return Ecosystem(sppl, habitat, nichefit, seed = seed)
end

@testset "CachedEcosystem save/load/resume" begin
    seed = 20250708
    times = (0.0year):(1.0year):(5.0year)

    # Reference run: never cached, five one-year steps
    ref = cache_test_eco(seed)
    for _ in 1:5
        update!(ref, 1.0year)
    end
    expected = copy(ref.abundances.matrix)

    # Cached run, only simulated (and saved to disk) as far as year 3
    dir = mktempdir()
    eco1 = cache_test_eco(seed)
    cache1 = CachedEcosystem(eco1, dir, times)
    abundances(cache1, 3.0year)
    @test !isempty(filter(f -> endswith(f, ".jld2"), readdir(dir)))

    # Fresh cache over the same folder: its only in-memory state is year 0, so
    # asking for year 5 must load the year-3 snapshot from disk (restoring the
    # per-species RNG streams) and simulate the remaining two years.
    eco2 = cache_test_eco(seed)
    cache2 = CachedEcosystem(eco2, dir, times)
    resumed = abundances(cache2, 5.0year)

    # Resuming from the cache reproduces the uncached run exactly
    @test resumed.matrix == expected

    # clearcache! removes the saved files — the `!` is owed, it deletes them
    @test_nowarn EcoSISTEM.clearcache!(cache2)
    @test isempty(filter(f -> endswith(f, ".jld2"), readdir(dir)))
end

@testset "CachedEcosystem matches Ecosystem" begin
    seed = 424242
    times = (0.0year):(1.0year):(5.0year)

    # Plain ecosystem, five one-year steps
    eco = cache_test_eco(seed)
    for _ in 1:5
        update!(eco, 1.0year)
    end

    # Same configuration driven as a CachedEcosystem, run to year 5. Its results
    # come from its own in-memory cache; the disk-reload path is covered above.
    cached = CachedEcosystem(cache_test_eco(seed), mktempdir(), times)
    result = abundances(cached, 5.0year)

    @test result.matrix == eco.abundances.matrix
end

@testset "CachedEcosystem caching interval independence" begin
    seed = 987654
    times = (0.0year):(1.0year):(5.0year)

    # Two caches, same seed and timestep but different disk-save cadence
    dirA = mktempdir()
    dirB = mktempdir()
    cacheA = CachedEcosystem(cache_test_eco(seed), dirA, times,
                             saveinterval = 1.0year)
    cacheB = CachedEcosystem(cache_test_eco(seed), dirB, times,
                             saveinterval = 2.0year)

    rA = abundances(cacheA, 4.0year)
    rB = abundances(cacheB, 4.0year)

    # The caching interval must not affect the simulation result
    @test rA.matrix == rB.matrix

    # ...but the two caches genuinely saved at different frequencies
    countfiles(d) = length(filter(f -> endswith(f, ".jld2"), readdir(d)))
    @test countfiles(dirA) != countfiles(dirB)
end

# **The cached ecosystem must answer diversity questions like any other**, and this is the only
# testset that asks it one: every other compares abundances or counts files. That gap is how
# `_getabundance(::CachedEcosystem, …)` can skip `_calcabundance` — the hook `SpeciesList` forwards
# to `sl.types` through — without anything noticing.
#
# **A `UniqueTypes` fixture cannot see it.** `_calcabundance` returns a `UniqueTypes` array
# unchanged, so the two paths agree exactly and the bug is invisible; it bites only where the types
# carry structure. `Test1Ecosystem` builds a phylogeny, so the non-raw abundances are over
# **branches**, of which there are more than there are species.
@testset "CachedEcosystem answers diversity like an Ecosystem" begin
    eco = Test1Ecosystem(seed = 1)

    # `_gettypes` answers the `SpeciesList` itself, not the phylogeny inside it: `SpeciesList <:
    # AbstractTypes` but *holds* the real types in `sl.types` and forwards every hook there. That
    # forwarding is the whole subject of this testset, so it is asserted rather than assumed.
    @test Diversity.API._gettypes(eco) === eco.spplist

    # The non-vacuity guard: with a `UniqueTypes` list every assertion below is trivially true.
    # Matched by name because `PhyloBranches` lives in `DiversityPhyloExt`, an extension of
    # Diversity — there is no `Diversity.PhyloBranches` to name.
    @test nameof(typeof(eco.spplist.types)) === :PhyloBranches

    times = (0.0month_mean_duration):(1.0month_mean_duration):(2.0month_mean_duration)
    cache = CachedEcosystem(eco, mktempdir(), times)

    raweco = Diversity.API._getabundance(eco, true)
    rawcache = Diversity.API._getabundance(cache, true)
    @test rawcache == raweco

    proceco = Diversity.API._getabundance(eco, false)
    proccache = Diversity.API._getabundance(cache, false)
    @test proccache ≈ proceco

    # The assertion the fix is really about: processing maps species onto branches, so the answer is
    # a *different shape* from the raw abundances. Returning the raw ones — the old behaviour —
    # would make this test fail on size alone.
    @test size(proccache, 1) > size(rawcache, 1)
end

end
