# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDynamics

using EcoSISTEM
# `update!` is `public` rather than exported, so it must be named.
using EcoSISTEM: update!
using Test
using Unitful
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using JLD2
using DataFrames: DataFrame

include("TestCases.jl")

@testset "Update functions" begin
    eco = Test1Ecosystem()
    @test_nowarn update!(eco, 1month_mean_duration)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.habitat.topology,
                                              1, 1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.habitat.topology,
                                              1, 1, 1, eco, 10)) ==
          Vector{Int64}
    @test_nowarn repopulate!(eco, verbosity = :silent)

    # Test Cylinder
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Cylinder(), 1, 1, 1, eco, 10)
    # Test Torus
    @test_nowarn EcoSISTEM.calc_lookup_moves!(Torus(), 1, 1, 1, eco, 10)

    # Test ecosystem with multiple supplies
    eco = TestMultiEcosystem()
    @test_nowarn update!(eco, 1month_mean_duration)
    @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.habitat.topology,
                                              1, 1, 1, eco, 10)
    @test typeof(EcoSISTEM.calc_lookup_moves!(eco.habitat.topology,
                                              1, 1, 1, eco, 10)) ==
          Vector{Int64}
    @test_nowarn repopulate!(eco, verbosity = :silent)
end

@testset "Dispersal from a cell with no reachable destination" begin
    # If every destination a species' kernel can reach is inactive, the dispersal weights sum to
    # zero. Normalising that gives `NaN` probabilities and a throw from `Multinomial`, so the case
    # is instead defined as "no move this step": `moves` is zeroed and the species' RNG stream is
    # left untouched (a stranded cell must not consume draws and re-phase every later one).
    #
    # `update!` cannot reach this today - it only disperses from an *active* cell, and a cell is
    # always its own kernel's highest-weighted destination - so the guard is exercised directly.
    eco = Test1Ecosystem(seed = 4)
    eco.habitat.active .= false
    lookup = EcoSISTEM.getlookup(eco, 1)
    before = copy(eco.rngs[1])
    for bound in (Island(), Cylinder(), Torus())
        fill!(lookup.moves, 99)
        @test_nowarn EcoSISTEM.calc_lookup_moves!(bound, 5, 5, 1, eco, 100)
        @test all(iszero, lookup.moves)
        @test !any(isnan, lookup.pnew)
    end
    @test eco.rngs[1] == before

    # An ordinary draw is unaffected: all `abun` individuals are still placed.
    eco = Test1Ecosystem(seed = 4)
    lookup = EcoSISTEM.getlookup(eco, 1)
    EcoSISTEM.calc_lookup_moves!(Torus(), 5, 5, 1, eco, 100)
    @test sum(lookup.moves) == 100
end

@testset "Simulation clock" begin
    # `update!` advances the ecosystem's own clock by one timestep. Nothing reads it yet - it
    # exists for the layer-change and intervention work built on top of it.
    eco = Test1Ecosystem(seed = 5)
    @test EcoSISTEM.simulationtime(eco) == 0.0s
    update!(eco, 1month_mean_duration)
    @test EcoSISTEM.simulationtime(eco) == uconvert(s, 1.0month_mean_duration)
    update!(eco, 1month_mean_duration)
    @test EcoSISTEM.simulationtime(eco) == uconvert(s, 2.0month_mean_duration)
    # The seed is retained rather than consumed and discarded.
    @test eco.seed === UInt64(5)
end

@testset "NoGrowth freezes with one or two supplies" begin
    # NoGrowth must zero the birth/death adjustment regardless of how many resource
    # supplies the environment has: the multi-supply path must not skip it.
    N = 8
    grid = (4, 4)
    area = 100.0km^2
    birth = 0.1 / month_mean_duration
    nogrowth = NoGrowth{typeof(unit(birth))}(fill(birth, N), fill(birth, N),
                                             1.0, 0.0)
    tolerance = NicheTolerance(Temperature, Normal, fill(274.0K, N),
                               fill(0.5K, N))
    movement = fill(InsectDispersal()(meandist = 1.0km,
                                      thresh = 1.0e-3),
                    N)
    nichefit = NicheSuitability{Temperature, typeof(1.0K)}()

    # single supply
    sppl1 = SpeciesList(N, tolerance,
                        Demand{SolarRadiation}(fill(10.0kJ / day, N)),
                        movement, nogrowth)
    # Grid decided first: 100 km^2 over 4 × 4 is 2.5 km cells.
    habitat1 = GridHabitat(regime = UniformSpec(274.0K,
                                                axis = Temperature),
                           supply = UniformSpec(10000.0kJ / km^2 / day,
                                                axis = SolarRadiation),
                           area = StudyArea(extent = (sqrt(area),
                                                      sqrt(area)),
                                            cellsize = sqrt(area) /
                                                       grid[1],
                                            verbosity = :silent))
    eco1 = build_ecosystem(BySupply(each = 10), sppl1, habitat1,
                           nichefit = nichefit, distributed = false)
    @test EcoSISTEM.resource_adjustment(eco1, eco1.habitat.supply, 1, 1) ==
          (0.0, 0.0)

    # two supplies - the path that must not skip the adjustment
    resource2 = SpeciesRequirementCollection((Demand{SolarRadiation}(fill(10.0kJ /
                                                                          day,
                                                                          N)),
                                              Demand{Precipitation}(fill(2.0Unitful.L /
                                                                         day,
                                                                         N))))
    sppl2 = SpeciesList(N, tolerance, resource2, movement, nogrowth)
    side = sqrt(area)
    habitat2 = GridHabitat(regime = UniformSpec(274.0K,
                                                axis = Temperature),
                           supply = (UniformSpec(10000.0kJ / km^2 / day,
                                                 axis = SolarRadiation),
                                     UniformSpec(10.0mm / day,
                                                 axis = Precipitation)),
                           area = StudyArea(extent = (side, side),
                                            cellsize = side / grid[1],
                                            verbosity = :silent))
    eco2 = build_ecosystem(BySupply(each = 10), sppl2, habitat2,
                           nichefit = nichefit, distributed = false)
    @test EcoSISTEM.resource_adjustment(eco2, eco2.habitat.supply, 1, 1) ==
          (0.0, 0.0)
end

# Build an ecosystem whose species have *differing* per-capita demands, on a supply that can never
# limit anything when a second resource is added. The whole shipped corpus uses *scalar* demands,
# so ϵ̄ = 1 for every species and none of the properties below is observable in it - which is exactly
# why a multi-resource arithmetic bug survived for years unnoticed.
#
# `demand` is the three species' demands on one resource, or a tuple of them, one per supply, on the
# axes `demandaxis` names; each is a column of the species' pool.
function _demandeco(demand, supply, demandaxis = SolarRadiation;
                    tolerance = ToleranceSpec(Temperature, μ = 280.0K,
                                              σ = 5.0K))
    area = StudyArea(extent = (3.0km, 5.0km), cellsize = 1.0km,
                     verbosity = :silent)
    env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                      supply = supply, area = area)
    members = demand isa Tuple ? demand : (demand,)
    resourceaxes = demandaxis isa Tuple ? demandaxis : (demandaxis,)
    columns = ntuple(i -> Symbol(:demand, i), length(members))
    table = DataFrame(species = string.(1:3))
    for (column, member) in zip(columns, members)
        table[!, column] = member
    end
    specs = map(resourceaxes, columns) do axis, column
        return DemandSpec(axis, TraitColumn(column))
    end
    spp = build_species(SpeciesPool(table, verbosity = :silent),
                        tolerance = tolerance,
                        demand = demand isa Tuple ? specs : only(specs),
                        seed = 1, verbosity = :silent)
    eco = build_ecosystem(BySupply(total = 300), spp, env, seed = 1)
    eco.abundances.matrix .= 10
    EcoSISTEM.update_resource_usage!(eco)
    return eco
end

const _SOLAR = UniformSpec(1.0e5kJ / (m^2 * day), axis = SolarRadiation)
# So abundant it can never be the limiting resource, whatever any species demands of it.
const _SPAREWATER = UniformSpec(1.0e12Unitful.L / (m^2 * day),
                                axis = Precipitation)

@testset "demand varying across species drives the demographic terms" begin
    solar = [1.0e9, 2.0e9, 3.0e9] .* kJ / day        # ϵ̄ = 0.5, 1.0, 1.5
    one = _demandeco(solar, _SOLAR)
    adj(eco, sp) = EcoSISTEM._resourceadjustment(eco, eco.habitat.supply, 1,
                                                 sp)

    # The demand term carries the SAME exponent sign on birth and death, so it cancels from their
    # ratio: demand sets the *tempo* of turnover, not the equilibrium. Only suitability (opposite
    # signs) moves the balance. Pinned because it is what tells you which term to look at when a
    # change moves the wrong thing.
    ratios = [((b, d) = adj(one, sp); b / d)
              for sp in 1:3]
    @test all(≈(first(ratios)), ratios)
    # ...and the rates themselves genuinely differ, or the test above would be vacuous.
    @test !≈(first(adj(one, 1)), first(adj(one, 3)))

    # REGRESSION GUARD. A second resource that can never limit anything, and that every species
    # demands equally (ϵ̄ = 1, the product's fixed point), must leave birth and death *exactly*
    # unchanged. Applying the suitability term once per resource would move birth by a factor of
    # 0.5457 and death by 1.8325 here - `suit^∓survival` - for a resource that cannot possibly
    # matter.
    inert = _demandeco((solar, fill(1.0Unitful.L / day, 3)),
                       (_SOLAR, _SPAREWATER),
                       (SolarRadiation, Precipitation))
    for sp in 1:3
        @test adj(inert, sp) == adj(one, sp)
    end

    # ...but two *genuine* demands compound, which is required: a species needing 2× solar and 2×
    # water really is more demanding than one needing 2× solar alone. Here the water demands are
    # proportional to the solar ones, so ϵ̄ is identical on both resources and `demanded` is ϵ̄^2.
    twice = _demandeco((solar, [1.0, 2.0, 3.0] .* Unitful.L / day),
                       (_SOLAR, _SPAREWATER),
                       (SolarRadiation, Precipitation))
    longevity = one.spplist.params.longevity
    for sp in 1:3
        ϵ = only(EcoSISTEM.values(one.spplist.demand)).resource[sp] *
            only(EcoSISTEM.values(one.spplist.demand)).exchange_rate
        # birth and death both carry `demanded^-longevity`, so both scale by ϵ̄^-longevity again.
        @test all(adj(twice, sp) .≈ adj(one, sp) .* ϵ^-longevity[sp])
    end
    # The average species is the fixed point of that compounding and is *not* a witness to it -
    # ϵ̄ = 1 there, so it alone would pass whatever the combining rule.
    @test adj(twice, 2) == adj(one, 2)
end

# Each species' own longevity is the exponent on its own demand, on birth and death alike, so changing
# one species' longevity scales that species' two terms by its demand to the change, leaves their
# ratio - the equilibrium - as it was, and touches no other species; one ecosystem holding several
# exponents does so for every species at once, single resource or several.
@testset "each species' longevity sets its own tempo" begin
    solar = [1.0e9, 2.0e9, 3.0e9] .* kJ / day        # ϵ̄ = 0.5, 1.0, 1.5
    adj(eco, sp) = EcoSISTEM._resourceadjustment(eco, eco.habitat.supply, 1,
                                                 sp)
    for (demand, supply, axes) in ((solar, _SOLAR, SolarRadiation),
                                   ((solar, fill(1.0Unitful.L / day, 3)),
                                    (_SOLAR, _SPAREWATER),
                                    (SolarRadiation, Precipitation)))
        even = _demandeco(demand, supply, axes)
        varied = _demandeco(demand, supply, axes)
        varied.spplist.params.longevity .= [0.5, 1.0, 2.0]
        for sp in 1:3
            ϵ = solar[sp] / mean(solar)
            factor = ϵ^-(varied.spplist.params.longevity[sp] - 1.0)
            @test all(adj(varied, sp) .≈ adj(even, sp) .* factor)
            (b, d) = adj(varied, sp)
            (b0, d0) = adj(even, sp)
            @test b / d ≈ b0 / d0
        end
        # ϵ̄ = 1 is the fixed point: the average species' terms do not move at all.
        @test adj(varied, 2) == adj(even, 2)
    end
end

# A supply that never runs short leaves a whole run unchanged, and a stock does so as a flow does. A
# soil-water volume every species demands equally has a demand term of exactly 1 and is never the
# scarcest resource, so neither `min(K/E)` nor `max(E/K)` picks it and every rate is bit-identical;
# the abundances after two years are compared exactly, since a last-bit difference in any rate would
# change a draw. The scarce control shows the comparison can see a difference.
@testset "a supply that never runs short changes nothing in a run" begin
    solar = [1.0e9, 2.0e9, 3.0e9] .* kJ / day
    water = fill(1.0m^3, 3)
    spare = UniformSpec(1.0e12mm, axis = SoilWaterVolume)
    scarce = UniformSpec(0.01mm, axis = SoilWaterVolume)
    function finalabundances(demand, supply, demandaxis)
        built = _demandeco(demand, supply, demandaxis)
        simulate!(built, 2year, 1month_mean_duration)
        return copy(built.abundances.matrix)
    end
    # The rates themselves, exactly: a resource that was not the scarcest but still moved a rate
    # would do so by about one part in 10^13 here, too little to change a draw in two years, so the
    # run comparison below cannot see it on its own.
    one = _demandeco(solar, _SOLAR, SolarRadiation)
    withspare = _demandeco((solar, water), (_SOLAR, spare),
                           (SolarRadiation, SoilWaterVolume))
    rates(eco, sc, sp) = EcoSISTEM._resourceadjustment(eco, eco.habitat.supply,
                                                       sc, sp)
    @test all(rates(withspare, sc, sp) == rates(one, sc, sp)
              for sc in axes(one.abundances.matrix, 2), sp in 1:3)

    alone = finalabundances(solar, _SOLAR, SolarRadiation)
    @test sum(alone) > 0
    @test finalabundances((solar, water), (_SOLAR, spare),
                          (SolarRadiation, SoilWaterVolume)) == alone
    @test finalabundances((solar, water), (_SOLAR, scarce),
                          (SolarRadiation, SoilWaterVolume)) != alone
end

@testset "non-square grid: y/x dimension order regression" begin
    # A square grid cannot detect an X/Y mixup (sizes and indices coincide by
    # accident). This grid is deliberately asymmetric, so confused code either
    # throws a BoundsError (indexing a 2-row array with a column-range value up
    # to 6) or silently reads the wrong cells.
    N = 4
    ny, nx = 2, 6
    grid = (ny, nx)
    area = 100.0km^2

    # A partial active mask - only column x=1 is active, every row - to check
    # the active-mask indexing itself corresponds to (y, x), not (x, y):
    # deliberately asymmetric (not a full row/column/diagonal) so a swap would
    # show up as reading the wrong cells, not coincidentally the same ones.
    partial_active = fill(false, ny, nx)
    partial_active[:, 1] .= true
    # **`within` takes the mask directly**, so a deliberately asymmetric one is still
    # expressible without the deprecated builder's positional `active` argument - which matters
    # here, because a disc or any other shape spec could not make *only column 1* active.
    partial_habitat = GridHabitat(regime = UniformSpec(274.0K,
                                                       axis = Temperature),
                                  supply = UniformSpec(10000.0kJ / km^2 /
                                                       day,
                                                       axis = SolarRadiation),
                                  area = StudyArea(extent = (ny * 1.0km,
                                                             nx * 1.0km),
                                                   cellsize = 1.0km,
                                                   within = partial_active,
                                                   verbosity = :silent))
    @test size(partial_habitat.regime.matrix) == grid
    @test size(partial_habitat.active) == grid
    @test all(partial_habitat.active[y, 1] for y in 1:ny)
    @test all(!partial_habitat.active[y, x] for y in 1:ny, x in 2:nx)

    tolerance = NicheTolerance(Temperature, Normal, fill(274.0K, N),
                               fill(0.5K, N))
    param = EqualPop(0.6 / month_mean_duration, 0.6 / month_mean_duration, 1.0,
                     0.0)
    kernel = fill(GaussianKernel(meandist = 1.0km, thresh = 1.0e-3), N)
    movement = InsectDispersal().(kernel)
    resource = Demand{SolarRadiation}(fill(2.0kJ / day, N))
    sppl = SpeciesList(N, tolerance, resource, movement, param)

    # A fully-active habitat for exercising calc_lookup_moves!/update! at
    # extreme y/x coordinates - every cell reachable, so no cell's move
    # probabilities can go to zero regardless of dimension order.
    habitat = GridHabitat(regime = UniformSpec(274.0K,
                                               axis = Temperature),
                          supply = UniformSpec(10000.0kJ / km^2 / day,
                                               axis = SolarRadiation),
                          area = StudyArea(extent = (ny * 1.0km,
                                                     nx * 1.0km),
                                           cellsize = 1.0km,
                                           verbosity = :silent))
    @test size(habitat.regime.matrix) == grid
    @test size(habitat.active) == grid

    nichefit = NicheSuitability{EcoSISTEM.axisof(habitat.regime),
                                eltype(habitat.regime)}()
    eco = build_ecosystem(BySupply(each = 10), sppl, habitat,
                          nichefit = nichefit, distributed = false)

    # getgridshape must report (ny, nx) - the actual array shape - not (nx, ny).
    @test EcoSISTEM.getgridshape(eco) == grid

    # Calling calc_lookup_moves! at every combination of extreme y/x coordinates
    # must not throw - a BoundsError here is exactly what an X/Y mixup produces
    # on a non-square grid (e.g. checking a y=2-valid offset against nx=6's
    # bound instead of ny=2's, then indexing a nonexistent row of `active`).
    for (y, x) in ((1, 1), (ny, 1), (1, nx), (ny, nx))
        @test_nowarn EcoSISTEM.calc_lookup_moves!(eco.habitat.topology,
                                                  y, x, 1, eco, 5)
        @test_nowarn EcoSISTEM.calc_lookup_moves!(Cylinder(), y, x, 1, eco, 5)
        @test_nowarn EcoSISTEM.calc_lookup_moves!(Torus(), y, x, 1, eco, 5)
    end

    # convert_coords must round-trip correctly through the (y,x) convention -
    # height (dimension 1) is ny, not nx.
    for i in 1:(ny * nx)
        (y, x) = EcoSISTEM.convert_coords(eco, i, ny)
        @test 1 <= y <= ny && 1 <= x <= nx
        @test EcoSISTEM.convert_coords(y, x, ny) == i
    end

    # A full timestep must run without error on this asymmetric grid.
    @test_nowarn update!(eco, 1month_mean_duration)
end

@testset "Multithreaded reproducibility" begin
    # A seeded run must give identical results regardless of the number of
    # threads. Run the same seeded simulation in child processes forced to use
    # 1, 2 and 4 threads (so the parallel `update!` path is genuinely exercised)
    # and check that all three final abundance matrices are identical. Each child
    # runs the simulation twice, which confirms same-thread-count repeatability,
    # so a failure distinguishes "not reproducible at all" from "not reproducible
    # across thread counts".
    #
    # **The repeat is inside the children rather than a fourth child**, because a
    # child spends about 14 of its 15 seconds starting Julia and compiling the
    # fixture build, and 0.3 s simulating: a fourth process bought the same
    # assurance at a whole process's price.
    script = pkgdir(EcoSISTEM, "test", "threading_reproducibility.jl")
    dir = mktempdir()
    function run_with(nthreads, tag)
        out = joinpath(dir, "repro_$tag.jld2")
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t $nthreads $script $out`
        @test success(pipeline(cmd, stdout = stdout, stderr = stderr))
        result = load(out)
        @test result["matrix"] == result["again"]
        return result["matrix"]
    end
    m1 = run_with(1, "t1")
    m2 = run_with(2, "t2")
    m4 = run_with(4, "t4")
    @test m1 == m2 == m4
end

# **`_getneighbours` takes `(y, x)` - row first**, as everything else in the package does, and
# its parameter *names* must say so: a caller that believes a name over the body passes the column
# first. On a **square** matrix that is invisible, so every assertion below is on a non-square one.
@testset "_getneighbours takes (y, x), row first" begin
    M = zeros(5, 3)         # 5 rows, 3 columns - a transposed call goes out of bounds

    # The four-neighbourhood of an interior cell, as (row, column) pairs.
    @test sort(collect(eachrow(EcoSISTEM._getneighbours(M, 3, 2)))) ==
          sort([[2, 2], [4, 2], [3, 1], [3, 3]])
    # The asymmetry is the whole point: row 5 exists and column 5 does not.
    @test EcoSISTEM._getneighbours(M, 5, 3) isa Matrix{Int}
    @test_throws ErrorException EcoSISTEM._getneighbours(M, 3, 5)

    # Cells off the edge are dropped, so a corner has two orthogonal neighbours and three diagonal.
    @test size(EcoSISTEM._getneighbours(M, 1, 1), 1) == 2
    @test size(EcoSISTEM._getneighbours(M, 1, 1, 8), 1) == 3
    # ...and every pair returned indexes `M` - which is what a transposed call would break.
    for chess in (4, 8), y in Base.axes(M, 1), x in Base.axes(M, 2)
        ns = EcoSISTEM._getneighbours(M, y, x, chess)
        @test all(n -> checkbounds(Bool, M, n[1], n[2]), eachrow(ns))
        # A neighbour is one step away in at least one direction and never the cell itself.
        @test all(n -> 0 < abs(n[1] - y) + abs(n[2] - x) <= 2, eachrow(ns))
    end

    # The vector form pairs them up positionally, in the same order.
    @test EcoSISTEM._getneighbours(M, [1, 5], [1, 3]) ==
          vcat(EcoSISTEM._getneighbours(M, 1, 1),
               EcoSISTEM._getneighbours(M, 5, 3))
    @test_throws ErrorException EcoSISTEM._getneighbours(M, 1, 1, 5)
end

# -- The hot loop must not allocate per cell -----------------------------------------------------
#
# **Why these exist, and why there are two.** `Ecosystem` declares `abundances::GridLandscape`,
# the bare name of a parametric type, so the field is abstract. Were a field the loop reads through
# it to name a type parameter, every `eco.abundances.matrix` reached inside a threaded body would be
# boxed, at about 176 bytes per cell, growing with the grid.
#
# The cause would be in `Landscape.jl` and the damage in `dynamics.jl`, whose text need not change,
# so reading the loop finds nothing.
#
# The two checks catch different halves and neither subsumes the other: the first names the boxing
# directly and fails the moment a type on the path gains a parameter, the second catches any
# per-cell allocation whatever its cause.

# Accessors at module scope rather than closures, so inference is asked about a stable signature.
_readmatrix(eco) = eco.abundances.matrix
_readgrid(eco) = eco.abundances.grid
_readtotaldemand(eco) = eco.cache.totaldemand
_readnetmigration(eco) = eco.cache.netmigration
_readtopology(eco) = eco.habitat.topology
_readmovement(eco) = eco.spplist.movement

# Behind a function barrier, as `test_collections.jl`'s allocation checks are: `eco` arrives as a
# typed argument, so the test's own necessarily-unstable local cannot land on the measurement, and
# the warming call compiles the same specialisation the measured one runs.
function _usagealloc(eco)
    eco.cache.valid = false
    EcoSISTEM.update_resource_usage!(eco)
    eco.cache.valid = false
    return @allocated EcoSISTEM.update_resource_usage!(eco)
end

@testset "the hot loop's reads infer concretely" begin
    eco = Test1Ecosystem(seed = 1)
    T = typeof(eco)
    # **Deliberately NOT `isconcretetype(fieldtype(T, :abundances))`.** That field is still the bare
    # `GridLandscape` `UnionAll` and is *meant* to be: the landscape holds concrete raw arrays
    # inside an abstract container, and Julia infers a field through it so long as the field's own
    # declared type names no type parameter. A field-concreteness walk would fail on exactly this
    # field, and would have to carve it out -- worse than no check.
    @test !isconcretetype(fieldtype(T, :abundances))
    @test Base.return_types(_readmatrix, (T,))[1] === Matrix{Int64}
    @test Base.return_types(_readgrid, (T,))[1] === Array{Int64, 3}
    @test Base.return_types(_readtotaldemand, (T,))[1] === Matrix{Float64}
    @test Base.return_types(_readnetmigration, (T,))[1] === Matrix{Int64}
    # `GridHabitat` is parameterised on its topology for this reason alone. Declared as the bare
    # `AbstractTopology` it did not infer, so the `calc_lookup_moves!` call in `_move!` dispatched
    # dynamically and boxed its `Int64` arguments once per species per cell per timestep.
    @test isconcretetype(Base.return_types(_readtopology, (T,))[1])
    # Each species' movement is read per species per cell: the vector must infer concretely, its
    # element a concrete movement or a union of at most four that the call to `move!` splits.
    for eco in (eco, Test1Ecosystem(seed = 1, mixed = true))
        M = Base.return_types(_readmovement, (typeof(eco),))[1]
        @test isconcretetype(M)
        @test length(Base.uniontypes(eltype(M))) <= 4
    end
    @test length(Base.uniontypes(eltype(Test1Ecosystem(seed = 1,
                                                       mixed = true).spplist.movement))) ==
          4
end

# **A sweep, because the checks above are a LIST and a list cannot catch a field it does not name.**
# That is not hypothetical: `topology` was abstract from the day it was written, cost 240 bytes a
# cell, and was absent from the list until it was measured - the four paths there were added when
# `abundances` was fixed, and the next such field simply was not among them.
#
# Every entry here must be looked at rather than counted: a new abstract field on any type the loop
# reaches fails this until it is either fixed or listed with the reason it is safe.
#
# This deliberately does NOT subsume the path checks above, and neither subsumes the other. A
# field walk asks whether the *container* is concrete; the path checks ask whether inference reaches
# *through* an abstract one. `Ecosystem.abundances` is the case that separates them - abstract by
# design, and its inner fields infer perfectly.
const _ABSTRACT_BY_DESIGN = Dict((:Ecosystem,
                                  :abundances) =>
                                     "abstract container by design - its " *
                                     "inner fields name no type parameter, " *
                                     "so inference reaches through it. The " *
                                     "path assertions above are what cover it.",
                                 (:Ecosystem,
                                  :epoch) =>
                                     "a `Union{Nothing, TimeType}`, because a " *
                                     "run may have no epoch. Not read in the loop.",
                                 (:Ecosystem,
                                  :calendar) =>
                                     "an `AbstractRunCalendar`, a marker set " *
                                     "by `build_ecosystem` and read only when " *
                                     "dates become elapsed time. Not read in the loop.",
                                 (:Ecosystem,
                                  :placement) =>
                                     "`Any`, because a placement may be any " *
                                     "function. Read only by `repopulate!` and " *
                                     "the run's provenance, never in the loop.",
                                 (:SpeciesList,
                                  :authority) =>
                                     "a `Union{Nothing, AbstractNameAuthority}`, " *
                                     "because species may have no name authority. " *
                                     "Not read in the loop.",
                                 (:GridHabitat,
                                  :active) =>
                                     "abstract, and measured to cost nothing: " *
                                     "indexing an abstract array dispatches " *
                                     "dynamically without boxing, where calling " *
                                     "a generic through one must box.")

@testset "no new abstract field on the types the hot loop reaches" begin
    eco = Test1Ecosystem(seed = 1)
    unexplained = Tuple{Symbol, Symbol}[]
    for obj in (eco, eco.habitat, eco.spplist, eco.cache)
        T = typeof(obj)
        for f in fieldnames(T)
            isconcretetype(fieldtype(T, f)) && continue
            haskey(_ABSTRACT_BY_DESIGN, (nameof(T), f)) && continue
            push!(unexplained, (nameof(T), f))
        end
    end
    @test unexplained == Tuple{Symbol, Symbol}[]

    # And the list may not rot: every entry must still name a field that is actually abstract.
    stale = [k
             for k in keys(_ABSTRACT_BY_DESIGN)
             if !any((eco, eco.habitat, eco.spplist, eco.cache)) do obj
        return nameof(typeof(obj)) == k[1] && k[2] in fieldnames(typeof(obj)) &&
               !isconcretetype(fieldtype(typeof(obj), k[2]))
    end]
    @test stale == []
end

@testset "the hot loop allocates nothing per cell" begin
    # Two grids, and the assertion is about the *slope*, not a byte count. `Threads.@threads` has
    # its own per-loop cost which is constant in the number of cells, so it cancels; a fixed
    # threshold would instead have to be retuned for every thread count and Julia version.
    small = Test1Ecosystem(seed = 1, grid = (5, 7))
    large = Test1Ecosystem(seed = 1, grid = (20, 28))
    ncells(eco) = size(eco.abundances.matrix, 2)
    extracells = ncells(large) - ncells(small)
    @test extracells > 500
    # Fewer than one byte per additional cell: a boxed field costs about 176 of them.
    @test _usagealloc(large) - _usagealloc(small) < extracells
end

# The measured allocation of one whole timestep, warmed first so the figure is the steady state
# rather than the compilation.
function _updatealloc(eco, timestep)
    EcoSISTEM.update!(eco, timestep)
    EcoSISTEM.update!(eco, timestep)
    return @allocated EcoSISTEM.update!(eco, timestep)
end

@testset "a whole timestep allocates nothing per cell or per species" begin
    # **Both directions, because a fixture that varies one can only prove half of it.** The check
    # above covers `update_resource_usage!` alone, which never touches the topology - so it stayed
    # green while a whole timestep cost 240 bytes a cell.
    #
    # Slopes rather than byte counts, as above: `Threads.@threads` costs a few hundred bytes per
    # thread, constant in the problem's size, so it cancels from a difference but would force a
    # fixed threshold to be retuned per thread count.
    timestep = 1.0month_mean_duration

    smallgrid = Test1Ecosystem(seed = 1, grid = (5, 7))
    largegrid = Test1Ecosystem(seed = 1, grid = (20, 28))
    extracells = size(largegrid.abundances.matrix, 2) -
                 size(smallgrid.abundances.matrix, 2)
    @test extracells > 500
    # Under a byte per added cell. Before the topology was parameterised this was about 240.
    @test _updatealloc(largegrid, timestep) -
          _updatealloc(smallgrid, timestep) < extracells

    fewspecies = Test1Ecosystem(seed = 1, numspecies = 10)
    manyspecies = Test1Ecosystem(seed = 1, numspecies = 80)
    extraspecies = 80 - 10
    # And under a byte per added species, where it was about 560.
    @test _updatealloc(manyspecies, timestep) -
          _updatealloc(fewspecies, timestep) < extraspecies
end

# The same with four kinds of movement, the most a species list may mix: each species' `move!` is
# chosen from a union the compiler splits, which must cost nothing per cell or species, as it would
# not were the union abstract or larger.
@testset "a timestep with four kinds of movement allocates nothing per cell or species" begin
    timestep = 1.0month_mean_duration
    smallgrid = Test1Ecosystem(seed = 1, grid = (5, 7), mixed = true)
    largegrid = Test1Ecosystem(seed = 1, grid = (20, 28), mixed = true)
    extracells = size(largegrid.abundances.matrix, 2) -
                 size(smallgrid.abundances.matrix, 2)
    @test _updatealloc(largegrid, timestep) -
          _updatealloc(smallgrid, timestep) < extracells
    fewspecies = Test1Ecosystem(seed = 1, numspecies = 12, mixed = true)
    manyspecies = Test1Ecosystem(seed = 1, numspecies = 80, mixed = true)
    @test _updatealloc(manyspecies, timestep) -
          _updatealloc(fewspecies, timestep) < 80 - 12
end

end
