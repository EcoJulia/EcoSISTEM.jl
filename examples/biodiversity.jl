# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Ecological **property** tests: does the simulation reproduce the qualitative results it is supposed
# to? Species should be commonest near their temperature optimum, narrow specialists should beat broad
# generalists in a constant regime, more resource should support more individuals, a bigger area should
# hold more of them, most species should persist, and a longer dispersal kernel should spread a
# population further across the grid.
#
# These are properties, not blessed numbers — the companion to `test/canonical/`, not a member of
# it. A canonical test pins *what the model produced*; these pin *what the model must never stop
# doing*, and no re-blessing can silence them.
#
# **It lives here, among the examples, rather than in `test/`** because what it demonstrates is of
# interest in its own right: these are the ecological claims the model is supposed to reproduce, and
# they should be easy to read and to re-run. It is an example that happens to assert.
#
# **It never runs under the test suite or in CI.** `test/extras_examples.jl` does execute every
# top-level `examples/*.jl`, but `runtests.jl` sets `ECOSISTEM_SCALE=small` and this file skips itself
# whenever it is set. It is ~20 whole simulations, several of 10–50 simulated years over 100–500
# species, measured at **3 minutes** — far too slow for a routine gate. Run it deliberately:
#
#     julia --project=examples examples/biodiversity.jl
#     ECOSISTEM_SCALE=large julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_examples.jl"])'
#
# The second works because `runtests.jl` sets the variable with `get!`, so an explicit value wins.
#
# It is seeded throughout and must stay that way. Every assertion below is a statement about a
# *stochastic* simulation, so an unseeded run fails intermittently — as this file did, in its previous
# life, with two failures per run and a different two each time. A property test that cries wolf gets
# ignored, which is worse than not having it.

module Biodiversity

using Test
using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Distributions
using Random

# The gate. `runtests.jl` sets `ECOSISTEM_SCALE=small` for every `Pkg.test` run, so naming this set
# there is not enough to run it — that is deliberate, because it is far too slow to sit in any routine
# gate. A bare script run leaves the variable unset (read as full scale), and an explicit
# `ECOSISTEM_SCALE=large` overrides the suite's default, so both deliberate routes work.
if get(ENV, "ECOSISTEM_SCALE", "large") == "small"
    # `println`, not `@info`: logging goes to **stderr**, and `test/extras_examples.jl` wraps every
    # example in `@test_nowarn`, which fails on *any* stderr output rather than only on warnings — so
    # announcing the skip through the logger fails the very gate the skip exists to keep fast.
    println("Skipping the biodiversity property tests: they are the most expensive thing in the " *
            "repository and never run under the test suite. Set ECOSISTEM_SCALE=large to run them.")
else
    # One seed for the whole file. `build_species`/`build_ecosystem` take it and build one
    # deterministic RNG stream per species, so results do not depend on the thread or MPI-rank count.
    const SEED = 20260804

    # The published configuration: supply is **per unit area**, demand is **per unit of plant**, so a
    # species' actual demand is `DEMAND .* size`. Keeping demand per m² rather than per individual is
    # what lets the large-pool test below give every species its own body size.
    const SUPPLY = (sunlight = 4.5e11kJ / (km^2 * day), water = 192.0mm / day)
    const DEMAND = (sunlight = 450000.0kJ / (m^2 * day),
                    water = 192.0Unitful.L / (m^2 * day))

    const NUMSPECIES = 100
    const CELLS = 10
    const AREA = 100.0km^2

    # A square study area of total `area`, divided into `cells × cells`. Cell size follows from the
    # two rather than being fixed — that is what lets the grid-resolution test carve one landscape
    # four ways and the area test grow the landscape at a fixed grid.
    function _studyarea(cells, area)
        side = sqrt(area)
        return StudyArea(extent = (side, side), cellsize = side / cells,
                         verbosity = :silent)
    end

    # `topology` is an **environment** keyword since step 19: whether the grid wraps is a property
    # of the map, not of the species dispersing over it. Defaults to `Torus()` here because these
    # experiments were written against a wrapping grid, and several of their assertions depend on
    # having no edge — the ones that want an edge say so.
    """
        uniform_environment(cells, area; temperature = 298.0K)

    Build a `cells × cells` environment over `area`, uniform in every respect: one temperature
    everywhere, with uniform solar and water supply.

    Returned on its own rather than folded into the community builder so that a test can impose a
    supply gradient on it before the species arrive — see "more resource supports more individuals".
    """
    function uniform_environment(cells, area; temperature = 298.0K,
                                 topology = Torus())
        return GridHabitat(regime = UniformSpec(temperature,
                                                axis = Temperature),
                           supply = (sunlight = UniformSpec(SUPPLY.sunlight,
                                                            axis = SolarRadiation),
                                     water = UniformSpec(SUPPLY.water,
                                                         axis = Precipitation)),
                           area = _studyarea(cells, area),
                           topology = topology)
    end

    """
        community(env, opts, widths; kw...)

    Put `length(opts)` species on `env`, with Gaussian temperature tolerances (`opts` the optima,
    `widths` the niche widths) and Gaussian dispersal of mean distance `dispersal`.

    Every per-species keyword takes a scalar or a length-`n` vector, so a test can vary body size,
    dispersal or mortality across the pool without a different builder. `individuals` is a **total**,
    split across species from the seeded stream — so the starting state is identical on every run,
    without which none of the assertions here are reproducible.
    """
    function community(env, opts, widths; dispersal = 2.4km,
                       individuals = 100_000_000,
                       birth = 0.15 / year, death = 0.15 / year,
                       sizes = 1.0m^2, seed = SEED)
        species = build_species(length(opts), tolerance = (opts, widths),
                                toleranceaxis = Temperature,
                                demand = map(d -> d .* sizes, DEMAND),
                                demandaxis = (sunlight = SolarRadiation,
                                              water = Precipitation),
                                dispersal = dispersal,
                                birth = birth, death = death, survival = 0.1,
                                abundance = individuals, seed = seed)
        return build_ecosystem(species, env, seed = seed)
    end

    # Run for `times` and hand back the ecosystem, so each assertion reads as one statement.
    function simulated(eco; times = 10year, timestep = 1month_mean_duration)
        simulate!(eco, times, timestep)
        return eco
    end

    # Total abundance, over the whole grid, of the species whose niche parameter falls in `[lo, hi)`.
    #
    # This replaces a histogram built by `repeat`ing each optimum `endabun[i]` times for
    # `i in eachindex(opts)`. `endabun` is a species × cells matrix, so those linear indices ran down
    # column 1: the old assertions were measured on a *single grid cell* out of 100, and materialised a
    # million-element vector of unit-carrying temperatures to do it. Summing per species over `dims = 2`
    # asks the same ecological question of the entire landscape, and allocates one vector of 100.
    function abundance_in(eco, values, lo, hi)
        byspecies = vec(sum(eco.abundances.matrix, dims = 2))
        return sum(byspecies[(lo .<= values) .& (values .< hi)])
    end

    # Total abundance per cell, laid back out on the grid.
    function percell(eco)
        return dropdims(sum(eco.abundances.grid, dims = 1), dims = 1)
    end

    @testset "Biodiversity properties" begin
        @testset "commonest at the temperature optimum" begin
            # Optima spread ±3 widths around the environment's own 298 K, so the species matched to
            # the regime should end up commoner than those at either extreme of the gradient.
            widths = fill(2.0K, NUMSPECIES)
            opts = 298.0K .+
                   widths .* range(-3, stop = 3, length = NUMSPECIES)
            eco = simulated(community(uniform_environment(CELLS, AREA), opts,
                                      widths))
            @test abundance_in(eco, opts, 298.0K, 299.0K) >
                  abundance_in(eco, opts, 292.0K, 293.0K)
            @test abundance_in(eco, opts, 298.0K, 299.0K) >
                  abundance_in(eco, opts, 303.0K, 304.0K)
        end

        @testset "specialists beat generalists in a matched regime" begin
            # Every species shares the environment's optimum and differs only in niche width. In a
            # constant regime there is nothing for a wide niche to buy, so the narrowest should win.
            widths = collect(range(0.0001K, stop = 5.0K, length = NUMSPECIES))
            opts = fill(298.0K, NUMSPECIES)
            eco = simulated(community(uniform_environment(CELLS, AREA), opts,
                                      widths))
            @test abundance_in(eco, widths, 0.1K, 0.3K) >
                  abundance_in(eco, widths, 4.7K, 4.9K)
        end

        @testset "a mismatched regime favours intermediate niche widths" begin
            # The same species, but the environment is 1 K off every optimum. Now the narrowest species
            # are excluded by the mismatch and the widest are still too diffuse, so the peak moves
            # inwards — the one case in this file where the answer is non-monotonic, and why it is here.
            widths = collect(range(0.0001K, stop = 5.0K, length = NUMSPECIES))
            opts = fill(298.0K, NUMSPECIES)
            env = uniform_environment(CELLS, AREA, temperature = 299.0K)
            eco = simulated(community(env, opts, widths))
            @test abundance_in(eco, widths, 1.1K, 1.3K) >
                  abundance_in(eco, widths, 4.7K, 4.9K)
            @test abundance_in(eco, widths, 1.1K, 1.3K) >
                  abundance_in(eco, widths, 0.1K, 0.3K)
        end

        @testset "more resource supports more individuals" begin
            # Each resource is graded in its **own** environment, with the other left uniform. That
            # is forced by the model, not a stylistic choice: `_resourceadjustment` scales birth by
            # `min(K/E)` over the resources (Liebig's law of the minimum, `src/Generate.jl`), so in a
            # cell where water binds, more solar buys exactly nothing. Grading both at once therefore
            # ties every cell along the scarce edge — `min(0.25, 0.25) == min(1.0, 0.25)` — and the
            # comparison tests the minimum operator rather than the ecology.
            #
            # The gradient is a **fraction of the environment's own uniform supply**, read back off
            # the matrix, rather than a hard-coded quantity. The original wrote literals — `192.0 L/day`
            # for water, against a true per-cell supply of ~1.92e8 L/day once supplies became rates —
            # which made water limiting in every cell, wiped out the entire grid, and left both
            # assertions comparing 0 with 0. Deriving the gradient cannot drift out of step again.
            #
            # It also runs 0.25 → 1, not 0 → 1: a zero-supply edge is uninhabitable, and `0 < 0` fails
            # an increasing-abundance test for a reason that has nothing to do with the property.
            fracs = range(0.25, stop = 1.0, length = CELLS)
            opts, widths = fill(298.0K, NUMSPECIES), fill(2.0K, NUMSPECIES)

            # Named supplies, so these read as the resource they grade rather than as `.one`/`.two`.
            byrow = uniform_environment(CELLS, AREA, topology = Island())
            for (i, f) in enumerate(fracs)
                byrow.supply.sunlight.matrix[i, :] .*= f
            end
            rowtotal = percell(simulated(community(byrow, opts, widths)))
            @test all(rowtotal[1, :] .< rowtotal[end, :])

            bycol = uniform_environment(CELLS, AREA, topology = Island())
            for (j, f) in enumerate(fracs)
                bycol.supply.water.matrix[:, j] .*= f
            end
            coltotal = percell(simulated(community(bycol, opts, widths)))
            @test all(coltotal[:, 1] .< coltotal[:, end])
        end

        @testset "invariant to grid resolution" begin
            # The same landscape carved into 1, 4, 25 or 100 cells: total abundance is a property of
            # the area and its resource, so cutting it up more finely must not change the answer.
            totals = map([1, 2, 5, 10]) do cells
                eco = community(uniform_environment(cells, AREA,
                                                    topology = Island()),
                                fill(298.0K, NUMSPECIES),
                                fill(2.0K, NUMSPECIES))
                return sum(simulated(eco).abundances.matrix)
            end
            @test all(isapprox.(totals[1], totals, rtol = 0.1))
        end

        @testset "abundance scales with area" begin
            # Resource is specified per unit area, so a larger landscape carries proportionally more of
            # it and must support strictly more individuals.
            totals = map([10.0, 20.0, 50.0, 100.0]) do a
                eco = community(uniform_environment(CELLS, a * km^2,
                                                    topology = Island()),
                                fill(298.0K, NUMSPECIES),
                                fill(2.0K, NUMSPECIES))
                return sum(simulated(eco).abundances.matrix)
            end
            @test totals[1] < totals[2] < totals[3] < totals[4]
        end

        @testset "sustains a large number of species" begin
            # Heterogeneous species — random widths, optima, dispersal, mortality and body size — to
            # check that coexistence is not an artefact of every species being identical, and that it
            # holds as the pool grows tenfold.
            #
            # Every one of these is a plain per-species vector handed to `community`: the builder
            # takes a scalar or a length-`n` vector for each, so heterogeneity needs no separate path.
            pool = [50, 100, 200, 500]
            surviving = map(pool) do n
                rng = Xoshiro(SEED)
                widths = rand(rng, Uniform(1.0, 5.0), n) .* K
                opts = 298.0K .+ widths .* range(-3, stop = 3, length = n)
                death = abs.(rand(rng, Normal(0.15, 0.135), n)) ./ year
                eco = community(uniform_environment(CELLS, AREA), opts, widths,
                                dispersal = rand(rng, Uniform(0.6, 2.4), n) .*
                                            km,
                                birth = death, death = death,
                                sizes = rand(rng, Normal(1.0, 0.05), n) .* m^2)
                return count(>(0),
                             sum(simulated(eco).abundances.matrix, dims = 2))
            end
            @test all(surviving ./ pool .>= 0.9)
        end

        @testset "longer dispersal spreads further" begin
            # Two species seeded on opposite edges of the grid, with no boundary to wrap around. Larger
            # mean dispersal distance must drain the edges and fill the middle, so the edge columns fall
            # monotonically with kernel width while the centre column rises.
            totals = map([0.5, 1.0, 2.0, 4.0] .* km) do d
                eco = community(uniform_environment(CELLS, AREA,
                                                    topology = Island()),
                                fill(298.0K, 2), fill(2.0K, 2),
                                dispersal = d, individuals = 0)
                # Start each species pinned to one edge, rather than spread over the grid, so that what
                # the assertions measure is dispersal and nothing else.
                eco.abundances.grid[1, :, 1] .= 100
                eco.abundances.grid[2, :, end] .= 100
                return percell(simulated(eco, times = 50year))
            end
            # `>= 9` of 10 rows, not all 10: this is a stochastic simulation and one row may tie or
            # invert by chance. Seeded, so the count is now reproducible rather than a coin toss.
            @test count(i -> totals[1][i, 1] > totals[2][i, 1] >
                             totals[3][i, 1] > totals[4][i, 1], 1:CELLS) >= 9
            @test count(i -> totals[1][i, end] > totals[2][i, end] >
                             totals[3][i, end] > totals[4][i, end], 1:CELLS) >=
                  9
            @test count(i -> totals[1][i, 5] < totals[2][i, 5] <
                             totals[3][i, 5] < totals[4][i, 5], 1:CELLS) >= 9
        end
    end
end

end
