# SPDX-License-Identifier: LGPL-3.0-or-later

module TestMovement

using EcoSISTEM
# `[C7-VIS]` C: these are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: Periodic, Bounded, EdgeTopology
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units

@testset "Movement" begin
    numSpecies = 10

    @test_nowarn GaussianKernel.(fill(0.2km, numSpecies), 10e-4)
    @test_nowarn LongTailKernel.(fill(0.2km, numSpecies), 1.0, 10e-4)
    kernel = GaussianKernel.(fill(0.2m, numSpecies), 10e-4)
    @test_nowarn AlwaysMovement(kernel)
    @test_nowarn AlwaysMovement(kernel)
    @test_nowarn AlwaysMovement(kernel)
    @test_nowarn BirthOnlyMovement(kernel)
    @test_nowarn NoMovement(kernel)

    kernel = GaussianKernel.(fill(0.2m, numSpecies), 10e-4)
    mov = AlwaysMovement(kernel)
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    @test EcoSISTEM.dispersesafely(mov, 1) == mov.disperse_safely[1]
    mov = BirthOnlyMovement(kernel)
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    @test EcoSISTEM.dispersesafely(mov, 1) == mov.disperse_safely[1]
    mov = NoMovement(kernel)
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    kernel = LongTailKernel.(fill(0.2m, numSpecies), 1.0, 10e-4)
    mov = AlwaysMovement(kernel)
    @test EcoSISTEM.getkernels(mov) == mov.kernels
    @test EcoSISTEM.dispersesafely(mov, 1) == mov.disperse_safely[1]
end

# **The test that would have caught the proposal's original axis order.** `EdgeTopology{BCY, BCX}`
# is positional, so which parameter is which is invisible at every call site - and the proposal as
# offered had `Cylinder = EdgeTopology{Periodic, Bounded}`, i.e. Y wrapping, when the code
# measurably wraps X. Asserting it *behaviourally* rather than structurally is what makes this
# worth keeping: a structural `Cylinder === EdgeTopology{Bounded, Periodic}` restates the definition,
# while this would fail if the wiring in `calc_lookup_moves!` were reversed too.
@testset "boundary conditions apply to the axis they name" begin
    # `_stepto` is the whole of the per-axis rule, so pin it directly first.
    @test isnothing(EcoSISTEM._stepto(EcoSISTEM.Bounded, 1, -1, 10))   # off the low edge
    @test isnothing(EcoSISTEM._stepto(EcoSISTEM.Bounded, 10, 1, 10))   # off the high edge
    @test EcoSISTEM._stepto(EcoSISTEM.Bounded, 4, 1, 10) == 5
    @test EcoSISTEM._stepto(EcoSISTEM.Periodic, 1, -1, 10) == 10        # wraps to the far side
    @test EcoSISTEM._stepto(EcoSISTEM.Periodic, 10, 1, 10) == 1

    # ...then behaviourally, on a **non-square** grid so a Y/X mix-up cannot hide. One species is
    # pinned to the first row and the first column; after dispersal, a wrapping axis has carried
    # individuals to the *opposite* end of that axis and a bounded one has not.
    function _reach(topology)
        area = StudyArea(extent = (11.0km, 7.0km), cellsize = 1.0km,
                         verbosity = :silent)
        env = GridHabitat(regime = UniformSpec(298.0K,
                                               axis = Temperature),
                          supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                               axis = SolarRadiation),
                          area = area, topology = topology)
        sp = build_species(1, tolerance = (298.0K, 3.0K),
                           toleranceaxis = Temperature,
                           demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                           dispersal = 1.2km, abundance = 0, seed = 1)
        eco = build_ecosystem(sp, env, seed = 1)
        ny, nx = size(eco.abundances.grid, 2), size(eco.abundances.grid, 3)
        eco.abundances.grid .= 0
        eco.abundances.grid[1, 1, 1] = 100_000        # a corner
        # `move!(eco, movement, sp, loc, netmigration, abun)` - `loc` is the flat cell index of
        # the corner the individuals were placed in.
        EcoSISTEM.move!(eco, eco.spplist.movement, 1, 1,
                        eco.cache.netmigration, 100_000)
        # `netmigration` is the **flat** `species × location` cache, not a grid - reshaped here
        # the same way `GridLandscape` views its own matrix, with Y varying fastest.
        moved = reshape(view(eco.cache.netmigration, 1, :), ny, nx)
        # did anything land at the far end of each axis?
        return (y = sum(moved[ny, :]) > 0, x = sum(moved[:, nx]) > 0)
    end

    # `Cylinder` wraps **x** and not **y** - the assertion the proposal's order would have failed.
    cyl = _reach(Cylinder())
    @test cyl.x
    @test !cyl.y
    # ...and the newly expressible mirror image wraps y and not x.
    flipped = _reach(EdgeTopology(y = Periodic, x = Bounded))
    @test flipped.y
    @test !flipped.x
    # `Island` wraps neither, `Torus` both.
    isl = _reach(Island())
    @test !isl.y && !isl.x
    tor = _reach(Torus())
    @test tor.y && tor.x
end

@testset "Dispersal into dead cells" begin
    # `disperse_safely` must not be read off the boundary object: that bundles a **grid** property
    # (does the world wrap) with a **species** one (is a disperser aimed at a dead cell
    # redistributed, or lost). They are separate - the topology is on the habitat, the flag is one
    # entry per species on the movement.
    for T in (Island, Cylinder, Torus)
        @test isempty(fieldnames(T))          # topologies carry no per-species state at all
    end
    # The flag defaults to the safe (redistributing) behaviour, one entry per species.
    mov = BirthOnlyMovement(GaussianKernel.(fill(0.2km, 3), 10e-4))
    @test mov.disperse_safely == fill(true, 3)
    @test all(sp -> EcoSISTEM.dispersesafely(mov, sp), 1:3)
    # One flag per species, or two vectors would silently describe different species.
    @test_throws "one entry per species" BirthOnlyMovement(GaussianKernel.(fill(0.2km,
                                                                                3),
                                                                           10e-4),
                                                           [true, false])

    # One landscape, one seed, one species set - only the topology and the flag vary.
    function _total(topology, safely = true; holes = false)
        area = StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                         verbosity = :silent)
        env = GridHabitat(regime = UniformSpec(298.0K,
                                               axis = Temperature),
                          supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                               axis = SolarRadiation),
                          area = area, topology = topology)
        if holes
            for y in 1:10, x in 1:10
                iseven(y + x) && (env.active[y, x] = false)
            end
        end
        species = build_species(20, tolerance = (298.0K, 3.0K),
                                toleranceaxis = Temperature,
                                demand = 1.0e9kJ / day,
                                demandaxis = SolarRadiation,
                                movement = BirthOnlyMovement,
                                disperse_safely = safely,
                                abundance = 10_000, seed = 1)
        eco = build_ecosystem(species, env, seed = 1)
        simulate!(eco, 2.0year, 1.0month_mean_duration)
        return sum(eco.abundances.matrix)
    end

    # The flag must be a **no-op when nothing is blocked**, not merely "close". A torus has no
    # edge, so with every cell active there is nothing to lose and the two settings must agree
    # exactly - which they only do because `_drawmoves!` skips the extra `Binomial` draw when
    # `lost` is zero. Without that guard the draw re-phases the species' RNG stream and the totals
    # diverge (measured: 5883 against 5805) despite identical ecology.
    @test _total(Torus(), true) == _total(Torus(), false)

    # A torus still has dead cells if any are *inactive*, and that is the case the flag exists
    # for as much as the map edge is.
    @test _total(Torus(), false, holes = true) <
          _total(Torus(), true, holes = true)

    # At the edge, losing dispersers costs abundance under both edged topologies.
    @test _total(Island(), false) < _total(Island(), true)
    @test _total(Cylinder(), false) < _total(Cylinder(), true)

    # And the ordering follows the geometry: a torus loses nothing, a cylinder loses its two
    # north-south edges, and `Island` loses all four.
    @test _total(Island(), false) < _total(Cylinder(), false) <
          _total(Torus(), false)

    # **The case that makes the flag per-species**: two species on ONE grid, differing only in
    # `disperse_safely`. A flag living on a single boundary object shared by the whole species list
    # cannot express it - a wind-dispersed species and an animal-dispersed one would have to agree
    # about what happens to a seed blown off the map.
    area = StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                     verbosity = :silent)
    env = GridHabitat(regime = UniformSpec(298.0K, axis = Temperature),
                      supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                           axis = SolarRadiation),
                      area = area)
    mixed = build_species(2, tolerance = (298.0K, 3.0K),
                          toleranceaxis = Temperature,
                          demand = 1.0e9kJ / day, demandaxis = SolarRadiation,
                          movement = BirthOnlyMovement,
                          disperse_safely = [true, false],
                          abundance = 10_000, seed = 1)
    @test mixed.movement.disperse_safely == [true, false]
    eco = build_ecosystem(mixed, env, seed = 1)
    simulate!(eco, 2.0year, 1.0month_mean_duration)
    # The lossy species must end up worse off than the safe one on the *same* grid, under the
    # same seed - which is only a meaningful comparison because they now genuinely differ.
    totals = sum(eco.abundances.matrix, dims = 2)[:, 1]
    @test totals[2] < totals[1]
end

end
