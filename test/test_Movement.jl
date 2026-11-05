# SPDX-License-Identifier: LGPL-3.0-or-later

module TestMovement

using EcoSISTEM
# These are `public` rather than exported - a spec is what a user writes,
# and these are what it materialises into.
using EcoSISTEM: Periodic, Bounded, EdgeTopology
using Test
using Distributions
using Unitful.DefaultSymbols
using EcoSISTEM.Units

@testset "Movement" begin
    numSpecies = 10

    @test_nowarn GaussianKernel(meandist = 0.2km, thresh = 10e-4)
    @test_nowarn LongTailKernel(scale = 0.2km, shape = 2.0, thresh = 10e-4)
    kernel = fill(GaussianKernel(meandist = 0.2m, thresh = 10e-4), numSpecies)
    for mode in (SeedDispersal(GaussianKernel, disperse_safely = true), InsectDispersal(), AlwaysDispersal())
        movement = mode.(kernel)
        @test EcoSISTEM.getkernels(movement) == kernel
        @test EcoSISTEM.dispersesafely(movement, 1) == mode.disperse_safely
    end
    none = fill(NoDispersal()(), numSpecies)
    @test all(isnothing, EcoSISTEM.getkernels(none))
    @test EcoSISTEM.dispersesafely(none, 1)
    kernel = fill(LongTailKernel(scale = 0.2m, shape = 2.0, thresh = 10e-4),
                  numSpecies)
    @test EcoSISTEM.getkernels(VertebrateDispersal().(kernel)) == kernel
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
        EcoSISTEM.move!(eco, eco.spplist.movement[1], 1, 1,
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
    # The flag is the mode's, and defaults to the safe (redistributing) behaviour.
    kernels = fill(GaussianKernel(meandist = 0.2km, thresh = 10e-4), 3)
    @test all(sp -> EcoSISTEM.dispersesafely(SeedDispersal(GaussianKernel,
                                                           disperse_safely = true).(kernels),
                                             sp), 1:3)
    @test !EcoSISTEM.dispersesafely(SeedDispersal(GaussianKernel,
                                                  disperse_safely = false).(kernels),
                                    1)

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
                                movement = InsectDispersal(),
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
                          movement = InsectDispersal(),
                          disperse_safely = [true, false],
                          abundance = 10_000, seed = 1)
    @test [EcoSISTEM.dispersesafely(mixed.movement, sp) for sp in 1:2] ==
          [true, false]
    eco = build_ecosystem(mixed, env, seed = 1)
    simulate!(eco, 2.0year, 1.0month_mean_duration)
    # The lossy species must end up worse off than the safe one on the *same* grid, under the
    # same seed - which is only a meaningful comparison because they now genuinely differ.
    totals = sum(eco.abundances.matrix, dims = 2)[:, 1]
    @test totals[2] < totals[1]
end

# A species list's movement prints in one phrase whatever the species count: the mode and how far
# its species reach where they share one, else each mode and how many species take it.
@testset "a species list's movement prints in one phrase" begin
    phrase = EcoSISTEM._movementphrase
    @test phrase(fill(SeedDispersal(GaussianKernel, disperse_safely = true)(meandist = 1.0km,
                                                                            thresh = 1e-3),
                      10)) ==
          "SeedDispersal(GaussianKernel, disperse_safely = true), GaussianKernel mean 1.0 km"
    @test phrase([SeedDispersal(GaussianKernel, disperse_safely = true)(meandist = d,
                                                                        thresh = 1e-3)
                  for d in range(0.5km, 4.0km, length = 1000)]) ==
          "SeedDispersal(GaussianKernel, disperse_safely = true), GaussianKernel mean 0.5 km " *
          "to 4.0 km"
    long = fill(VertebrateDispersal()(scale = 2.0km, shape = 2.0,
                                      thresh = 1e-3),
                1000)
    @test phrase(long) ==
          "VertebrateDispersal(LongTailKernel, disperse_safely = true), LongTailKernel scale 2.0 km"
    @test phrase(fill(NoDispersal()(), 1000)) == "NoDispersal()"
    mixed = [AlwaysDispersal(disperse_safely = sp > 3)(meandist = 1.0km)
             for sp in 1:1000]
    @test phrase(mixed) ==
          "AlwaysDispersal(GaussianKernel, disperse_safely = false) for 3, " *
          "AlwaysDispersal(GaussianKernel, disperse_safely = true) for 997"
    @test phrase(SpeciesDispersal{NoDispersal, Nothing}[]) == "none"
end

# A dispersal mode fixes which individuals disperse, its kernel's family and the dead-cell fate, so
# species sharing a mode differ only in their kernels' parameters.
@testset "dispersal modes" begin
    # The generic mode states its family and fate; the named ones default them for their agent.
    @test_throws MethodError SeedDispersal()
    @test_throws UndefKeywordError SeedDispersal(GaussianKernel)
    @test !SeedDispersal(WaldKernel, disperse_safely = false).disperse_safely
    for (mode, family, safely) in ((WindDispersal(), WaldKernel, false),
        (InsectDispersal(), GaussianKernel, true),
        (VertebrateDispersal(), LongTailKernel, true),
        (AlwaysDispersal(), GaussianKernel, true))
        @test EcoSISTEM._family(mode) === family
        @test mode.disperse_safely == safely
    end
    insect = InsectDispersal()
    @test insect(meandist = 2.0km) ==
          SpeciesDispersal(insect, GaussianKernel(meandist = 2.0km))
    @test insect(GaussianKernel(meandist = 2.0km)) == insect(meandist = 2.0km)
    @test NoDispersal()() == SpeciesDispersal(NoDispersal(), nothing)
    @test SeedDispersal(LongTailKernel, disperse_safely = true)(scale = 1.0km,
                                                                shape = 2.0).kernel isa
          LongTailKernel
    @test_throws "of its family" SpeciesDispersal(insect,
                                                  WaldKernel(meandist = 1.0km,
                                                             scale = 1.0km))
    @test_throws "of its family" SpeciesDispersal(NoDispersal(),
                                                  GaussianKernel(meandist = 1.0km))
    @test_throws "concrete kernel type" SeedDispersal(EcoSISTEM.AbstractKernel,
                                                      disperse_safely = true)
    @test repr(WindDispersal()) ==
          "WindDispersal(WaldKernel, disperse_safely = false)"
    @test repr(insect(meandist = 2.0km)) ==
          "SpeciesDispersal(InsectDispersal(GaussianKernel, disperse_safely = true), " *
          "GaussianKernel(meandist = 2.0 km, thresh = 1.0e-9))"

    # A species list narrows mixed movements to the union of their own types, at most four, so
    # the step loop splits them rather than dispatching dynamically.
    narrow = EcoSISTEM._dispersalvector
    kinds = [
        SeedDispersal(GaussianKernel, disperse_safely = true)(meandist = 1.0km),
        AlwaysDispersal()(meandist = 1.0km),
        InsectDispersal()(meandist = 1.0km), NoDispersal()()]
    @test !isconcretetype(eltype(kinds))
    @test Set(Base.uniontypes(eltype(narrow(kinds)))) == Set(typeof.(kinds))
    @test narrow(kinds) == kinds
    fifth = [kinds; WindDispersal(GaussianKernel)(meandist = 1.0km)]
    @test_throws "at most four kinds" narrow(fifth)
    # Two species in one mode with different fates are one kind: the fate is a value, not a type.
    @test eltype(narrow([
                            SeedDispersal(GaussianKernel,
                                          disperse_safely = true)(meandist = 1.0km),
                            SeedDispersal(GaussianKernel,
                                          disperse_safely = false)(meandist = 2.0km)
                        ])) ==
          SpeciesDispersal{SeedDispersal{GaussianKernel}, GaussianKernel}
end

# A kernel's distances are named for what they are. A Gaussian kernel's mean distance and per-axis
# standard deviation describe it alike; a long-tailed kernel's scale and mean distance do not, and
# its mean exists only for a shape above 3/2.
@testset "kernels take their distances by name" begin
    g = GaussianKernel(meandist = 2.4km)
    @test g.thresh == 1.0e-9
    @test EcoSISTEM._meandist(g) ≈ 2.4km
    @test GaussianKernel(sd = 2.4km / sqrt(pi / 2)).scale ≈ g.scale
    @test GaussianKernel(meandist = 2400.0m) == g
    # The 2Dt's mean is `scale * sqrt(pi) / 2 * gamma(shape - 3/2) / gamma(shape - 1)`: at shape 5/2
    # the gamma ratio is 2 / sqrt(pi), so mean and scale coincide, and at shape 3 it is sqrt(pi) / 2.
    @test LongTailKernel(meandist = 3.0km, shape = 2.5).scale ≈ 3.0km
    @test LongTailKernel(meandist = 3.0km, shape = 3.0).scale ≈ 12 / pi * km
    @test LongTailKernel(scale = 1.0km, shape = 1.5).shape == 1.5

    @test_throws "exactly one of `meandist` and `sd`; got neither" GaussianKernel()
    @test_throws "got both" GaussianKernel(meandist = 1.0km, sd = 1.0km)
    @test_throws "must be positive" GaussianKernel(meandist = 0.0km)
    @test_throws "must exceed 1" LongTailKernel(scale = 1.0km, shape = 1.0)
    @test_throws "no finite mean distance" LongTailKernel(meandist = 1.0km,
                                                          shape = 1.5)
    @test_throws "exactly one of `scale` and `meandist`" LongTailKernel(shape = 2.0)
    @test_throws UndefKeywordError LongTailKernel(scale = 1.0km)

    @test repr(g) == "GaussianKernel(meandist = 2.4 km, thresh = 1.0e-9)"
    @test repr(LongTailKernel(scale = 2.0km, shape = 2.0, thresh = 1e-3)) ==
          "LongTailKernel(shape = 2.0, scale = 2.0 km, thresh = 0.001)"

    w = WaldKernel(meandist = 5000.0m, scale = 10.0km)
    @test w.meandist == 5.0km && w.scale == 10.0km && w.thresh == 1.0e-9
    @test_throws UndefKeywordError WaldKernel(meandist = 5.0km)
    @test_throws "must be positive" WaldKernel(meandist = 5.0km, scale = 0.0km)
    @test repr(w) ==
          "WaldKernel(meandist = 5.0 km, scale = 10.0 km, thresh = 1.0e-9)"
    @test EcoSISTEM._movementphrase(WindDispersal().(fill(w, 3))) ==
          "WindDispersal(WaldKernel, disperse_safely = false), WaldKernel mean 5.0 km"
end

# A lookup must put individuals at the distances its kernel was given: the mean displacement from
# cell centres, on 1 km cells and a grid wide enough not to truncate, is the kernel's mean to within
# the cells' discreteness, about half a percent for every kernel alike. A WALD lookup built from its
# distance density without dividing by `2pi d` puts its mean at 7.5 km instead, and fails.
@testset "a lookup puts individuals at the distances its kernel was given" begin
    area = StudyArea(extent = (81.0km, 81.0km), cellsize = 1.0km,
                     verbosity = :silent)
    regime = build_habitat(area = area, verbosity = :silent).regime
    lookup(k) = EcoSISTEM.genlookups(regime, k)
    meanof(l) = sum(p * sqrt(y^2 + x^2) for (y, x, p) in zip(l.y, l.x, l.p))
    spreadof(l) = sqrt(sum(p * (sqrt(y^2 + x^2) - meanof(l))^2
                           for (y, x, p) in zip(l.y, l.x, l.p)))
    loose = lookup(WaldKernel(meandist = 5.0km, scale = 10.0km))
    tight = lookup(WaldKernel(meandist = 5.0km, scale = 50.0km))
    @test meanof(lookup(GaussianKernel(meandist = 5.0km))) ≈ 5.0 rtol = 0.02
    @test meanof(loose) ≈ 5.0 rtol = 0.02
    @test meanof(tight) ≈ 5.0 rtol = 0.02
    # The variance is `meandist^3 / scale`, so the larger scale gathers the distances closer.
    @test spreadof(loose) ≈ sqrt(5.0^3 / 10.0) rtol = 0.05
    @test spreadof(tight) < spreadof(loose)
end

end
