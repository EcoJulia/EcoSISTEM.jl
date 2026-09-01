# SPDX-License-Identifier: LGPL-3.0-or-later

module TestNicheAxisMacro

using EcoSISTEM
# `[C7-VIS]` B1/B2/B3: these are `public` rather than exported, so they must be named.
using EcoSISTEM: getdist
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using InteractiveUtils
using Distributions
using Test

# `@nicheaxis` is the *public* extension point - the `canonicalunit`/`supplytype`/`demandtype`/
# `bounds` methods behind it are all internal - so it is tested the way a user meets it: declared
# from **this** module, not from inside `EcoSISTEM`. That is load-bearing. The macro emits its
# methods fully qualified, and the first version did not; declared from inside the package that is
# invisible, while from out here it fails with "function EcoSISTEM.canonicalunit must be explicitly
# imported to be extended". Only an outside caller can catch it.
@testset "@nicheaxis declares axes from another module" begin
    @nicheaxis(TestGroupAxis<:EcoSISTEM.NicheAxis, condition=K,
               bounds=(0.0K, nothing))
    @nicheaxis(TestLeaf <: TestGroupAxis)
    # A subtype on the *same* unit is the case D24 exists for (`SoilTemperature <: Temperature`),
    # and it stays legal - only a **contradicting** unit is refused, below.
    @nicheaxis(TestSubLeaf <: TestLeaf)
    @nicheaxis(TestResource<:EcoSISTEM.NicheAxis, condition=nothing,
               resource=g/day,
               supply=Supply{CarbonFlux}, demand=Demand{CarbonFlux})
    @nicheaxis(TestBoth<:EcoSISTEM.NicheAxis, condition=mm/day,
               resource=u"L"/day,
               supply=Supply{Precipitation}, demand=Demand{Precipitation})
    @nicheaxis(TestReference <: EcoSISTEM.NicheAxis, reference)

    # Every axis is abstract, so a leaf can gain children later without a breaking change.
    @test all(isabstracttype,
              (TestGroupAxis, TestLeaf, TestSubLeaf, TestReference))

    # A leaf with no declarations inherits the group's - the case that needs `Type{<:A}`, not
    # `Type{A}`, and the reason the reachability testset below exists.
    @test EcoSISTEM.canonicalunit(TestLeaf) == K
    @test EcoSISTEM.bounds(TestLeaf) == (0.0K, nothing)
    # ...and inheritance reaches an arbitrary depth, not just one level.
    @test EcoSISTEM.canonicalunit(TestSubLeaf) == K

    # A resource does **not** imply "no condition unit" - only `condition = nothing` says that.
    @test isnothing(EcoSISTEM.canonicalunit(TestResource))
    @test EcoSISTEM.canonicalunit(EcoSISTEM.Resource, TestResource) == g / day
    @test EcoSISTEM.supplytype(TestResource) === Supply{CarbonFlux}
    @test EcoSISTEM.demandtype(TestResource) === Demand{CarbonFlux}
    # ...and an axis can be both, which is what would break if `resource` forced `nothing`.
    @test EcoSISTEM.canonicalunit(TestBoth) == mm / day
    @test EcoSISTEM.canonicalunit(EcoSISTEM.Resource, TestBoth) == u"L" / day

    # `reference` states "neither", which is different from having said nothing.
    @test isnothing(EcoSISTEM.canonicalunit(TestReference))
    @test isnothing(EcoSISTEM.supplytype(TestReference))

    # **An axis may not CONTRADICT the unit it inherits**, and this is a by-construction guard
    # rather than a runtime check. Why it matters: `typeintersect` promotion rewraps a layer onto a
    # more specific axis **without converting**, so a group declaring `K` with a `K*day` leaf beneath
    # it would let degree-day data be relabelled as a temperature. And there is no meaningful
    # conversion to do instead - `K` and `K*day` are different *quantities*, not one quantity in two
    # units, exactly as a ratio and a count are.
    # This is why `TemperatureAxis` is now a pure grouping node: it was the only axis in the
    # shipped hierarchy declaring a unit its own leaves overrode.
    # `@eval`, because the macro emits an `abstract type` and the guard fires while that block
    # runs - a bare call inside the testset would be an expression, not a top-level declaration.
    msg = try
        @eval @nicheaxis(TestBadChild<:TestGroupAxis, condition=K * day)
        ""
    catch e
        sprint(showerror, e)
    end
    @test !isempty(msg)
    # The message must name both units and point at the remedy, not merely refuse.
    @test occursin("TestGroupAxis", msg)
    @test occursin("grouping node", msg)

    # Declarations are inherited **one at a time**, not all-or-nothing: declaring `bounds` does not
    # stop `condition` being inherited. This is what makes a group worth having.
    @nicheaxis(TestOnlyBounds <: TestGroupAxis, bounds = (5.0K, 9.0K))
    @test EcoSISTEM.canonicalunit(TestOnlyBounds) == K       # still the group's
    @test EcoSISTEM.bounds(TestOnlyBounds) == (5.0K, 9.0K)   # its own

    # ...but `reference` must **shed** the group's bounds, not just be refused its own: it has no
    # canonical unit left to state them in, and `_enforcebounds!` would compare a reference layer
    # against the group's `0.0K`. Inheriting here was a real bug.
    @nicheaxis(TestRefUnderGroup <: TestGroupAxis, reference)
    @test isnothing(EcoSISTEM.canonicalunit(TestRefUnderGroup))
    @test EcoSISTEM.bounds(TestRefUnderGroup) == (nothing, nothing)
end

# The expansion-time validations. A macro sees *expressions*, not values, so these are the only
# checks it can make - anything needing the catalogue or a dimension lives in a test or at first use.
@testset "@nicheaxis refuses malformed declarations" begin
    # The one that matters most: the three resource declarations are all-or-nothing, which is what
    # makes the "resource unit with no supply type" drift unrepresentable rather than merely absent.
    @test_throws LoadError @eval @nicheaxis(BadPartial <: EcoSISTEM.NicheAxis,
                                            resource = g / day)
    @test_throws LoadError @eval @nicheaxis(BadPartial2<:EcoSISTEM.NicheAxis,
                                            resource=g/day,
                                            supply=Supply{CarbonFlux})

    # `reference` is a positive statement of "neither", so it contradicts both roles...
    @test_throws LoadError @eval @nicheaxis(BadRef <: EcoSISTEM.NicheAxis,
                                            reference, condition = K)
    # ...and has no unit for a bound to be expressed in.
    @test_throws LoadError @eval @nicheaxis(BadRefBounds <: EcoSISTEM.NicheAxis,
                                            reference, bounds = (0.0, nothing))

    # Typos and malformed heads, which is what a user actually hits.
    @test_throws LoadError @eval @nicheaxis(BadOpt <: EcoSISTEM.NicheAxis,
                                            units = K)
    @test_throws LoadError @eval @nicheaxis(BadDup<:EcoSISTEM.NicheAxis,
                                            condition=K,
                                            condition=mm)
    @test_throws LoadError @eval @nicheaxis BadHead
end

@testset "canonicalunit" begin
    @testset "SolarRadiationAxis" begin
        # WorldClim's srad is kJ*m^-2*day^-1, CHELSA's rsds*/BioClimPlus's rsds_* are MJ*m^-2*day^-1 - the
        # override reconciles both to one scale (see `_tocanon`/Layer.jl).
        @test EcoSISTEM.canonicalunit(SolarRadiation) == kJ / (m^2 * day)
        @test EcoSISTEM.canonicalunit(SolarRadiationRange) == kJ / (m^2 * day)
    end

    @testset "CarbonFlux is resource-only" begin
        # `npp` is potential productivity - a resource species compete for, not a condition they
        # are matched against - so `CarbonFlux` declares a `Resource` unit and **deliberately no
        # `Condition` one**. Pinned here so the axis-wide canonical-unit audit (which does have
        # genuine omissions to fill in - `SnowWaterEquivalent`, `WindSpeed`, the `*Range` family)
        # cannot "fix" this one by reflex: it is a decision, not a gap.
        @test EcoSISTEM.canonicalunit(EcoSISTEM.Resource, CarbonFlux) ==
              g / day
        @test isnothing(EcoSISTEM.canonicalunit(CarbonFlux))
        @test isnothing(EcoSISTEM.canonicalunit(EcoSISTEM.Condition,
                                                CarbonFlux))
    end
end

# Canonicalising a regime value is the one place a wrongly-united number could enter a layer, so
# every route in is checked. A unitful value converts (and a dimension mismatch throws); a **bare
# number carries nothing to check**, so it is refused rather than silently stamped.
#
# This closed the last hole in that net: `_checkaxisunit` already covered every shipped catalogue
# row and every user `axis =` override, but a bare `Real` bypassed both. Measured before the
# change: no test, canonical run, docs example or example script reached the branch, so refusing
# cost nothing.
@testset "a regime value must carry a unit" begin
    # A unitful value converts, affine units included.
    @test EcoSISTEM._canonical(285.0K, Temperature) == 285.0K
    @test EcoSISTEM._canonical(12.0u"°C", Temperature) ≈ 285.15K
    # A dimension mismatch is caught by the conversion itself.
    @test_throws Unitful.DimensionError EcoSISTEM._canonical(285.0mm,
                                                             Temperature)
    # The point of the testset: a bare number on a united axis is an error, not 285 K. It could
    # equally have meant 285 °C, and nothing in the value says which.
    @test_throws ErrorException EcoSISTEM._canonical(285.0, Temperature)
    # But an axis with *no* canonical unit still takes a bare number - there is nothing to
    # disagree with, and `EcoSISTEM.NicheAxis` layers are dimensionless by construction.
    @test EcoSISTEM._canonical(0.5, EcoSISTEM.NicheAxis) == 0.5
end

# The three `_tocanon` branches are three different *statements*, pinned here because they share an
# implementation and so could be collapsed back together by someone tidying.
@testset "the three canonicalisation cases stay distinct" begin
    # 1. No condition unit at all: passes through. **Must stay permissive** - reference layers are
    #    exactly the layers with no canonical unit. `CarbonFlux` states `condition = nothing`.
    @test EcoSISTEM._canonical(0.5, CarbonFlux) == 0.5
    @test EcoSISTEM._canonical(285.0K, CarbonFlux) == 285.0K

    # 2. Declared dimensionless - and now **strict**, which is the mirror of the bare-number hole
    #   and was the last piece of it. It could not be closed while `EcoSISTEM.NicheAxis` also said
    #    `NoUnits`: the default axis for a layer built without one had to keep accepting `285.0K`,
    #    since unit-bearing unclassified layers are ordinary. Now `EcoSISTEM.NicheAxis` says `nothing`
    #    instead, so the two statements are separable and only the *declared* one is tightened.
    @test EcoSISTEM._canonical(0.5, Isothermality) == 0.5
    @test EcoSISTEM._canonical(64.85u"percent", Isothermality) ≈ 0.6485
    @test_throws "must be dimensionless" EcoSISTEM._canonical(285.0K,
                                                              Isothermality)
    # ...while the axis-less case is untouched, which is the whole point of separating them: no
    #    axis means nothing to disagree with the unit, so any unit is legitimate.
    @test EcoSISTEM._canonical(285.0K, EcoSISTEM.NicheAxis) == 285.0K
    @test EcoSISTEM._canonical(0.5, EcoSISTEM.NicheAxis) == 0.5

    # 3. A real unit converts, and an affine one is absolutised on the way.
    @test EcoSISTEM._canonical(12.0u"°C", Temperature) ≈ 285.15K
end

# A tolerance cannot be built on an axis with no canonical unit - but the *reason* differs, and so
# does the remedy, so neither may surface as `MethodError: no method matching dimension(::Nothing)`
# from Unitful, which names neither the axis nor the tolerance.
@testset "a tolerance on an axis with no condition unit is refused by name" begin
    err = try
        EcoSISTEM.NicheTolerance(CarbonFlux, Normal, [1.0], [0.5])
        nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test !occursin("dimension(::Nothing)", err)   # the failure it replaces
    @test occursin("CarbonFlux", err)
    # `CarbonFlux` *states* `condition = nothing`, so the message must say it is a resource rather
    # than suggesting a missing declaration - told apart by which method answers.
    @test occursin("condition = nothing", err)
    @test occursin("resource", err)
end

@testset "the carbon resource family" begin
    # The axis knows its supply/demand types...
    @test EcoSISTEM.supplytype(CarbonFlux) === Supply{CarbonFlux}
    @test EcoSISTEM.demandtype(CarbonFlux) === Demand{CarbonFlux}
    # ...and the **supply** side must never answer from a value's dimension: `m/s` × cell area and
    # `mm/day` × cell area share one dimension, so no such table can tell wind from water. A
    # supply's type comes from its axis or nowhere.
    @test_throws MethodError EcoSISTEM.supplytype(typeof(1.0g / day))

    # **...and neither does the demand side.** `demandtype(typeof(1.0g/day))` must not answer
    # `Demand{CarbonFlux}`: a dimension cannot say which quantity it belongs to on either side.
    # Both sides answer from the axis or not at all, which is the whole design in one
    # pair of assertions.
    @test_throws MethodError EcoSISTEM.demandtype(typeof(1.0g / day))
    @test_throws MethodError EcoSISTEM.demandtype(typeof(1.0kg /
                                                         month_mean_duration))
    # the rate-dimension identity the Demand guard in `test_Demand.jl` checks
    @test EcoSISTEM._basedimension(typeof(1.0g / day)) == dimension(g)

    # `cancel` - a per-area carbon flux against a cell area - is what a carbon supply is built
    # through: 𝐌𝐋^-2𝐓^-1 × 𝐋^2 -> 𝐌𝐓^-1, stated in the axis's own canonical resource unit.
    # **Asks the three-argument axis form, which is the one on the live path.** Asking the
    # two-argument dimension-dispatched form instead would test a method the build path never
    # reaches: that form serves v0.4.0 compatibility only and lives in
    # `deprecations.jl`; it keeps one assertion below so the move does not leave it uncovered.
    @test EcoSISTEM.cancel(2.0g / (m^2 * day), 3.0m^2, CarbonFlux) == 6.0g / day
    @test EcoSISTEM.cancel(1.0g / (km^2 * day), 1.0km^2, CarbonFlux) ==
          1.0g / day
    @test EcoSISTEM.cancel(2.0g / (m^2 * day), 3.0m^2) == 6.0g / day
end

@testset "CircleMaskSpec" begin
    # `centre` in the wrong unit is rejected at `CircleMaskSpec` construction (not silently
    # wrong, and not a confusing error deep inside `_circle` later) - includes the same unit
    # as `radius`, since centre and radius are genuinely different physical quantities
    # (an angle vs. a length), not interchangeable
    @test_throws TypeError CircleMaskSpec(radius = 100.0km,
                                          centre = (2.0km, 2.0km))
    @test_throws TypeError CircleMaskSpec(radius = 100.0km, centre = (2.0, 2.0))
end

@testset "every axis with disagreeing shipped units has a canonicalunit override" begin
    # Catalogue-level audit, generalising the `SolarRadiation` fix: for every concrete niche axis, if the
    # shipped sources (across all tables) disagree on unit, `canonicalunit` must reconcile them - otherwise
    # a regime built from one source and a tolerance defaulted from another can silently diverge in scale
    # (exactly the bug `SolarRadiation` had). Reuses the catalogue helpers that already back
    # `layerinfo`/`layerunit` - no new registry.
    for A in EcoSISTEM._leafaxes()
        recs = EcoSISTEM.layersbyaxis(A)
        isempty(recs) && continue
        units = unique(r.unit for r in recs)
        if length(units) > 1
            @test !isnothing(EcoSISTEM.canonicalunit(A))
        end
    end
end

@testset "bounds" begin
    b = EcoSISTEM.bounds

    # Stated in the axis' CANONICAL unit - temperature is ≥ 0 in K, not °C.
    @test b(Temperature) == (0.0K, nothing)
    @test b(CumulativeHeat) == (0.0K * day, nothing)
    @test b(Precipitation) == (0.0mm / day, nothing)
    # A `*Range` is a max - min, so it inherits the group's floor.
    @test b(TemperatureRange) == (0.0K, nothing)
    @test b(SolarRadiationRange) == b(SolarRadiation)

    # The bounded-above axes: percentages on the 0-100 scale the shipped layers use.
    @test b(RelativeHumidity) == (0.0, 1.0)
    @test b(CloudCover) == (0.0, 1.0)
    @test b(SurfaceArea) == (0.0, 1.0)
    # ...and `Isothermality` is dimensionless, overriding its `TemperatureAxis` group.
    @test b(Isothermality) == (0.0, nothing)

    # `CarbonFlux` is stated in its *Resource* unit, having no Condition one - the shipped `npp`
    # is potential NPP (Miami model), non-negative by construction.
    @test b(CarbonFlux) == (0.0EcoSISTEM.g / day, nothing)

    # The deliberately unbounded axes.
    @test b(ClimateMoisture) == (nothing, nothing)     # a balance: P - PET
    @test b(SiteWaterBalance) == (nothing, nothing)    # likewise
    # `Altitude` is below sea level in places and is not a balance either - it simply declares no
    # bound, like every other axis that has none. It is *not* a special case: bounds are declared,
    # never derived from the catalogue, so there is nothing for it to be an exception to.
    @test b(Altitude) == (nothing, nothing)
end

@testset "bounds agree with the catalogue's `balance` marker" begin
    # The floor is the same fact the catalogue records as `Category = balance`: an axis with a
    # balance row is sign-indefinite, and is therefore not supply-eligible either. Asserted here so
    # the two statements of one fact cannot drift apart - in the style of `_checkrangerows`.
    CP = EcoSISTEM
    catalogue = CP._catalogue()
    balanceaxes = Set(rec.axis for rec in catalogue if rec.category == :balance)
    @test !isempty(balanceaxes)

    for axis in balanceaxes
        # A sign-indefinite axis must declare no floor.
        @test isnothing(first(EcoSISTEM.bounds(axis)))
    end

    # ...and conversely, every axis that DOES declare a floor must not be a balance - with `Altitude`
    # the one axis unbounded for a reason the column cannot express.
    for rec in catalogue
        isnothing(rec.axis) && continue
        lo = first(EcoSISTEM.bounds(rec.axis))
        isnothing(lo) && continue
        @test rec.category != :balance
    end
end

# The guard against a missed `<:` in an axis declaration. Every axis is an abstract type and
# every interface function bottoms out in a `::Type{<:EcoSISTEM.NicheAxis}` fallback, so writing
# `canonicalunit(::Type{TemperatureAxis})` instead of `::Type{<:TemperatureAxis}` does NOT raise
# "no method matching" - it silently stops applying to `Temperature`, which then answers the root
# default instead. What that breaks is *inheritance*: a declaration on a parent must reach every
# leaf beneath it.
#
# The invariant is about **which method answers**, not what it returns. An earlier version compared
# *values* against the default, and could not tell a deliberate `nothing` (`reference`, or
# `condition = nothing` on a supply-only axis) from a declaration that failed to reach the leaf - both
# read as "the default". Asking `which` separates them exactly: a deliberate refusal has its own
# method; a missed `<:` falls all the way through to the root fallback.
@testset "axis declarations are reachable from their leaves" begin
    # Every `EcoSISTEM.NicheAxis` subtype, however deep, found by walking down from the root.
    function allaxes(T = EcoSISTEM.NicheAxis, acc = Type[])
        push!(acc, T)
        for S in subtypes(T)
            allaxes(S, acc)
        end
        return acc
    end

    axes = allaxes()
    # A guard on the guard: a broken walk must not pass by finding nothing.
    @test length(axes) > 40

    interface = (EcoSISTEM.canonicalunit, EcoSISTEM.bounds,
                 EcoSISTEM.supplytype, EcoSISTEM.demandtype)
    # The method serving the root *is* "nobody has declared anything for this axis".
    isfallback(f, A) = which(f, Tuple{Type{A}}) ===
                       which(f, Tuple{Type{EcoSISTEM.NicheAxis}})

    for f in interface, A in axes
        # No axis may fail to dispatch at all.
        @test (f(A); true)
        # If any ancestor is served by a real declaration, so must this axis be.
        P = supertype(A)
        while P <: EcoSISTEM.NicheAxis
            isfallback(f, P) || @test !isfallback(f, A)
            P = supertype(P)
        end
    end

    # ...and the check is only meaningful if declarations are genuinely *inherited* rather than
    # restated on every leaf.
    #
    # **Inheritance is by an explicit delegating method, not by method-table fall-through**, so this
    # asserts delegation rather than requiring a leaf to be answered by its group's *own* method.
    # `@nicheaxis` emits a delegating method for an axis that declares nothing, because that is what
    # lets the root fallback **error** for a hand-rolled axis: with fall-through, a pure grouping
    # node and
    # an undeclared type were indistinguishable, both having no method of their own.
    # So every axis has its own method now, and what inheritance means is that the leaf **answers
    # its group's value** without restating it - which is what these assert.
    # **Inheritance still happens, and this is the shipped case that proves it**: `SolarRadiation`
    # declares no `condition`, and answers the unit its group declares.
    @test EcoSISTEM.canonicalunit(SolarRadiation) ==
          EcoSISTEM.canonicalunit(EcoSISTEM.SolarRadiationAxis)
    @test !isnothing(EcoSISTEM.canonicalunit(SolarRadiation))
    # A group that declares nothing delegates to *its* parent, all the way to the root.
    @test isnothing(EcoSISTEM.canonicalunit(EcoSISTEM.TemperatureAxis))
    @test isnothing(EcoSISTEM.canonicalunit(EcoSISTEM.NicheAxis))
    # ...and a leaf that declares its own is unaffected by its group saying nothing.
    @test EcoSISTEM.canonicalunit(Temperature) === EcoSISTEM.K
end

# **The regression test for `[TF-DENSITY-SCALE]`: an axis's canonical unit must NOT change the
# dynamics.** Before this, a stripped `pdf` was "probability mass in a window one canonical unit
# wide", so re-declaring an axis from `mm/d` to `cm/d` multiplied every suitability by 10 and moved
# every equilibrium (`b/d` by `10^(2*survival)` - 100× at `survival = 1.0`).
# Two axes declared here at **different canonical units for the same physical quantity**, each with
# the *same* fixed physical `densitywidth`, must give the **same** suitability for the same niche.
@testset "an axis's canonical unit does not change suitability" begin
    @nicheaxis(TestMM<:EcoSISTEM.NicheAxis, condition=mm/day,
               densitywidth=1.0mm/day)
    @nicheaxis(TestCM<:EcoSISTEM.NicheAxis, condition=u"cm/d",
               densitywidth=1.0mm/day)

    # the SAME physical niche, expressed on each axis
    tmm = NicheTolerance(TestMM, Normal, [50.0mm / day], [5.0mm / day])
    tcm = NicheTolerance(TestCM, Normal, [50.0mm / day], [5.0mm / day])
    # the distributions genuinely differ, which is what must not leak into the answer
    @test params(getdist(tmm, 1)) != params(getdist(tcm, 1))

    smm = NicheSuitability(tmm)(getdist(tmm, 1), 50.0mm / day)
    scm = NicheSuitability(tcm)(getdist(tcm, 1), 50.0mm / day)
    @test smm ≈ scm                       # the assertion the defect broke
    # ...and the scale is what does it: 1.0 in the mm frame, 0.1 in the cm frame
    @test EcoSISTEM._densityscale(TestMM, typeof(1.0mm / day)) == 1.0
    @test EcoSISTEM._densityscale(TestCM, typeof(1.0u"cm/d")) ≈ 0.1

    # **Specialist advantage survives**, which is the requirement that kills peak-normalising:
    # the width is species-independent, so a narrower niche still peaks higher.
    narrow = NicheTolerance(TestMM, Normal, [50.0mm / day], [1.0mm / day])
    @test NicheSuitability(narrow)(getdist(narrow, 1), 50.0mm / day) > smm

    # An axis declaring no width is unscaled, exactly as before.
    @test isnothing(EcoSISTEM.densitywidth(EcoSISTEM.NicheAxis))
    @test EcoSISTEM._densityscale(EcoSISTEM.NicheAxis, Float64) == 1.0

    # **A DIMENSIONLESS axis's width is a bare `Float64`, not a `Quantity`.** `1.0NoUnits`
    # collapses to `1.0`, so `_asscale` needs a `Real` method; without one every continuous niche on
    # one of the eight dimensionless axes died with a `MethodError` in the hot loop. The tests
    # above cannot catch it - they declare unitful widths only, so they pass either way.
    @test EcoSISTEM.densitywidth(EcoSISTEM.SurfaceArea) isa Real
    @test !(EcoSISTEM.densitywidth(EcoSISTEM.SurfaceArea) isa Unitful.Quantity)
    @test EcoSISTEM._densityscale(EcoSISTEM.SurfaceArea, Float64) == 1.0
    # ...and it has to survive the route that actually broke: evaluating a niche on such an axis.
    surf = NicheTolerance(EcoSISTEM.SurfaceArea, Normal, [0.5], [0.1])
    @test NicheSuitability(surf)(getdist(surf, 1), 0.5) > 0.0
end

end
