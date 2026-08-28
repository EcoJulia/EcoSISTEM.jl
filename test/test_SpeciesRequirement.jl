# SPDX-License-Identifier: LGPL-3.0-or-later

module TestSpeciesRequirement

using EcoSISTEM
using EcoSISTEM: AbstractSpeciesRequirement, AbstractTolerance, AbstractDemand,
                 Condition, Resource
using EcoSISTEM.Units
using Distributions: Normal
using Unitful, Unitful.DefaultSymbols
using Test

const NSP = 5
_tol() = NicheTolerance(Temperature, Normal, fill(298.0K, NSP), fill(2.0K, NSP))
function _rain()
    return NicheTolerance(Precipitation, Normal, fill(3.0mm / day, NSP),
                          fill(1.0mm / day, NSP))
end
_dem() = Demand{SolarRadiation}(fill(450000.0kJ / m^2 / day * 1.0m^2, NSP))

# **The guarantee the merge had to preserve, and it is the whole reason the aliases are
# parametric.** Before `AbstractTolerance` and `AbstractDemand` were merged they were *unrelated*
# types, so neither could stand in for the other. A single family keyed on `Role` keeps that only if
# each alias pins its role — a non-parametric alias, or one that forgot the role, would let a demand
# satisfy a tolerance annotation with nothing to notice.
@testset "the role-pinned aliases still refuse each other" begin
    tol, dem = _tol(), _dem()
    @test tol isa AbstractTolerance
    @test dem isa AbstractDemand
    @test !(dem isa AbstractTolerance)
    @test !(tol isa AbstractDemand)
    # …at the level of types, which is what a `where` clause and a struct field see.
    @test typeof(tol) <: AbstractTolerance
    @test !(typeof(dem) <: AbstractTolerance)
    @test !(typeof(tol) <: AbstractDemand)
    # Both are the one family, which is what lets anything generic be written over the pair.
    @test typeof(tol) <: AbstractSpeciesRequirement
    @test typeof(dem) <: AbstractSpeciesRequirement
    # And the alias carries the **axis and** the eltype as well as the role, so it is usable in
    # a `where`. The axis slot came first (D3(a)) and is what the three families now share.
    @test eltype(tol) == typeof(1.0K)
    @test typeof(tol) <: AbstractTolerance{Temperature, typeof(1.0K)}
    @test !(typeof(tol) <: AbstractTolerance{Temperature, Float64})
    # …and a different axis is refused even at the same eltype, which is the whole point
    @test !(typeof(tol) <: AbstractTolerance{Precipitation, typeof(1.0K)})
end

# The role is *read off* the members, exactly as `LayerCollection` does it, so there is no role to
# pass and none to get wrong.
@testset "a collection infers its own role" begin
    tc = SpeciesRequirementCollection((temperature = _tol(),
                                       rainfall = _rain()))
    dc = SpeciesRequirementCollection((sun = _dem(),))
    @test tc isa SpeciesRequirementCollection{Condition}
    @test dc isa SpeciesRequirementCollection{Resource}
    # …and a collection is itself a member of the family whose role it inferred, so the accessors
    # and the pairing checks reach it through the same alias a single member uses.
    @test tc isa AbstractTolerance
    @test dc isa AbstractDemand
    @test !(tc isa AbstractDemand)
    @test !(dc isa AbstractTolerance)

    @test length(values(tc)) == 2
    @test keys(tc) == (:temperature, :rainfall)
    @test length(values(dc)) == 1
    @test keys(dc) == (:sun,)
    # Members by name, never by index — see `src/collections.jl`.
    @test tc.rainfall === values(tc)[2]
    # A collection built from a plain `Tuple` is named by its members' **axes** where those are
    # distinguishable, so this pair reads the same as the named one above without the caller
    # writing anything — which is what lets a tolerance and its regime check out unaided.
    positional = SpeciesRequirementCollection((_tol(), _rain()))
    @test keys(positional) == (:Temperature, :Precipitation)
    @test positional.Precipitation === values(positional)[2]
end

# **The guarantee the merge GAINS**, and it must be asserted or it would be easy to ship without
# it firing at all. Before the merge each collection's `_checkbacking` knew only its own family, so a
# tolerance among demands was caught as *"not a demand"* and the cross-role case could not arise —
# the two collections were different types. Now they are one, so a mixed collection is representable
# and has to be refused explicitly.
@testset "a collection mixing roles is refused" begin
    @test_throws ErrorException SpeciesRequirementCollection((t = _tol(),
                                                              d = _dem()))
    @test_throws ErrorException SpeciesRequirementCollection((_dem(), _tol()))
    # And the message names both roles, so it says *what* disagreed rather than only that
    # something did.
    msg = try
        SpeciesRequirementCollection((t = _tol(), d = _dem()))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("same role", msg)
    @test occursin("Condition", msg)
    @test occursin("Resource", msg)

    # The pre-existing guards are unchanged: a non-member, and an empty collection.
    @test_throws ErrorException SpeciesRequirementCollection((_tol(), 1))
    @test_throws ErrorException SpeciesRequirementCollection(())
end

# `SpeciesList`'s two type parameters are the place the swap would actually happen, and they are
# written exactly as they were before the merge — `TL <: AbstractTolerance, DM <: AbstractDemand`.
# This asserts that spelling still rejects the swap, which is what makes the merge free.
@testset "SpeciesList still refuses the two sides swapped" begin
    tol, dem = _tol(), _dem()
    abun = fill(10, NSP)
    movement = BirthOnlyMovement(GaussianKernel.(fill(1.0km, NSP), 10e-4))
    param = EqualPop(0.1 / year, 0.1 / year, 1.0, 0.1, 1.0)
    native = fill(true, NSP)
    @test_nowarn SpeciesList(NSP, tol, abun, dem, movement, param, native)
    # The swap is a `MethodError`, not a runtime check — the signature refuses it.
    @test_throws MethodError SpeciesList(NSP, dem, abun, tol, movement, param,
                                         native)
end

end
