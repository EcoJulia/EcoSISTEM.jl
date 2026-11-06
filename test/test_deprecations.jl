# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDeprecations

using EcoSISTEM
using RasterDataSources
using Test
using Unitful

# Coverage for `src/deprecations.jl` and the dataset-typed shims in the RasterDataSources
# extension, which holds none at present. A shim is checked on both halves: it warns and its
# result matches the current API it forwards to. A method's warning is a log record
# (`@test_deprecated`); a deprecated binding's is printed by the runtime to stderr, so it is caught
# with `@test_warn`, and the access is evaluated inside it rather than written in the test's own
# code, which on Julia 1.11 would warn when that code is compiled, before anything is captured.

# A kernel's positional form builds the very kernel its keyword form does, and so the same lookup.
# That both match the released code to the bit is for the canonical tests, which still build their
# kernels positionally.
@testset "kernels' positional forms forward to their keyword forms" begin
    regime = build_habitat(DefaultEcosystem(), verbosity = :silent).regime
    function samelookup(a, b)
        la, lb = EcoSISTEM.genlookups(regime, a),
                 EcoSISTEM.genlookups(regime, b)
        return (la.y, la.x, la.p) == (lb.y, lb.x, lb.p)
    end
    old = @test_deprecated GaussianKernel(2.4u"km", 1.0e-9)
    new = GaussianKernel(meandist = 2.4u"km", thresh = 1.0e-9)
    @test old == new && samelookup(old, new)
    old = @test_deprecated LongTailKernel(2.4u"km", 2.0, 1.0e-6)
    new = LongTailKernel(scale = (2 * 2.4u"km") / sqrt(pi), shape = 2.0,
                         thresh = 1.0e-6)
    @test old == new && samelookup(old, new)
end

# A list-wide movement is the vector of its species' movements, each in the mode the old type
# named, with its own kernel, the kernel's type as the family and its own flag as the dead-cell
# fate.
@testset "list-wide movements forward to each species' movement" begin
    kernels = [GaussianKernel(meandist = 1.0u"km"),
        GaussianKernel(meandist = 2.0u"km")]
    @test (@test_deprecated BirthOnlyMovement(kernels)) ==
          SeedDispersal(GaussianKernel, disperse_safely = true).(kernels)
    @test (@test_deprecated BirthOnlyMovement(kernels, [true, false])) ==
          [SeedDispersal(GaussianKernel, disperse_safely = true)(kernels[1]),
        SeedDispersal(GaussianKernel, disperse_safely = false)(kernels[2])]
    @test (@test_deprecated AlwaysMovement(kernels)) ==
          AlwaysDispersal().(kernels)
    @test (@test_deprecated NoMovement(kernels)) == fill(NoDispersal()(), 2)
    long = [LongTailKernel(scale = 1.0u"km", shape = 2.0)]
    @test (@test_deprecated BirthOnlyMovement(long)) ==
          SeedDispersal(LongTailKernel, disperse_safely = true).(long)
    @test_throws "one entry per species" EcoSISTEM._listmovement(SeedDispersal,
                                                                 kernels,
                                                                 [true])
end

# A form that took a species list's starting abundances, or read them from one, raises an error saying
# how to rewrite the call; `build_species(numspecies; ...)`'s is written from the keywords it was given.
@testset "forms reading starting abundances are removed, each saying how to rewrite it" begin
    rewrite(f) =
        try
            f()
            ""
        catch err
            sprint(showerror, err)
        end
    message = rewrite() do
        return build_species(5, tolerance = (290.0u"K", 2.0u"K"),
                             toleranceaxis = Temperature,
                             demand = 1.0e9u"kJ/d", demandaxis = SolarRadiation,
                             dispersal = 5.0u"km", abundance = 500, seed = 1)
    end
    @test occursin("`build_species(5; ...)` was removed in v0.11.0", message)
    @test occursin("build_species(SpeciesPool(5); tolerance = ToleranceSpec(Temperature, " *
                   "μ = 290.0 K, σ = 2.0 K)", message)
    @test occursin("movement = MovementSpec(InsectDispersal(), meandist = 5.0 km), seed = 1)",
                   message)
    @test occursin("build_ecosystem(BySupply(total = 500), species, habitat)",
                   message)
    message = rewrite() do
        return build_species(2, tolerance = ([290.0, 291.0]u"K", 2.0u"K"),
                             toleranceaxis = Temperature, demand = 1.0u"kJ/d",
                             demandaxis = SolarRadiation, abundance = [3, 4],
                             native = [true, false])
    end
    @test occursin("μ = TraitColumn(:μ)", message)
    @test occursin("BySupply(each = row -> abundance[row.key])", message)
    @test occursin("deactivate_species!", message)
    # The `DefaultEcosystem()` shorthand builds from a pool, refusing what only the removed form took.
    @test_throws "builds species from a `SpeciesPool` since v0.11.0" build_species(DefaultEcosystem(),
                                                                                   abundance = 10,
                                                                                   verbosity = :silent)
    @test_throws "builds species from a `SpeciesPool` since v0.11.0" build_species(DefaultEcosystem(),
                                                                                   movement = BirthOnlyMovement,
                                                                                   verbosity = :silent)
    species = build_species(DefaultEcosystem(), numspecies = 3,
                            verbosity = :silent)
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    eco = build_ecosystem(NoAbundance(), species, habitat, verbosity = :silent)
    @test_throws "`abun` and `native` arguments were removed" SpeciesList(3,
                                                                          species.tolerance,
                                                                          [
                                                                              1,
                                                                              2,
                                                                              3
                                                                          ],
                                                                          species.demand,
                                                                          species.movement,
                                                                          species.params,
                                                                          fill(true,
                                                                               3))
    @test_throws "`Ecosystem(spplist, habitat, nichefit)` was removed" Ecosystem(species,
                                                                                 habitat,
                                                                                 eco.nichefit)
    @test_throws "`Ecosystem(popfun, spplist, habitat, nichefit)` was removed" Ecosystem(identity,
                                                                                         species,
                                                                                         habitat,
                                                                                         eco.nichefit)
    @test_throws "`MPIEcosystem(spplist, habitat, nichefit)` was removed" MPIEcosystem(species,
                                                                                       habitat,
                                                                                       eco.nichefit)
    for (populate, placement) in ((populate!, "BySupply()"),
        (populate_by_tolerance!, "BySuitability()"),
        (emptypopulate!, "NoAbundance()"))
        @test_throws "`$placement` places as" populate(eco.abundances, species,
                                                       habitat, eco.nichefit,
                                                       eco.rngs)
    end
end

# Placing again by suitability is `repopulate!` with the suitability placement.
@testset "repopulate_by_tolerance! places again by suitability" begin
    species = build_species(DefaultEcosystem(), numspecies = 3,
                            verbosity = :silent)
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    eco = build_ecosystem(BySupply(total = 300), species, habitat, seed = 1,
                          verbosity = :silent)
    @test (@test_deprecated repopulate_by_tolerance!(eco)) === eco
    @test sum(eco.abundances.matrix) == 3000
    @test eco.placement == BySuitability()
end

@testset "InputRecord is ProvenanceRecord" begin
    @test Base.isdeprecated(EcoSISTEM, :InputRecord)
    @test (@test_warn "deprecated" Core.eval(@__MODULE__,
                                             :(EcoSISTEM.InputRecord))) ===
          EcoSISTEM.ProvenanceRecord
    @test Base.ispublic(EcoSISTEM, :InputRecord) &&
          !Base.isexported(EcoSISTEM, :InputRecord)
end

# Shims removed at a breaking release must leave nothing behind: a stale `export` or `public` of a
# deleted binding is only a precompile warning, and a method left on a function that stays is
# silent.
@testset "removed shims leave nothing behind" begin
    for n in (:simulate_action!, :simulate_record!, :simulate_record_diversity!,
        :retrieve_era5, :readfile, :ConstructedSpec,
        :emptyMPIgridlandscape,
        :empty_mpi_gridlandscape)
        @test !isdefined(EcoSISTEM, n)
        @test !Base.isexported(EcoSISTEM, n) && !Base.ispublic(EcoSISTEM, n)
    end
    # Methods taken from functions that stay: the six-argument `simulate!` that cached to files,
    # the five-argument demographic constructors, `build_species`'s `boost`, and `read` on the
    # dataset types.
    @test !hasmethod(simulate!,
                     Tuple{Ecosystem, Unitful.Time, Unitful.Time, Unitful.Time,
                           String, String})
    @test !hasmethod(EqualPop, NTuple{5, Any})
    @test_throws MethodError build_species(DefaultEcosystem(), numspecies = 2,
                                           boost = 1.0, verbosity = :silent)
    # Another package's generic method may still accept these (FlatBuffers has `read(::Type{T},
    # bytes)`), so the question is whether any method is the package's own.
    ext = Base.get_extension(EcoSISTEM, :EcoSISTEMRasterDataSourcesExt)
    ours(sig) = [m for m in methods(read, sig) if m.module in (EcoSISTEM, ext)]
    @test isempty(ours(Tuple{Type{WorldClim{BioClim}}, Any}))
    @test isempty(ours(Tuple{Type{CHELSA{Climate}}, String, String}))
    @test isempty(ours(Tuple{Type{EcoSISTEM.ERA}, String, String}))
    @test isempty(ours(Tuple{Type{EcoSISTEM.CRUTS}, String, String}))
end

end
