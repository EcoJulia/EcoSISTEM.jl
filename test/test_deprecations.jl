# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDeprecations

using EcoSISTEM
using RasterDataSources
using Test
using Unitful

# Coverage for `src/deprecations.jl` and the dataset-typed shims in the RasterDataSources
# extension, which hold none at present. A shim is checked on both halves: it warns
# (`@test_deprecated`) and its result matches the current API it forwards to.

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
