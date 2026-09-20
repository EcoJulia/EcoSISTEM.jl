# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Type piracy: a method on a generic this package does not own, on types it does not own either,
# changes that generic for every package loading both. Aqua's check runs over the parent and over
# each loaded extension separately, since it sees only the methods a module itself defines.
#
# Runs through `test/extras_clean.jl`; needs the test environment for `Aqua` and the extensions:
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_clean.jl"])'

module CleanAqua

using Test
using Aqua
using EcoSISTEM
using RasterDataSources
using CDSAPI
using MPI
using Phylo

# Aqua's checks, each decided once. Two duplicate gates the suite already has and are not run:
# `test_ambiguities` (`clean_Ambiguities.jl`, which names the two known cases and has a control) and
# `test_undocumented_names` (`test_EcoSISTEM.jl`, whose named list exempts the deprecation shims and
# the generated unit constants that Aqua reports). `test_persistent_tasks` passes but is deferred:
# it precompiles a package depending on this one, about 26 s, to catch a `Task` left running on
# load.

@testset "Project and module hygiene" begin
    Aqua.test_unbound_args(EcoSISTEM)
    Aqua.test_undefined_exports(EcoSISTEM)
    Aqua.test_project_extras(EcoSISTEM)
    # The parent never loads these three: `BlockArrays` serves the MPI extension and `Calculus` and
    # `Optim` the Phylo one, and an extension can load only its trigger packages and the parent's
    # dependencies, so they stay dependencies of the parent.
    Aqua.test_stale_deps(EcoSISTEM, ignore = [:BlockArrays, :Calculus, :Optim])
    Aqua.test_deps_compat(EcoSISTEM)
end

@testset "Type piracy" begin
    # The parent pirates nothing.
    Aqua.test_piracies(EcoSISTEM)

    # Every extension is checked on its own, since Aqua sees only the methods a module itself
    # defines - and from an extension's point of view the parent's functions and types are foreign
    # too, so they are declared its own: an extension exists to add methods to them. What is left is
    # a method on someone else's generic over someone else's types.
    own = Union{Function, Type}[getfield(m, n)
                                for m in (EcoSISTEM, EcoSISTEM.Units)
                                for n in names(m, all = true)
                                if isdefined(m, n) &&
        getfield(m, n) isa Union{Function, Type}]
    # The piracies an extension still carries, named: the two `Base.read` methods on
    # RasterDataSources' types are deprecated shims, kept for one release, and go with the
    # deprecations. Named rather than counted, so a new one has to be looked at and their removal
    # shrinks the list rather than passing silently. A method and its keyword wrapper count once.
    known = Dict(:EcoSISTEMRasterDataSourcesExt =>
                     Set([Type{<:RasterDataSources.RasterDataSource},
                             Type{RasterDataSources.CHELSA{RasterDataSources.Climate}}
                         ]))
    for name in (:EcoSISTEMRasterDataSourcesExt, :EcoSISTEMERAExt,
        :EcoSISTEMMPIExt, :EcoSISTEMPhyloExt)
        ext = Base.get_extension(EcoSISTEM, name)
        @test !isnothing(ext)
        found = Aqua.Piracy.hunt(ext, treat_as_own = own)
        @test all(m -> m.name === :read, found)
        # The type the pirated method dispatches on: the first positional argument, after the
        # function type (and the keyword wrapper's two leading parameters, where it has them).
        dispatched(m) = Base.unwrap_unionall(m.sig).parameters[Aqua.is_kwcall(m.sig) ?
                                                               4 : 2]
        @test Set(dispatched(m) for m in found) == get(known, name, Set())
    end
end

end
