# SPDX-License-Identifier: LGPL-3.0-or-later

module ExtEcoSISTEMDataPipelineExt

using EcoSISTEM
using Test

# **This extension cannot be exercised, and this file says so rather than pretending otherwise.**
#
# `DataPipeline` is a weak dependency but is **not in the test target**, so `EcoSISTEMDataPipelineExt`
# never loads under any gate — it is neither precompiled nor run. Adding it is not currently possible
# either: `Pkg.add("DataPipeline")` into a scratch environment fails outright with unsatisfiable
# `DiskArrays` against this dependency set. See **A48** in the master plan.
#
# What that gap already cost: the extension called `assetdir(DataPipeline)` **without importing
# `DataPipeline`**, so `unziptemp` threw `UndefVarError` on every call. A trigger package is *loaded*
# when an extension activates, but its name is only in scope if the module asks — and the module
# still precompiles, because a function body resolves at call time, so nothing surfaced until the
# function was finally read.
#
# An earlier version of this file wrapped a single commented-out `@test` in a `@testset` behind two
# `if`s. That is worse than no test: it reports as a passing testset while asserting nothing.
@testset "the unziptemp hook exists in the parent, awaiting a working DataPipeline" begin
    # The parent declares the hook so the extension has something to add a method to. That much is
    # testable without the trigger package, and it is what a missing extension would break.
    @test isdefined(EcoSISTEM, :unziptemp)
    @test EcoSISTEM.unziptemp isa Function

    # Method-less until `DataPipeline` is loaded: the extension supplies the sole method, so an
    # empty method table here is the *correct* state for a run without the trigger, and a non-empty
    # one would mean the parent had grown an implementation it should not have.
    @test isempty(methods(EcoSISTEM.unziptemp))
end

end
