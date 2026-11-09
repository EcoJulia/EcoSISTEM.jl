# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Is every type declared after everything its declaration needs? Run with the other hygiene checks:
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["extras_clean.jl"])'
#
# Must go through `Pkg.test`: without the weak deps the extension types are invisible and the
# audit covers less than it claims - the same constraint `clean_Architecture.jl` has.
#
# **Why a test.** The include order in `src/EcoSISTEM.jl` is a real invariant - a struct's
# supertype, field types and parameter bounds must all exist when it is declared - and a comment
# saying which file sits above which can only *describe* the order. This asserts it.

module TestCleanTypeOrder

using EcoSISTEM
using Test

include(joinpath(pkgdir(EcoSISTEM), "data", "src", "typeorder.jl"))
using .TypeOrderAudit

@testset "Type order" begin
    report = TypeOrderAudit.typeorder_report()

    # **The invariant.** A type declared before something its own declaration needs would not
    # load at all, so this cannot fail on a working package - it fails on a *proposed* order, which
    # is the point: it is what makes moving a declaration between files a checkable operation
    # rather than a hopeful one.
    @test isempty(report.violations)

    # **No exception list, deliberately.** A type whose declaration cannot be located is not a
    # tolerable gap - it is a type the check silently says nothing about. This is not
    # hypothetical: 51 types (every niche axis, via `@nicheaxis`, and both SimpleTraits markers, via
    # `@traitdef`) were invisible until the parser learned to read macro declarations, and the audit
    # reported a clean bill of health the whole time.
    @test isempty(report.unplaced)

    # **Every type and every type alias the package owns lives in a concept-named file**: none is
    # declared in a `src/_*.jl` (`src/deprecations.jl` is unprefixed and never counts). Value
    # constants are not counted - `report.legacyvalues` holds those separately - since a `const`
    # that specialises no type belongs with the functions that read it.
    @test isempty(report.legacy)

    # **And the detector is not vacuous** - the failure it guards against is an audit that checks
    # nothing and passes. A parser change that stopped finding declarations, or a module list that
    # stopped finding types, would otherwise read as success.
    @test report.ntypes > 200
    @test report.nedges > 200
    @test report.nfiles > 50
end

end
