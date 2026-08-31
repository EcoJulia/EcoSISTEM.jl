# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Is every type declared after everything its declaration needs? Run with the other hygiene checks:
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["extras_clean.jl"])'
#
# Must go through `Pkg.test`: without the weak deps the extension types are invisible and the
# audit covers less than it claims - the same constraint `clean_Architecture.jl` has.
#
# **What this replaces.** The include order in `src/EcoSISTEM.jl` is a real invariant - a struct's
# supertype, field types and parameter bounds must all exist when it is declared - and it was
# maintained by hand, a comment explaining which file sits above which could only ever *describe*
# the order. This asserts it.

module TestCleanTypeOrder

using EcoSISTEM
using Test

include(joinpath(pkgdir(EcoSISTEM), "docs", "typeorder.jl"))
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

    # **The Phase A progress meter, as a RATCHET.** Every `src/_*.jl` is a file whose
    # declarations have not yet been redistributed into the concept-named files
    # (`docs/type_partition.md`). **Lower this number as each file moves** - that is what stops
    # the work silently going backwards, and it is the same discipline as the named exception lists
    # in `test_EcoSISTEM.jl`: a bound that has to shrink, not a count that can drift.
    # **When it reaches 0, replace this with `@test isempty(report.legacy)`** - that is the
    # completion test for Phase A.
    # `src/deprecations.jl` is deliberately unprefixed and never counts: it is the one file the
    # reorganisation does not touch.
    # **Value constants are NOT counted** - `report.legacyvalues` holds those separately. A
    # `const` that specialises no type (`_WINDOW_PAD`, `px`, `_SUPPLY_SIZE`) belongs to the functions
    # that read it and travels in Phase B, so counting it here would give this ratchet a floor it
    # could never reach.
    # **PHASE A IS COMPLETE (2026-08-21): this is now the completion test, not a ratchet.**
    # Every type and every type alias the package owns lives in a concept-named file.
    @test isempty(report.legacy)

    # **And the detector is not vacuous** - the failure it guards against is an audit that checks
    # nothing and passes. A parser change that stopped finding declarations, or a module list that
    # stopped finding types, would otherwise read as success.
    @test report.ntypes > 200
    @test report.nedges > 200
    @test report.nfiles > 50
end

end
