# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Does `data/architecture.md` still describe the types that exist? Run with the other hygiene checks:
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["extras_clean.jl"])'
#
# Must go through `Pkg.test`: the extension types are invisible without the weak deps, and the
# audit would then report false gaps.
#
# **When this fails, run the fixer AND READ THE PROSE** - see the protocol in `CLAUDE.md`
# ("Keeping `data/architecture.md` honest"). The diagrams are the cheap half.

module TestCleanArchitecture

using EcoSISTEM
using Test

include(joinpath(pkgdir(EcoSISTEM), "data", "src", "architecture.jl"))
using .ArchitectureAudit

@testset "Architecture" begin
    report = ArchitectureAudit.architecture_report()

    # **A name in a diagram that resolves nowhere is always wrong** - the type has been deleted or
    # renamed and the doc was not followed through. There is no legitimate case, so no exception list.
    # This is not hypothetical: it caught `Reference`, a `ClimatePref` type deleted while both the
    # diagram *and* the sentence "`ERA`/`CERA`/`CRUTS`/`Reference` remain for their data sources"
    # stayed behind.
    @test report.stale == Symbol[]

    # **A named list, not a count** - the pattern `test_EcoSISTEM.jl` already uses. A hierarchy may
    # be deliberately undiagrammed, but a *new* one has to be looked at, and documenting one shrinks
    # this list rather than passing silently.
    # Empty, and that is the point: an entry here excuses a real gap, so the list shrinking is what
    # a newly documented hierarchy looks like. `EcoSISTEMSource` was the last one, excused as
    # "internal plumbing" - until `ERA`/`CERA`/`CRUTS` joined it and made it the answer to *where did
    # this data come from*, which is a question a user does ask.
    known = Symbol[]
    @test isempty(setdiff(report.undocumented, known))

    # ...and the exception list must not rot: every entry still has to be a real, undocumented
    # hierarchy, or it is silently excusing nothing.
    @test isempty(setdiff(known, report.undocumented))

    # The inheritance edges must already be what the code says, or `--fix` has work to do. Written
    # against a throwaway `IOBuffer` so a failing run reports the count rather than printing a repair
    # plan into the test log.
    @test ArchitectureAudit.regenerate(write = false, io = IOBuffer()) == 0

    # **And the detector is not vacuous.** A check that silently stopped seeing the doc - a moved
    # file, a changed fence - would pass every assertion above while proving nothing.
    @test length(report.documented) > 150
    @test length(report.roots) > 20
    @test :AbstractLayer in report.documented
    @test :AbstractSpeciesRequirement in report.documented
    # Aliases are collected, not derived: no `supertype` walk produces `AbstractTolerance`, so if
    # this ever empties, `--fix` would start rewriting hand-written alias edges into their underlying
    # types and quietly flatten the diagrams.
    @test :AbstractTolerance in report.aliases
    @test :AbstractRegime in report.aliases
end

end
