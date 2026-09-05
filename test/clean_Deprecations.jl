# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Does every deprecation section say which release it belongs to? Run with the other hygiene checks:
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["extras_clean.jl"])'
#
# Why this exists. `src/deprecations.jl` is grouped by concept rather than by version, because that
# is how a reader arrives at it - from a name they were using. Removing a release's shims is then a
# matter of deleting the sections carrying that release's label, which only works while every
# section has one. A new section without a label is invisible until someone tries to do that.
#
# It is fast, needs no network, and reads the files as text rather than loading anything.

module TestCleanDeprecations

using EcoSISTEM
using Test

# Every `# ----` ... `# ----` comment block that opens a section, as (line number, text).
function sections(path)
    lines = readlines(path)
    isrule(l) = occursin(r"^# -{20,}$", l)
    found = Tuple{Int, String}[]
    i = 1
    while i <= length(lines)
        if isrule(lines[i])
            j = i + 1
            while j <= length(lines) && startswith(lines[j], "#") &&
                  !isrule(lines[j])
                j += 1
            end
            if j <= length(lines) && isrule(lines[j]) && j > i + 1
                push!(found, (i, join(lines[i:j], "\n")))
                i = j + 1
                continue
            end
        end
        i += 1
    end
    return found
end

@testset "Deprecation sections are labelled by release" begin
    path = pkgdir(EcoSISTEM, "src", "deprecations.jl")
    @test isfile(path)
    secs = sections(path)
    # Not vacuous: the file has many sections, and a detector finding none would pass everything.
    @test length(secs) > 20

    unlabelled = [n
                  for (n, text) in secs
                  if !occursin(r"Deprecated in v\d+\.\d+\.\d+", text)]
    @test isempty(unlabelled)

    # The labels must name releases that exist, so a typo cannot pass as a new one. Named rather
    # than counted: a release appearing here has to be looked at, and one disappearing shrinks the
    # list instead of passing silently.
    labels = Set(m.match
                 for (_, text) in secs
                 for m in eachmatch(r"Deprecated in v\d+\.\d+\.\d+", text))
    @test labels == Set(["Deprecated in v0.5.0", "Deprecated in v0.6.0",
                  "Deprecated in v0.7.0"])
end

end
