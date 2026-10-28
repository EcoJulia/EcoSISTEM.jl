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
    # The releases whose shims the files carry, named rather than counted, so a release appearing
    # here has to be looked at and one disappearing shrinks the list instead of passing silently.
    # Empty, as there are no shims at present: a new deprecation adds its release here.
    known = String[]
    version(label) = VersionNumber(label[(length("Deprecated in v") + 1):end])
    labelsin(secs) = [m.match
                      for (_, text) in secs
                      for m in eachmatch(r"Deprecated in v\d+\.\d+\.\d+", text)]
    main = pkgdir(EcoSISTEM, "src", "deprecations.jl")
    for path in [main,
            pkgdir(EcoSISTEM, "ext", "EcoSISTEMRasterDataSourcesExt",
                   "deprecations.jl")]
        @test isfile(path)
        secs = sections(path)
        unlabelled = [n
                      for (n, text) in secs
                      if !occursin(r"Deprecated in v\d+\.\d+\.\d+", text)]
        @test isempty(unlabelled)
        # The labels must name releases it knows, so a typo cannot pass as a new one.
        labels = labelsin(secs)
        @test issubset(labels, known)
        # Newest release first, so a release's worth of shims is deleted from the foot of the file:
        # the labels must never get newer going down.
        @test issorted(version.(labels), rev = true)
    end
    # The main file carries every release; the extension's may carry a subset.
    @test Set(labelsin(sections(main))) == Set(known)

    # The detector is not vacuous: of three sections, two labelled, it finds all three and the
    # unlabelled one, and reads the labels in order.
    mktemp() do path, io
        rule = "# " * "-"^75
        write(io,
              join([rule, "# One", "#", "# Deprecated in v0.9.0.", rule,
                       "f() = 1", "",
                       rule, "# Two", "#", "# Deprecated in v0.8.0.", rule,
                       "g() = 2", "",
                       rule, "# Three, unlabelled", rule, "h() = 3", ""], "\n"))
        close(io)
        secs = sections(path)
        @test length(secs) == 3
        @test count(t -> !occursin("Deprecated in v", last(t)), secs) == 1
        @test labelsin(secs) == ["Deprecated in v0.9.0", "Deprecated in v0.8.0"]
    end
end
end
