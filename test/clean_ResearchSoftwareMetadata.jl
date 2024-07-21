# SPDX-License-Identifier: LGPL-3.0-or-later

module CleanRSMD
using Test
using Phylo
using Git
using Logging
using Pkg
using ResearchSoftwareMetadata

function is_repo_clean(repo_path; ignore_untracked = true)
    # Get the status of the repository
    statuses = readlines(`$(Git.git()) status -s $repo_path`)

    if ignore_untracked
        # Repo must be clean except for untracked files
        statuses = filter((!) ∘ contains("??"), statuses)
    end

    is_clean = isempty(statuses)

    # If not clean then report on dirty files
    is_clean || @error "\n" * join(statuses, "\n")

    return is_clean
end

# Metadata crosswalk testing only works on Julia v1.8 and after due to Project.toml changes
# Also does not currently work on Windows runners on GitHub due to file writing issues
if VERSION ≥ VersionNumber("1.8.0") &&
   (!haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows")
    @testset "RSMD" begin
        git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
        @test isnothing(ResearchSoftwareMetadata.crosswalk())
        global_logger(SimpleLogger(stderr, Logging.Warn))
        @test_nowarn ResearchSoftwareMetadata.crosswalk()
        @test is_repo_clean(git_dir)
    end
else
    @test_broken VERSION ≥ VersionNumber("1.8.0") &&
                 (!haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows")
end

end
