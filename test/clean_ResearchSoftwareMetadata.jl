# SPDX-License-Identifier: LGPL-3.0-or-later

module CleanRSMD
using Test
using Phylo
using Git
using Logging
using Pkg
using ResearchSoftwareMetadata

include("GitUtils.jl")
using .GitUtils

# Does not currently work on Windows runners on GitHub due to file writing issues
if !haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows"
    @testset "RSMD" begin
        git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
        @test isnothing(ResearchSoftwareMetadata.crosswalk())
        global_logger(SimpleLogger(stderr, Logging.Warn))
        @test_nowarn ResearchSoftwareMetadata.crosswalk()
        global_logger(SimpleLogger(stderr, Logging.Info))
        @test is_repo_clean(git_dir; strict = haskey(ENV, "RUNNER_OS"))
    end
else
    @test_broken !haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows"
end

end
