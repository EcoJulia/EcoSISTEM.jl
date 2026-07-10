# SPDX-License-Identifier: LGPL-3.0-or-later

module CleanJuliaFormatter
using Test
using EcoSISTEM
using Git
using Pkg
using JuliaFormatter

include("GitUtils.jl")
using .GitUtils

# Does not currently work on Windows runners on GitHub due to file writing issues
if !haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows"
    @testset "JuliaFormatter" begin
        git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
        @test_nowarn format(EcoSISTEM)
        @test is_repo_clean(git_dir; strict = haskey(ENV, "RUNNER_OS"))
    end
else
    @test_broken !haskey(ENV, "RUNNER_OS") || ENV["RUNNER_OS"] ≠ "Windows"
end

end
