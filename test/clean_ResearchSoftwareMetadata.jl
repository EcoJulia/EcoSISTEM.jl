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

# No Windows guard is needed here any more: `runtests.jl` skips the whole `extras_*` group on a
# Windows CI runner, so this file is never reached there. It still runs on Windows *locally*, where
# the file-writing problem that motivated the old guard does not arise.
@testset "RSMD" begin
    git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
    @test isnothing(ResearchSoftwareMetadata.crosswalk())
    global_logger(SimpleLogger(stderr, Logging.Warn))
    @test_nowarn ResearchSoftwareMetadata.crosswalk()
    global_logger(SimpleLogger(stderr, Logging.Info))
    @test is_repo_clean(git_dir, strict = haskey(ENV, "RUNNER_OS"))
end

end
