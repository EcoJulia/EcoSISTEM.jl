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

# No Windows guard: `runtests.jl` skips the whole `extras_*` group on a Windows CI runner, so this
# file is never reached there, and it runs on Windows locally.
@testset "RSMD" begin
    git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
    @test isnothing(ResearchSoftwareMetadata.crosswalk())
    # `@test_logs`, not `@test_nowarn`: the crosswalk's warnings are log records, which
    # `@test_nowarn` sees only as text on a redirected `stderr`, and so not at all under a logger
    # holding a stream of its own, the default one included.
    @test_logs min_level = Logging.Warn ResearchSoftwareMetadata.crosswalk()
    @test is_repo_clean(git_dir, strict = haskey(ENV, "RUNNER_OS"))
end

end
