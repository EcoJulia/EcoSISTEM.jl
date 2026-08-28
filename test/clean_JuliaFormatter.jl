# SPDX-License-Identifier: LGPL-3.0-or-later

module CleanJuliaFormatter
using Test
using EcoSISTEM
using Git
using Pkg
using JuliaFormatter

include("GitUtils.jl")
using .GitUtils

# No Windows guard is needed here any more: `runtests.jl` skips the whole `extras_*` group on a
# Windows CI runner, so this file is never reached there. It still runs on Windows *locally*, where
# the file-writing problem that motivated the old guard does not arise.
@testset "JuliaFormatter" begin
    git_dir = readchomp(`$(Git.git()) rev-parse --show-toplevel`)
    @test_nowarn format(EcoSISTEM)
    @test is_repo_clean(git_dir, strict = haskey(ENV, "RUNNER_OS"))
end

end
