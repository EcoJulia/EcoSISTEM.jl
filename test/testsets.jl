# SPDX-License-Identifier: LGPL-3.0-or-later

# What the test sets share about how their files are run. `include`d by `core_test.jl`,
# `core_rasters.jl`, `core_ext.jl` and `runtests.jl`; not named `test_*.jl`, `core_*.jl` or
# `extras_*.jl`, so nothing runs it as a set of its own.

using ParallelTestRunner: parse_args

# The unit test files that read real rasters, which `core_rasters.jl` runs and `core_test.jl` leaves
# out, so that on a CI runner they have a job of their own. Measured peak RSS on a GitHub runner:
# test_rasters 4140 MB, test_datasetread 4127 MB, test_StudyArea 3933 MB, test_deprecations 3753 MB,
# against a few hundred for most of the rest. ParallelTestRunner sizes its worker pool at one per
# 2 GiB of available memory, roughly half what any of these needs, so four of them at once exhausted a
# 16 GB Linux runner and the job was killed outright.
#
# The set is those that read real rasters, not those that merely allocate: what makes them expensive
# is holding a downloaded layer, so a file joins this list when it starts naming a dataset, not when
# it gets slower. Regenerate the candidates with
# `grep -l "WorldClim{\|EarthEnv{\|CHELSA{" test/*.jl`, which over-reports - most of those touch
# only a fixture - and keep the ones that actually materialise a layer.
const RASTERTESTS = ["test_datasetread", "test_rasters", "test_StudyArea",
    "test_deprecations", "test_GridHabitat"]

# The files to run one at a time rather than beside each other: all of `names` on a CI runner, where
# memory and cores are scarce, and none on a development machine, where serialising costs real time
# (measured on one, 2 minutes became 7.5) and buys nothing. The environment variable `override`,
# when it is set, forces the question either way.
function serialonrunner(names; override = nothing)
    !isnothing(override) && haskey(ENV, override) &&
        return ENV[override] == "true" ? names : String[]
    return haskey(ENV, "RUNNER_OS") ? names : String[]
end

# ParallelTestRunner's arguments for a set: `ECOSISTEM_TEST_JOBS` workers when it is set, and its own
# choice from the cores and free memory otherwise. The macOS runner reports too little free memory
# for a second worker, which would run every file in sequence, so the workflow asks for two there.
#
# Never `parse_args(ARGS)`: `ARGS` is the `test_args` that selected the set, and ParallelTestRunner
# keeps only the tests whose names start with one of them, so forwarding them runs nothing.
function setargs()
    jobs = get(ENV, "ECOSISTEM_TEST_JOBS", "")
    return parse_args(isempty(jobs) ? String[] : ["--jobs=$jobs"])
end
