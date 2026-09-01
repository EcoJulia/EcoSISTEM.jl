# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Relax Rasters' pre-read memory guard when the tests are running on a CI runner.
#
# Rasters refuses a read whose result would not fit in memory, and decides that with
# `Sys.free_memory() > bytes`. On Linux that is a reasonable proxy; on macOS it is not, because
# `Sys.free_memory()` reports only genuinely *free* pages and macOS keeps almost none - the rest is
# file cache, purgeable and compressed pages, all of which are reclaimable on demand. A GitHub macOS
# runner reported 54 MB free while refusing a 60 MB read, with ParallelTestRunner having already
# dropped to a single worker, so this is not a symptom of the suite asking for too much at once.
#
# ParallelTestRunner faces the same question when it decides how many workers to start, and answers
# it correctly: its `available_memory()` adds the inactive, purgeable and compressor page counts on
# Darwin rather than trusting the free count. Rasters has no such correction, so its guard fires on
# reads that are three orders of magnitude smaller than the memory actually available.
#
# What is given up: a read that genuinely does not fit is now killed by the operating system rather
# than refused with a clear message. That trade is only sound because the suite bounds its own
# footprint by other means -- `core_test.jl` runs the raster-heavy files one at a time, and
# `test/canonical/canonical.jl`'s `heavydata()` keeps the largest CHELSA layers off runners
# entirely. Weakening either of those without revisiting this is what would make it dangerous.
#
# Deliberately keyed on `RUNNER_OS`, matching `heavydata()` and `runtests.jl`'s Windows skip: a
# developer running the suite by hand keeps the guard, because the reason to relax it is that a
# runner's memory accounting misreports, not that the check is unwanted.
#
# Held as an expression as well as being evaluated, so that ParallelTestRunner can run it in each
# worker through `init_worker_code`. The guard against redefinition is needed because more than one
# test set includes this file, and each is also documented as runnable on its own.
if !@isdefined(RELAXRASTERMEMCHECK)
    const RELAXRASTERMEMCHECK = quote
        if haskey(ENV, "RUNNER_OS")
            import Rasters
            Rasters.checkmem!(false)
        end
    end
end

eval(RELAXRASTERMEMCHECK)
