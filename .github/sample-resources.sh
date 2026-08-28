#!/usr/bin/env bash
#
# Start a disk and memory sampler, run the command given, then stop it.
#
# Passed to julia-actions/julia-runtest as its `prefix`, so it wraps the julia invocation and the
# sampler's output goes to THE TEST STEP'S OWN STDOUT. That is the whole point: GitHub streams a
# step's output to the log as it is produced, so every reading is already in the log before the
# machine is gone.
#
# An earlier attempt wrote the readings to a file under RUNNER_TEMP and printed them from a later
# step with `if: always()`. That does not work, and the reason is worth keeping: when the runner
# receives a shutdown signal the job is aborted rather than failed, so NO further step runs -- not an
# `always()` one, not a post step -- and the file goes with the machine. Measured: the sampler step
# ran, the reporting step never did, and nothing was recovered.
#
# Lines are prefixed `RES` so they can be picked out of the interleaved test output with grep.
set -u

(
    while true; do
        disk=$(df -Pk / | awk 'NR == 2 {printf "disk_avail=%dM", $4 / 1024}')
        # `free` is GNU-only; the fallback keeps this harmless if it is ever run on macOS.
        if command -v free > /dev/null 2>&1; then
            mem=$(free -m | awk '/^Mem:/ {printf "mem_used=%sM mem_avail=%sM", $3, $7}')
        else
            mem=$(vm_stat | awk '/Pages free/ {printf "pages_free=%s", $3}')
        fi
        echo "RES $(date -u +%H:%M:%S) ${disk} ${mem}"
        sleep 5
    done
) &
sampler=$!

# Not `exec`, because the sampler has to be stopped afterwards and the command's exit status has to
# reach the action -- an `exec` would replace this shell and lose both.
"$@"
status=$?
kill "${sampler}" 2> /dev/null || true
exit ${status}
