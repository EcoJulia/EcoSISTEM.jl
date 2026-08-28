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
# receives a shutdown signal the job is ABORTED rather than failed, so NO further step runs -- the
# runner's own debug log says "Skip evaluate condition on runner shutdown" for every remaining step,
# post steps included -- and the file goes with the machine. Measured: the sampler step ran, the
# reporting step never did, and nothing was recovered.
#
# It earned its keep: it identified an out-of-memory kill on the first run that carried it, after
# several runs had shown only `signal 15` and no cause.
#
# QUIET UNTIL IT MATTERS. Sampling every five seconds and printing every reading added around 180
# lines to each job, which is a lot of noise for something that is usually flat. It now prints only
# when available memory is under the threshold below, plus the first reading and one a minute as a
# heartbeat -- so a healthy job carries a dozen lines and a job in trouble carries a dense trace of
# exactly the part that matters. Grep the log for `RES`.
set -u

# Print every reading below this many MB of available memory. Chosen from the measured failure: the
# runner died going from 5.1 GB used to 15.7 GB in twenty seconds, so a threshold near a fifth of the
# machine catches the descent with several readings to spare rather than only the last one.
readonly RES_THRESHOLD_MB=3072
readonly RES_INTERVAL_S=5
readonly RES_HEARTBEAT=12  # one in twelve, so once a minute at the interval above

(
    i=0
    while true; do
        i=$((i + 1))
        disk=$(df -Pk / | awk 'NR == 2 {printf "disk_avail=%dM", $4 / 1024}')
        # `free` is GNU-only; the fallback keeps this harmless if it is ever run on macOS, where
        # there is no cheap equivalent of "available", so those runs get the heartbeat alone.
        if command -v free > /dev/null 2>&1; then
            # One call, both answers: `avail` drives the threshold and `mem` is what gets printed.
            read -r avail mem <<< "$(free -m |
                awk '/^Mem:/ {printf "%s mem_used=%sM mem_avail=%sM", $7, $3, $7}')"
        else
            avail=""
            mem=$(vm_stat | awk '/Pages free/ {printf "pages_free=%s", $3}')
        fi
        if [ "${i}" -eq 1 ] ||
            [ $((i % RES_HEARTBEAT)) -eq 0 ] ||
            { [ -n "${avail}" ] && [ "${avail}" -lt "${RES_THRESHOLD_MB}" ]; }; then
            echo "RES $(date -u +%H:%M:%S) ${disk} ${mem}"
        fi
        sleep "${RES_INTERVAL_S}"
    done
) &
sampler=$!

# Not `exec`, because the sampler has to be stopped afterwards and the command's exit status has to
# reach the action -- an `exec` would replace this shell and lose both.
"$@"
status=$?
kill "${sampler}" 2> /dev/null || true
exit ${status}
