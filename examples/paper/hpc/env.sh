# The settings every step of the paper's HPC run shares, for `setup.bash` and `run.bash`, which
# source it from the root of the clone. Source it yourself to use the same Julia and depot by hand:
#
#     source examples/paper/hpc/env.sh
#
# The first three can be set beforehand, in `~/.bashrc` say, rather than edited here: an edited
# file in the clone makes `fair run` refuse to run.

# Where everything large lives - the Julia depot, with EcoSISTEM's asset cache and DataPipeline's
# copy of `fair` in its scratch spaces, and the local registry and its data store. On MARS, project
# storage: home's quota would not hold the data AfricaPlants reads.
export ECOSISTEM_HPC_ROOT="${ECOSISTEM_HPC_ROOT:-/mnt/data/project0000/$USER}"
# The juliaup channel to run, as `juliaup add` named it: `juliaup add 1.13.3` makes `1.13.3`.
export JULIAUP_CHANNEL="${JULIAUP_CHANNEL:-1.13.3}"
# The port the local registry listens on. Not 8000: a node is shared, and another user's server
# there would stop the registry starting.
export ECOSISTEM_HPC_PORT="${ECOSISTEM_HPC_PORT:-8030}"

# MARS provides juliaup as a module; elsewhere it is expected on the path already.
if type module > /dev/null 2>&1; then
    module load apps/juliaup/1.18.9
fi

# The trailing colon keeps Julia's own depots after this one, or every standard library recompiles.
export JULIA_DEPOT_PATH="$ECOSISTEM_HPC_ROOT/julia:"
# Where `fair registry install` puts the registry; once `fair init` has run, its configuration says.
export FAIR_REGISTRY_DIR="$ECOSISTEM_HPC_ROOT/fair/registry"
# The local registry is reached directly, never through a web proxy: MARS sets `http_proxy` on every
# node, and `fair` then takes the proxy's answer for another registry's and will not start its own.
export no_proxy="127.0.0.1,localhost${no_proxy:+,$no_proxy}"
export NO_PROXY="127.0.0.1,localhost${NO_PROXY:+,$NO_PROXY}"
# A compute node has no display: GR draws the figures without one.
export GKSwstype=100
# Keep threaded system libraries to one thread each, under Julia's own threads.
export OMP_NUM_THREADS=1
# Inside a job, as many Julia threads, and precompilation tasks, as cores were allocated.
if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
    export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
    export JULIA_NUM_PRECOMPILE_TASKS="$SLURM_CPUS_PER_TASK"
fi
