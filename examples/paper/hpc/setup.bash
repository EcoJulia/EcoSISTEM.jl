#!/bin/bash -l
#
# Set up the paper's HPC runs, or bring a set-up clone to the newest releases. It updates the
# juliaup channel, resolves `examples/paper`'s environment at the newest versions and precompiles
# it on the kind of node that runs it. Then `setup.jl` installs DataPipeline's copy of `fair` and
# the local registry, and runs the clone's `fair init`. Submit from the root of the clone, after
# filling in `$ECOSISTEM_HPC_ROOT/fair/user.yaml` (see `fair-user.yaml`):
#
#     sbatch examples/paper/hpc/setup.bash
#
# The `#SBATCH` lines are for MARS's `smp` nodes and project `project0000`: edit them for yours, or
# override them on the command line. See `examples/paper/README.md`.

############# SLURM SETTINGS #############
#SBATCH --account=project0000     # the project the job is charged to
#SBATCH --job-name=paper-setup    # the logs are <job name>-<job ID>.out and .err
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=smp           # the partition the runs use, so the compiled code suits its CPUs
#SBATCH --time=0-02:00:00         # d-hh:mm:ss
#SBATCH --mem=32G                 # per node
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16        # precompilation tasks

############# ENVIRONMENT #############
if [ ! -f examples/paper/hpc/env.sh ]; then
    echo "Submit this from the root of the EcoSISTEM clone, not $PWD" >&2
    exit 1
fi
source examples/paper/hpc/env.sh
set -euo pipefail
mkdir -p "$ECOSISTEM_HPC_ROOT/fair"

############# SET-UP #############
# The channel's newest Julia, then the newest versions of the environment's packages.
juliaup update "$JULIAUP_CHANNEL"
julia --project=examples/paper -e 'using Pkg; Pkg.update()'
julia --project=examples/paper examples/paper/hpc/setup.jl
