#!/bin/bash -l
#
# The paper's figures and tables as a `fair run` on one node: `fair pull` registers the parameters,
# then `fair run examples/paper/paper.yaml` draws every figure and registers it, in the local
# registry that `setup.bash` made, which this job starts and stops. Submit from the root of the
# clone, once `setup.bash` has run:
#
#     sbatch examples/paper/hpc/run.bash
#
# The small-scale run first, as a test of everything above, in minutes:
#
#     ECOSISTEM_SCALE=small sbatch --time=01:00:00 examples/paper/hpc/run.bash
#
# The clone must be clean: `fair run` refuses a modified checkout. The `#SBATCH` lines are for
# MARS's `smp` nodes and project `project0000`: edit them for yours, or override them on the command
# line. See `examples/paper/README.md`.

############# SLURM SETTINGS #############
#SBATCH --account=project0000     # the project the job is charged to
#SBATCH --job-name=paper          # the logs are <job name>-<job ID>.out and .err
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=smp           # one node of 128 cores and 1 TB
#SBATCH --time=3-00:00:00         # d-hh:mm:ss
#SBATCH --mem=64G                 # per node; the full run peaks at about 25 GB on twelve threads
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128       # Julia threads
#SBATCH --threads-per-core=1

############# ENVIRONMENT #############
if [ ! -f examples/paper/hpc/env.sh ]; then
    echo "Submit this from the root of the EcoSISTEM clone, not $PWD" >&2
    exit 1
fi
source examples/paper/hpc/env.sh
set -euo pipefail

############# THE RUN #############
# Brings the environment up to date with the clone, should it have moved since `setup.bash`, and
# keeps the versions it resolved beside this job's logs: the manifest is not part of the repository.
julia --project=examples/paper -e 'using Pkg; Pkg.instantiate()'
cp examples/paper/Manifest.toml "${SLURM_JOB_NAME:-paper}-${SLURM_JOB_ID:-local}.Manifest.toml"

# `run` without `--local` checks that the commit is on the clone's remote, so the registry records a
# run that anyone can reproduce. `pull` takes `--local`, or it would look the parameters up in the
# remote registry.
julia --project=examples/paper -e '
    using DataPipeline
    fair"registry start"
    try
        fair"pull --local examples/paper/paper.yaml"
        fair"run examples/paper/paper.yaml"
    finally
        fair"registry stop"
    end'
