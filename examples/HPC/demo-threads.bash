#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0000     # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --output=%x-%j.out        # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err         # error file name will contain job name + job ID
#SBATCH --partition=smp           # which partition to use, default on MARS is “nodes" (or "smp", "gpu", "gpuplus")
#SBATCH --time=0-12:00:00         # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=512GB               # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1                 # number of nodes to allocate, default is 1
#SBATCH --ntasks=1                # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=128       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1       # number of tasks to be launched on each allocated node
#SBATCH --threads-per-core=1      # Threads per core

############# LOADING MODULES (optional) #############
module load apps/julia

############# ENVIRONMENT #############
# Set the number of OpenMP threads to 1 to prevent
# any threaded system libraries from automatically
# using threading. Then manually set Julia threads
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=128

############# MY CODE #############
julia -t 128 --project=examples examples/Africa_run.jl
