#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0000     # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --output=%x-%j.out        # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err         # error file name will contain job name + job ID
#SBATCH --partition=nodes         # which partition to use, default on MARS is “nodes" (or "smp", "gpu", "gpuplus")
#SBATCH --time=0-12:00:00         # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=256G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=4                 # number of nodes to allocate, default is 1
#SBATCH --ntasks=32               # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=8         # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=8       # number of tasks to be launched on each allocated node
#SBATCH --threads-per-core=1      # Threads per core

############# LOADING MODULES (optional) #############
module load apps/julia
julia --project=examples -e 'using Pkg; Pkg.instantiate(); Pkg.build("MPI"); using MPI; MPI.install_mpiexecjl(destdir = "bin", force = true)'

############# ENVIRONMENT #############
# Set the number of OpenMP threads to 1 to prevent
# any threaded system libraries from automatically
# using threading. Then manually set Julia threads
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=8

############# MY CODE #############
bin/mpiexecjl --project=examples -n 32 julia -t 8 --project=examples examples/HPC/MPIRun.jl
