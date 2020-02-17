#!/bin/bash --login

# PBS job options (name, compute nodes, job time)
#PBS -N Simulation4-4-144
# Select 4 full nodes
#PBS -l select=4:ncpus=1
# Parallel jobs should always specify exclusive node access
#PBS -l place=scatter:excl
#PBS -l walltime=00:20:00

# Replace [budget code] below with your project code - ec108
#PBS -A ec108

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

# Load any required modules
module load intel-mpi-17
module load intel-compilers-17

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=36

# Launch the parallel job
#   Using 4 MPI processes and 1 MPI processes per node
mpirun -ppn 1 -n 4 julia --project=../git/Simulation/examples \
  ../git/Simulation/examples/CirrusMPIRun.jl
