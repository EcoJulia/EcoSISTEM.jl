#!/bin/bash --login

# PBS job options (name, compute nodes, job time)
#PBS -N Simulation20-40-1360
# Select 20 full nodes
#PBS -l select=20:ncpus=2
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
export JULIA_NUM_THREADS=34

# Launch the parallel job
#   Using 40 MPI processes and 2 MPI processes per node
mpirun -ppn 2 -n 40 julia --project=../git/Simulation/examples \
  ../git/Simulation/examples/CirrusRunUK.jl
