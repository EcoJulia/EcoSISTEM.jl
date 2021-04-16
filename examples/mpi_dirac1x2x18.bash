#!/bin/bash
#SBATCH --job-name run1
#SBATCH --account DIRAC-DC003-CPU
#SBATCH --ntasks 36
#SBATCH --time 1:00:00
#SBATCH --mail-type ALL
#SBATCH --no-requeue
#SBATCH --partition skylake
#SBATCH --output log.txt
. /etc/profile.d/modules.sh
# module purge
# module load rhel7/default-peta4
module load julia/1.4
module load R/3.6
module load intel-mpi-2017.4-gcc-5.4.0-rjernby

#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=18

# stdbuf -oL -eL julia --project=@. ../git/Simulation/examples/CirrusMPIRun.jl

# Launch the parallel job
#   Using 2 MPI processes and 2 MPI processes per node
mpirun -ppn 2 -n 2 julia --project=../git/Simulation/examples \
  ../git/Simulation/examples/CirrusMPIRun.jl
