# MPI testing

These tests use a simple, but relatively large, example of 256 x 256 grid
and 128k species. The MPI testing run from EcoSISTEM package directory using
Julia's built-in MPI libraries. Note that this uses the MPIRun.jl code, and
a folder is set in there (SAVEDIR), which may not be appropriate.

## Comparison of different process vs thread counts on a single node with 2 processors x 32 cores

sbatch -J MPIRun-1x1x64 examples/HPC/MARS-demo-MPI-1x1x64.bash
sbatch -J MPIRun-1x2x32 examples/HPC/MARS-demo-MPI-1x2x32.bash
sbatch -J MPIRun-1x8x8 examples/HPC/MARS-demo-MPI-1x8x8.bash
sbatch -J MPIRun-1x64x1 examples/HPC/MARS-demo-MPI-1x64x1.bash

## Comparison with multiple nodes

sbatch -J MPIRun-4x64x1 examples/HPC/MARS-demo-MPI-4x64x1.bash
