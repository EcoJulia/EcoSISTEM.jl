# Run tests for SmallMPItest.jl
using MPI
using JLD
using Test

cd("examples")

# Compare 1 thread 4 processes vs. 4 threads 1 process vs. 2 threads 2 processes
ENV["JULIA_NUM_THREADS"] = 1
mpiexec(cmd -> run(`$cmd -n 4 julia SmallMPItest.jl`));

ENV["JULIA_NUM_THREADS"] = 4
mpiexec(cmd -> run(`$cmd -n 1 julia SmallMPItest.jl`));

## All answers should be the same
abuns1thread = load("Test_abuns1.jld", "abuns")
abuns4thread = load("Test_abuns4.jld", "abuns")

@test abuns1thread == abuns4thread
