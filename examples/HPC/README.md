# HPC testing

All of the testing here can be run directly from the root directory of a
checked out version of the [EcoSISTEM package][ecosistem-git], with the
scripts and the project environment being found in the `examples` folder and
its `examples/HPC` subfolder.

> You will probably need to install the relevant packages on the login node or
> some other node with internet access, as by default the package will attempt
> to install the package on startup.

## Multithreaded testing

An example of running a standard (multithreaded) job, in this case on a node
with 2 processors of 64 cores each. By default any run is multithreaded,
and will run in parallel on all available threads, so here the job runs with
128 threads.

```sh
sbatch examples/HPC/demo-threads.bash
```

## MPI testing

To run the code using MPI, you may need to configure it correctly. Here we use
Julia's MPI libraries (installed with the MPI package), but on HPC it is
usually the case that the HPC's own MPI libraries will be (potentially much)
faster, as they will be configured to take advantage of the exact topology and
hardware of the system.

These tests use a simple, but relatively large, example of 256 x 256 grid
containing 64k species. The MPI testing run from EcoSISTEM package directory
using Julia's built-in MPI libraries. Note that this uses the MPIRun.jl code,
and a folder is set in there (SAVEDIR) for outputs that may not be appropriate.

### Comparison of different process vs thread counts on a single node with 2 processors x 32 cores

The first example runs one task on each processor, with 32 threads per task
(one thread per core). The second runs one task per core, with each process
running single-threaded.

```sh
sbatch examples/HPC/demo-MPI-threads.bash
sbatch examples/HPC/demo-MPI-processes.bash
```

### Comparison with four nodes and a mixture of multi-threading and multi-process

```sh
sbatch examples/HPC/demo-MPI-nodes.bash
```

[ecosistem-git]: https://github.com/EcoJulia/EcoSISTEM.jl.git
