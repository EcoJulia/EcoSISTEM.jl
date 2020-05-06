## Benchmarking with PkgBenchmark.jl

### To benchmark example scripts locally:

#### Sequential
```
julia --project=benchmark -e '
    using Pkg; Pkg.instantiate();
    include("benchmark/runbenchmarks.jl");'
```

#### Multi-thread
```
julia --project=benchmark -e '
    using Pkg; Pkg.instantiate();
    include("benchmark/runbenchmarks_multithread.jl");'
```
