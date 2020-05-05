## Benchmarking with PkgBenchmark.jl

### To benchmark example scripts locally:
```
julia --project=benchmark -e '
    using Pkg; Pkg.instantiate();
    include("benchmark/runbenchmarks.jl");'
```
