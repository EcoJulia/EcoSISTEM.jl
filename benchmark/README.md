## Benchmark simulate in Canonical tests

The benchmark is done on the simulate call (`simulate!` or `simulate_record!`) in [canonical test files](../test/canonical/). If the precomputed EpiSystem and other inputs already exist (under the same num_thread setting), they will be loaded to run benchmark on. If not or `gen_epi` in [benchmarks.jl](benchmarks.jl) is set to `true`, the EpiSystem will be regenerated and saved for benchmarking.

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

### To get the benchmark results after running it locally:

#### Sequential
```
julia --project=benchmark -e '
    using Pkg; Pkg.instantiate();
    include("benchmark/display_results.jl");
    get_result("result.json")'
```

#### Multi-thread
```
julia --project=benchmark -e '
    using Pkg; Pkg.instantiate();
    include("benchmark/display_results.jl");
    get_result("result_multithread.json")'
```
