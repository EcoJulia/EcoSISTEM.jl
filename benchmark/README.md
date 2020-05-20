## Benchmark simulate in canonical tests

The benchmark is done on the simulate call (`simulate!` or `simulate_record!`) in [canonical test files](../test/canonical/). If the precomputed EpiSystem and other inputs already exist (under the same num_thread setting), they will be loaded to run benchmark on. If not or `gen_epi` in [benchmarks.jl](benchmarks.jl) is set to `true`, the EpiSystem will be regenerated and saved for benchmarking.

Please note, if you are running a newer version of the code. It's safer to regenerate a EpiSystem by setting `gen_epi` in [benchmarks.jl](benchmarks.jl) to be `true` (or simply delete any locally saved EpiSystems under `benchmark/episystem/` if it exists) for the first time. This will update any locally pre-saved EpiSystems in case they are outdated.

### To benchmark canonical scripts locally:

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

### To profile canonical scripts

If you just want to profile the `simulate` part in any canonical file, run
```
julia --project=benchmark benchmark/profile.jl <canonical_file>
```
`<canonical_file>` is the canonical test file name you want to do profiling, e.g. `test_SEIR.jl`.

The precomputed `EpiSystem` is handled in the same way as in benchmark.

We use [StatProfilerHTML.jl](https://github.com/tkluck/StatProfilerHTML.jl) to visualise the output from Julia's Profile module. It's a lightweight debugging tool and you can call `@profilehtml` to profile any function you are interested in.

For example, to profile the function that generates the `EpiSystem` in the canonical file, you can do
```
using StatProfilerHTML

# ...run the relevant code in the canonical file...

# add `@profilehtml` in front of the function you want to profile, e.g.
epi = @profilehtml EpiSystem(epilist, epienv, rel);
```
