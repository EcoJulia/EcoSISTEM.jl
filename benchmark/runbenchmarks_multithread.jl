using PkgBenchmark
benchmarkpkg(dirname(@__DIR__),
             BenchmarkConfig(env = Dict("JULIA_NUM_THREADS" => "2")),
             resultfile = joinpath(@__DIR__, "result_multithread.json"))
