# Identify files in test/ that are testing matching files in src/
#  - src/Source.jl will be matched by test/test_Source.jl
using Compat
using Random
using Test

filebase = map(file -> replace(file, r"(.*).jl" => s"\1"),
                filter(file -> occursin(r".*\.jl", file), readdir("../src")))
testbase = map(file -> replace(file, r"test_(.*).jl" => s"\1"),
                filter(str -> occursin(r"^test_.*\.jl$", str), readdir()))

# Seed RNG to make tests reproducible
Random.seed!(1234)

@testset "Simulation.jl" begin
    @info "Running tests for files:"
    for t in testbase
        println("    = $t.jl")
    end
    println()

    @info "Running tests..."
    @testset for t in testbase
        fn = "test_$t.jl"
        println("    * Testing $t.jl ...")
        include(fn)
        println()
    end
end

@testset "Canonical tests" begin
    @info "Running canonical tests ..."
    canonical_testbase = map(file -> replace(file, r"test_(.*).jl" => s"\1"), filter(str -> occursin(r"^test_.*\.jl$", str), readdir("canonical")))
    for t in canonical_testbase
        fn = "canonical/test_$t.jl"
        println("    * Testing $t.jl ...")
        include(fn)
        println()
    end
end
