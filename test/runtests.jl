# Identify files in test/ that are testing matching files in src/
#  - src/Source.jl will be matched by test/test_Source.jl

filebase = map(file -> replace(file, r"(.*).jl", s"\1"),
                filter(file -> ismatch(r".*\.jl", file), readdir("../src")))
testbase = map(file -> replace(file, r"test_(.*).jl", s"\1"),
                filter(str -> ismatch(r"^test_.*\.jl$", str), readdir()))

info("Running tests for files:")
for t in testbase
    println("    = $t.jl")
end
println()

info("Running tests...")
for t in testbase
    fn = "test_$t.jl"
    println("    * Testing $t.jl ...")
    include(fn)
    println()
end
