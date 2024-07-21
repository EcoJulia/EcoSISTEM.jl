# SPDX-License-Identifier: LGPL-3.0-or-later

using Random
using Test
using EcoSISTEM
using Pkg

rsmd = get(ENV, "RSMD_CROSSWALK", "FALSE")

if rsmd == "FALSE"
    # Normal testing

    # Identify files in test/ that are testing matching files in src/
    #  - src/Source.jl will be matched by test/test_Source.jl
    filebase = String[]
    for (root, dirs, files) in walkdir("../src")
        append!(filebase,
                map(file -> replace(file, r"(.*).jl" => s"\1"),
                    filter(file -> occursin(r".*\.jl", file), files)))
    end

    testbase = map(file -> replace(file, r"test_(.*).jl" => s"\1"),
                   filter(str -> occursin(r"^test_.*\.jl$", str), readdir()))

    # Identify tests with no matching file
    superfluous = filter(f -> f ∉ filebase, testbase)
    if length(superfluous) > 0
        println()
        @info "Potentially superfluous tests:"
        for f in superfluous
            println("    + $f.jl")
        end
        println()
    end

    # Identify files with no matching test
    notest = filter(f -> f ∉ testbase, filebase)
    if length(notest) > 0
        println()
        @info "Potentially missing tests:"
        for f in notest
            println("    - $f.jl")
        end
        println()
    end

    # Identify files in test/ that are testing matching files in ext/
    #  - ext/SourceExt.jl will be matched by test/ext_SourceExt.jl
    filebase = String[]
    for (root, dirs, files) in walkdir("../ext")
        append!(filebase,
                map(file -> replace(file, r"(.*).jl" => s"\1"),
                    filter(file -> occursin(r".*\.jl", file), files)))
    end

    extbase = map(file -> replace(file, r"ext_(.*).jl" => s"\1"),
                  filter(str -> occursin(r"^ext_.*\.jl$", str), readdir()))

    # Identify tests with no matching file
    superfluous = filter(f -> f ∉ filebase, extbase)
    if length(superfluous) > 0
        println()
        @info "Potentially superfluous extension tests:"
        for f in superfluous
            println("    + $f.jl")
        end
        println()
    end

    # Identify files with no matching test
    notest = filter(f -> f ∉ extbase, filebase)
    if length(notest) > 0
        println()
        @info "Potentially missing extension tests:"
        for f in notest
            println("    - $f.jl")
        end
        println()
    end

    # Seed RNG to make tests reproducible
    Random.seed!(1234)

    @testset "EcoSISTEM.jl" begin
        @test isfile(EcoSISTEM.path("runtests.jl"))
        println()
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
        end

        println()
        @info "Running tests for extensions:"
        for t in extbase
            println("    = $t.jl")
        end
        println()

        @info "Running extension tests..."
        @testset for t in extbase
            fn = "ext_$t.jl"
            println("    * Testing $t.jl extension...")
            include(fn)
        end
    end

    # Identify files that are cross-validating results against other packages
    # test/pkg_Package.jl should validate results against the Package package

    pkgbase = map(file -> replace(file, r"pkg_(.*).jl$" => s"\1"),
                  filter(str -> occursin(r"^pkg_.*\.jl$", str),
                         readdir()))

    if length(pkgbase) > 0
        @info "Cross validation packages:"
        @testset begin
            for p in pkgbase
                println("    = $p")
            end
            println()

            @testset for p in pkgbase
                fn = "pkg_$p.jl"
                println("    * Validating $p.jl ...")
                include(fn)
            end
        end
    end

    # Identify files that are cross-validating results against other packages
    # test/pkg_Package.jl should validate results against the Package package

    @testset "Examples folder" begin
        println()
        @info "Running from examples folder ..."
        Pkg.activate(EcoSISTEM.path(dir = "examples"))
        Pkg.rm("EcoSISTEM")
        Pkg.develop(PackageSpec(path = EcoSISTEM.path(dir = "")))
        Pkg.instantiate()
        Pkg.update()
        example_testbase = map(file -> replace(file, r"test_(.*).jl" => s"\1"),
                               filter(str -> occursin(r"^test_.*\.jl$", str),
                                      readdir("../examples/")))
        for t in example_testbase
            fn = "../examples/test_$t.jl"
            println("    * Testing $t.jl ...")
            @test_nowarn include(fn)
        end
    end
end

if rsmd == "TRUE" || !haskey(ENV, "RUNNER_OS") # Crosswalk runner or local testing
    # Test RSMD crosswalk and other hygene issues

    # Identify files that are checking package hygene
    cleanbase = map(file -> replace(file, r"clean_(.*).jl$" => s"\1"),
                    filter(str -> occursin(r"^clean_.*\.jl$", str),
                           readdir()))

    if length(cleanbase) > 0
        @info "Crosswalk and clean testing:"
        @testset begin
            for c in cleanbase
                println("    = $c")
            end
            println()

            @testset for c in cleanbase
                fn = "clean_$c.jl"
                println("    * Verifying $c.jl ...")
                include(fn)
            end
        end
    end
end
