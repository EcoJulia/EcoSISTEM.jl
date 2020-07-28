using Simulation
using Distributions
using Random
using Test

Random.seed!(0)
@testset "MedianGenerator" begin

@testset "median of Multinomial" begin
    @testset "issorted" begin
        for _ in 1:20
            @test_broken issorted(Multinomial(rand(1:5), rand(1:5)))
        end
    end
    @testset "n = 1" begin
        for i in 0:4
            @test median(Multinomial(i, 1)) == [i]
        end
    end

    @testset "n = 2" begin
        @test median(Multinomial(0, 2)) == [0, 0]
        @test_throws BoundsError median(Multinomial(1, 2)) == [1, 1]
        @test median(Multinomial(2, 2)) == [1, 1]
        @test median(Multinomial(3, 2)) == [1, 2]
        @test median(Multinomial(4, 2)) == [2, 2]
        @test_throws BoundsError median(Multinomial(5, 2))
    end

    @testset "n = 3" begin
        @test median(Multinomial(0, 3)) == [0, 0, 0]
        @test median(Multinomial(1, 3)) == [1, 1, 1]
        @test median(Multinomial(2, 3)) == [0, 1, 1]
        @test median(Multinomial(3, 3)) == [1, 1, 1]
        @test median(Multinomial(4, 3)) == [2, 2, 2]
        @test median(Multinomial(5, 3)) == [1, 2, 2]
        @test median(Multinomial(6, 3)) == [2, 2, 2]
        @test median(Multinomial(7, 3)) == [3, 3, 3]
        @test median(Multinomial(8, 3)) == [2, 3, 3]
        @test median(Multinomial(9, 3)) == [3, 3, 3]
        @test median(Multinomial(10, 3)) == [4, 4, 4]
        @test median(Multinomial(11, 3)) == [3, 4, 4]
        @test median(Multinomial(12, 3)) == [4, 4, 4]
        @test median(Multinomial(13, 3)) == [5, 5, 5]
        @test median(Multinomial(14, 3)) == [4, 5, 5]
        @test median(Multinomial(15, 3)) == [5, 5, 5]
    end
end

mg = MedianGenerator()

@testset "rand(::MedianGenerator, ::Binomial)" begin
    for _ in 1:20
        n = rand(1:5)
        b = Binomial(n, rand())
        @test rand(mg, b) == median(b)
    end
end

@testset "rand!(::MedianGenerator, ::Multivariate)" begin
    for _ in 1:20
        n = rand(1:5)
        p = rand(1:5)
        m = Multinomial(n, p)
        @test rand(mg, m) == median(m)
    end
end

@testset "rand!(::MedianGenerator, <:Distribution, output)" begin
    for _ in 1:20, d in (Binomial(rand(1:5), rand()), Multinomial(rand(1:10), rand(1:3)))
        n = rand(1:5)
        output = zeros(Int64, rand(1:10))
        rand!(mg, d, output)
        @test all(rand(mg, d) .== output)
    end
end

end

