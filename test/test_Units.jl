module TestUnits

using EcoSISTEM, EcoSISTEM.Units
using Test

@testset "Checking Units" begin
    @test 12months == 12month == 1year
end

end
