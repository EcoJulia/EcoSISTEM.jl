# SPDX-License-Identifier: LGPL-3.0-or-later

module TestHabitats

using EcoSISTEM
using Unitful.DefaultSymbols
using Test
using EcoSISTEM.Units
using EcoSISTEM.ClimatePref
using AxisArrays
using RasterDataSources

@testset "Habitats" begin
    grid = (10, 10)
    area = 25.0km^2
    totalK = 10000.0kJ / km^2
    numNiches = 4
    active = fill(true, grid)

    # TEST simplehabitat
    fillval = 0.0
    habitat = simplehabitat(fillval, grid, totalK, area)
    @test EcoSISTEM.iscontinuous(habitat.regime) == true
    @test eltype(habitat.regime) == typeof(habitat.regime.matrix[1])
    @test size(habitat.regime, 1) == grid[1]
    @test EcoSISTEM.xmin(habitat.regime) == 0
    @test EcoSISTEM.ymin(habitat.regime) == 0
    @test EcoSISTEM.xcellsize(habitat.regime) == sqrt(area / prod(grid)) / km
    @test EcoSISTEM.ycellsize(habitat.regime) == sqrt(area / prod(grid)) / km
    @test EcoSISTEM.xcells(habitat.regime) == size(habitat.regime, 1)
    @test EcoSISTEM.ycells(habitat.regime) == size(habitat.regime, 2)
    @test EcoSISTEM.indices(habitat.regime) ==
          EcoSISTEM.coordinates(habitat.regime)
    @test EcoSISTEM.indices(habitat.regime, 1) ==
          repeat(collect(1:grid[1]), grid[2])

    # TEST tempgradhabitat
    habitat = tempgradhabitat(-10.0K, 10.0K, grid, totalK, area, 0.01K / month)
    @test EcoSISTEM.iscontinuous(habitat.regime) == true
    @test eltype(habitat.regime) == typeof(habitat.regime.matrix[1])
    @test size(habitat.regime, 1) == grid[1]

    # TEST simplenichehabitat
    habitat = simplenichehabitat(numNiches, grid, totalK, area)
    @test EcoSISTEM.iscontinuous(habitat.regime) == false
    @test eltype(habitat.regime) == typeof(habitat.regime.matrix[1])
    @test size(habitat.regime, 1) == grid[1]

    # TEST erahabitat
    temp = AxisArray(fill(1.0K, 10, 10, 3),
                     Axis{:latitude}((1:10) .* °),
                     Axis{:longitude}((1:10) .* °),
                     Axis{:time}(collect(1:3) .* s))
    eratemp = ERA(temp)
    active = fill(true, 10, 10)
    solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
    ea = erahabitat(eratemp, solar, active)
    @test EcoSISTEM.iscontinuous(ea.regime) == true
    @test eltype(ea.regime) == typeof(ea.regime.matrix[1])
    ea.regime.time = 2
    EcoSISTEM._resettime!(ea.regime)
    @test ea.regime.time == 1
    @test size(ea.regime, 1) == grid[1]

    # TEST worldclimhabitat
    temp = AxisArray(fill(1.0K, 10, 10, 12),
                     Axis{:latitude}(collect(1:10) .* m),
                     Axis{:longitude}(collect(1:10) .* m),
                     Axis{:time}(collect(1:12) .* month))
    worldclimtemp = ClimateRaster(WorldClim{Climate}, temp)
    active = fill(true, 10, 10)
    solar = SolarTimeSupply(fill(10.0kJ, 10, 10, 3), 1)
    worldclim = worldclimhabitat(worldclimtemp, solar, active)
    @test EcoSISTEM.iscontinuous(worldclim.regime) == true
    @test eltype(worldclim.regime) == typeof(worldclim.regime.matrix[1])
    worldclim.regime.time = 2
    EcoSISTEM._resettime!(worldclim.regime)
    @test worldclim.regime.time == 1
    @test size(worldclim.regime, 1) == grid[1]

    # Test multi regimes
    regime = RegimeCollection2(worldclim.regime, ea.regime)
    @test EcoSISTEM.iscontinuous(regime) == [true, true]
    @test eltype(regime) ==
          [typeof(regime.one.matrix[1]), typeof(regime.two.matrix[1])]
    @test size(regime, 1) == grid[1]
    regime.one.time = 2
    EcoSISTEM._resettime!(regime)
    @test regime.one.time == 1
    @test EcoSISTEM._getgridsize(regime) == regime.one.size
    @test isapprox(EcoSISTEM._getsize(regime), regime.one.size^2 * prod(grid),
                   rtol = 1e-5)
    @test EcoSISTEM._getdimension(regime) == grid
    @test EcoSISTEM._countsubcommunities(regime) == prod(grid)

    regime = RegimeCollection3(habitat.regime, worldclim.regime, ea.regime)
    @test EcoSISTEM.iscontinuous(regime) == [false, true, true]
    @test eltype(regime) == [
        typeof(regime.one.matrix[1]),
        typeof(regime.two.matrix[1]),
        typeof(regime.three.matrix[1])
    ]
    @test size(regime, 1) == grid[1]
    @test EcoSISTEM._getgridsize(regime) == regime.one.size
    @test isapprox(EcoSISTEM._getsize(regime), regime.one.size^2 * prod(grid),
                   rtol = 1e-5)
    @test EcoSISTEM._getdimension(regime) == grid
    @test EcoSISTEM._countsubcommunities(regime) == prod(grid)

    regime = RegimeCollection3(worldclim.regime, worldclim.regime, ea.regime)
    regime.two.time = 2
    EcoSISTEM._resettime!(regime)
    @test regime.two.time == 1
end

end
