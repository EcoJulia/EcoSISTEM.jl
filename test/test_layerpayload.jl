# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `src/layerpayload.jl` - a built layer taken apart into plain data and rebuilt.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_layerpayload.jl"])'

module TestLayerPayload

using Test
using EcoSISTEM
using EcoSISTEM: materialise, in_memory_raster
using EcoSISTEM: _layerpayload, _rebuildlayer, _allocatepayload,
                 _checkdescriptor, _reportpayload, _rebuildreport
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using Dates: Date
using DimensionalData: DimensionalData, DimArray, Y, dims, refdims
using Rasters
using ArchGDAL
include("rasterfixtures.jl")
include("layercompare.jl")

const SD = EcoSISTEM.SyntheticData
const COND = EcoSISTEM.Condition
const RES = EcoSISTEM.Resource

# What a receiver does: allocate from the descriptors, fill, rebuild.
function roundtrip(layer)
    payload = _layerpayload(layer)
    received = map(_allocatepayload, payload.descriptors)
    foreach(copyto!, received, payload.arrays)
    return _rebuildlayer(payload.skeleton, received)
end

const EAST = (245000.0:2500.0:262500.0) .* m
const NORTH = (640000.0:2500.0:655000.0) .* m

function bng(values; kw...)
    return _testraster(SD, values; lat = NORTH, long = EAST,
                       crs = Rasters.EPSG(27700), kw...)
end

flat() = bng([291.0K + (i + 2j) * 0.1K for i in 1:7, j in 1:8])

function cube(time)
    values = [4.0mm / day * (i + 2j + 3k) for i in 1:7, j in 1:8, k in 1:3]
    values[2, 3, 2] = NaN * mm / day
    return bng(values, time = time)
end

function projected()
    return StudyArea(regime = _reg(flat(), axis = Temperature),
                     verbosity = :silent)
end

@testset "every layer a spec builds survives the round trip" begin
    area = projected()
    synthetic = StudyArea(extent = (70.0km, 80.0km), cellsize = 10.0km,
                          verbosity = :silent)
    regime = GradientSpec(280.0K, 300.0K, axis = Temperature)
    supply = UniformSpec(10.0kJ / m^2 / day, axis = SolarRadiation)
    niche = NicheSpec(4, axis = EcoSISTEM.NicheAxis, seed = 1)
    data = in_memory_raster(flat(), axis = Temperature)
    celsius = in_memory_raster(bng(fill(18.0u"°C", 7, 8)), axis = Temperature)
    undated = in_memory_raster(cube((1:3) .* month_mean_duration),
                               axis = Precipitation)
    dates = [Date(2000, 1, 15), Date(2000, 2, 15), Date(2000, 3, 15)]
    dated = in_memory_raster(cube(dates), axis = Precipitation,
                             atend = HoldAtEnd())
    named = (summer = UniformSpec(290.0K, axis = Temperature), winter = data)
    layers = [
        "a synthetic regime" => materialise(regime, synthetic, role = COND),
        "a synthetic supply" => materialise(supply, synthetic, role = RES),
        "a niche layout" => materialise(niche, synthetic),
        "a 2-D data layer" => materialise(data, area, role = COND),
        "a layer in °C" => materialise(celsius, area),
        "an undated series" => materialise(undated, area, role = RES),
        "a dated series" => materialise(dated, area, role = COND),
        "a named collection" => materialise(named, area, role = COND)]
    for (label, layer) in layers
        @testset "$label" begin
            @test samelayer(layer, roundtrip(layer))
        end
    end
    # The cases the fixtures exist for, so a change in what `materialise` builds cannot quietly
    # leave them untested.
    kinds = Dict(layers)
    @test kinds["a niche layout"] isa EcoSISTEM.CategoricalLayer
    @test eltype(kinds["a niche layout"].matrix) === Int64
    @test unit(eltype(kinds["a layer in °C"].matrix)) == u"°C"
    series = kinds["an undated series"]
    @test series.change isa EcoSISTEM.SeriesLayerChange
    @test !isempty(refdims(series.matrix))
    @test any(isnan, series.change.slices)
    @test kinds["a dated series"].change.calendar isa EcoSISTEM.DatedSeries
    @test Rasters.crs(dims(series.matrix, Y)) == Rasters.EPSG(27700)
    @test keys(kinds["a named collection"]) == (:summer, :winter)
    # One descriptor per array, describing it exactly: the matrix, the slices and the baseline.
    payload = _layerpayload(series)
    @test length(payload.arrays) == length(payload.descriptors) == 3
    @test all(d.eltype === eltype(a) && d.size == size(a)
              for (d, a) in zip(payload.descriptors, payload.arrays))
end

@testset "a mask's bits travel as Bools" begin
    active = projected().report.active
    @test active isa Raster
    @test parent(active) isa BitMatrix
    arrays, descriptors = Array[], NamedTuple[]
    skeleton = EcoSISTEM._packdimarray!(arrays, descriptors, active, "`active`")
    @test only(arrays) isa Matrix{Bool}
    @test only(descriptors).bits
    rebuilt = EcoSISTEM._unpackdimarray(skeleton, Array[copy(only(arrays))])
    @test samedimarray(active, rebuilt)
end

@testset "a study area's report survives the round trip" begin
    reports = ["a data area" => projected().report,
        "a synthetic area" => StudyArea(extent = (70.0km, 80.0km),
                  cellsize = 10.0km,
                  verbosity = :silent).report,
        "an investigation" => investigate_study_area(regime = _reg(flat(),
                                             axis = Temperature))]
    for (label, report) in reports
        @testset "$label" begin
            payload = _reportpayload(report)
            received = map(_allocatepayload, payload.descriptors)
            foreach(copyto!, received, payload.arrays)
            rebuilt = _rebuildreport(payload.skeleton, received, report.specs,
                                     report.constraints, report.cache)
            @test samereport(report, rebuilt)
            # The three a receiver supplies are its own, never copies.
            @test rebuilt.specs === report.specs
            @test rebuilt.constraints === report.constraints
            @test rebuilt.cache === report.cache
            @test length(payload.arrays) == 1
        end
    end
    # A mismatched array is refused, so a report cannot be rebuilt from a layer's arrays.
    report = projected().report
    payload = _reportpayload(report)
    @test_throws "where its skeleton describes" _rebuildreport(payload.skeleton,
                                                               [zeros(2, 2)],
                                                               report.specs,
                                                               report.constraints,
                                                               report.cache)
end

@testset "what cannot be sent is refused by name" begin
    layer = materialise(in_memory_raster(flat(), axis = Temperature),
                        projected(), role = COND)
    steady = typeof(layer)(layer.matrix, layer.size,
                           EcoSISTEM.SteadyLayerChange(1.0K / year))
    @test_throws "SteadyLayerChange, which cannot yet be sent" _layerpayload(steady)

    # The element limit is checked on the descriptor, so nothing is allocated.
    big = (eltype = Float64, size = (2^16, 2^15), bits = false)
    @test_throws "has 2147483648 elements" _checkdescriptor(big,
                                                            "the layer's `matrix`")
    @test _checkdescriptor((eltype = Float64, size = (2^16, 2^15 - 1),
                            bits = false), "fits") isa NamedTuple

    payload = _layerpayload(layer)
    @test_throws "ran out of arrays" _rebuildlayer(payload.skeleton, Array[])
    @test_throws "1 more array(s)" _rebuildlayer(payload.skeleton,
                                                 [payload.arrays;
                                                  payload.arrays])
    @test_throws "where its skeleton describes" _rebuildlayer(payload.skeleton,
                                                              [Float64.(ustrip.(only(payload.arrays)))])
end

end
