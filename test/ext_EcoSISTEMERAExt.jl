# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `ext/EcoSISTEMERAExt/`, the Climate Data Store fetch. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["ext_EcoSISTEMERAExt.jl"])'
# or as part of `core_ext.jl`.

module TestEcoSISTEMERAExt

using Test
using EcoSISTEM
using EcoSISTEM: CDSRequest, assetpath
using CDSAPI
using Unitful: ustrip

# Whether the one test that talks to the Climate Data Store runs: never on a CI runner, only
# where a CDS key is present, and not when switched off. The request it sends is one month of one
# variable over a two-degree box, so the queue rather than the transfer is the cost.
function cdsfetch()
    haskey(ENV, "RUNNER_OS") && return false
    get(ENV, "ECOSISTEM_CDS_FETCH", "true") == "false" && return false
    return isfile(expanduser("~/.cdsapirc"))
end

@testset "the extension loads on CDSAPI and supplies the fetch hook" begin
    @test !isnothing(Base.get_extension(EcoSISTEM, :EcoSISTEMERAExt))
    @test hasmethod(EcoSISTEM._fetchcds, Tuple{CDSRequest})
end

@testset "a request whose file exists is never sent" begin
    dir = mktempdir()
    path = joinpath(dir, "have.nc")
    write(path, "not a netCDF file, and never opened")
    r = CDSRequest(EcoSISTEM._ERA5_MONTHLY,
                   Dict("variable" => ["2m_temperature"]), path)
    # A fetch would need a key and the network; presence is the whole check.
    @test assetpath(r) == path
    @test r.request == Dict("variable" => ["2m_temperature"])
end

@testset "an ERA5 monthly-means request is built as the CDS takes it" begin
    # The request the Africa_plants downloads were made with, as their provenance records it.
    r = EcoSISTEM._era5request("2m_temperature", 1990:1999,
                               area = [40, -25, -36, 56])
    @test r == Dict("product_type" => ["monthly_averaged_reanalysis"],
               "variable" => ["2m_temperature"],
               "year" => string.(1990:1999),
               "month" => ["01", "02", "03", "04", "05", "06", "07", "08",
                   "09", "10", "11", "12"],
               "time" => ["00:00"],
               "area" => [40, -25, -36, 56],
               "data_format" => "netcdf",
               "download_format" => "unarchived")
    @test !haskey(EcoSISTEM._era5request("total_precipitation", 2000:2000),
                  "area")
    @test EcoSISTEM._era5request("total_precipitation", 2000:2000,
                                 months = 6:7)["month"] == ["06", "07"]
end

@testset "a live fetch, where a CDS key is present" begin
    if cdsfetch()
        dir = mktempdir()
        r = CDSRequest(EcoSISTEM._ERA5_MONTHLY,
                       EcoSISTEM._era5request("2m_temperature", 2020:2020,
                                              months = 1:1,
                                              area = [56, -4, 55, -3]),
                       joinpath(dir, "t2m_2020"))
        path = assetpath(r)
        @test isfile(path)
        @test !isfile(path * ".part")
        # The file reads through the same row as any other ERA5 file - no extension on its name.
        cr = read(SourceSpec(EcoSISTEM.ERA, "t2m", file = path))
        @test size(cr.array, 3) == 1
        @test 250 < ustrip(cr.array[1, 1, 1]) < 300
    else
        @info "Climate Data Store fetch not exercised: no ~/.cdsapirc, on a CI runner, or ECOSISTEM_CDS_FETCH=false"
    end
end

end
