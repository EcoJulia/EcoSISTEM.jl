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
using Unitful: ustrip, °

include("tomlstrings.jl")

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

@testset "a request is built from a catalogued code, and split by decade" begin
    r = CDSRequest(EcoSISTEM.ERA, "t2m", years = 1990:1999,
                   area = [40, -25, -36, 56], path = "era5_t2m_1990s.nc")
    @test r.dataset == EcoSISTEM._ERA5_MONTHLY && r.path == "era5_t2m_1990s.nc"
    @test r.request ==
          EcoSISTEM._era5request("2m_temperature", 1990:1999,
                                 area = [40, -25, -36, 56])
    # A code the table does not know is refused by the catalogue lookup itself.
    @test_throws ErrorException CDSRequest(EcoSISTEM.ERA, "t2m_no_such",
                                           years = 1990:1990, path = "x")
    # A whole decade is named as one; a partial block by its years, so that asking later for
    # one more year names a new file instead of finding the old one and stopping.
    rs = EcoSISTEM.era5requests("tp", 1985:2004, dir = "d")
    @test length(rs) == 3
    @test [basename(q.path) for q in rs] ==
          ["era5_tp_1985-1989.nc", "era5_tp_1990s.nc", "era5_tp_2000-2004.nc"]
    @test rs[1].request["year"] == string.(1985:1989)
    @test rs[3].request["year"] == string.(2000:2004)
    @test basename(only(EcoSISTEM.era5requests("tp", 2020:2025, dir = "d")).path) ==
          "era5_tp_2020-2025.nc"
    @test basename(only(EcoSISTEM.era5requests("tp", 2020:2026, dir = "d")).path) ==
          "era5_tp_2020-2026.nc"
    # The requests are a spec's files as they stand.
    spec = SourceSpec(EcoSISTEM.ERA, "tp", files = rs)
    @test length(spec.files) == 3
end

@testset "every field a request varies by is in its file's name" begin
    E = EcoSISTEM
    filenames(; kw...) = [basename(q.path)
                          for q in E.era5requests("t2m", 1985:1999; dir = "d",
                                                  kw...)]
    # The globe and the whole year keep the plain name, so files fetched that way are found.
    @test filenames() == ["era5_t2m_1985-1989.nc", "era5_t2m_1990s.nc"]
    # A box is its edges with their hemispheres, so no name holds a minus sign; months are their
    # runs; the box comes first.
    box = (40, -25, -36, 56)
    @test filenames(area = box) ==
          [
        "era5_t2m_1985-1989_40N25W36S56E.nc",
        "era5_t2m_1990s_40N25W36S56E.nc"
    ]
    @test filenames(months = 6:8) ==
          ["era5_t2m_1985-1989_m06-08.nc", "era5_t2m_1990s_m06-08.nc"]
    @test last(filenames(area = box, months = [12, 1, 2])) ==
          "era5_t2m_1990s_40N25W36S56E_m01-02.12.nc"
    @test last(filenames(months = [3])) == "era5_t2m_1990s_m03.nc"
    @test last(filenames(area = (0.5, 0, -10.25, 30))) ==
          "era5_t2m_1990s_0.5N0E10.25S30E.nc"
    # Two boxes, or two sets of months, never share a file.
    @test last(filenames(area = box)) !=
          last(filenames(area = (40, -25, -35, 56)))
    @test last(filenames(months = 6:8)) != last(filenames(months = 6:9))
    # Every name is a valid data product name, and portable.
    for name in vcat(filenames(area = box, months = [12, 1, 2]),
                     filenames(area = (0.5, 0, -10.25, 30)))
        @test occursin(r"^[A-Za-z0-9_.-]+$", name)
    end
end

@testset "an area is a box in any of its forms, or a region by name" begin
    E = EcoSISTEM
    want = [40.0, -25.0, -36.0, 56.0]
    request(area) = CDSRequest(E.ERA, "t2m", years = 1990:1990, area = area,
                               path = "x.nc").request["area"]
    @test request((40, -25, -36, 56)) == want
    @test request([40, -25, -36, 56]) == want
    # A named tuple is read by its names, whatever order they are written in.
    @test request((west = -25, south = -36, east = 56, north = 40)) == want
    @test request(E.Extents.Extent(X = (-25.0, 56.0), Y = (-36.0, 40.0))) ==
          want
    @test request(E.Extents.Extent(X = (-25.0°, 56.0°), Y = (-36.0°, 40.0°))) ==
          want
    @test !haskey(CDSRequest(E.ERA, "t2m", years = 1990:1990, path = "x.nc").request,
                  "area")
    # A region's box comes from the shipped table, with no download, rounded outwards to whole
    # degrees; a name meaning different ground at different levels needs its level.
    exact = E.boundingbox("Scotland")
    scotland = request("Scotland")
    @test all(isinteger, scotland)
    @test scotland[1] >= ustrip(exact.Y[2]) && scotland[3] <= ustrip(exact.Y[1])
    @test scotland[2] <= ustrip(exact.X[1]) && scotland[4] >= ustrip(exact.X[2])
    @test_throws "NaturalEarthSpec(\"Africa\", level = ...)" request("Africa")
    africa = request(NaturalEarthSpec("Africa", level = "CONTINENT"))
    @test africa ==
          request(E.boundingbox("Africa", level = "CONTINENT", round = 1°))
    # What is not a box is refused, naming the form.
    @test_throws "four numbers" request((40, -25, -36))
    @test_throws "north edge" request((-36, -25, 40, 56))
    @test_throws "names `north`" request((top = 40, left = -25, bottom = -36,
                                          right = 56))
    # Months are months of the year, and a request names them in order whatever order they
    # were asked in, so two requests for the same data compare equal.
    @test_throws "1 to 12" CDSRequest(E.ERA, "t2m", years = 1990:1990,
                                      months = 0:2,
                                      path = "x.nc")
    @test CDSRequest(E.ERA, "t2m", years = [1991, 1990], months = [8, 6, 7],
                     path = "x.nc").request ==
          CDSRequest(E.ERA, "t2m", years = 1990:1991, months = 6:8,
                     path = "x.nc").request
end

@testset "an ERA5 spec given years fetches into the asset cache" begin
    E = EcoSISTEM
    spec = SourceSpec(E.ERA, "t2m", years = 1990:2009, area = "Scotland",
                      months = 6:8)
    @test length(spec.files) == 2 && all(f -> f isa CDSRequest, spec.files)
    @test all(f -> dirname(f.path) == E.assetdir(owner = E.ERA), spec.files)
    @test [basename(f.path) for f in spec.files] ==
          [basename(q.path)
           for q in E.era5requests("t2m", 1990:2009, dir = "d",
                                   area = "Scotland", months = 6:8)]
    # Nothing is fetched at construction.
    @test !any(f -> isfile(f.path), spec.files)
    # `years` is required, goes with no location, asks for one layer, and only of a source the
    # Climate Data Store serves.
    @test_throws "give `years`" SourceSpec(E.ERA, "t2m")
    @test_throws "give `years` too" SourceSpec(E.ERA, "t2m", area = "Scotland")
    @test_throws "with none of" SourceSpec(E.ERA, "t2m", years = 1990:1999,
                                           file = "x.nc")
    @test_throws "name the layer's code" SourceSpec(E.ERA, years = 1990:1999)
    @test_throws "not fetched from it" SourceSpec(E.TwentyCR, "air",
                                                  years = 1990:1999)
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
        # The fetch wrote the file's provenance: the request as sent, the CDS job, the row's DOI.
        rec = provenance(path)
        @test rec.dataset == "ERA" && rec.code == "t2m"
        @test rec.request["variable"] == ["2m_temperature"] &&
              !isnothing(rec.job)
        @test rec.doi == "10.24381/cds.f17050d7" && !isnothing(rec.sha256)
        # The CDS's record of converting the file, which differs between fetches of the same data,
        # is replaced in the file and kept in the record.
        EcoSISTEM.NCDatasets.Dataset(path) do ds
            return @test ds.attrib["history"] ==
                         rpad(first(EcoSISTEM._HISTORYNOTE,
                                    length(rec.history)),
                              length(rec.history))
        end
        @test occursin("GRIB to CDM", rec.history) && rec.primary === false
        @test !any(s -> occursin(homedir(), s),
                   tomlstrings(EcoSISTEM._sidecarpath(path)))
        # The file reads through the same row as any other ERA5 file - no extension on its name.
        cr = read(SourceSpec(EcoSISTEM.ERA, "t2m", file = path))
        @test size(cr.array, 3) == 1
        @test 250 < ustrip(cr.array[1, 1, 1]) < 300
    else
        @info "Climate Data Store fetch not exercised: no ~/.cdsapirc, on a CI runner, or ECOSISTEM_CDS_FETCH=false"
    end
end

end
