# SPDX-License-Identifier: LGPL-3.0-or-later

module TestProvenance

using EcoSISTEM
using EcoSISTEM: InputRecord
using Dates: Dates
using ArchGDAL
using Test

# Write a zip holding the version file Natural Earth ships inside each of its own, through GDAL's
# zip filesystem.
function _versionedzip(path, entry, text)
    handle = ArchGDAL.GDAL.vsifopenl("/vsizip/" * path * "/" * entry, "wb")
    bytes = Vector{UInt8}(text)
    ArchGDAL.GDAL.vsifwritel(bytes, 1, length(bytes), handle)
    ArchGDAL.GDAL.vsifclosel(handle)
    return path
end

@testset "an input record holds where a published input came from, and nothing personal" begin
    r = InputRecord(role = :species, dataset = "GBIF occurrence download",
                    doi = "10.15468/dl.abc123", licence = "CC BY-NC 4.0",
                    path = "records.csv", fetched = "2026-09-13T07:21:07Z",
                    bytes = 12, sha256 = "00")
    @test r.role === :species && r.dataset == "GBIF occurrence download"
    @test r.fetched == Dates.DateTime(2026, 9, 13, 7, 21, 7)
    @test isnothing(r.url) && isnothing(r.request) && isempty(r.citation)
    @test string(r) ==
          "InputRecord(species, \"GBIF occurrence download\", \"records.csv\", doi = \"10.15468/dl.abc123\")"
    shown = sprint(show, MIME("text/plain"), r)
    @test occursin("InputRecord (species)", shown) && occursin("doi", shown) &&
          !occursin("url", shown)
    # A role outside the closed set, or an absolute path, is refused where it was written; the
    # set covers what a run is seeded from and restarted from as well as what it is built from.
    @test_throws ArgumentError InputRecord(role = :weather, dataset = "x")
    @test InputRecord(role = :abundance, dataset = "GBIF").role === :abundance
    @test InputRecord(role = :state, dataset = "burn-in.jld2").role === :state
    @test_throws ArgumentError InputRecord(role = :habitat, dataset = "x",
                                           path = joinpath(homedir(), "x.nc"))
    # A request body is kept as sent, keys as strings.
    q = InputRecord(role = :habitat, dataset = "ERA",
                    request = Dict(:variable => ["2m_temperature"]))
    @test q.request == Dict{String, Any}("variable" => ["2m_temperature"])
end

@testset "an input record carries what its provider asks, filled in when shown" begin
    E = EcoSISTEM
    r = InputRecord(role = :habitat, dataset = "ERA",
                    title = "ERA5, t2m (era5_t2m_1990s)",
                    primary = false, releasedate = "2019-04-18",
                    references = [(doi = "10.24381/cds.f17050d7",
                                   text = "C3S (2023), ERA5. (Accessed on DD-MMM-YYYY)")],
                    notices = [Dict("strength" => "required",
                                    "text" => "Contains Copernicus information YYYY.")],
                    requestpage = "https://example.org/requests/#era5",
                    region = "NaturalEarthSpec(\"Scotland\", level = \"ADMIN\")",
                    fetched = "2026-10-25T11:00:00Z")
    @test r.primary === false && r.releasedate == "2019-04-18"
    @test only(r.notices) ==
          (strength = :required, text = "Contains Copernicus information YYYY.")
    # The texts are kept as the provider wrote them, and filled from the fetch time when shown.
    @test only(r.references).text ==
          "C3S (2023), ERA5. (Accessed on DD-MMM-YYYY)"
    shown = sprint(show, MIME("text/plain"), r)
    @test occursin("cite    C3S (2023), ERA5. (Accessed on 25-Oct-2026)", shown)
    @test occursin("notice  (required) Contains Copernicus information 2026.",
                   shown)
    @test occursin("asks    https://example.org/requests/#era5", shown)
    @test occursin("Scotland", shown) && occursin("primary  false", shown)
    # With no fetch time the placeholders stay, so nothing is claimed that is not known.
    @test E._fillplaceholders("on DD-MMM-YYYY in YYYY", nothing) ==
          "on DD-MMM-YYYY in YYYY"
    @test E._fillplaceholders("on DD-MMM-YYYY in YYYY",
                              Dates.DateTime(2025, 7, 2)) ==
          "on 02-Jul-2025 in 2025"
    # A strength a provider cannot ask with is refused.
    @test_throws ArgumentError InputRecord(role = :habitat, dataset = "x",
                                           notices = [(strength = :urgent,
                                                       text = "x")])
    # The new fields survive a sidecar: written by the catalogue's record, read back the same.
    dir = mktempdir()
    path = joinpath(dir, "era5_t2m_1990s.nc")
    write(path, "x")
    record = E._baserecord(path, :habitat)
    merge!(record, E._cataloguerecord(E.ERA, E.layerinfo(E.ERA, "t2m"), path))
    record["primary"] = false
    record["region"] = r.region
    E._writesidecar(path, record)
    back = provenance(path)
    rec = E.datasetinfo(E.ERA)
    @test back.title == "$(rec.title), 2 metre temperature (era5_t2m_1990s)"
    @test back.primary === false && back.region == r.region
    @test back.references == rec.references && back.notices == rec.notices
    @test back.releasedate == rec.releasedate &&
          back.requestpage == E._requestpage(rec)
    # Two records of one file are equal, and hash alike, though their lists are separate objects.
    again = provenance(path)
    @test again == back && isequal(again, back) && hash(again) == hash(back)
    @test again != InputRecord(role = :habitat, dataset = "ERA")
end

@testset "a file with no record answers nothing; a sidecar reads back as a record" begin
    dir = mktempdir()
    path = joinpath(dir, "plain.nc")
    write(path, "x")
    @test isnothing(provenance(path))
    write(EcoSISTEM._sidecarpath(path),
          """
          writer = "EcoSISTEM 0.8.0"
          role = "habitat"
          file = "plain.nc"
          source = "TwentyCR"
          code = "air"
          url = "https://example.org/air.nc"
          final_url = "https://mirror.example.org/air.nc"
          fetched = "2026-09-13T08:13:59Z"
          bytes = 1
          sha256 = "ab"
          doi = "10.1002/qj.3598"
          licence = "free"
          version = "3"
          citation = "Slivinski et al. (2019)"
          """)
    r = provenance(path)
    @test r isa InputRecord && r.role === :habitat
    @test r.dataset == "TwentyCR" && r.code == "air" && r.path == "plain.nc"
    @test r.url == "https://mirror.example.org/air.nc"      # the URL answered from
    @test r.fetched == Dates.DateTime(2026, 9, 13, 8, 13, 59) && r.bytes == 1
    @test r.doi == "10.1002/qj.3598" && r.citation == "Slivinski et al. (2019)"
    # A CDS record: the request and the job come back too, the dataset being the CDS's name.
    cds = joinpath(dir, "era5_t2m_1990s.nc")
    write(cds, "x")
    write(EcoSISTEM._sidecarpath(cds),
          """
          writer = "EcoSISTEM 0.8.0"
          role = "habitat"
          file = "era5_t2m_1990s.nc"
          dataset = "reanalysis-era5-single-levels-monthly-means"
          source = "ERA"
          code = "t2m"
          job = "50fb750a-172d-4af8-ab53-9f7dff9d73d2"
          [request]
          variable = ["2m_temperature"]
          year = ["1990"]
          """)
    c = provenance(cds)
    @test c.dataset == "ERA" && c.job == "50fb750a-172d-4af8-ab53-9f7dff9d73d2"
    @test c.request["variable"] == ["2m_temperature"]
    # A sidecar another program wrote under the same name - prose in `role`, a sentence in
    # `source` - is neither read nor overwritten: it reports nothing, once with a warning.
    foreign = joinpath(dir, "air.2m.mon.mean.nc")
    write(foreign, "x")
    text = """
           source = "NOAA PSL, 20th Century Reanalysis V3 monthly means"
           role = "2 m air temperature (K)"
           writer = "data/src/twentycr_download.jl"
           """
    write(EcoSISTEM._sidecarpath(foreign), text)
    @test_logs (:warn, r"was not written by EcoSISTEM") match_mode=:any (@test isnothing(provenance(foreign)))
    @test read(EcoSISTEM._sidecarpath(foreign), String) == text
    # ...and a present file with no sidecar at all is recorded on first use, without a fetch time,
    # so a copy fetched by hand still ends up with a record.
    held = joinpath(dir, "held.bin")
    write(held, "held")
    asset = EcoSISTEM.CachedAsset(TwentyCR, "https://example.org/held.bin",
                                  path = held)
    @test EcoSISTEM.assetpath(asset) == held
    h = provenance(held)
    @test h.url == "https://example.org/held.bin" && isnothing(h.fetched)
    @test h.sha256 == EcoSISTEM._sha256(held) && h.bytes == 4
    # The foreign sidecar stays foreign through the same route.
    @test EcoSISTEM.assetpath(EcoSISTEM.CachedAsset(TwentyCR,
                                                    "https://example.org/air.nc",
                                                    path = foreign)) == foreign
    @test read(EcoSISTEM._sidecarpath(foreign), String) == text
end

@testset "a Natural Earth zip records its source's facts and the version it states" begin
    E = EcoSISTEM
    dir = mktempdir()
    zip = _versionedzip(joinpath(dir, "ne_10m_fixture.zip"),
                        "ne_10m_fixture.VERSION.txt", "5.1.1\r\n")
    @test E._zipversion(zip) == "5.1.1"
    # A zip without the entry, and no zip at all, state nothing.
    other = _versionedzip(joinpath(dir, "other.zip"), "readme.txt", "x")
    @test isnothing(E._zipversion(other))
    @test isnothing(E._zipversion(joinpath(dir, "absent.zip")))
    url = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_fixture.zip"
    rec = E._assetrecord(E.CachedAsset(E.NaturalEarthLevel, url, path = zip),
                         zip,
                         nothing)
    @test rec["role"] == "region" && rec["source"] == "NaturalEarth"
    @test rec["version"] == "5.1.1" && rec["licence"] == "CC-PDM-1.0"
    # The role follows what owns the download: a raster file fetched from a URL is a layer, a
    # vector file an outline.
    tif = joinpath(dir, "a.tif")
    write(tif, "x")
    @test E._assetrecord(E.CachedAsset(E.RasterSpec,
                                       "https://example.org/a.tif",
                                       path = tif), tif, nothing)["role"] ==
          "habitat"
    @test E._assetrecord(E.CachedAsset(E.ShapeSpec, "https://example.org/a.zip",
                                       path = tif), tif, nothing)["role"] ==
          "region"
    # A record holding only the file and its URL reads back with the row's facts and the zip's
    # own version filled in.
    write(E._sidecarpath(zip),
          """
          writer = "EcoSISTEM 0.8.0"
          role = "region"
          file = "ne_10m_fixture.zip"
          url = "$url"
          """)
    r = E._regionrecord(zip)
    @test r.dataset == "NaturalEarth" && r.version == "5.1.1" && r.url == url
    @test r.licence == E.datasetinfo(E.NaturalEarthLevel).licence
    @test isnothing(E._regionrecord(other))
end

@testset "a shape spec reports the records of the files it outlines" begin
    E = EcoSISTEM
    dir = mktempdir()
    a, b = joinpath(dir, "a.zip"), joinpath(dir, "b.gpkg")
    write(a, "a")
    write(b, "b")
    write(E._sidecarpath(a),
          """
          writer = "EcoSISTEM 0.8.0"
          role = "region"
          file = "a.zip"
          url = "https://example.org/a.zip"
          """)
    sa, sb = E.ShapeSpec(a), E.ShapeSpec(b)
    @test only(provenance(sa)).url == "https://example.org/a.zip"
    @test provenance(sb) == [nothing]
    # A combination lists its members' records in order, and reads nothing to do it.
    both = provenance(E.ConstructedShapeSpec(E.ShapeUnion(), sa, sb))
    @test length(both) == 2 && both[1].path == "a.zip" && isnothing(both[2])
    # A named region answers with its level's zip: its record where it has been fetched, `nothing`
    # where not, and never a download.
    named = provenance(E.NaturalEarthSpec("Scotland"))
    @test length(named) == 1 &&
          (isnothing(only(named)) || only(named).dataset == "NaturalEarth")
end

end
