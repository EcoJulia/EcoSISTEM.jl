# SPDX-License-Identifier: LGPL-3.0-or-later

module TestProvenance

using EcoSISTEM
using EcoSISTEM: InputRecord
using Dates: Dates
using Test

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
    # A role outside the closed set, or an absolute path, is refused where it was written.
    @test_throws ArgumentError InputRecord(role = :weather, dataset = "x")
    @test_throws ArgumentError InputRecord(role = :habitat, dataset = "x",
                                           path = joinpath(homedir(), "x.nc"))
    # A request body is kept as sent, keys as strings.
    q = InputRecord(role = :habitat, dataset = "ERA",
                    request = Dict(:variable => ["2m_temperature"]))
    @test q.request == Dict{String, Any}("variable" => ["2m_temperature"])
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
end

end
