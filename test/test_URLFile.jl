# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/Provenancer/URLFile.jl` and the fetch verb in `src/Provenancer/assets.jl`: a file
# fetched from a URL once, recorded with what is known of its source, and left alone after. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_URLFile.jl"])'

module TestURLFile

using EcoSISTEM
using EcoSISTEM.Provenancer
using Test

const P = EcoSISTEM.Provenancer

@testset "a URL file is fetched once and recorded with its source's facts" begin
    dir = mktempdir()
    served = joinpath(dir, "served.csv")
    write(served, "a,b\n1,2\n")
    url = "file://" * served
    facts = SourceFacts(source = "My table", doi = "10.1234/abc",
                        title = "A table",
                        licence = "CC-BY-4.0", identifies = :data,
                        releasedate = "2024-05")
    f = URLFile(url, joinpath(dir, "cache", "table.csv"), facts = facts)
    @test P.route(f) == VersionOfRecord() &&
          P.istrait(P.IsFetchable{URLFile})
    # Compared against `repr`, as `show` writes the two strings: a Windows path's backslashes are
    # escaped.
    @test string(f) ==
          "URLFile($(repr(url)), $(repr(f.path)), facts = " *
          "SourceFacts(\"My table\", doi 10.1234/abc))"
    # The record is returned, saying where the file and its record are.
    fetched = P.fetchfile(f, role = :species)
    @test fetched isa ProvenanceRecord && fetched.location == f.path &&
          fetched.sidecar == P.sidecarpath(f.path) &&
          fetched.route == VersionOfRecord()
    path = fetched.location
    @test isfile(path) && read(path, String) == "a,b\n1,2\n"
    r = provenance(path)
    @test r == fetched
    @test r.role === :species && r.source == "My table" &&
          r.doi == "10.1234/abc" &&
          r.url == url && r.primary === true && r.licence == "CC-BY-4.0" &&
          r.releasedate == "2024-05" && r.title == "A table" &&
          r.path == "table.csv"
    @test !isnothing(r.fetched) && r.bytes == 8 &&
          r.sha256 == P._filesha256(path)
    # A second call fetches nothing and leaves the record alone.
    sidecar = P.sidecarpath(path)
    mtimes = (stat(path).mtime, stat(sidecar).mtime)
    @test P.fetchfile(f, role = :species, verify = true) == fetched
    @test (stat(path).mtime, stat(sidecar).mtime) == mtimes
    # A file present with no record is recorded without a fetch time, the facts still written.
    rm(sidecar)
    P.fetchfile(f, role = :species)
    r2 = provenance(path)
    @test isnothing(r2.fetched) && r2.source == "My table" &&
          r2.sha256 == r.sha256
    # Without facts, the download alone is recorded.
    g = URLFile(url, joinpath(dir, "cache", "plain.csv"))
    P.fetchfile(g, role = :habitat)
    r3 = provenance(g.path)
    @test r3.source == "file" && isnothing(r3.doi) && isnothing(r3.primary) &&
          r3.url == url &&
          r3.role === :habitat
    @test string(g) == "URLFile($(repr(url)), $(repr(g.path)))"
    # A file altered since it was fetched is refused when verified.
    write(path, "a,b\n1,3\n")
    @test_throws "does not match the checksum" P.fetchfile(f, role = :species,
                                                           verify = true)
end

@testset "a URL file given a catalogued source's facts is not the catalogue's" begin
    # The catalogue's name in a record comes from the catalogue's own fetch, so a download of your
    # own recorded with a row's facts is a file of yours, which `record_provenance` returns.
    dir = mktempdir()
    served = joinpath(dir, "served.csv")
    write(served, "a,b\n1,2\n")
    row = SourceFacts(EcoSISTEM.datasetinfo(EcoSISTEM.NaturalEarthLevel))
    f = URLFile("file://" * served, joinpath(dir, "cache", "ne.csv"),
                facts = row)
    r = P.fetchfile(f, role = :region)
    path = r.location
    @test !isnothing(r.fetched) && r.source == row.source &&
          isnothing(r.catalogue) && !P.iscatalogued(r)
    @test record_provenance(path, row, role = :region) == r
end

@testset "a URL file's record is what a spec of it takes, and its track a run's" begin
    dir = mktempdir()
    served = joinpath(dir, "served")
    mkpath(served)
    write(joinpath(served, "region.geojson"), "{}")
    write(joinpath(served, "layer.tif"), "x")
    facts = SourceFacts(source = "My region", doi = "10.1234/abc",
                        identifies = :data)
    shape = URLFile("file://" * joinpath(served, "region.geojson"),
                    joinpath(dir, "cache", "region.geojson"), facts = facts)
    layer = URLFile("file://" * joinpath(served, "layer.tif"),
                    joinpath(dir, "cache", "layer.tif"), facts = facts)
    # A spec of a download takes its record in place of a path, as a spec of a recorded file does.
    r = P.fetchfile(shape, role = :region)
    spec = ShapeSpec(r)
    @test spec.path === r && EcoSISTEM._specpaths(spec) == [shape.path]
    raster = RasterFileSpec(P.fetchfile(layer, role = :habitat),
                            axis = EcoSISTEM.NicheAxis)
    @test EcoSISTEM._specpaths(raster) == [layer.path]
    # A track with no handle is the untracked call; any other handle is no pipeline run's.
    @test P.fetchfile(shape, role = :region,
                      track = Track(nothing, "tables")) == r
    @test_throws "FAIR data pipeline run" P.fetchfile(shape, role = :region,
                                                      track = Track(Ref(0),
                                                                    "tables"))
    # A download has nothing to do untracked, since it is fetched on first use, there or not.
    absent = URLFile(shape.url, joinpath(dir, "cache", "absent.geojson"),
                     facts = facts)
    @test isempty(EcoSISTEM.missingfiles(absent)) &&
          isempty(EcoSISTEM.missingfiles(absent,
                                         track = Track(nothing, "tables")))
    @test_throws "FAIR data pipeline run" EcoSISTEM.missingfiles(absent,
                                                                 track = Track(Ref(0),
                                                                               "tables"))
    # A file beside a same-named sidecar that is not a record of ours has no record to return.
    foreign = joinpath(dir, "cache", "foreign.geojson")
    write(foreign, "{}")
    write(P.sidecarpath(foreign), "writer = \"another program\"\n")
    @test_throws "has no provenance record of ours" P.fetchfile(URLFile(shape.url,
                                                                        foreign),
                                                                role = :region)
    # Nor is such a sidecar overwritten by a download: the file is not fetched beside it.
    unfetched = joinpath(dir, "cache", "unfetched.geojson")
    write(P.sidecarpath(unfetched), "writer = \"another program\"\n")
    @test_throws "neither read nor overwritten" P.fetchfile(URLFile(shape.url,
                                                                    unfetched),
                                                            role = :region)
    @test !isfile(unfetched) &&
          read(P.sidecarpath(unfetched), String) ==
          "writer = \"another program\"\n"
    # A stale record of ours beside a file not there is replaced by the download's.
    stale = joinpath(dir, "cache", "stale.geojson")
    P.writesidecar(stale,
                   Dict{String, Any}("writer" => P._writer(),
                                     "role" => "region",
                                     "source" => "stale"))
    @test P.fetchfile(URLFile(shape.url, stale, facts = facts),
                      role = :region).source ==
          "My region"
end

end
