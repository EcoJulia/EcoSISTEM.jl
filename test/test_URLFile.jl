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
    @test P.route(f) == PublishedByProvider() &&
          P.istrait(P.IsFetchable{URLFile})
    # Compared against `repr`, as `show` writes the two strings: a Windows path's backslashes are
    # escaped.
    @test string(f) ==
          "URLFile($(repr(url)), $(repr(f.path)), facts = " *
          "SourceFacts(\"My table\", doi 10.1234/abc))"
    path = P.fetchfile(f, role = :species)
    @test path == f.path && isfile(path) && read(path, String) == "a,b\n1,2\n"
    r = provenance(path)
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
    @test P.fetchfile(f, role = :species, verify = true) == path
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
    path = P.fetchfile(f, role = :region)
    r = provenance(path)
    @test !isnothing(r.fetched) && r.source == row.source &&
          isnothing(r.catalogue) && !P.iscatalogued(r)
    @test record_provenance(path, row, role = :region) == r
end

end
