# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests the recording of a present file in `src/Provenancer/assets.jl`: `record_input` writing a
# file's record from its source's facts or from its DOI, once, beside the file or where `sidecar`
# says, and the sidecar functions reading from there. The two records AfricaPlants writes by hand
# - its climate envelopes table and its GBIF occurrence download - are the cases, their facts from
# the saved replies in `test/doi/`. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_assets.jl"])'

module TestAssets

using EcoSISTEM
using EcoSISTEM.Provenancer
using TOML: TOML
using Test

const P = EcoSISTEM.Provenancer
include("doi/replies.jl")

const ENVELOPES_DOI = "10.5285/ca339c86-3674-4030-b891-35326e71141e"
const PAPER_DOI = "10.1111/geb.13564"
const GBIF_DOI = "10.15468/dl.n4a2sa"

# A small file standing for a table the user downloaded.
function tablefile(dir, name)
    path = joinpath(dir, name)
    write(path, "species,tmin,tmax\nAcacia tortilis,10.0,35.0\n")
    return path
end

@testset "a published table is recorded from its DOI's facts, once" begin
    dir = mktempdir()
    path = tablefile(dir, "envelopes.csv")
    facts = lookup_doi(ENVELOPES_DOI, fetch = saved)
    r = record_input(path, facts, role = :species)
    @test r isa InputRecord && r.role === :species &&
          r.path == "envelopes.csv"
    @test r.source == "Bioclimatic envelopes for global plant species" &&
          r.doi == ENVELOPES_DOI
    @test r.primary === true && r.licence == "ogl-uk-3.0" &&
          r.releasedate == "2022-06-22"
    @test isnothing(r.fetched) && r.bytes == filesize(path) &&
          r.sha256 == P._filesha256(path) && isnothing(r.url)
    @test startswith(r.citation, "C. Harris, ") &&
          occursin("(2022). Bioclimatic envelopes for global plant species. NERC",
                   r.citation) &&
          endswith(r.citation, "https://doi.org/" * ENVELOPES_DOI)
    # The sidecar is beside the file, ours, and reads back as the same record.
    sidecar = P.sidecarpath(path)
    @test isfile(sidecar) && provenance(path) == r
    t = TOML.parsefile(sidecar)
    @test startswith(t["writer"], "EcoSISTEM.Provenancer") &&
          !haskey(t, "fetched") && t["primary"] === true &&
          t["role"] == "species"
    # Recorded once: a second call, by a DOI that would not resolve, looks nothing up and returns
    # the record as it stands.
    mtime = stat(sidecar).mtime
    @test record_input(path, "10.9999/never-looked-up", role = :habitat) == r
    @test stat(sidecar).mtime == mtime
    # The facts as a user amends them: the paper describing the data, as a work to cite.
    amended = SourceFacts(facts,
                          references = [(doi = PAPER_DOI,
                                         text = "Harris et al. (2022), Global Ecology and " *
                                                "Biogeography 31, 2191-2203")])
    @test amended.doi == facts.doi && amended.title == facts.title &&
          amended.authors == facts.authors &&
          only(amended.references).doi == PAPER_DOI && facts.references == []
    @test_throws ArgumentError SourceFacts(facts, identifies = :dataset)
    rm(sidecar)
    r2 = record_input(path, amended, role = :species)
    @test only(r2.references).doi == PAPER_DOI && r2.citation == r.citation &&
          r2.primary === true
    shown = sprint(show, MIME("text/plain"), r2)
    @test occursin("cite    Harris et al. (2022)", shown)
end

@testset "an occurrence download is an extract, primary by its own DOI" begin
    dir = mktempdir()
    path = tablefile(dir, "records_africa.csv")
    facts = lookup_doi(GBIF_DOI, fetch = saved)
    # The citation GBIF asks for, in place of the one composed from the record's creator.
    cited = SourceFacts(facts,
                        citation = "GBIF.org (12 September 2026) GBIF Occurrence Download " *
                                   "https://doi.org/" * GBIF_DOI)
    url = "https://api.gbif.org/v1/occurrence/download/request/0012345-250000000000000.zip"
    r = record_input(path, cited, role = :abundance,
                     route = ExtractedByService(GBIF_DOI), url = url)
    @test r.role === :abundance && r.primary === true && r.doi == GBIF_DOI &&
          r.url == url
    @test r.licence == "cc-by-nc-4.0" && r.releasedate == "2026-09-12" &&
          r.source == "Occurrence Download"
    @test startswith(r.citation, "GBIF.org (12 September 2026)")
    # Without the service's identifier on the route the extract is not primary, nor is a file the
    # user altered; the citation composed from the looked-up facts names the record's creator.
    rm(P.sidecarpath(path))
    @test record_input(path, facts, role = :abundance,
                       route = ExtractedByService()).primary === false
    rm(P.sidecarpath(path))
    r3 = record_input(path, facts, role = :abundance, route = SuppliedByUser())
    @test r3.primary === false && isnothing(r3.url)
    @test r3.citation ==
          "GBIF.org User (2026). Occurrence Download. The Global Biodiversity " *
          "Information Facility. https://doi.org/" * GBIF_DOI
end

@testset "the sidecar goes where `sidecar` says, and is read from there" begin
    dir = mktempdir()
    path = tablefile(dir, "readonly.csv")
    elsewhere = joinpath(mktempdir(), "records", "readonly.csv.provenance.toml")
    facts = SourceFacts(source = "A table", title = "A table",
                        licence = "CC-BY-4.0")
    r = record_input(path, facts, role = :species, sidecar = elsewhere,
                     title = "Chosen rows", code = "rows")
    @test isfile(elsewhere) && !isfile(P.sidecarpath(path))
    @test r.path == "readonly.csv" && r.title == "Chosen rows" &&
          r.code == "rows" &&
          r.primary === false && r.doi == "" && r.licence == "CC-BY-4.0"
    @test isnothing(provenance(path)) &&
          provenance(path, sidecar = elsewhere) == r
    @test P.recordat(path, sidecar = elsewhere) == r &&
          isnothing(P.recordat(path))
    @test P.verifyfile(path, sidecar = elsewhere) && !P.verifyfile(path)
    # A second call finds the record where it was put.
    @test record_input(path, facts, role = :habitat, sidecar = elsewhere) == r
    write(path, "species,tmin,tmax\nAcacia tortilis,11.0,35.0\n")
    @test_throws "does not match the checksum" P.verifyfile(path,
                                                            sidecar = elsewhere)
    # The writer and the present-file recorder take the same keyword, making the directory.
    other = tablefile(dir, "other.csv")
    side = joinpath(dir, "side", "other.toml")
    @test P.writesidecar(other, P.baserecord(other, :habitat),
                         sidecar = side) ==
          side
    @test P.recordat(other, sidecar = side).role === :habitat
    rm(side)
    P.recordpresent(other, sidecar = side) do
        return P.baserecord(other, :region)
    end
    @test P.recordat(other, sidecar = side).role === :region &&
          isnothing(P.recordat(other, sidecar = side).fetched)
end

@testset "an absent file, and another program's sidecar, are refused" begin
    dir = mktempdir()
    facts = SourceFacts(source = "x")
    @test_throws ArgumentError record_input(joinpath(dir, "missing.csv"), facts,
                                            role = :species)
    path = tablefile(dir, "theirs.csv")
    theirs = "writer = \"SomeoneElse 1.0\"\nfile = \"theirs.csv\"\n"
    write(P.sidecarpath(path), theirs)
    @test_throws "neither read from it nor recorded over it" record_input(path,
                                                                          facts,
                                                                          role = :species)
    @test read(P.sidecarpath(path), String) == theirs
    # Named another place, the record is written there.
    r = record_input(path, facts, role = :species,
                     sidecar = joinpath(dir, "ours.toml"))
    @test r.source == "x" && r.primary === false
end

# The live lookup runs off a CI runner, and not with ECOSISTEM_DOI_LOOKUP=false.
function livelookup()
    haskey(ENV, "RUNNER_OS") && return false
    return get(ENV, "ECOSISTEM_DOI_LOOKUP", "true") != "false"
end

@testset "a file is recorded from its DOI, looked up live" begin
    if livelookup()
        dir = mktempdir()
        path = tablefile(dir, "envelopes.csv")
        r = record_input(path, ENVELOPES_DOI, role = :species)
        @test r.doi == ENVELOPES_DOI && r.primary === true &&
              r.licence == "ogl-uk-3.0" &&
              r.source == "Bioclimatic envelopes for global plant species"
    else
        @info "DOI lookup not exercised live: on a CI runner, or ECOSISTEM_DOI_LOOKUP=false"
        @test true
    end
end

end
