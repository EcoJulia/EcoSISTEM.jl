# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests the recording of a present file in `src/Provenancer/assets.jl`: `record_provenance` writing
# a file's record from its source's facts or from its DOI, once, beside the file or where `sidecar`
# says, refusing a record that says otherwise, and the sidecar functions reading from there. The
# two records AfricaPlants writes by hand
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
    r = record_provenance(path, facts, role = :species)
    @test r isa ProvenanceRecord && r.role === :species &&
          r.path == "envelopes.csv"
    # The record says where the file and its sidecar are, how the file came to be and what its
    # identifier names, none of which is written but the last two.
    @test r.location == path && r.sidecar == P.sidecarpath(path) &&
          r.route == PublishedByProvider() && r.identifies === :data &&
          isnothing(r.namespace)
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
          t["role"] == "species" && t["route"] == "PublishedByProvider" &&
          !haskey(t, "route_identifier") && !haskey(t, "location") &&
          !haskey(t, "sidecar")
    # Recorded once: a second call from the same source under the same role, by a DOI that would
    # not resolve, looks nothing up and returns the record as it stands; one naming another
    # source, or another role, is refused, and the record is left as it was.
    mtime = stat(sidecar).mtime
    @test record_provenance(path, uppercase(ENVELOPES_DOI), role = :species) ==
          r
    @test record_provenance(path, facts, role = :species) == r
    @test_throws "rather than the DOI `10.9999/never-looked-up`" record_provenance(path,
                                                                                   "10.9999/never-looked-up",
                                                                                   role = :species)
    @test_throws "under the role `species` rather than `habitat`" record_provenance(path,
                                                                                    facts,
                                                                                    role = :habitat)
    @test stat(sidecar).mtime == mtime && provenance(path) == r
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
    r2 = record_provenance(path, amended, role = :species)
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
    r = record_provenance(path, cited, role = :abundance,
                          route = ExtractedByService(GBIF_DOI), url = url)
    @test r.role === :abundance && r.primary === true && r.doi == GBIF_DOI &&
          r.url == url
    # The route and the service's identifier are kept, and read back.
    @test r.route == ExtractedByService(GBIF_DOI)
    t = TOML.parsefile(P.sidecarpath(path))
    @test t["route"] == "ExtractedByService" &&
          t["route_identifier"] == GBIF_DOI
    @test r.licence == "cc-by-nc-4.0" && r.releasedate == "2026-09-12" &&
          r.source == "Occurrence Download"
    @test startswith(r.citation, "GBIF.org (12 September 2026)")
    # Without the service's identifier on the route the extract is not primary, nor is a file the
    # user altered; the citation composed from the looked-up facts names the record's creator.
    rm(P.sidecarpath(path))
    @test record_provenance(path, facts, role = :abundance,
                            route = ExtractedByService()).primary === false
    rm(P.sidecarpath(path))
    r3 = record_provenance(path, facts, role = :abundance,
                           route = SuppliedByUser())
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
    r = record_provenance(path, facts, role = :species, sidecar = elsewhere,
                          title = "Chosen rows", code = "rows")
    @test isfile(elsewhere) && !isfile(P.sidecarpath(path))
    @test r.path == "readonly.csv" && r.title == "Chosen rows" &&
          r.code == "rows" &&
          r.primary === false && isnothing(r.doi) && r.licence == "CC-BY-4.0"
    @test isnothing(provenance(path)) &&
          provenance(path, sidecar = elsewhere) == r
    @test P.recordat(path, sidecar = elsewhere) == r &&
          isnothing(P.recordat(path))
    @test P.verifyfile(path, sidecar = elsewhere) && !P.verifyfile(path)
    # A second call finds the record where it was put.
    @test record_provenance(path, facts, role = :species,
                            sidecar = elsewhere) ==
          r
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
    @test_throws ArgumentError record_provenance(joinpath(dir, "missing.csv"),
                                                 facts,
                                                 role = :species)
    path = tablefile(dir, "theirs.csv")
    theirs = "writer = \"SomeoneElse 1.0\"\nfile = \"theirs.csv\"\n"
    write(P.sidecarpath(path), theirs)
    @test_throws "neither read from it nor recorded over it" record_provenance(path,
                                                                               facts,
                                                                               role = :species)
    @test read(P.sidecarpath(path), String) == theirs
    # Named another place, the record is written there.
    r = record_provenance(path, facts, role = :species,
                          sidecar = joinpath(dir, "ours.toml"))
    @test r.source == "x" && r.primary === false &&
          r.sidecar == joinpath(dir, "ours.toml")
    # A record a downloader of your own wrote, fetch time and all, is returned as any record of the
    # same source under the same role is.
    fetched = tablefile(dir, "fetched.csv")
    P.writesidecar(fetched,
                   merge!(P.baserecord(fetched, :species),
                          P.sourceentries(facts, PublishedByProvider())))
    @test !isnothing(provenance(fetched).fetched) &&
          !P.iscatalogued(provenance(fetched))
    @test record_provenance(fetched, facts, role = :species) ==
          provenance(fetched)
    # A file EcoSISTEM fetched from its catalogue - its record naming a catalogued dataset, with a
    # fetch time - is refused whatever it is given; a file cut from the same source and recorded
    # with its facts, which has no fetch time, is not.
    catalogued = SourceFacts(EcoSISTEM.datasetinfo(EcoSISTEM.NaturalEarthLevel))
    held = tablefile(dir, "ne.csv")
    P.writesidecar(held,
                   merge!(P.baserecord(held, :region),
                          P.sourceentries(catalogued, PublishedByProvider())))
    @test P.iscatalogued(provenance(held))
    @test_throws "fetched from the `EcoSISTEM` catalogue" record_provenance(held,
                                                                            catalogued,
                                                                            role = :region)
    @test_throws "fetched from the `EcoSISTEM` catalogue" record_provenance(held,
                                                                            "10.1/x",
                                                                            role = :region)
    cut = tablefile(dir, "cut.csv")
    r = record_provenance(cut, catalogued, role = :region,
                          route = SuppliedByUser())
    @test !P.iscatalogued(r) && isnothing(r.fetched) &&
          record_provenance(cut, catalogued, role = :region) == r
    # The mark is the catalogue's name with a fetch time, read from the record alone: a sidecar
    # from before carries no name, and is not refused.
    @test r.catalogue == "EcoSISTEM" &&
          provenance(held).catalogue == "EcoSISTEM"
    @test !P.iscatalogued(ProvenanceRecord(role = :habitat, source = "x",
                                           catalogue = "EcoSISTEM")) &&
          P.iscatalogued(ProvenanceRecord(role = :habitat, source = "x",
                                          catalogue = "EcoSISTEM",
                                          fetched = "2026-01-01T00:00:00Z")) &&
          !P.iscatalogued(ProvenanceRecord(role = :habitat,
                                           source = "WorldClim{BioClim}",
                                           fetched = "2026-01-01T00:00:00Z"))
    @test length(methods(P.iscatalogued)) == 1 &&
          only(methods(P.iscatalogued)).module === P
    # A record naming its source by a citation key, or by name alone, is compared by that.
    keyed = tablefile(dir, "keyed.csv")
    bykey = SourceFacts(source = "Keyed", key = "keyed", identifies = :data)
    @test record_provenance(keyed, bykey, role = :species).key == "keyed"
    @test record_provenance(keyed, bykey, role = :species).key == "keyed"
    @test_throws "from the citation `keyed` rather than the citation `other`" record_provenance(keyed,
                                                                                                SourceFacts(bykey,
                                                                                                            key = "other"),
                                                                                                role = :species)
    @test_throws "rather than `Named`" record_provenance(keyed,
                                                         SourceFacts(source = "Named"),
                                                         role = :species)
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
        r = record_provenance(path, ENVELOPES_DOI, role = :species)
        @test r.doi == ENVELOPES_DOI && r.primary === true &&
              r.licence == "ogl-uk-3.0" &&
              r.source == "Bioclimatic envelopes for global plant species"
    else
        @info "DOI lookup not exercised live: on a CI runner, or ECOSISTEM_DOI_LOOKUP=false"
        @test true
    end
end

@testset "a track with no run records as an untracked call does, and another handle is refused" begin
    dir = mktempdir()
    path = tablefile(dir, "t.csv")
    facts = SourceFacts(source = "A table", doi = "10.1/x", title = "A table",
                        identifies = :data, namespace = "Tables")
    r = record_provenance(path, facts, role = :species,
                          track = Track(nothing, "tables"))
    @test r isa ProvenanceRecord && r.doi == "10.1/x" &&
          isfile(P.sidecarpath(path))
    @test r.location == path && r.sidecar == P.sidecarpath(path) &&
          r.namespace == "Tables" && r.identifies === :data
    # The record keeps what a later registration needs beyond its entries.
    t = TOML.parsefile(P.sidecarpath(path))
    @test t["identifies"] == "data" && t["namespace"] == "Tables" &&
          !haskey(t, "key")
    # A handle that is no pipeline run is refused by name, before anything is written.
    other = tablefile(dir, "u.csv")
    @test_throws "FAIR data pipeline run" record_provenance(other, facts,
                                                            role = :species,
                                                            track = Track(Ref(0),
                                                                          "tables"))
    @test_throws "FAIR data pipeline run" record_provenance(other, "10.1/x",
                                                            role = :species,
                                                            track = Track(Ref(0),
                                                                          "tables"))
    @test !isfile(P.sidecarpath(other))
    # What is still to do for a file of your own: nothing while it is there, untracked.
    @test isempty(EcoSISTEM.missingfiles(path)) &&
          isempty(EcoSISTEM.missingfiles(path,
                                         track = Track(nothing, "tables")))
    gone = joinpath(dir, "gone.csv")
    todo = EcoSISTEM.missingfiles(gone)
    @test only(todo).file == gone && occursin("not there", only(todo).reason)
    @test_throws "FAIR data pipeline run" EcoSISTEM.missingfiles(path,
                                                                 track = Track(Ref(0),
                                                                               "tables"))
end
end
