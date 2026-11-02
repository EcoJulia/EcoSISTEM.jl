# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/Provenancer/SourceFacts.jl` and `src/Provenancer/doi.jl`: what is known of a source,
# whether a file of it is primary, the entries of a record written from it, and the facts the
# metadata of a DOI gives, read from the saved replies in `test/doi/`. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_SourceFacts.jl"])'

module TestSourceFacts

using EcoSISTEM
using EcoSISTEM.Provenancer
using RasterDataSources
using Test

const P = EcoSISTEM.Provenancer
include("doi/replies.jl")

# The facts of a table, as a user writes them.
function tablefacts()
    return SourceFacts(source = "My table", doi = "10.1234/abc",
                       title = "A table",
                       authors = ["A. Person"],
                       identifiers = Dict("A. Person" => "https://orcid.org/0000-0000-0000-0000"),
                       releasedate = "2024-05", licence = "CC-BY-4.0",
                       identifies = :data,
                       notices = [(strength = :required, text = "Say so")])
end

@testset "the facts of a source are kept as given, blank for unknown" begin
    f = tablefacts()
    @test f.source == "My table" && f.doi == "10.1234/abc" && f.url == "" &&
          f.version == "" && f.publisher == ""
    @test f.identifies === :data && only(f.notices).strength === :required
    @test f.identifiers["A. Person"] == "https://orcid.org/0000-0000-0000-0000"
    @test f == tablefacts() && hash(f) == hash(tablefacts())
    @test f != SourceFacts(source = "Other")
    @test_throws ArgumentError SourceFacts(source = "x", identifies = :dataset)
    @test_throws "strength" SourceFacts(source = "x",
                                        notices = [(strength = :must,
                                                    text = "t")])
    @test string(f) == "SourceFacts(\"My table\", doi 10.1234/abc)"
    shown = sprint(show, MIME("text/plain"), f)
    @test occursin("(data)", shown) && occursin("A. Person", shown) &&
          occursin("CC-BY-4.0", shown)
    @test P.IDENTIFIES == (:data, :paper, :software, :other)
end

@testset "a file is primary by what the identifier names and the route it came by" begin
    data = SourceFacts(source = "d", doi = "10.1/d", identifies = :data)
    paper = SourceFacts(source = "p", doi = "10.1/p", identifies = :paper)
    @test P.isprimary(data, PublishedByProvider())
    @test !P.isprimary(paper, PublishedByProvider())
    # An extract is the cited item only where the service gave it an identifier of its own.
    @test !P.isprimary(data, ExtractedByService())
    @test P.isprimary(data, ExtractedByService("10.1/d"))
    @test !P.isprimary(paper, ExtractedByService("10.1/p"))
    @test !P.isprimary(data, SuppliedByUser())
    @test string(ExtractedByService()) == "ExtractedByService()"
    @test string(ExtractedByService("10.1/d")) ==
          "ExtractedByService(\"10.1/d\")"
    @test ExtractedByService("10.1/d") == ExtractedByService("10.1/d") &&
          isnothing(ExtractedByService().identifier)
end

@testset "a source's entries in a file's record" begin
    facts = SourceFacts(source = "WorldClim{BioClim}", doi = "10.1002/joc.5086",
                        url = "https://www.worldclim.org/",
                        title = "WorldClim 2.1",
                        releasedate = "2020-01", version = "2.1",
                        licence = "https://www.worldclim.org/about.html",
                        identifies = :paper,
                        citation = "Fick and Hijmans (2017)",
                        references = [(doi = "10.1002/joc.5086",
                                       text = "Fick, S.E.")],
                        notices = [(strength = :requested, text = "Cite it")],
                        requestpage = "https://docs/requests/#worldclim")
    e = P.sourceentries(facts, PublishedByProvider(),
                        title = "WorldClim 2.1, bio1 (wc2.1_10m_bio_1)",
                        code = "bio1")
    @test e["source"] == "WorldClim{BioClim}" &&
          e["title"] == "WorldClim 2.1, bio1 (wc2.1_10m_bio_1)" &&
          e["primary"] === false
    @test e["doi"] == "10.1002/joc.5086" && e["licence"] == facts.licence &&
          e["version"] == "2.1" && e["release_date"] == "2020-01" &&
          e["citation"] == "Fick and Hijmans (2017)" &&
          e["request_page"] == facts.requestpage && e["code"] == "bio1"
    @test e["references"] ==
          [Dict("doi" => "10.1002/joc.5086", "text" => "Fick, S.E.")] &&
          e["notices"] == [Dict("strength" => "requested", "text" => "Cite it")]
    # Blank facts write no key, and the title defaults to the source's.
    bare = P.sourceentries(SourceFacts(source = "s", title = "T"),
                           SuppliedByUser())
    @test Set(keys(bare)) == Set(["source", "title", "primary"]) &&
          bare["title"] == "T" &&
          bare["primary"] === false
    # A citation is composed where the source gives none and has a DOI and authors.
    composed = P.sourceentries(SourceFacts(source = "s", doi = "10.1/x",
                                           title = "A set",
                                           authors = ["A. One", "B. Two"],
                                           releasedate = "2022-06-22",
                                           publisher = "Pub",
                                           version = "2", identifies = :data),
                               PublishedByProvider())
    @test composed["citation"] ==
          "A. One, B. Two (2022). A set (version 2). Pub. https://doi.org/10.1/x"
    @test composed["primary"] === true
    @test !haskey(P.sourceentries(SourceFacts(source = "s", doi = "10.1/x"),
                                  PublishedByProvider()), "citation")
    # The catalogue's record of a file is these entries, from the row's facts, with the file's
    # own title and layer; a type with no row gives none.
    layer = EcoSISTEM.layerinfo(WorldClim{BioClim}, 1)
    rec = EcoSISTEM._cataloguerecord(WorldClim{BioClim}, layer,
                                     "/x/wc2.1_10m_bio_1.tif",
                                     route = PublishedByProvider())
    @test rec["source"] == "WorldClim{BioClim}" &&
          rec["code"] == first(layer.aliases) &&
          rec["primary"] === false && rec["doi"] == "10.1002/joc.5086" &&
          endswith(rec["title"], "(wc2.1_10m_bio_1)")
    @test isempty(EcoSISTEM._cataloguerecord(Int, nothing, "/x/y.tif",
                                             route = PublishedByProvider()))
    era = EcoSISTEM._sourcefacts(EcoSISTEM.datasetinfo(ERA))
    @test era.source == "ERA" && era.identifies === :data &&
          era.key == "era5" &&
          era.doi == "10.24381/cds.f17050d7" && era.licence == "CC-BY-4.0" &&
          startswith(era.requestpage, "https://")
end

@testset "a DOI's DataCite record gives the facts of a dataset" begin
    f = lookup_doi("10.5285/ca339c86-3674-4030-b891-35326e71141e",
                   fetch = saved)
    @test f.title == "Bioclimatic envelopes for global plant species" &&
          f.source == f.title
    @test f.doi == "10.5285/ca339c86-3674-4030-b891-35326e71141e"
    @test f.identifies === :data && f.releasedate == "2022-06-22" &&
          f.version == ""
    @test f.licence == "ogl-uk-3.0"
    @test length(f.authors) == 4 && occursin("Harris", first(f.authors)) &&
          occursin("Reeve", last(f.authors))
    @test f.identifiers[first(f.authors)] ==
          "https://orcid.org/0000-0003-0852-2340"
    @test startswith(f.publisher, "NERC") &&
          startswith(f.url, "https://catalogue.ceh.ac.uk/")
    @test startswith(f.description,
                     "Bioclimatic envelopes for over 200,000 plant species")
    # The same facts whichever way the DOI is written.
    @test lookup_doi("doi:10.5285/ca339c86-3674-4030-b891-35326e71141e",
                     fetch = saved) == f
    @test lookup_doi("https://doi.org/10.5285/ca339c86-3674-4030-b891-35326e71141e",
                     fetch = saved) == f
    # An occurrence download: an extract with a DOI of its own, primary by it.
    g = lookup_doi("10.15468/dl.n4a2sa", fetch = saved)
    @test g.title == "Occurrence Download" && g.identifies === :data &&
          g.licence == "cc-by-nc-4.0" && g.releasedate == "2026-09-12"
    @test startswith(g.description,
                     "A dataset containing 17208657 species occurrences")
    @test g.publisher == "The Global Biodiversity Information Facility"
    @test P.isprimary(g, ExtractedByService(g.doi)) &&
          !P.isprimary(g, ExtractedByService())
    # A creator with given and family names, the stray space in the given name dropped; the year
    # alone.
    m = lookup_doi("10.5067/MODIS/MOD13C2.061", fetch = saved)
    @test m.authors == ["Kamel Didan"] && m.releasedate == "2021" &&
          m.identifies === :data && m.licence == ""
    @test startswith(m.publisher, "NASA")
    # An organisation as creator, and no rights at all.
    e = lookup_doi("10.24381/cds.f17050d7", fetch = saved)
    @test e.authors == ["Copernicus Climate Change Service"] &&
          e.releasedate == "2019" &&
          e.licence == "" && e.identifies === :data
    @test startswith(e.title, "ERA5 monthly averaged data on single levels")
end

@testset "a paper's DOI has no DataCite record, and its CSL record is read" begin
    p = lookup_doi("10.1111/geb.13564", fetch = saved)
    @test p.identifies === :paper &&
          p.title ==
          "Strong phylogenetic signals in global plant bioclimatic envelopes"
    @test first(p.authors) == "Claire Harris" &&
          p.identifiers["Claire Harris"] ==
          "https://orcid.org/0000-0003-0852-2340"
    @test p.releasedate == "2022-07-20" && p.publisher == "Wiley" &&
          p.licence == "" &&
          p.doi == "10.1111/geb.13564" &&
          p.url == "http://dx.doi.org/10.1111/geb.13564"
    @test !P.isprimary(p, PublishedByProvider())
    e = P.sourceentries(p, PublishedByProvider())
    @test startswith(e["citation"], "Claire Harris, ") &&
          endswith(e["citation"], "https://doi.org/10.1111/geb.13564") &&
          occursin("(2022). Strong phylogenetic", e["citation"])
end

@testset "what is no DOI, or has no metadata, is refused by name" begin
    @test_throws "is not a DOI" lookup_doi("not a doi", fetch = saved)
    @test_throws "no metadata is to be had" lookup_doi("10.9999/nothing",
                                                       fetch = (u,
                                                                t) -> nothing)
    @test P._trydoi(" 10.1234/x ") == "10.1234/x" &&
          isnothing(P._trydoi("10.12/short"))
    @test P._trydoi("http://dx.doi.org/10.1111/geb.13564") ==
          "10.1111/geb.13564"
    @test P._isodate("2019-04-18T12:00:00Z") == "2019-04-18" &&
          P._isodate("2019/2020") == "2019" && isnothing(P._isodate("April"))
    @test P._dateasfaras(Any[2022, 7]) == "2022-07" &&
          P._dateasfaras(Any[2022]) == "2022" &&
          P._dateasfaras(Any[]) == ""
    @test P._markupless("<jats:title>Abstract</jats:title><p>Text &amp; more</p>") ==
          "Text & more"
    @test P._identifies("JournalArticle") === :paper &&
          P._identifies("Software") === :software &&
          P._identifies("Event") === :other && P._identifies(nothing) === :other
    @test P._cslidentifies("article-journal") === :paper &&
          P._cslidentifies("dataset") === :data &&
          P._cslidentifies("journal-article") === :paper &&
          P._cslidentifies("other") === :other
    @test P._creatorname(Dict("name" => "Reeve, R.", "nameType" => "Personal")) ==
          "R. Reeve"
    @test P._creatorname(Dict("givenName" => "Richard",
                              "familyName" => "Reeve")) ==
          "Richard Reeve"
    @test P._creatorname(Dict("name" => "GBIF.org User",
                              "nameType" => "Organizational")) ==
          "GBIF.org User"
    @test isnothing(P._creatorname(Dict("nameType" => "Personal")))
    @test P._creatoridentifier(Dict("nameIdentifiers" =>
                                        [Dict("nameIdentifierScheme" => "ORCID",
                                              "nameIdentifier" => "0000-0003-2589-8091")])) ==
          "https://orcid.org/0000-0003-2589-8091"
    @test isnothing(P._creatoridentifier(Dict("nameIdentifiers" =>
                                                  [Dict("nameIdentifierScheme" => "Other",
                                                        "nameIdentifier" => "")])))
    @test P._datacitelicence([Dict("rightsUri" => "https://x.org/l")]) ==
          "https://x.org/l" &&
          isnothing(P._datacitelicence([]))
end

# The live lookup runs off a CI runner, and not with ECOSISTEM_DOI_LOOKUP=false.
function livelookup()
    haskey(ENV, "RUNNER_OS") && return false
    return get(ENV, "ECOSISTEM_DOI_LOOKUP", "true") != "false"
end

@testset "the resolver answers today as the saved replies do" begin
    if livelookup()
        doi = "10.5285/ca339c86-3674-4030-b891-35326e71141e"
        live = lookup_doi(doi)
        kept = lookup_doi(doi, fetch = saved)
        @test live.title == kept.title && live.doi == kept.doi &&
              live.identifies === :data
        @test lookup_doi("10.1111/geb.13564").identifies === :paper
        # An address that cannot be reached errors naming the way round.
        @test_throws "cannot reach" P._fetchmetadata("https://doi.invalid/10.1/x",
                                                     "application/vnd.datacite.datacite+json")
    else
        @info "DOI lookup not exercised live: on a CI runner, or ECOSISTEM_DOI_LOOKUP=false"
        @test true
    end
end

end
