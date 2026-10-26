# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Do the catalogue's provider requests, citations and authors agree with one another, and with the
# documentation's page quoting the requests? Run with the other hygiene checks:
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_clean.jl"])'
#
# `data/catalogue/requests/` quotes what each provider asks of anyone using its data, and each row of
# `datasets.csv` names its request; `citations.bib` and `authors.toml`, made by
# `data/src/citations.jl`, are what a registration of a catalogued file cites. Reads files only: no
# network. It cannot see whether a quoted text is still the provider's own, nor whether doi.org
# still says what `citations.bib` holds; running `data/src/citations.jl` answers the second.

module TestCleanCatalogueRequests

using BibParser: BibParser
using EcoSISTEM
using Test
using TOML

const CATALOGUE = joinpath(pkgdir(EcoSISTEM), "data", "catalogue")

# `requestspage`, which `docs/make.jl` writes the documentation's page of requests with.
include(joinpath(pkgdir(EcoSISTEM), "docs", "requests.jl"))

# Every entry of `citations.bib`. A month written as a bare word, as Crossref writes `month=May`,
# is braced first, as DataPipeline braces it before reading: BibParser takes a bare word for an
# abbreviation, and none is defined.
function bibentries()
    text = read(joinpath(CATALOGUE, "citations.bib"), String)
    text = replace(text,
                   r"(\bmonth\s*=\s*)([A-Za-z]+)(\s*[,}])"i => s"\1{\2}\3")
    return BibParser.parse_bibliography(IOBuffer(text), format = :BibTeX,
                                        check = :none).entries
end

# An author as a citation names one: its given names, particle, family name and suffix, as
# DataPipeline reads a BibTeX author, so `authors.toml` must spell a name this way to be matched.
function personname(name)
    parts = (name.first, name.middle, name.particle, name.last, name.junior)
    return replace(join(filter(!isempty, collect(parts)), " "), r"[{}]" => "")
end

@testset "the catalogue's requests, citations and authors agree" begin
    # Loading the table reads every row's request, refusing a missing file or a misshapen one.
    rows = EcoSISTEM._datasets()
    files = filter(endswith(".toml"), readdir(joinpath(CATALOGUE, "requests")))
    names = [first(splitext(f)) for f in files]
    entries = bibentries()
    keys_ = [e.raw.key for e in entries]
    bibdois = Set(lowercase(e.canonical.access.doi)
                  for e in entries if !isempty(e.canonical.access.doi))

    @testset "each request file is a row's, named as a documentation anchor" begin
        # Lower-case ASCII and distinct ignoring case, so no file system can merge two.
        for n in names
            @test occursin(r"^[a-z0-9]+(-[a-z0-9]+)*$", n)
        end
        @test allunique(lowercase.(names))
        @test Set(r.request for r in rows) == Set(names)
    end

    @testset "each DOI has its citation, and a written DOI is the reference's own" begin
        for r in rows
            # A row with no DOI is cited by the entry keyed by its request's name.
            if isempty(r.doi)
                @test r.request in keys_
            else
                @test lowercase(r.doi) in bibdois
            end
            for ref in r.references
                @test lowercase(ref.doi) in bibdois
                # A provider's text need not write its DOI, but one it writes must be this one.
                for m in eachmatch(r"10\.\d{4,9}/[^\s\"<>]+", ref.text)
                    @test lowercase(rstrip(m.match, ['.', ',', ')'])) ==
                          lowercase(ref.doi)
                end
            end
        end
        @test allunique(keys_)
    end

    @testset "only the known placeholders are left in the providers' texts" begin
        for r in rows
            for text in vcat([x.text for x in r.references],
                             [x.text for x in r.notices])
                for m in eachmatch(r"\b[DMY]{2,}(?:-[DMY]{2,})*\b", text)
                    @test m.match in EcoSISTEM._REQUEST_PLACEHOLDERS
                end
                # The providers' own forms of a year are written as `YYYY`.
                @test !occursin(r"<[A-Z0-9_]+>|\[year\]"i, text)
            end
        end
    end

    @testset "each licence is an SPDX identifier or the provider's URL" begin
        for r in rows
            @test occursin(r"^https?://\S+$", r.licence) ||
                  occursin(r"^[A-Za-z0-9][A-Za-z0-9.+-]*$", r.licence)
        end
    end

    @testset "each identified author is an author of a citation" begin
        authors = TOML.parsefile(joinpath(CATALOGUE, "authors.toml"))
        @test !isempty(authors)
        cited = Set(personname(n) for e in entries for n in e.canonical.authors)
        for (name, id) in authors
            @test name in cited
            @test occursin(r"^https://(orcid\.org/\d{4}-\d{4}-\d{4}-\d{3}[\dX]|ror\.org/0[a-z0-9]{8})$",
                           id)
        end
    end

    # The Markdown is checked, so whether Documenter renders each anchor is the docs build's to show.
    @testset "the documentation quotes each request where a record's link points" begin
        page = requestspage()
        for n in names
            @test count("(@id $n)", page) == 1
        end
        for r in rows
            @test EcoSISTEM._requestpage(r) ==
                  EcoSISTEM._REQUESTS_PAGE * "#" * r.request
            for ref in r.references
                @test occursin(mdescape(ref.text), page)
            end
            for notice in r.notices
                @test occursin("*$(notice.strength)*: " * mdescape(notice.text),
                               page)
            end
        end
        # The link is the page `docs/make.jl` builds, under the site's canonical address.
        make = read(joinpath(pkgdir(EcoSISTEM), "docs", "make.jl"), String)
        @test occursin("\"requests.md\"", make)
        canonical = match(r"canonical = \"([^\"]+)\"", make)
        @test !isnothing(canonical) &&
              EcoSISTEM._REQUESTS_PAGE ==
              rstrip(canonical[1], '/') * "/requests/"
    end
end

end
