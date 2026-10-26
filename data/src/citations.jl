# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Make the catalogue's citations and authors from what doi.org says of each DOI the catalogue names:
# `data/catalogue/citations.bib`, the BibTeX of every DOI in `datasets.csv` and in the request files
# under `data/catalogue/requests/`, and `data/catalogue/authors.toml`, the identifier (ORCID or ROR)
# of each author that the DOIs' metadata identifies. A registration of a catalogued file cites its
# source from these, so it needs no network, and every user registers the same facts.
#
# Compare what doi.org says now with the files, writing nothing:
#
#     julia --project=data/src data/src/citations.jl
#
# Write the files:
#
#     julia --project=data/src data/src/citations.jl --write
#
# Each entry is kept as doi.org returns it. Anything written here rather than fetched is in this
# script: the one entry with no DOI, Natural Earth's, keyed by its request file's name, which is how
# a catalogue row with no DOI names its citation. Needs the network.

module CitationsBuild

using CSV: CSV
using Downloads: Downloads
using JSON: JSON
using TOML: TOML

const CATALOGUE = joinpath(dirname(dirname(@__DIR__)), "data", "catalogue")
const BIBFILE = joinpath(CATALOGUE, "citations.bib")
const AUTHORSFILE = joinpath(CATALOGUE, "authors.toml")

# The citation of the catalogue's one source with no DOI, written by hand from its terms of use
# (`requests/naturalearth.toml`), its authors those the site names as the primary ones.
const NATURALEARTH = """
@misc{naturalearth,
  author = {Patterson, Tom and Kelso, Nathaniel Vaughn},
  title = {Natural Earth},
  publisher = {Natural Earth},
  url = {https://www.naturalearthdata.com/}
}"""

# Every DOI the catalogue names, in lower case, sorted: each row's `DOI`, and each reference of each
# request file.
function cataloguedois()
    dois = Set{String}()
    for row in CSV.File(joinpath(CATALOGUE, "datasets.csv"))
        ismissing(row.DOI) || push!(dois, lowercase(strip(row.DOI)))
    end
    for file in readdir(joinpath(CATALOGUE, "requests"), join = true)
        endswith(file, ".toml") || continue
        for reference in get(TOML.parsefile(file), "reference", [])
            push!(dois, lowercase(reference["doi"]))
        end
    end
    return sort!(collect(dois))
end

# The body of a doi.org answer for `url` in the media type `accept`, as text.
function fetchtext(url::AbstractString, accept::AbstractString)
    io = IOBuffer()
    Downloads.request(url, output = io, headers = ["Accept" => accept],
                      throw = true)
    return String(take!(io))
end

# The registration agency of a DOI, `"Crossref"` or `"DataCite"`, as doi.org reports it.
function agency(doi::AbstractString)
    return only(JSON.parse(fetchtext("https://doi.org/ra/$doi",
                                     "application/json")))["RA"]
end

# A DOI's BibTeX as doi.org returns it, refused where the entry names another DOI.
function bibtex(doi::AbstractString)
    entry = strip(fetchtext("https://doi.org/$doi", "application/x-bibtex"))
    startswith(entry, '@') ||
        error("doi.org returned no BibTeX entry for $doi: $(first(entry, 200))")
    named = match(r"\bdoi\s*=\s*\{([^}]*)\}"i, entry)
    (isnothing(named) || lowercase(named[1]) == doi) ||
        error("the BibTeX entry for $doi names another DOI, $(named[1])")
    return String(entry)
end

# The authors a DOI's metadata identifies, as `name => identifier`: from Crossref's CSL record, each
# author with an ORCID; from DataCite's own record, each creator with an ORCID or a ROR id. A person
# is named by given names then family name, as a citation's author is read.
function identifiedauthors(doi::AbstractString)
    found = Pair{String, String}[]
    if agency(doi) == "DataCite"
        record = JSON.parse(fetchtext("https://doi.org/$doi",
                                      "application/vnd.datacite.datacite+json"))
        for creator in get(record, "creators", [])
            ids = [String(id["nameIdentifier"])
                   for id in get(creator, "nameIdentifiers", [])
                   if occursin(r"^https://(orcid|ror)\.org/",
                               id["nameIdentifier"])]
            isempty(ids) && continue
            name = haskey(creator, "familyName") ?
                   join(strip.([get(creator, "givenName", ""),
                                   creator["familyName"]]), " ") :
                   creator["name"]
            push!(found, normalisename(name) => first(ids))
        end
    else
        record = JSON.parse(fetchtext("https://doi.org/$doi",
                                      "application/vnd.citationstyles.csl+json"))
        for author in get(record, "author", [])
            haskey(author, "ORCID") || continue
            name = join(strip.([get(author, "given", ""), author["family"]]),
                        " ")
            orcid = replace(author["ORCID"],
                            r"^https?://orcid\.org/" => "https://orcid.org/")
            push!(found, normalisename(name) => orcid)
        end
    end
    return found
end

# A name with its whitespace collapsed to single spaces.
normalisename(name::AbstractString) = join(split(name), " ")

# What doi.org says now: the BibTeX file's text and the authors' identifiers, an author identified
# by two DOIs with different identifiers being an error naming both.
function build()
    entries = String[]
    authors = Dict{String, String}()
    for doi in cataloguedois()
        println("fetching ", doi)
        push!(entries, bibtex(doi))
        for (name, id) in identifiedauthors(doi)
            known = get!(authors, name, id)
            known == id ||
                error("$name is identified as both $known and $id; one of them is wrong")
        end
    end
    push!(entries, NATURALEARTH)
    header = "% Generated by data/src/citations.jl from doi.org; regenerate rather than edit.\n\n"
    return (bib = header * join(entries, "\n\n") * "\n", authors = authors)
end

# The authors' file's text: a header, then one line per author in order of name.
function authorstext(authors::AbstractDict)
    io = IOBuffer()
    println(io,
            "# Generated by data/src/citations.jl from doi.org's metadata; regenerate rather than edit.")
    println(io)
    TOML.print(io, authors, sorted = true)
    return String(take!(io))
end

# Script entry point: compare by default, `--write` to write the two files.
function main(args)
    made = build()
    files = (BIBFILE => made.bib, AUTHORSFILE => authorstext(made.authors))
    if "--write" in args
        for (path, text) in files
            write(path, text)
            println("wrote ", path)
        end
        return 0
    end
    stale = [path
             for (path, text) in files
             if !isfile(path) || read(path, String) != text]
    isempty(stale) && (println("both files agree with doi.org"); return 0)
    println("differ from doi.org now: ", join(stale, ", "),
            " - re-run with `--write` to replace them")
    return 1
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(CitationsBuild.main(ARGS))
end
