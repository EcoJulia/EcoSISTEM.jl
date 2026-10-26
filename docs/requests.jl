# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The documentation page quoting what each data provider asks of anyone using its data, made from the
# request files in `data/catalogue/requests/`: one section per file, anchored by the file's name,
# which is where `EcoSISTEM._requestpage` points a fetched file's record. `docs/make.jl` writes it to
# `docs/src/requests.md` before building, and `test/clean_CatalogueRequests.jl` checks it against
# the request files and the package's link to it.

using EcoSISTEM

# Where the page's "View source" link goes, the page itself being generated.
const REQUESTS_SOURCE = "https://github.com/EcoJulia/EcoSISTEM.jl/tree/main/data/catalogue/requests"

# Text quoted from a provider, with the characters Markdown would read as markup escaped.
mdescape(s::AbstractString) = replace(s, r"([\\`*_\[\]<>$])" => s"\\\1")

# The name a request's section goes by: its rows' title where they share one, else the words their
# titles begin with, as EarthEnv's two land-cover versions do.
function sectiontitle(titles)
    allequal(titles) && return first(titles)
    words = split.(titles, ' ')
    n = 0
    while all(w -> length(w) > n, words) && allequal(w[n + 1] for w in words)
        n += 1
    end
    n == 0 && return join(titles, "; ")
    return rstrip(join(first(first(words), n), ' '), [',', ';', ':'])
end

# What the page opens with: what it quotes, how strongly a notice can be asked for, and where it is
# made from.
const REQUESTS_INTRO = """
```@meta
EditURL = "$REQUESTS_SOURCE"
```

# Providers' requests

What each data provider asks of anyone using its data, quoted from the provider's own pages: the
works to cite, and the notices to carry with the data. Every file the package fetches records the
same in its provenance (see [Provenance](provenance.md)), with a link to its section here. Where a
text holds `DD-MMM-YYYY` or `YYYY`, a record fills in the date or the year the file was fetched;
here they are shown as the provider writes them.

A notice is *required* where the provider or the data's licence says it must accompany the data,
*requested* where the provider asks for it, and *optional* where the provider offers it for those
who wish to use it.

The page is made from the files in `data/catalogue/requests/`, a section for each.

"""

# A licence as the page shows it: an SPDX identifier as code, a provider's terms as a link.
showlicence(licence) = startswith(licence, "http") ? "<$licence>" : "`$licence`"

# The page's Markdown: a section per request file, in the order `datasets.csv` first names them.
function requestspage()
    rows = EcoSISTEM._datasets()
    io = IOBuffer()
    print(io, REQUESTS_INTRO)
    names = unique(r.request for r in rows)
    sections = [filter(r -> r.request == name, rows) for name in names]
    titles = [sectiontitle([r.title for r in named]) for named in sections]
    for (name, named, title) in zip(names, sections, titles)
        request = EcoSISTEM._readrequest(name, first(named).dataset)
        # A title two sections share, as CHELSA V2.1's two requests do, is told apart by the
        # datasets each covers.
        types = join(("`$(r.dataset)`" for r in named), ", ")
        heading = count(==(title), titles) > 1 ? "$(mdescape(title)) ($types)" :
                  mdescape(title)
        println(io, "## [$heading](@id $name)\n")
        for r in named
            licence = showlicence(r.licence)
            println(io,
                    "- `$(r.dataset)`: $(mdescape(r.title)), licence $licence")
        end
        pages = join(("<$p>" for p in request.pages), " and ")
        when = isnothing(request.read) ? "" : " on $(request.read)"
        println(io, "\nRead from $pages$when.\n")
        if !isempty(request.references)
            println(io, "**Cite**\n")
            for ref in request.references
                doi = occursin(ref.doi, ref.text) ? "" :
                      " ([doi:$(ref.doi)](https://doi.org/$(ref.doi)))"
                println(io, "- ", mdescape(ref.text), doi)
            end
            println(io)
        end
        if !isempty(request.notices)
            println(io, "**Notices**\n")
            for notice in request.notices
                println(io, "- *$(notice.strength)*: $(mdescape(notice.text))")
            end
            println(io)
        end
    end
    return String(take!(io))
end
