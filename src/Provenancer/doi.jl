# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What a DOI resolves to, read into the facts of a source: DataCite's record of it where the DOI
# is DataCite's (a dataset's, usually), else the Citation Style Language record every registration
# agency serves (a paper's).

"""
    lookup_doi(doi::AbstractString; fetch = <asks the network>)

Return the [`SourceFacts`](@ref) the metadata of `doi` gives, the DOI written bare
(`10.5285/...`), as `doi:10.5285/...` or as a `doi.org` address. The resolver is asked for
DataCite's record first, which a dataset's DOI has: its main title; its creators, a person as
given names then family name, an organisation as it is written, with the ORCID or ROR identifier
each is given; its publisher; the first of its Issued, Available, Created and Submitted dates,
else its year; its version; its abstract; its licence, as the SPDX identifier a right gives,
else the address of the first right that has one; and what the DOI identifies, from its resource
type. Where there is none, the Citation Style Language record every DOI has is read, a paper's
usually: title, authors with their ORCIDs, publisher, the date issued, version and abstract; a
CSL record's licence is the publisher's terms for the work, not for data, and is not taken. The
facts' `source` is the title.

Errors where neither record is to be had, and where the resolver cannot be reached, naming facts
written by hand as the way round: nothing is looked up without a network.

# Arguments

  - `doi`: the DOI, in any of the three forms.
  - `fetch`: what fetches an address with an `Accept` type, returning the reply as text, or
    `nothing` where the server has no such form; by default the network is asked, and a test
    gives saved replies.
"""
function lookup_doi(doi::AbstractString; fetch::Function = _fetchmetadata)
    bare = _baredoi(doi)
    url = "https://doi.org/" * bare
    for (type, reader) in _METADATA_READERS
        reply = fetch(url, type)
        isnothing(reply) || return reader(JSON.parse(reply), bare)
    end
    return error("no metadata is to be had for $bare, as DataCite's record or as CSL: write " *
                 "the source's `SourceFacts` by hand.")
end

# A DOI as written: bare, as `doi:` text, or as a `doi.org` address, with or without `dx.`.
const _DOI = r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)?(10\.\d{4,9}/\S+)$"i

# A DOI written in any of its forms as the bare name, or `nothing` for text that is no DOI.
function _trydoi(text::AbstractString)
    found = match(_DOI, strip(text))
    return isnothing(found) ? nothing : String(only(found.captures))
end

# The bare DOI, or a refusal naming what was given.
function _baredoi(text::AbstractString)
    doi = _trydoi(text)
    isnothing(doi) &&
        throw(ArgumentError("`$text` is not a DOI: write it `10.xxxx/...`, `doi:10.xxxx/...` " *
                            "or as a `doi.org` address."))
    return doi
end

# The reply of an address asked for one `type` of metadata, as text: `nothing` where the server
# has no such form (406) or answers with something other than JSON; an error where it cannot be
# reached or refuses, naming facts written by hand as the way round.
function _fetchmetadata(url::AbstractString, type::AbstractString)
    buffer = IOBuffer()
    response = try
        Downloads.request(url, headers = ["Accept" => type], output = buffer,
                          downloader = _http11downloader(), timeout = 60,
                          throw = true)
    catch e
        e isa Downloads.RequestError || rethrow()
        error("cannot reach $url to look its metadata up ($(e.message)): with no network, " *
              "write the source's `SourceFacts` by hand.")
    end
    response.status == 406 && return nothing
    response.status == 200 ||
        error("no metadata is to be had from $url (HTTP status $(response.status)): write " *
              "the source's `SourceFacts` by hand.")
    occursin("json",
             lowercase(something(_headervalue(response.headers, "content-type"),
                                 ""))) ||
        return nothing
    return String(take!(buffer))
end

# The entries of a list in parsed JSON, none where it is missing or null.
function _items(record::AbstractDict, key::AbstractString)
    items = get(record, key, nothing)
    return items isa AbstractVector ? items : Any[]
end

# Text from metadata with its white space single, or `nothing` where that leaves nothing; a value
# that is not text - a year written as a number - as its text.
_plain(::Nothing) = nothing

function _plain(text::AbstractString)
    plain = strip(replace(text, r"\s+" => " "))
    return isempty(plain) ? nothing : String(plain)
end

_plain(value) = _plain(string(value))

# Text from metadata with its markup taken out - the HTML of a description, the JATS of an
# abstract, a heading of one dropped - and its white space single.
_markupless(::Nothing) = nothing

function _markupless(text::AbstractString)
    text = replace(text, r"<(jats:)?title>.*?</(jats:)?title>"is => " ")
    text = replace(text, r"<[^<>]*>" => " ")
    for (entity, character) in ("&lt;" => "<", "&gt;" => ">", "&quot;" => "\"",
        "&#39;" => "'", "&amp;" => "&")
        text = replace(text, entity => character)
    end
    return _plain(text)
end

# A text or `nothing` as a facts field takes it: text, blank for nothing.
_text(x) = something(_plain(x), "")

# The main title of a DataCite record: the first with no title type, else the first.
function _maintitle(record::AbstractDict)
    titles = _items(record, "titles")
    isempty(titles) && return nothing
    main = findfirst(t -> t isa AbstractDict &&
                          isnothing(get(t, "titleType", nothing)), titles)
    title = titles[something(main, 1)]
    return title isa AbstractDict ? _plain(get(title, "title", nothing)) :
           _plain(title)
end

# The abstract of a DataCite record, or `nothing`.
function _abstract(record::AbstractDict)
    for description in _items(record, "descriptions")
        description isa AbstractDict || continue
        get(description, "descriptionType", nothing) == "Abstract" || continue
        return _markupless(_plain(get(description, "description", nothing)))
    end
    return nothing
end

# A publisher as DataCite writes one, a table with a `name` in its newer records and text in its
# older ones.
function _publishername(publisher::AbstractDict)
    return _plain(get(publisher, "name", nothing))
end

_publishername(publisher) = _plain(publisher)

# A creator of a DataCite record as a reference list names an author: a person's given names and
# then family name, an organisation as it is written.
function _creatorname(creator::AbstractDict)
    given = _plain(get(creator, "givenName", nothing))
    family = _plain(get(creator, "familyName", nothing))
    !isnothing(given) && !isnothing(family) && return "$given $family"
    name = _plain(get(creator, "name", nothing))
    isnothing(name) && return nothing
    parts = split(name, ',', limit = 2)
    (get(creator, "nameType", nothing) == "Personal" && length(parts) == 2) ||
        return name
    return _plain(string(parts[2], " ", parts[1]))
end

# The ORCID or ROR identifier a DataCite record gives a creator, as an `https` address, or
# `nothing`.
function _creatoridentifier(creator::AbstractDict)
    for given in _items(creator, "nameIdentifiers")
        given isa AbstractDict || continue
        scheme = lowercase(something(_plain(get(given, "nameIdentifierScheme",
                                                nothing)), ""))
        scheme in ("orcid", "ror") || continue
        value = something(_plain(get(given, "nameIdentifier", nothing)), "")
        isempty(value) && continue
        if !occursin(r"^https?://", value)
            base = something(_plain(get(given, "schemeUri", nothing)),
                             "https://$scheme.org/")
            value = string(rstrip(base, '/'), "/", value)
        end
        occursin(r"^https?://\S+$", value) &&
            return replace(value, r"^http:" => "https:")
    end
    return nothing
end

# The kinds of date in a DataCite record that say when a source was released, in the order they
# are looked for.
const _RELEASE_DATE_TYPES = ("Issued", "Available", "Created", "Submitted")

# A date as far as metadata gives it, `YYYY`, `YYYY-MM` or `YYYY-MM-DD`, from text beginning so -
# a time of day or a range after it dropped - or `nothing`.
function _isodate(text)
    found = match(r"^(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?",
                  something(_plain(text), ""))
    return isnothing(found) ? nothing : String(found.match)
end

# The date a DataCite record gives for a source's release: the first of the kinds that say so, as
# far as it goes, else its year.
function _datacitedate(record::AbstractDict)
    for kind in _RELEASE_DATE_TYPES, entry in _items(record, "dates")
        entry isa AbstractDict || continue
        get(entry, "dateType", nothing) == kind || continue
        date = _isodate(get(entry, "date", nothing))
        isnothing(date) || return date
    end
    return _isodate(get(record, "publicationYear", nothing))
end

# The licence a DataCite record's rights give: the identifier of the first right in the SPDX
# scheme, as the record writes it, else the address of the first right that has one, else
# `nothing`.
function _datacitelicence(rights::AbstractVector)
    for right in rights
        right isa AbstractDict || continue
        scheme = something(_plain(get(right, "rightsIdentifierScheme", nothing)),
                           "")
        identifier = something(_plain(get(right, "rightsIdentifier", nothing)),
                               "")
        (lowercase(scheme) == "spdx" && !isempty(identifier)) &&
            return identifier
    end
    for right in rights
        right isa AbstractDict || continue
        address = something(_plain(get(right, "rightsUri", nothing)), "")
        occursin(r"^https?://\S+$", address) && return address
    end
    return nothing
end

# What a DataCite resource type says a DOI identifies: data, a work about it, software, or
# something else.
function _identifies(kind)
    kind in ("Dataset", "Collection") && return :data
    kind in
    ("Book", "BookChapter", "ConferencePaper", "ConferenceProceeding",
     "DataPaper",
     "Dissertation", "Journal", "JournalArticle", "PeerReview", "Preprint",
     "Report",
     "Standard", "Text") && return :paper
    kind in ("Software", "ComputationalNotebook", "Workflow") &&
        return :software
    return :other
end

# The same of a CSL record's type, in CSL's own names and in Crossref's, which writes
# `journal-article` for CSL's `article-journal`.
function _cslidentifies(kind)
    kind == "dataset" && return :data
    kind in
    ("article", "article-journal", "article-magazine", "article-newspaper",
     "book",
     "chapter", "document", "manuscript", "paper-conference", "periodical",
     "post",
     "post-weblog", "report", "review", "review-book", "standard", "thesis",
     "journal-article", "book-chapter", "proceedings-article", "posted-content",
     "monograph", "reference-book", "edited-book", "reference-entry",
     "dissertation",
     "journal", "peer-review") &&
        return :paper
    kind == "software" && return :software
    return :other
end

# The facts DataCite's record of a DOI gives of its source.
function _datacitefacts(record::AbstractDict, doi::AbstractString)
    authors = String[]
    identifiers = Dict{String, String}()
    for creator in _items(record, "creators")
        creator isa AbstractDict || continue
        name = _creatorname(creator)
        isnothing(name) && continue
        push!(authors, name)
        identifier = _creatoridentifier(creator)
        isnothing(identifier) || get!(identifiers, name, identifier)
    end
    title = _text(_maintitle(record))
    types = get(record, "types", nothing)
    kind = types isa AbstractDict ?
           _plain(get(types, "resourceTypeGeneral", nothing)) : nothing
    return SourceFacts(source = isempty(title) ? doi : title, doi = doi,
                       url = _text(get(record, "url", nothing)), title = title,
                       publisher = _text(_publishername(get(record, "publisher",
                                                            nothing))),
                       authors = authors, identifiers = identifiers,
                       releasedate = something(_datacitedate(record), ""),
                       version = _text(get(record, "version", nothing)),
                       description = something(_abstract(record), ""),
                       licence = something(_datacitelicence(_items(record,
                                                                   "rightsList")),
                                           ""),
                       identifies = _identifies(kind))
end

# The facts the Citation Style Language record of a DOI gives of its source.
function _cslfacts(record::AbstractDict, doi::AbstractString)
    authors = String[]
    identifiers = Dict{String, String}()
    for author in _items(record, "author")
        author isa AbstractDict || continue
        family = _plain(get(author, "family", nothing))
        name = isnothing(family) ? _plain(get(author, "literal", nothing)) :
               _plain(string(something(_plain(get(author, "given", nothing)),
                                       ""), " ", family))
        isnothing(name) && continue
        push!(authors, name)
        orcid = something(_plain(get(author, "ORCID", nothing)), "")
        occursin(r"^https?://\S+$", orcid) &&
            get!(identifiers, name, replace(orcid, r"^http:" => "https:"))
    end
    title = get(record, "title", nothing)
    title isa AbstractVector &&
        (title = isempty(title) ? nothing : first(title))
    title = _text(title)
    issued = get(record, "issued", nothing)
    parts = issued isa AbstractDict ? _items(issued, "date-parts") : Any[]
    date = isempty(parts) || !(first(parts) isa AbstractVector) ? "" :
           _dateasfaras(first(parts))
    return SourceFacts(source = isempty(title) ? doi : title, doi = doi,
                       url = _text(get(record, "URL", nothing)), title = title,
                       publisher = _text(get(record, "publisher", nothing)),
                       authors = authors, identifiers = identifiers,
                       releasedate = date,
                       version = _text(get(record, "version", nothing)),
                       description = something(_markupless(_plain(get(record,
                                                                      "abstract",
                                                                      nothing))),
                                               ""),
                       identifies = _cslidentifies(_plain(get(record, "type",
                                                              nothing))))
end

# A CSL date, its year, month and day as given, as `YYYY`, `YYYY-MM` or `YYYY-MM-DD`.
function _dateasfaras(parts::AbstractVector)
    numbers = [tryparse(Int, something(_plain(p), "")) for p in parts]
    (isempty(numbers) || isnothing(first(numbers))) && return ""
    text = lpad(first(numbers), 4, '0')
    for n in numbers[2:min(3, end)]
        isnothing(n) && break
        text *= "-" * lpad(n, 2, '0')
    end
    return text
end

# The types of metadata a DOI is asked for, in order, and what reads each.
const _METADATA_READERS = (("application/vnd.datacite.datacite+json",
                            _datacitefacts),
                           ("application/vnd.citationstyles.csl+json",
                            _cslfacts))
