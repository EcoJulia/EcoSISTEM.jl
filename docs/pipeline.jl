# SPDX-License-Identifier: LGPL-3.0-or-later
#
# What the data pipeline page, `docs/src/pipeline.md`, shows of the example's runs besides their
# output: a registered image, and the provenance graph of a data product, drawn from the registry's
# own PROV record with every node named, back through the runs to the items the providers' data is
# cited by. The page includes this when it runs the example.

using DataPipeline
using Downloads
using JSON

# A PNG file, which Documenter shows as the image it is.
struct PNGFile
    path::String
end
Base.show(io::IO, ::MIME"image/png", file::PNGFile) = write(io, read(file.path))

# An SVG drawing, which Documenter shows as the image it is.
struct SVGDrawing
    svg::String
end
function Base.show(io::IO, ::MIME"image/svg+xml", drawing::SVGDrawing)
    return print(io, drawing.svg)
end

# The image a data product of a trace is, from the store where the registry says it is.
function registeredimage(product::DataPipeline.TraceProduct)
    product.format == "png" ||
        throw(ArgumentError("$(product.name) is a $(product.format) file, not a PNG"))
    return PNGFile(chopprefix(product.location, "file://"))
end

# The registry's PROV record of the data product at `url`, as far back as `depth` code runs.
function provrecord(registry::DataPipeline.RegistryEndpoint,
                    url::AbstractString,
                    depth::Integer)
    id = match(r"/data_product/(\d+)/", url)[1]
    report = string(registry.url, "prov-report/", id, "/?format=json&depth=",
                    depth)
    headers = isnothing(registry.token) ? Pair{String, String}[] :
              ["Authorization" => "token $(registry.token)"]
    return JSON.parse(String(take!(Downloads.download(report, IOBuffer(),
                                                      headers = headers))))
end

# The relations of one kind in a PROV record, which keys each by a blank node.
relations(record, kind) = collect(values(get(record, kind, Dict())))

# A PROV value, which is either the value itself or a typed value holding it.
provvalue(x) = x isa AbstractDict ? x["\$"] : x

# Text broken at spaces into lines of at most `width` characters.
function wrapped(text::AbstractString, width::Integer = 24)
    lines = String[]
    for word in split(text)
        if !isempty(lines) && length(lines[end]) + 1 + length(word) <= width
            lines[end] *= " " * word
        else
            push!(lines, word)
        end
    end
    return join(lines, "\n")
end

# Text in GraphViz's quoted form.
function dotquoted(text)
    return string('"',
                  replace(text, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n"),
                  '"')
end

# An identifier as a label shows it: a DOI as `doi:...`, else the link without its scheme.
function shortlink(link::AbstractString)
    startswith(link, r"https?://(dx\.)?doi\.org/") &&
        return replace(link, r"^https?://(dx\.)?doi\.org/" => "doi:")
    return rstrip(replace(link, r"^[a-z]+://" => ""), '/')
end

# "Ann", "Ann and Bob", "Ann, Bob and Cy".
function listed(names)
    length(names) <= 2 && return join(names, " and ")
    return string(join(names[1:(end - 1)], ", "), " and ", names[end])
end

# The provenance of the data product `trace` starts from, as the registry it was traced in records
# it `depth` code runs back, drawn by GraphViz from the providers' items on the left to the data
# product on the right.
function provenancegraph(trace::DataPipeline.Trace; depth::Integer)
    root = only(trace.roots)
    record = provrecord(trace.registry, root, depth)
    prefix = record["prefix"]["lreg"]
    return drawn(provenancedot(record,
                               "lreg:" * rstrip(chopprefix(root, prefix), '/')))
end

# A PROV record as a GraphViz graph: each item a provider's data is cited by (a folded corner) with
# its title, version, year, authors and identifier; each data product (an ellipse) with its name,
# namespace and version, `root` outlined; each code run (a box) with its description, when and by
# whom it was started, and the code it ran, with that code's authors. An arrow from a cited item is
# a data product derived from it, dashed where the data product's file is the item itself; the
# others join each run to what it read and what it wrote. The configs and scripts the registry also
# holds for each run are left out.
function provenancedot(record::AbstractDict, root::AbstractString)
    entities = get(record, "entity", Dict())
    activities = get(record, "activity", Dict())
    agents = get(record, "agent", Dict())
    person(id) = agents[id]["foaf:name"]
    attributed = Dict{String, Vector{String}}()
    for relation in relations(record, "wasAttributedTo")
        push!(get!(attributed, relation["prov:entity"], String[]),
              relation["prov:agent"])
    end
    started = Dict(relation["prov:activity"] => relation["prov:trigger"]
                   for relation in relations(record, "wasStartedBy"))
    runners = Set(values(started))
    used = relations(record, "used")
    software = Dict(relation["prov:activity"] => relation["prov:entity"]
                    for relation in used
                    if provvalue(relation["prov:role"]) == "fair:software")
    generated = relations(record, "wasGeneratedBy")
    made = Set(relation["prov:entity"] for relation in generated)
    # Each cited item's data products: derived from it, or specialising it, its file being the item.
    derived = [(relation["prov:usedEntity"], relation["prov:generatedEntity"])
               for relation in relations(record, "wasDerivedFrom")
               if occursin("external_object", relation["prov:usedEntity"])]
    specialised = [(relation["prov:generalEntity"],
                    relation["prov:specificEntity"])
                   for relation in relations(record, "specializationOf")]
    # A cited item's authors are those of the data products made from it, other than whoever
    # started a run.
    authors(item) = unique(person(agent)
                           for (cited, product) in vcat(derived, specialised)
                           if cited == item
                           for agent in get(attributed, product, String[])
                           if agent ∉ runners)
    lines = ["digraph provenance {",
        "  rankdir=LR; bgcolor=white; nodesep=0.12; ranksep=0.25;",
        "  node [fontname=Helvetica, fontsize=12, style=filled,",
        "        penwidth=0.8, margin=\"0.06,0.03\"];",
        "  edge [color=\"#555555\", arrowsize=0.6];"]
    function node(id, style, label)
        return push!(lines,
                     string("  ", dotquoted(id), " [", style, ", label=",
                            dotquoted(label), "];"))
    end
    function edge(from, to, style = "")
        return push!(lines,
                     string("  ", dotquoted(from), " -> ", dotquoted(to),
                            isempty(style) ? "" : " [$style]", ";"))
    end
    for (id, entity) in entities
        if occursin("external_object", id)
            version = get(entity, "dcat:hasVersion", "")
            issued = provvalue(get(entity, "dcterms:issued", ""))
            when = join(filter(!isempty,
                               [isempty(version) ? "" : "v" * version,
                                   first(issued, 4)]), ", ")
            identifier = shortlink(get(entity, "dcterms:identifier", ""))
            facts = [wrapped(entity["dcterms:title"]), when,
                wrapped(listed(authors(id))), identifier]
            node(id, "shape=note, fillcolor=\"#e8e8e8\"",
                 join(filter(!isempty, facts), "\n"))
        elseif occursin("data_product", id)
            # A long name's last segment on a line of its own.
            name = entity["dcterms:title"]
            count('/', name) > 1 &&
                (name = replace(name, r"/([^/]*)$" => s"/\n\1"))
            outline = id == root ? ", penwidth=2" : ""
            node(id, "shape=ellipse, fillcolor=\"#fffc87\"" * outline,
                 string(name, "\n", entity["fair:namespace"], ", v",
                        entity["dcat:hasVersion"]))
        end
    end
    for (id, activity) in activities
        facts = [wrapped(activity["dcterms:description"])]
        haskey(started, id) && push!(facts, "run by " * person(started[id]))
        push!(facts, "on " * first(activity["prov:startTime"], 10))
        if haskey(software, id)
            code = entities[software[id]]
            push!(facts,
                  replace(code["prov:atLocation"], r"^[a-z]+://[^/]+/" => "",
                          r"\.git$" => ""))
            commit = get(code, "fair:commit", "")
            isempty(commit) ||
                push!(facts,
                      string("at ", first(commit, 8),
                             endswith(commit, "-dirty") ? "-dirty" : ""))
            coders = person.(get(attributed, software[id], String[]))
            isempty(coders) || push!(facts, wrapped("by " * listed(coders)))
        end
        node(id, "shape=box, fillcolor=\"#9fb1fc\"", join(facts, "\n"))
    end
    for (cited, product) in derived
        edge(cited, product)
    end
    for (cited, product) in specialised
        edge(cited, product, "style=dashed")
    end
    # A cited item whose data product no run made, a deposit, is stacked with it rather than given
    # a column of its own.
    for (cited, product) in vcat(derived, specialised)
        product in made ||
            push!(lines,
                  string("  { rank=same; ", dotquoted(cited), "; ",
                         dotquoted(product), "; }"))
    end
    for relation in used
        provvalue(relation["prov:role"]) == "fair:input_data" &&
            edge(relation["prov:entity"], relation["prov:activity"])
    end
    for relation in generated
        edge(relation["prov:activity"], relation["prov:entity"])
    end
    push!(lines, "}")
    return join(lines, "\n")
end

# A GraphViz graph drawn as SVG by `dot`.
function drawn(graph::AbstractString)
    dot = Sys.which("dot")
    isnothing(dot) &&
        error("drawing the provenance graph needs GraphViz's `dot` on the path")
    return SVGDrawing(read(pipeline(`$dot -Tsvg`, stdin = IOBuffer(graph)),
                           String))
end
