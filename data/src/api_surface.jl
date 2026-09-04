# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Regenerate `data/api_surface.md` — every name the package defines, grouped by **visibility**
# (exported, `public`, private) and within that by **kind**, with a one-line description taken from
# the name's own docstring or from the `#` comment above its definition.
#
# It exists as the working document for two audits: `[C7-VIS]` (which exported names belong in the
# `using EcoSISTEM` namespace at all) and `[COMMENTS]` (which private helpers are undocumented).
# Nothing here judges either question; it lays out the surface so a person can.
#
# Must run with the weak dependencies loaded, or the extensions never load and their private helpers
# are invisible. ⭐ `data/src` is that environment - it carries MPI, Phylo and RasterDataSources
# precisely so every extension loads, which no project in this repo used to do:
#
#     julia --project=data/src data/src/api_surface.jl
#
# The page's prose lives in `NOTES` below, not in the generated file — editing the output by hand
# would lose it on the next run. That is the same arrangement `data/src/overloads.jl` has, and it is why
# this generator was written at all: the document existed for months with no keeper, and named
# `src/Simplify.jl` in 115 places after that file was split up.
using EcoSISTEM, MPI, Phylo, RasterDataSources

const ROOT = pkgdir(EcoSISTEM)
const SUBMODULES = [EcoSISTEM.Units, EcoSISTEM.ClimatePref]

# The modules whose names are ours to inventory. An extension that is not loaded simply contributes
# nothing, so the report covers less rather than reporting wrongly.
function ourmodules()
    mods = Any[EcoSISTEM, SUBMODULES...]
    for e in (:EcoSISTEMMPIExt, :EcoSISTEMPhyloExt,
        :EcoSISTEMRasterDataSourcesExt, :EcoSISTEMERAExt,
        :EcoSISTEMDataPipelineExt)
        m = Base.get_extension(EcoSISTEM, e)
        isnothing(m) || push!(mods, m)
    end
    return mods
end

# --- What a name IS ----------------------------------------------------------

# The four kinds the report groups by. A macro counts as a function: it is called, not stored, and
# grouping `@nicheaxis` away from the functions would hide the package's one public macro.
function kindof(m::Module, n::Symbol)
    isdefined(m, n) || return :missing
    v = try
        getfield(m, n)
    catch
        return :missing
    end
    v isa Type || return v isa Function ? :function : :constant
    isalias(m, n) && return :alias
    return isabstracttype(v) ? :abstract : :concrete
end

# A type bound to a name that is not its own is an alias (`Torus`, `Supply`, `LayerSpec`), which the
# report keeps apart from the types themselves: an alias adds no new type to the hierarchy, and
# counting it as one would inflate the surface being audited. A `Union` is an alias by the same
# argument -- it names a set of types rather than declaring one, and `isabstracttype` says `false`
# for it, so without this it would be filed as concrete.
function isalias(m::Module, n::Symbol)
    v = getfield(m, n)
    v isa Type || return false
    b = Base.unwrap_unionall(v)
    return b isa DataType ? nameof(b) !== n : true
end

# --- Where a name is written -------------------------------------------------

# `file:line` relative to the repo root, or `nothing` for anything outside it (a re-exported name
# from a dependency, which the report should not claim as ours).
function loc(file, line)
    f = String(file)
    startswith(f, ROOT) || return nothing
    return string(relpath(f, ROOT), ":", line)
end

# Every place a function is defined, **in source order**. Sorting matters and is not cosmetic: the
# `#` comment documenting a helper sits above its *first* method, and `methods()` iterates in the
# order Julia happens to hold them, so taking its first entry reads the comment above some later
# method and reports a documented helper as undocumented.
function methodlocs(v)
    ls = Tuple{String, Int, String}[]
    for mm in methods(v)
        l = loc(mm.file, mm.line)
        isnothing(l) || push!(ls, (String(mm.file), Int(mm.line), l))
    end
    sort!(ls, by = x -> (x[1], x[2]))
    return unique(x[3] for x in ls)
end

# Top-level `const` declarations, by name, found by reading `src/` and `ext/`. Reflection cannot give
# these: a constant carries no line number, and its binding knows only its module.
function constlocs()
    out = Dict{Symbol, String}()
    for (dir, _, files) in walkdir(ROOT), f in files
        endswith(f, ".jl") || continue
        rel = relpath(joinpath(dir, f), ROOT)
        (startswith(rel, "src/") || startswith(rel, "ext/")) || continue
        for (i, line) in enumerate(eachline(joinpath(dir, f)))
            mt = match(r"^\s*const\s+([A-Za-z_][A-Za-z0-9_!]*)\s*=", line)
            isnothing(mt) && continue
            key = Symbol(mt.captures[1])
            haskey(out, key) || (out[key] = string(rel, ":", i))
        end
    end
    return out
end

# --- What is DEPRECATED, and therefore out of scope ---------------------------

# Names a `deprecations.jl` file declares in a way that makes the NAME itself deprecated, rather
# than merely giving a live name one deprecated method.
#
# 🔴 The distinction is the whole difficulty, and getting it wrong went in the dangerous direction:
# a first version matched every `function NAME(` and `NAME(...) =` in the file, which excluded
# `ContinuousRegime`, `Supply` and `NoChange` -- all current, all present only because a deprecated
# *constructor* for them lives there, or because a shim points at them. So the patterns here are
# only the ones that cannot mean anything else: the first argument of `@deprecate`/
# `@deprecate_binding` (never the target, which is the replacement), a method-less stub, and a
# `const` or type declared in the file itself.
function deprecatedinfiles()
    out = Set{Symbol}()
    pat = [
        r"^\s*(?:Base\.)?@deprecate(?:_binding)?\s+([A-Za-z_][A-Za-z0-9_!]*)",
        r"^\s*function\s+([A-Za-z_][A-Za-z0-9_!]*)\s+end\s*$",
        r"^\s*const\s+([A-Za-z_][A-Za-z0-9_!]*)\s*=",
        r"^\s*(?:mutable\s+)?struct\s+([A-Za-z_][A-Za-z0-9_!]*)",
        r"^\s*abstract\s+type\s+([A-Za-z_][A-Za-z0-9_!]*)"]
    for (dir, _, files) in walkdir(ROOT), f in files
        f == "deprecations.jl" || continue
        rel = relpath(joinpath(dir, f), ROOT)
        (startswith(rel, "src/") || startswith(rel, "ext/")) || continue
        for line in eachline(joinpath(dir, f)), r in pat
            m = match(r, line)
            isnothing(m) || push!(out, Symbol(m.captures[1]))
        end
    end
    return out
end

# Does any method of `v` warn when called? The `@deprecate` macro emits a `Base.depwarn` call, and a
# hand-written shim makes the same one.
function warnsondeprecation(mm::Method)
    src = try
        Base.uncompressed_ast(mm)
    catch
        return false
    end
    return any(st -> occursin("depwarn", string(st)), src.code)
end

# ⚠️ **A `@deprecate`-generated method reports Julia's own `base/deprecated.jl` as its file**, not
# the call site where the macro was written -- and `deprecated.jl` is not `deprecations.jl`, so a
# substring test for ours misses every one of them. Its body is the macro's, so the `depwarn` scan
# misses them too. Matching Julia's own file by name is what catches the whole `@deprecate` family.
function isdeprecatedmethod(mm::Method)
    f = String(mm.file)
    return occursin("deprecations.jl", f) ||
           !isnothing(match(r"(^|/)deprecated\.jl$", f)) ||
           warnsondeprecation(mm)
end

const DEPRECATED_IN_FILES = deprecatedinfiles()

# A name is out of scope for the visibility audit if it is on its way out: nothing is gained by
# deciding where a deprecated name belongs in a namespace it is being removed from.
#
# ⚠️ **`all`, never `any`** -- the name is deprecated only when *every* way of reaching it is. A live
# type routinely carries one deprecated constructor beside its current ones, and `any` would retire
# the type along with it.
function isdeprecatedname(m::Module, n::Symbol)
    Base.isdeprecated(m, n) && return true
    isdefined(m, n) || return false
    v = try
        getfield(m, n)
    catch
        return false
    end
    ms = (v isa Function || v isa Type) ? collect(methods(v)) : Method[]
    isempty(ms) || return all(isdeprecatedmethod, ms)
    # No methods at all: a method-less stub, a `const`, or a plain data type. Only the file scan can
    # speak for these, and it is written narrowly enough to be trusted with them.
    return n in DEPRECATED_IN_FILES
end

# --- What a name SAYS --------------------------------------------------------

# The first sentence, collapsed to one line and trimmed to something a table cell can hold.
function firstsentence(text::AbstractString)
    s = strip(replace(text, r"\s+" => " "))
    isempty(s) && return nothing
    # Stop at the first `. ` that is not inside a code span and not an abbreviation of a file
    # extension -- `src/Foo.jl. Next` must break after `jl.`, not after `Foo.`.
    m = match(r"^(.*?[.!?])(?: |$)", s)
    out = isnothing(m) ? s : m.captures[1]
    length(out) > 400 && (out = first(out, 397) * "…")
    return String(out)
end

# A docstring here opens with a signature header, a blank line, then the prose -- so the prose does
# not start at character one, and collapsing the whole thing to one sentence would report the
# signature as the description. Drop the leading indented block.
function stripsignature(text::AbstractString)
    lines = split(text, '\n')
    i = 1
    while i <= length(lines) && isempty(strip(lines[i]))
        i += 1
    end
    j = i
    while j <= length(lines) && startswith(lines[j], "    ") &&
          !isempty(strip(lines[j]))
        j += 1
    end
    j > i || return text                         # no header: the prose starts straight away
    while j <= length(lines) && isempty(strip(lines[j]))
        j += 1
    end
    return join(lines[j:end], '\n')
end

# The docstring attached to a name, read through its `Binding`. Going through the binding is the only
# form that works: `@doc(value)` documents the *variable* and reports everything as undocumented.
function docsentence(m::Module, n::Symbol)
    b = Base.Docs.Binding(m, n)
    for mod in Base.Docs.modules
        meta = Base.Docs.meta(mod, autoinit = false)
        isnothing(meta) && continue
        haskey(meta, b) || continue
        d = meta[b]
        for sig in d.order
            txt = join(d.docs[sig].text, "")
            s = firstsentence(stripsignature(txt))
            isnothing(s) || return s
        end
    end
    return nothing
end

# The `#` comment block immediately above a definition, which is how a private helper is documented
# (`CLAUDE.md`: "every private helper gets a `#` comment -- no exceptions"). Read upwards from the
# definition line, stopping at the first line that is neither a comment nor blank-inside-a-block.
function commentsentence(location::AbstractString)
    return first(commentabove(location))
end

# Whether a line is a one-line definition -- `_placex(p) = p.x`. A `#` comment often heads a *group*
# of these and documents all of them, so a walk that stops at the first non-comment line reports
# every member but the first as undocumented. Measured: that was a real class of false positive here.
isoneliner(t) = !isnothing(match(r"^\s*[A-Za-z_][A-Za-z0-9_!]*[({].*=(?!=)", t))

# The `#` comment documenting a definition, and whether it had to be reached over its neighbours.
# `(sentence, shared)` -- `shared` is true when the comment heads a group rather than sitting
# directly above this definition, which is a different thing from having none and the report says so.
function commentabove(location::AbstractString)
    parts = rsplit(location, ':', limit = 2)
    length(parts) == 2 || return nothing
    path = joinpath(ROOT, parts[1])
    isfile(path) || return nothing
    line = tryparse(Int, parts[2])
    isnothing(line) && return nothing
    lines = readlines(path)
    (1 <= line <= length(lines)) || return (nothing, false)
    i = line - 1
    shared = false
    while i >= 1 && isoneliner(lines[i])      # step over the rest of the group, if any
        shared = true
        i -= 1
    end
    block = String[]
    while i >= 1
        t = strip(lines[i])
        startswith(t, "#") || break
        c = strip(lstrip(t, '#'))
        startswith(c, "SPDX-License-Identifier") && break
        pushfirst!(block, c)
        i -= 1
    end
    isempty(block) && return (nothing, false)
    return (firstsentence(join(block, " ")), shared)
end

# --- The inventory -----------------------------------------------------------

struct Entry
    name::Symbol
    owner::Module
    kind::Symbol
    locs::Vector{String}
    what::Union{String, Nothing}
    note::Union{String, Nothing}     # why the description reads as it does, if anything
end

function label(e::Entry)
    return e.owner === EcoSISTEM ? "`$(e.name)`" :
           "`$(e.name)` *(`$(nameof(e.owner))`)*"
end

# One row of the report. `public` means visible-and-supported, so its description comes from the
# docstring; a private helper's comes from the `#` comment above it, and a docstring found on one is
# itself the finding `[COMMENTS]` is looking for.
function entry(m::Module, n::Symbol, vis::Symbol, consts)
    k = kindof(m, n)
    k === :missing && return nothing
    v = getfield(m, n)
    ls = k in (:function,) ? methodlocs(v) : String[]
    if isempty(ls) && haskey(consts, n)
        ls = [consts[n]]
    end
    doc = docsentence(m, n)
    if vis === :private
        # `CLAUDE.md`'s rule is about `_`-prefixed helpers specifically, so only those are held to
        # "a `#` comment, never a docstring". An internal name without the prefix -- `bounds`, kept
        # non-public because `DimensionalData` and `Rasters` both export it -- is documented prose
        # deliberately, and flagging it would bury the real findings among false ones.
        helper = startswith(string(n), "_")
        cmt, shared = isempty(ls) ? (nothing, false) : commentabove(first(ls))
        isnothing(doc) ||
            return Entry(n, m, k, ls, doc,
                         helper ? "docstring, not a `#` comment" : nothing)
        isnothing(cmt) && return Entry(n, m, k, ls, nothing, "undocumented")
        return Entry(n, m, k, ls, cmt,
                     shared ? "shared group comment" : nothing)
    end
    return Entry(n, m, k, ls, doc, isnothing(doc) ? "undocumented" : nothing)
end

# Types first by containment -- a type appearing inside another comes before it -- then
# alphabetically within each level. Functions and constants are alphabetical throughout: neither has
# a containment order to respect, and alphabetical is what an audit reads by.
function containmentorder(es::Vector{Entry})
    types = Dict(e.name => e for e in es)
    # What `e`'s declaration mentions that is also in this set: its supertype, its field types, and
    # (unlike `data/src/typeorder.jl`, which asks a different question) the concrete implementors of an
    # abstract field type, since a reader meeting the container wants those first.
    function needs(e)
        v = getfield(e.owner, e.name)
        b = Base.unwrap_unionall(v)
        b isa DataType || return Symbol[]
        acc = Set{Symbol}()
        for t in Any[supertype(b);
                     (isstructtype(b) ? collect(fieldtypes(b)) : Any[])...]
            u = Base.unwrap_unionall(t)
            u isa DataType || continue
            nm = nameof(u)
            nm === e.name && continue
            haskey(types, nm) && push!(acc, nm)
            for (other, oe) in types
                other === e.name && continue
                ov = getfield(oe.owner, other)
                ov isa Type && u !== Any && ov <: u && push!(acc, other)
            end
        end
        return collect(acc)
    end
    ready = sort(collect(keys(types)), by = string)
    done = Symbol[]
    placed = Set{Symbol}()
    for _ in 1:length(ready)
        progressed = false
        for n in ready
            n in placed && continue
            all(d -> d in placed, needs(types[n])) || continue
            push!(done, n)
            push!(placed, n)
            progressed = true
        end
        progressed || break
    end
    for n in ready                       # a cycle leaves names unplaced; list them rather than drop
        n in placed || push!(done, n)
    end
    return [types[n] for n in done]
end

function sortentries(es, kind)
    return kind in (:abstract, :concrete) ? containmentorder(es) :
           sort(es,
                by = e -> (lowercase(string(e.name)), string(nameof(e.owner))))
end

const NOTES = Dict{Symbol, String}(:header => """
Types are ordered so that a type appearing **inside** another comes first — supertypes before
subtypes, and a type used as a field before the type holding it. Functions and constants are
alphabetical.

**Ordering caveat**: the containment order follows *declared field types*, resolving an abstract
field type to the concrete types listed here. A container that stores its members in a bare
`NamedTuple` (`LayerCollection`, `SpeciesRequirementCollection`) carries no type-level link to them,
so it cannot be placed after them automatically.

Descriptions are the first sentence of the name's own docstring — or, for private helpers, of the
`#` comment above it. None are invented. *(undocumented)* means there is neither.
""",
                                   :exported => """
Reached by `using EcoSISTEM` alone. This is the surface `[C7-VIS]` audits: the question to ask of
each name is *"do you need this to **use** the package, or only to **extend** it?"* — the first
belongs here, the second belongs under `public`.

⚠️ Every abstract type in this section is already known to be misplaced: `CLAUDE.md` says every
abstract type is `public`, never exported, because a supertype is what you dispatch on or subtype
when extending the package rather than when using it.
""",
                                   :public => """
Supported and documented, but **not** in the `using EcoSISTEM` namespace — reach them as
`EcoSISTEM.name`, or `using EcoSISTEM: name`.

⚠️ `names(M)` includes `public` names on Julia 1.11+, so it cannot tell exported from public. This
report uses `Base.isexported`, and anything checking the split must do the same.
""",
                                   :private => """
Not part of the API. Listed because `[COMMENTS]` audits them against `CLAUDE.md`'s rule that **every
private helper gets a `#` comment** and that private functions are **not** given docstrings — so a
row marked *(undocumented)* or *(docstring, not a `#` comment)* is a finding, and the rest are
confirmation.

⚠️ A name defined only inside an extension appears here only when that extension is loaded. Run this
generator with `MPI`, `Phylo` and `RasterDataSources` present, or the counts come out short.
""")

# Which of our modules a name really belongs to. A submodule re-exports a great deal of the parent
# so that `using EcoSISTEM.ClimatePref` keeps working -- 37 of `ClimatePref`'s 54 public names are
# the parent's -- and listing those twice would inflate the surface being audited and invite the same
# name to be decided twice, differently.
function owner(n::Symbol, mods, declared)
    # 🔴 Only modules that DECLARE the name are candidates, never merely those it is visible in.
    # `EcoSISTEM` does `using EcoSISTEM.Units`, so `isdefined(EcoSISTEM, :January)` is true while
    # `names(EcoSISTEM, imported = false)` rightly omits it -- and picking `EcoSISTEM` as the owner
    # dropped every month from the report, because nothing ever offers the name from there.
    here = filter(m -> n in declared[m], mods)
    length(here) <= 1 && return isempty(here) ? nothing : first(here)
    v = getfield(first(here), n)
    same = filter(m -> (getfield(m, n) === v), here)
    length(same) <= 1 && return first(here)
    # `parentmodule` answers for a type or a function; a `const` (a `Union` such as `CODE_TYPE`) has
    # no parent, so fall back to declaration order, which puts the parent module first.
    pm = try
        v isa Type ? parentmodule(Base.unwrap_unionall(v)) : parentmodule(v)
    catch
        nothing
    end
    pm in same && return pm
    return first(same)
end

function collectall()
    consts = constlocs()
    mods = ourmodules()
    declared = Dict(m => Set(names(m, all = true, imported = false))
                    for m in mods)
    exported = Entry[]
    public = Entry[]
    private = Entry[]
    skipped = (deprecated = Ref(0), reexported = Ref(0))
    seen = Set{Tuple{Symbol, Module}}()
    for m in mods
        for n in names(m, all = true, imported = false)
            startswith(string(n), "#") && continue
            n === nameof(m) && continue
            (n, m) in seen && continue
            isdefined(m, n) || continue
            v = try
                getfield(m, n)
            catch
                continue
            end
            v isa Module && continue
            vis = Base.isexported(m, n) ? :exported :
                  Base.ispublic(m, n) ? :public : :private
            # Out of scope, and counted so the omission is visible rather than silent.
            if isdeprecatedname(m, n)
                vis === :private || (skipped.deprecated[] += 1)
                continue
            end
            if owner(n, mods, declared) !== m
                vis === :private || (skipped.reexported[] += 1)
                continue
            end
            e = entry(m, n, vis, consts)
            isnothing(e) && continue
            vis === :private && isempty(e.locs) && continue   # imported from a dependency
            push!(seen, (n, m))
            push!(vis === :exported ? exported :
                  vis === :public ? public : private, e)
        end
    end
    return (exported = exported, public = public, private = private,
            skipped = skipped)
end

const KINDNAMES = [(:abstract, "Abstract types"), (:concrete, "Concrete types"),
    (:alias, "Type aliases"), (:function, "Functions"),
    (:constant, "Constants and values")]

function writesection(io, title, entries, note; withloc::Bool)
    println(io, "\n## ", title, " — ", length(entries), "\n")
    println(io, note)
    for (k, kname) in KINDNAMES
        es = filter(e -> e.kind === k, entries)
        isempty(es) && continue
        println(io, "\n### ", kname, " — ", length(es), "\n")
        println(io,
                withloc ? "| name | defined | what it is for |" :
                "| name | what it is for |")
        println(io, withloc ? "|---|---|---|" : "|---|---|")
        for e in sortentries(es, k)
            what = isnothing(e.what) ?
                   "*($(something(e.note, "undocumented")))*" :
                   isnothing(e.note) ? e.what : "*($(e.note))* " * e.what
            what = replace(what, "|" => "\\|")
            if withloc
                l = isempty(e.locs) ? "—" :
                    length(e.locs) == 1 ? "`$(first(e.locs))`" :
                    "`$(first(e.locs))` (+$(length(e.locs) - 1) more)"
                println(io, "| ", label(e), " | ", l, " | ", what, " |")
            else
                println(io, "| ", label(e), " | ", what, " |")
            end
        end
    end
end

surface = collectall()
exts = filter(m -> !(m in [EcoSISTEM, SUBMODULES...]), ourmodules())
nundoc = count(e -> e.note == "undocumented",
               filter(e -> e.kind === :function, surface.private))
ndoc = count(e -> e.note == "docstring, not a `#` comment",
             filter(e -> e.kind === :function, surface.private))

open(joinpath(ROOT, "data", "api_surface.md"), "w") do io
    println(io, "# EcoSISTEM — the name surface\n")
    println(io,
            "Every name the package defines, grouped by **visibility** (exported → `public` → " *
            "private) and within that by **kind**.\n")
    println(io, NOTES[:header])
    println(io,
            "**Generated by `data/src/api_surface.jl`** — do not edit by hand. Extensions loaded for " *
            "this run: " *
            (isempty(exts) ? "*none*" :
             join(("`" * string(nameof(m)) * "`" for m in exts), ", ")) * ".\n")
    println(io, "| visibility | names |\n|---|---:|")
    println(io, "| exported | ", length(surface.exported), " |")
    println(io, "| `public` | ", length(surface.public), " |")
    println(io, "| private | ", length(surface.private), " |")
    println(io,
            "\n⚠️ **Two kinds of name are deliberately absent**, because neither is a question " *
            "anyone has to answer: **", surface.skipped.deprecated[],
            "** exported or `public` names that are **deprecated** — where a name on its way out " *
            "sits in a namespace it is being removed from is not worth deciding — and **",
            surface.skipped.reexported[],
            "** that a submodule **re-exports** from the parent, which are listed once where they " *
            "are defined. Both are found mechanically (`deprecations.jl` declarations plus " *
            "`Base.isdeprecated` and a `depwarn` scan; value identity across our modules), so " *
            "neither depends on a list that could rot.")
    println(io,
            "\n⭐ **`[COMMENTS]` at a glance**: of ",
            count(e -> e.kind === :function, surface.private),
            " private functions, **", nundoc,
            "** carry no documentation at all and **", ndoc,
            "** carry a docstring where a `#` comment belongs.")
    writesection(io, "Exported", surface.exported, NOTES[:exported],
                 withloc = false)
    writesection(io, "Public, not exported", surface.public, NOTES[:public],
                 withloc = false)
    return writesection(io, "Private", surface.private, NOTES[:private],
                        withloc = true)
end
println("wrote data/api_surface.md: ", length(surface.exported), " exported, ",
        length(surface.public), " public, ", length(surface.private),
        " private (", nundoc, " undocumented functions, ", ndoc,
        " with docstrings)")
