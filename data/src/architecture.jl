# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Keep `data/architecture.md` honest about the code's type hierarchies.
#
# Audit - what the doc gets wrong, and what it has never covered:
#
#     julia --project=data/src data/src/architecture.jl
#
# Regenerate the inheritance edges of every mermaid block, in place:
#
#     julia --project=data/src data/src/architecture.jl --fix
#
# 🔴 **`--fix` repairs the DIAGRAMS ONLY. The prose is the valuable half and it cannot see it.**
# ✅ Measured the day this was written: `data/architecture.md` declared `class Reference~A~` for a type
# deleted from `ClimatePref` - *and* said in prose that "`ERA`/`CERA`/`CRUTS`/`Reference` remain for
# their data sources". `--fix` removes the first and leaves the second silently wrong. Always read the
# text around any block it changes, and treat its output as a diff to review rather than an answer.
#
# ⚠️ Run this under the **test** environment (or with the weak deps loaded), or the extension types are
# invisible and it reports false gaps - `MPIEcosystem`'s concrete subtypes live in `EcoSISTEMMPIExt`.

module ArchitectureAudit

using EcoSISTEM

export audit, architecture_report, regenerate

const DOC = joinpath(pkgdir(EcoSISTEM), "data", "architecture.md")

# The modules whose own declarations the doc is expected to cover. ⭐ Extensions are included when
# loaded and skipped when not, so a run without the weak deps under-reports rather than lying.
function ourmodules()
    mods = Module[EcoSISTEM, EcoSISTEM.ClimatePref, EcoSISTEM.Units]
    for ext in (:EcoSISTEMMPIExt, :EcoSISTEMPhyloExt, :EcoSISTEMRasterDataSourcesExt,
        :EcoSISTEMDataPipelineExt, :EcoSISTEMERAExt)
        m = Base.get_extension(EcoSISTEM, ext)
        isnothing(m) || push!(mods, m)
    end
    return mods
end

# Where a name may legitimately resolve without being ours - the external supertypes we subtype from.
# ⭐ Resolving against real modules rather than keeping a hardcoded whitelist: a name is only "stale"
# if *nothing* reachable defines it, so an upstream rename shows up as staleness exactly as it should.
# 🔴🔴 `Base.identify_package` only sees DIRECT dependencies of the active project, so every module
# named below must be one, or its types silently drop out of `othermodules()` and the edges they
# carry are reported as stale. ✅ Measured: running under a project where Diversity and EcoBase were
# merely transitive produced 7 false differences - and `--fix` would have applied them, deleting real
# inheritance edges from a document the audit exists to protect. Run this under `data/src`, whose
# Project.toml lists all five, or under the test environment.
function othermodules()
    mods = Module[Base, Core]
    for name in (:Diversity, :EcoBase, :Distributions, :Unitful, :DimensionalData)
        id = Base.identify_package(String(name))
        isnothing(id) && continue
        m = get(Base.loaded_modules, id, nothing)
        isnothing(m) || push!(mods, m)
    end
    return mods
end

# Is `n` this module's own type *declaration*? Excludes aliases (whose value carries another name),
# imported names, and the gensym'd types Julia makes for closures.
function isowntype(m::Module, n::Symbol)
    startswith(String(n), "#") && return false
    isdefined(m, n) || return false
    v = getfield(m, n)
    v isa Type || return false
    u = Base.unwrap_unionall(v)
    u isa DataType || return false
    return u.name.name === n && parentmodule(v) === m
end

# Is `n` an alias - a binding to a type declared under a different name? ⭐ `const AbstractRegime =
# AbstractLayer{Condition}` and `const Torus = EdgeTopology{Periodic, Periodic}` both answer true.
# The doc names these deliberately, and no `supertype` walk can produce them, so they are collected
# rather than derived and `--fix` leaves their edges alone.
function isaliastype(m::Module, n::Symbol)
    startswith(String(n), "#") && return false
    isdefined(m, n) || return false
    v = getfield(m, n)
    v isa Type || return false
    u = Base.unwrap_unionall(v)
    return u isa DataType && u.name.name !== n
end

"""
    architecture_report()

Gather what the code says about its own type hierarchies, and what `data/architecture.md` says about
them. Returns a named tuple: `owned` (our type declarations), `aliases`, `children` (parent name ->
child names), `roots` (hierarchy roots with at least one child), `documented` (names appearing in a
mermaid block), `stale` (documented but resolvable nowhere) and `undocumented` (a root with children
that no block names).
"""
function architecture_report()
    owned = Dict{Symbol, Type}()
    aliases = Set{Symbol}()
    for m in ourmodules(), n in names(m, all = true, imported = false)
        isowntype(m, n) && (owned[n] = getfield(m, n))
        isaliastype(m, n) && push!(aliases, n)
    end

    # parent -> children in ONE pass over the owned types. ⚠️ Deliberately not `subtypes()`, which
    # rescans every loaded module per call (~0.6 s) and would make this a minute-long job.
    children = Dict{Symbol, Vector{Symbol}}()
    for (n, T) in owned
        S = supertype(Base.unwrap_unionall(T))
        S === Any && continue
        push!(get!(children, nameof(S), Symbol[]), n)
    end
    foreach(sort!, values(children))

    # A root is an abstract type of ours whose own supertype is not one of ours - so the hierarchy
    # starts here as far as this package is concerned, whether or not it continues upstream.
    roots = sort([n
                  for (n, T) in owned
                  if isabstracttype(Base.unwrap_unionall(T)) &&
                         !haskey(owned,
                                 nameof(supertype(Base.unwrap_unionall(T)))) &&
                         !isempty(get(children, n, Symbol[]))])

    documented = docnames(read(DOC, String))
    resolvable(n) = haskey(owned, n) || n in aliases ||
                    any(m -> isdefined(m, n), othermodules())
    stale = sort([n for n in documented if !resolvable(n)])
    undocumented = sort([r for r in roots if !(r in documented)])
    return (owned = owned, aliases = aliases, children = children,
            roots = roots,
            documented = documented, stale = stale, undocumented = undocumented)
end

# Every type name a mermaid block mentions - from `class X~P,Q~` declarations and from both sides of
# an `A <|-- B` or `A *-- B` edge. ⚠️ Only inside ```mermaid fences: prose names types in backticks
# constantly, and treating those as declarations would make the audit meaningless.
function docnames(text::AbstractString)
    found = Set{Symbol}()
    for block in mermaidblocks(text), line in split(block.body, '\n')
        for mm in eachmatch(r"\bclass\s+([A-Z][A-Za-z0-9_]*)", line)
            push!(found, Symbol(mm.captures[1]))
        end
        for mm in eachmatch(r"([A-Z][A-Za-z0-9_]*)\s*(?:\"[^\"]*\")?\s*(?:<\|--|\*--)\s*(?:\"[^\"]*\")?\s*([A-Z][A-Za-z0-9_]*)",
                            line)
            push!(found, Symbol(mm.captures[1]))
            push!(found, Symbol(mm.captures[2]))
        end
    end
    return found
end

# Each ```mermaid fence, as (body, first byte of the body, last byte of the body).
function mermaidblocks(text::AbstractString)
    blocks = NamedTuple{(:body, :start, :stop), Tuple{String, Int, Int}}[]
    for mm in eachmatch(r"```mermaid\n(.*?)```"s, text)
        push!(blocks,
              (body = mm.captures[1], start = mm.offsets[1],
               stop = mm.offsets[1] + ncodeunits(mm.captures[1]) - 1))
    end
    return blocks
end

"""
    audit(; io = stdout)

Print what `data/architecture.md` gets wrong and what it has never covered, and return `true` when
there is nothing to report. ⚠️ A clean audit means the *diagrams* agree with the code; it says nothing
about whether the prose around them is still true.
"""
function audit(; io = stdout)
    r = architecture_report()
    println(io, "data/architecture.md vs the code")
    println(io, "  type declarations found : ", length(r.owned))
    println(io, "  aliases found           : ", length(r.aliases))
    println(io, "  hierarchy roots         : ", length(r.roots))
    println(io, "  names in mermaid blocks : ", length(r.documented))
    if !isempty(r.stale)
        println(io,
                "\n🔴 STALE - named in the doc, defined nowhere (deleted or renamed):")
        for n in r.stale
            println(io, "     ", n)
        end
    end
    if !isempty(r.undocumented)
        println(io,
                "\n⚠️ UNDOCUMENTED - a hierarchy with children that no diagram names:")
        for n in r.undocumented
            println(io, "     ", rpad(n, 30), "-> ", join(r.children[n], ", "))
        end
    end
    isempty(r.stale) && isempty(r.undocumented) &&
        println(io,
                "\n✅ every diagram name resolves, and every hierarchy is named somewhere.")
    return isempty(r.stale) && isempty(r.undocumented)
end

"""
    regenerate(; write = false, io = stdout)

Repair the **inheritance edges** of every mermaid block in `data/architecture.md`, writing the file
only when `write = true`. Returns the number of lines added plus removed.

🔴 **Deliberately minimal, not a re-render.** It only *removes* what is provably wrong and *adds*
what is provably missing, leaving every other byte alone - so the hand-curated ordering, alignment,
`*--` composition edges and alias edges survive, and the diff is small enough to read. A full
re-render would be easier to write and impossible to review.

What it changes, per block, within the set of types that block already names:

  - removes a `class` line for a type that is defined nowhere;
  - removes an `A <|-- B` edge where either end is defined nowhere, or where `B` is one of our own
    declarations whose real supertype is a *different* type the block also names;
  - adds an `A <|-- B` edge that follows from `supertype` and is missing.

⚠️ It never touches an edge whose child is an **alias** (`AbstractTolerance`, `Torus`, ...): those are
written by hand because no `supertype` walk can produce them.
"""
function regenerate(; write::Bool = false, io = stdout)
    r = architecture_report()
    text = read(DOC, String)
    isstale(n) = n in r.stale
    # The real parent of one of our declarations, or `nothing` for an alias or a foreign name.
    function realparent(n)
        haskey(r.owned, n) || return nothing
        S = supertype(Base.unwrap_unionall(r.owned[n]))
        return S === Any ? nothing : nameof(S)
    end
    # 🔴 **An alias standing in for the real parent is a BETTER spelling, not a wrong one**, and this
    # is what stops the fixer destroying the documentation. The source says
    # `ContinuousTolerance <: AbstractTolerance{A, V}`, but `AbstractTolerance` is an alias, so
    # `supertype` resolves straight through it to `AbstractSpeciesRequirement` - rewriting the edge
    # that way would erase the very branch the Tolerances diagram exists to show.
    function namesthesame(parent, real)
        parent in r.aliases || return false
        for m in ourmodules()
            isdefined(m, parent) || continue
            u = Base.unwrap_unionall(getfield(m, parent))
            u isa DataType && u.name.name === real && return true
        end
        return false
    end

    out = IOBuffer()
    last = 0
    changes = 0
    for block in mermaidblocks(text)
        Base.write(out, SubString(text, last + 1, block.start - 1))
        scope = docnames("```mermaid\n" * block.body * "```")
        kept = String[]
        present = Set{Tuple{Symbol, Symbol}}()
        for line in split(block.body, '\n')
            mm = match(r"^(\s*)([A-Z][A-Za-z0-9_]*)\s+<\|--\s+([A-Z][A-Za-z0-9_]*)\s*$",
                       line)
            cls = match(r"^\s*class\s+([A-Z][A-Za-z0-9_]*)", line)
            if !isnothing(cls) && isstale(Symbol(cls.captures[1]))
                println(io, "  - ", strip(line), "   (defined nowhere)")
                changes += 1
                continue
            end
            if !isnothing(mm)
                parent, child = Symbol(mm.captures[2]), Symbol(mm.captures[3])
                real = realparent(child)
                if isstale(parent) || isstale(child)
                    println(io, "  - ", strip(line), "   (defined nowhere)")
                    changes += 1
                    continue
                elseif !isnothing(real) && real !== parent && real in scope &&
                       !namesthesame(parent, real)
                    println(io, "  - ", strip(line),
                            "   (its real supertype is ", real, ")")
                    changes += 1
                    continue
                end
                push!(present, (parent, child))
            end
            push!(kept, line)
        end
        # ...then anything the code says and the block does not. ⚠️ Only for a block that already
        # draws inheritance: one with none is a **composition** diagram (the `*--` edges), and adding
        # `<|--` lines to it would bury the picture it was drawn to give.
        isempty(present) &&
            (Base.write(out, join(kept, '\n')); last = block.stop; continue)
        for child in sort(collect(scope))
            real = realparent(child)
            (isnothing(real) || !(real in scope) || (real, child) in present) &&
                continue
            # ...nor one already drawn through an alias of that same parent.
            any(e -> e[2] === child && namesthesame(e[1], real), present) &&
                continue
            line = "    $real <|-- $child"
            println(io, "  + ", strip(line), "   (missing)")
            # ⚠️ **Before the trailing blank**, not after it. A fenced block's body ends in a newline,
            # so `split` leaves an empty final element; appending past it puts the new edge on the
            # same line as the closing fence and silently breaks the diagram. ✅ It did, once.
            insert!(kept, something(findlast(!isempty, kept), length(kept)) + 1,
                    line)
            changes += 1
        end
        Base.write(out, join(kept, '\n'))
        last = block.stop
    end
    Base.write(out, SubString(text, last + 1, lastindex(text)))
    new = String(take!(out))
    if write && new != text
        Base.write(DOC, new)
        println(io, "\n✅ rewrote data/architecture.md - ", changes,
                " line(s) changed.")
        println(io,
                "🔴 Now READ the prose around each change: `--fix` cannot see it, and a deleted type " *
                "is usually named in a sentence too.")
    elseif changes == 0
        println(io,
                "\n✅ every mermaid block's inheritance edges already match the code.")
    else
        println(io, "\n⚠️ ", changes,
                " line(s) would change - re-run with `--fix` to apply.")
    end
    return changes
end

end # module

# Script entry point: audit by default, `--fix` to repair the diagrams.
if abspath(PROGRAM_FILE) == @__FILE__
    using .ArchitectureAudit
    if "--fix" in ARGS
        ArchitectureAudit.regenerate(write = true)
    else
        ArchitectureAudit.audit()
        ArchitectureAudit.regenerate(write = false)
    end
end
