# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Is every type declared *after* everything its declaration needs? Run the audit with:
#
#     julia --project -e 'using Pkg; Pkg.test(; test_args = ["extras_clean.jl"])'
#
# or directly, for the report:
#
#     julia --project docs/typeorder.jl
#
# ⚠️ Must go through `Pkg.test` (or have the weak deps loaded) or the extension types are invisible
# and the audit reports false gaps — the same constraint `docs/architecture.jl` has.
#
# ⭐ **Why this exists.** The include order in `src/EcoSISTEM.jl` is a real invariant — a struct's
# field types, its supertype and its type-parameter bounds must all be defined when it is declared —
# and it used to be maintained by hand, defended by comments explaining which file had been hoisted
# above which. This checks it instead.

module TypeOrderAudit

using EcoSISTEM
using Base.JuliaSyntax
const JS = Base.JuliaSyntax

export typeorder_report

const ROOT = pkgdir(EcoSISTEM)

# The modules whose types are ours to order. The extensions are included when loaded; when they are
# not, their types simply do not appear and the audit covers less rather than reporting wrongly.
function ourmodules()
    mods = Any[EcoSISTEM, EcoSISTEM.ClimatePref, EcoSISTEM.Units]
    for e in (:EcoSISTEMMPIExt, :EcoSISTEMPhyloExt, :EcoSISTEMRasterDataSourcesExt,
        :EcoSISTEMERAExt, :EcoSISTEMDataPipelineExt)
        m = Base.get_extension(EcoSISTEM, e)
        m === nothing || push!(mods, m)
    end
    return mods
end

# The name a declaration head binds, unwrapping `where`, `<:`, `{}` and qualification.
function _headname(n)
    k = JS.kind(n)
    if k in (JS.K"where", JS.K"::", JS.K"call", JS.K"curly", JS.K"<:")
        cs = JS.children(n)
        return (cs === nothing || isempty(cs)) ? nothing : _headname(cs[1])
    elseif k == JS.K"Identifier"
        return string(n)
    elseif k == JS.K"."
        cs = JS.children(n)
        return (cs === nothing || isempty(cs)) ? nothing : _headname(cs[end])
    end
    return nothing
end

# Every `.jl` file under `src/` and `ext/`, and for each the type and const names it declares.
function _declarations()
    types = Dict{String, Tuple{String, Int}}()    # type name  => (file, line)
    consts = Dict{String, Tuple{String, Int}}()   # const name => (file, line)
    for d in ("src", "ext")
        isdir(joinpath(ROOT, d)) || continue
        for (root, _, files) in walkdir(joinpath(ROOT, d)), f in files
            endswith(f, ".jl") || continue
            path = joinpath(root, f)
            rel = relpath(path, ROOT)
            tree = try
                JS.parseall(JS.SyntaxNode, read(path, String), filename = path)
            catch
                continue
            end
            _collect!(tree, rel, types, consts)
        end
    end
    return (types = types, consts = consts)
end
function _collect!(n, rel, types, consts)
    k = JS.kind(n)
    cs = JS.children(n)
    if k in (JS.K"struct", JS.K"abstract") && cs !== nothing && !isempty(cs)
        nm = _headname(cs[1])
        nm === nothing || get!(types, nm, (rel, JS.source_line(n)))
    elseif k == JS.K"macrocall" && cs !== nothing && length(cs) >= 2
        # 🔴 Types declared BY A MACRO are invisible to a `struct`/`abstract type` walk, and they are
        # not a rare case — every niche axis (`@nicheaxis`) and both SimpleTraits markers
        # (`@traitdef`) arrive this way, 51 types in all. Each names its type in its first argument.
        mac = JS.kind(cs[1]) == JS.K"MacroName" ? string(cs[1]) :
              _headname(cs[1])
        if mac in ("@nicheaxis", "@traitdef")
            nm = _headname(cs[2])
            nm === nothing || get!(types, nm, (rel, JS.source_line(n)))
        end
    elseif k == JS.K"const" && cs !== nothing && !isempty(cs) &&
           JS.kind(cs[1]) == JS.K"="
        a = JS.children(cs[1])
        if a !== nothing && !isempty(a)
            nm = _headname(a[1])
            nm === nothing || get!(consts, nm, (rel, JS.source_line(n)))
        end
    end
    return cs === nothing || for c in cs
        _collect!(c, rel, types, consts)
    end
end

# The source as Julia actually sees it: a list of `(file, firstline, lastline)` segments in
# expansion order. 🔴 A whole-file ranking is NOT enough — `src/EcoSISTEM.jl` declares types *after*
# its own includes (`MPIEcosystem`), so its content is interleaved with every file it pulls in.
function _segments()
    segs = Tuple{String, Int, Int}[]
    visited = Set{String}()
    function walk(rel)
        rel in visited && return
        push!(visited, rel)
        path = joinpath(ROOT, rel)
        isfile(path) || return
        lines = readlines(path)
        dir = dirname(rel)
        from = 1
        for (i, l) in enumerate(lines)
            m = match(r"^\s*include\(\"([^\"]+)\"\)", l)
            m === nothing && continue
            i > from && push!(segs, (rel, from, i - 1))
            walk(normpath(joinpath(dir, m[1])))
            from = i + 1
        end
        return from <= length(lines) && push!(segs, (rel, from, length(lines)))
    end
    walk("src/EcoSISTEM.jl")
    # Extensions load after the parent in every case, so they simply follow, in a stable order.
    for (root, _, files) in walkdir(joinpath(ROOT, "ext")), f in sort(files)
        endswith(f, ".jl") || continue
        rel = relpath(joinpath(root, f), ROOT)
        rel in visited && continue
        push!(visited, rel)
        push!(segs, (rel, 1, countlines(joinpath(ROOT, rel))))
    end
    return segs
end

# Where a declaration sits in that expanded stream — the number a comparison actually needs.
function _position(segs, rel, line)
    for (i, (f, lo, hi)) in enumerate(segs)
        f == rel && lo <= line <= hi && return i
    end
    return nothing
end

# Every named type mentioned inside a type expression, following a TypeVar to its BOUND — a field
# declared as a bare parameter (`habitat::Part where Part <: AbstractHabitat`) carries its real
# constraint there, and reading only `fieldtypes` would see `Any`.
function _mentioned(t, acc = Set{Any}())
    t isa TypeVar && return _mentioned(t.ub, acc)
    t isa UnionAll && return _mentioned(t.body, acc)
    if t isa Union
        for u in Base.uniontypes(t)
            _mentioned(u, acc)
        end
        return acc
    end
    if t isa DataType
        push!(acc, t.name.wrapper)
        for p in t.parameters
            (p isa Type || p isa TypeVar) && _mentioned(p, acc)
        end
    end
    return acc
end

# 🔴 A struct often declares its own parameters UNBOUNDED and inherits the constraint from its
# supertype — `struct GridHabitat{H, B, L} <: AbstractHabitat{H, B, L}`, where it is `AbstractHabitat`
# that says `H <: AbstractRegime`. Recover those, or every such field reads as depending on nothing.
function _inheritedbounds(w)
    out = Dict{TypeVar, Any}()
    S = Base.unwrap_unionall(w)
    S isa DataType || return out
    P = supertype(S)
    while P isa DataType && P !== Any
        declared = TypeVar[]
        D = P.name.wrapper
        while D isa UnionAll
            push!(declared, D.var)
            D = D.body
        end
        for (i, arg) in enumerate(P.parameters)
            arg isa TypeVar || continue
            i <= length(declared) || continue
            b = declared[i].ub
            b === Any && continue
            (!haskey(out, arg) || out[arg] === Any) && (out[arg] = b)
        end
        P = supertype(P)
    end
    return out
end

# What must already be defined for `w`'s DECLARATION to compile: its literal supertype, its declared
# field types and its parameter bounds. ⚠️ Deliberately NOT the concrete implementors of an abstract
# field type — `Ecosystem` needs `AbstractHabitat`, not `GridHabitat`.
function _needs(w, universe)
    out = Set{Any}()
    S = Base.unwrap_unionall(w)
    sup = supertype(S)
    sup isa DataType && sup !== Any && union!(out, _mentioned(sup))
    D = w
    while D isa UnionAll
        D.var.ub === Any || union!(out, _mentioned(D.var.ub))
        D = D.body
    end
    if !isabstracttype(S)
        bounds = _inheritedbounds(w)
        try
            for f in fieldtypes(S)
                union!(out, _mentioned(f))
                f isa TypeVar && f.ub === Any && haskey(bounds, f) &&
                    union!(out, _mentioned(bounds[f]))
            end
        catch
        end
    end
    intersect!(out, universe)
    delete!(out, w)
    return out
end

# Does `nm` name a **type alias** rather than a value constant? Resolved through the binding, so a
# parametric alias (a `UnionAll`) counts exactly as a plain one does. ⚠️ `isdeprecated` first, or the
# lookup fires a `depwarn` for every shim.
function _isaliasof(mods, nm)
    for m in mods
        isdefined(m, Symbol(nm)) || continue
        Base.isdeprecated(m, Symbol(nm)) && continue
        v = try
            getfield(m, Symbol(nm))
        catch
            continue
        end
        return v isa Type
    end
    return false
end

"""
    typeorder_report()

Audit the include order against what each type's declaration actually needs.

Returns a named tuple:

  - `violations`: a vector of named tuples `(type, file, needs, needsfile)` — a type declared in a
    file that comes *before* the file declaring something it depends on. Should be empty.
  - `unplaced`: names of types whose declaring file could not be found in the source, so nothing
    could be checked about them.
  - `legacy`: types and **type aliases** still declared in a not-yet-redistributed `src/_*.jl`.
    Counts down to zero
    across the reorganisation, and is its completion test.
  - `ntypes`, `nedges`, `nfiles`: sizes, so a check can assert the audit is not vacuous.

⚠️ **Covers types and `const` type aliases.** A type that names an alias in its *declaration* is
checked through the alias's own target, since at run time an alias is indistinguishable from what it
aliases — so the alias's file is checked against its members, and the type's against the real
underlying types.
"""
function typeorder_report()
    mods = ourmodules()
    decls = _declarations()
    segs = _segments()

    universe = Set{Any}()
    wrappers = Dict{Any, String}()          # wrapper => the name it is declared under
    for m in mods, n in names(m, all = true)
        startswith(string(n), "#") && continue
        isdefined(m, n) || continue
        # ⚠️ Skip deprecated bindings BEFORE touching them: reaching one with `getfield` fires its
        # `depwarn`, and this walk would otherwise print ~20 warnings into every `extras_clean` run.
        # They are aliases onto types already reached under their real names, so nothing is lost.
        Base.isdeprecated(m, n) && continue
        v = try
            getfield(m, n)
        catch
            continue
        end
        v isa Type || continue
        S = Base.unwrap_unionall(v)
        S isa DataType || continue
        parentmodule(S) in mods || continue
        w = S.name.wrapper
        push!(universe, w)
        get!(wrappers, w, string(S.name.name))
    end

    # where each type's declaration sits in the expanded source
    at = Dict{Any, Int}()
    unplaced = String[]
    for w in universe
        nm = wrappers[w]
        d = get(decls.types, nm, nothing)
        p = d === nothing ? nothing : _position(segs, d[1], d[2])
        p === nothing ? push!(unplaced, nm) : (at[w] = p)
    end

    violations = NamedTuple[]
    nedges = 0
    where_(w) = (d = decls.types[wrappers[w]]; string(d[1], ":", d[2]))
    for w in universe
        haskey(at, w) || continue
        for d in _needs(w, universe)
            haskey(at, d) || continue
            nedges += 1
            # A name declared as both a parent abstract and an extension concrete resolves, for a
            # dependency in `src/`, to the parent's — one wrapper per name would otherwise report a
            # false violation against the extension copy.
            wrappers[d] == wrappers[w] && continue
            if at[d] > at[w]
                push!(violations,
                      (type = wrappers[w], file = where_(w),
                       needs = wrappers[d], needsfile = where_(d)))
            end
        end
    end

    # `const` type aliases: the alias must follow everything it names.
    for (nm, d) in decls.consts
        p = _position(segs, d[1], d[2])
        p === nothing && continue
        v = nothing
        for m in mods
            isdefined(m, Symbol(nm)) || continue
            Base.isdeprecated(m, Symbol(nm)) && continue
            v = try
                getfield(m, Symbol(nm))
            catch
                nothing
            end
            v === nothing || break
        end
        (v isa Type) || continue
        for t in intersect(_mentioned(v), universe)
            haskey(at, t) || continue
            wrappers[t] == nm && continue          # the alias naming its own target
            if at[t] > p
                push!(violations,
                      (type = nm, file = string(d[1], ":", d[2]),
                       needs = wrappers[t], needsfile = where_(t)))
            end
        end
    end

    # ⭐ **The Phase A progress meter.** Every `src/_*.jl` is a file whose declarations have not yet
    # been redistributed, so this counts down to zero as the reorganisation proceeds and becomes the
    # completion test. ⚠️ `src/deprecations.jl` is deliberately *not* prefixed and so never counts.
    legacy = String[]
    legacyvalues = String[]
    for (nm, d) in decls.types
        startswith(basename(d[1]), "_") && push!(legacy, nm)
    end
    # 🔴 **A value `const` is NOT a Phase A declaration.** `_WINDOW_PAD`, `px`, `_SUPPLY_SIZE` and the
    # rest specialise no type — they are tuning constants belonging to the *functions* that read them,
    # so they travel in Phase B. Counting them here would give the Phase A ratchet a floor it could
    # never reach. ⭐ A **type alias** does count, and the two are told apart exactly: an alias's
    # binding `isa Type`.
    for (nm, d) in decls.consts
        startswith(basename(d[1]), "_") || continue
        push!(_isaliasof(mods, nm) ? legacy : legacyvalues, nm)
    end

    return (violations = violations, unplaced = sort(unplaced),
            legacy = sort(legacy), legacyvalues = sort(legacyvalues),
            ntypes = length(at), nedges = nedges, nfiles = length(segs))
end

if abspath(PROGRAM_FILE) == @__FILE__
    r = typeorder_report()
    println("types placed: ", r.ntypes, "   dependency edges: ", r.nedges,
            "   source segments: ", r.nfiles)
    isempty(r.unplaced) || println("unplaced: ", join(r.unplaced, ", "))
    println("still declared in a legacy `src/_*.jl`: ", length(r.legacy),
            " types/aliases, plus ", length(r.legacyvalues),
            " value constants (Phase B)")
    if isempty(r.violations)
        println("✅ include order is valid — every type follows what its declaration needs")
    else
        println("🔴 ", length(r.violations), " violations:")
        for v in r.violations
            println("   ", v.type, " (", v.file, ") needs ", v.needs, " (",
                    v.needsfile, ")")
        end
    end
end

end
