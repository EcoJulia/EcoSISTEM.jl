# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Is every constructor written *after* the file that declares its type? Run the audit with:
#
#     julia --project=data/src -e 'using Pkg; Pkg.test(; test_args = ["extras_clean.jl"])'
#
# or directly, for the report:
#
#     julia --project=data/src data/src/constructororder.jl
#
# 🔴 **Why this exists.** A constructor written in a file included BEFORE its type's own file is
# SILENTLY LOST: the name is undefined at that point, so Julia makes a fresh generic function, and
# the later `struct` replaces the binding without a word. The loader is happy, every test that does
# not call it passes, and only a call fails - with a `MethodError` naming a constructor that is
# plainly there in the source.
#
# ⭐ This reads the source text rather than the loaded package, so it sees a constructor that has
# already been discarded - which reflection cannot, the binding being gone by then. That is why it
# needs no weak dependencies: nothing here asks the runtime what exists.
module ConstructorOrderAudit

using EcoSISTEM
using Base.JuliaSyntax

export constructororder_report

const JS = Base.JuliaSyntax

# The name at the head of a signature or declaration, unwrapping the layers a signature can carry
# (`where`, `::`, `curly` and friends). Returns `nothing` for anything qualified (`Base.foo`), which
# is someone else's name and cannot be one of our constructors.
function _headname(n)
    k = JS.kind(n)
    cs = JS.children(n)
    if k in
       (JS.K"where", JS.K"::", JS.K"call", JS.K"curly", JS.K"<:", JS.K"quote",
        JS.K"parens")
        return (cs === nothing || isempty(cs)) ? nothing : _headname(cs[1])
    elseif k == JS.K"."
        return nothing
    end
    return (cs === nothing || isempty(cs)) ? strip(string(n)) : nothing
end

# The `include(...)` order declared in `src/EcoSISTEM.jl` - the order Julia actually loads the
# files in, which is the only thing that decides whether a name exists yet.
function _includeorder(root)
    order = String[]
    for l in readlines(joinpath(root, "src", "EcoSISTEM.jl"))
        m = match(r"^\s*include\(\"([^\"]+)\"\)", l)
        m === nothing || push!(order, String(m.captures[1]))
    end
    return order
end

# A documented definition is wrapped in a `doc` node with the definition itself second, so both
# passes below must unwrap before asking what kind of thing they are looking at.
# 🔴 Missing this on the DECLARATION pass made the audit nearly vacuous: it saw only the handful of
# types with no docstring - 8 of 248 - so almost every type was unknown to it and every constructor
# of one was skipped as "not ours" rather than checked.
function _undoc(c)
    return JS.kind(c) == JS.K"doc" && length(JS.children(c)) >= 2 ?
           JS.children(c)[2] : c
end

# Parse one source file, returning `nothing` rather than throwing on a file that will not parse -
# a broken file is the syntax check's problem, not this audit's.
function _parsefile(path)
    isfile(path) || return nothing
    return try
        JS.parseall(JS.SyntaxNode, read(path, String), filename = path)
    catch
        nothing
    end
end

"""
    constructororder_report(root = pkgdir(EcoSISTEM))

Return what the source says about where each type is declared and where it is constructed.

`root` is the package directory to read, and exists so the audit can be pointed at a fabricated
tree: a check for something that has never happened cannot be believed until it has been shown to
fire, and only a synthetic tree can arrange the violation.

The result is a named tuple:

  - `violations`: a vector of named tuples `(type, file, declaredin)` - a constructor written in
    `file`, which is included *before* `declaredin` where its type is declared, and which Julia
    therefore discards silently.
  - `ntypes`: how many type declarations were located.
  - `nconstructors`: how many definitions were matched to a declared type.
  - `nfiles`: how many included source files were read.
"""
function constructororder_report(root = pkgdir(EcoSISTEM))
    order = _includeorder(root)
    pos = Dict(f => i for (i, f) in enumerate(order))

    # where each type is DECLARED
    decl = Dict{String, String}()
    nfiles = 0
    for f in order
        t = _parsefile(joinpath(root, "src", f))
        t === nothing && continue
        nfiles += 1
        for c0 in JS.children(t)
            c = _undoc(c0)
            JS.kind(c) in (JS.K"struct", JS.K"abstract") || continue
            nm = _headname(JS.children(c)[1])
            nm === nothing || (decl[nm] = f)
        end
    end

    violations = NamedTuple[]
    nconstructors = 0
    for f in order
        t = _parsefile(joinpath(root, "src", f))
        t === nothing && continue
        for c0 in JS.children(t)
            cc = _undoc(c0)
            JS.kind(cc) in (JS.K"function", JS.K"=") || continue
            nm = _headname(JS.children(cc)[1])
            nm === nothing && continue
            haskey(decl, nm) || continue
            nconstructors += 1
            if get(pos, f, 0) < get(pos, decl[nm], 0)
                push!(violations, (type = nm, file = f, declaredin = decl[nm]))
            end
        end
    end

    return (violations = violations, ntypes = length(decl),
            nconstructors = nconstructors, nfiles = nfiles)
end

if abspath(PROGRAM_FILE) == @__FILE__
    r = constructororder_report()
    println("constructor order vs the include order")
    println("  type declarations found : ", r.ntypes)
    println("  constructors matched    : ", r.nconstructors)
    println("  source files read       : ", r.nfiles)
    println()
    if isempty(r.violations)
        println("✅ every constructor follows its type's own file")
    else
        println("🔴 ", length(r.violations),
                " constructor(s) silently discarded:")
        for v in r.violations
            println("   ", v.type, " constructed in ", v.file,
                    " but declared in ",
                    v.declaredin)
        end
    end
end

end
