# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Regenerate `data/overloads.md` — every method this package defines on a generic it does **not**
# own, found by reflection rather than by grepping for `import` lines (which would miss every
# qualified `Base.keys(...)` definition, and there are dozens).
#
# ⚠️ **Must run with the weak dependencies loaded**, or the extensions never load and the Diversity,
# MPI and RasterDataSources methods are invisible. ⭐ `data/src` is that environment - it carries
# MPI, Phylo and RasterDataSources precisely so every extension loads:
#
#     julia --project=data/src data/src/overloads.jl
#
# ⭐ **The page's prose lives in `NOTES` below, not in the generated file** — editing the output by
# hand would lose it on the next run.
using EcoSISTEM, MPI, Phylo, RasterDataSources

const ROOT = pkgdir(EcoSISTEM)
pkgmod(name, uuid) = Base.loaded_modules[Base.PkgId(Base.UUID(uuid), name)]

function ourmodules()
    mods = Any[EcoSISTEM, EcoSISTEM.ClimatePref, EcoSISTEM.Units]
    for e in (:EcoSISTEMMPIExt, :EcoSISTEMPhyloExt, :EcoSISTEMRasterDataSourcesExt,
        :EcoSISTEMERAExt, :EcoSISTEMDataPipelineExt)
        m = Base.get_extension(EcoSISTEM, e)
        m === nothing || push!(mods, m)
    end
    return mods
end
const OURS = Set{Module}(ourmodules())
isours(m::Module) = m in OURS || Base.moduleroot(m) in OURS

function allmodules()
    seen = Set{Module}()
    todo = collect(values(Base.loaded_modules))
    while !isempty(todo)
        m = pop!(todo)
        m in seen && continue
        push!(seen, m)
        for n in names(m, all = true, imported = false)
            isdefined(m, n) || continue
            Base.isdeprecated(m, n) && continue
            v = try
                getfield(m, n)
            catch
                continue
            end
            v isa Module && v !== m && !(v in seen) && push!(todo, v)
        end
    end
    return seen
end

function inrepo(m)
    return startswith(String(m.file), ROOT) &&
           (occursin("/src/", String(m.file)) ||
            occursin("/ext/", String(m.file)))
end

# Does any argument type of `m` belong to us? If none does, the method is type piracy.
# ⚠️ follow BOTH wrappers to the bottom: `Type{<:Foo{S, C}}` is a `UnionAll` over a `TypeVar` whose
# bound is what actually names the type, and stopping at the first `DataType` reports it as foreign.
function _ours(t, depth = 0)
    depth > 8 && return false
    t isa Union && return any(u -> _ours(u, depth + 1), Base.uniontypes(t))
    t isa TypeVar && return _ours(t.ub, depth + 1)
    t = Base.unwrap_unionall(t)
    t isa DataType || return false
    t.name === Type.body.name && length(t.parameters) == 1 &&
        return _ours(t.parameters[1], depth + 1)
    isours(Base.moduleroot(parentmodule(t))) && return true
    return any(p -> (p isa Type || p isa TypeVar) && _ours(p, depth + 1),
               t.parameters)
end
function ispirate(m)
    ts = try
        Base.unwrap_unionall(m.sig).parameters[2:end]
    catch
        return false
    end
    return !any(_ours, ts)
end

function collectall()
    records = Dict{Any, Set{Method}}()
    for M in allmodules(), n in names(M, all = true, imported = false)
        isdefined(M, n) || continue
        Base.isdeprecated(M, n) && continue
        v = try
            getfield(M, n)
        catch
            continue
        end
        (v isa Function || v isa Type) || continue
        ml = try
            methods(v)
        catch
            continue
        end
        for m in ml
            inrepo(m) && push!(get!(records, v, Set{Method}()), m)
        end
    end
    return records
end

loc(m) = string(relpath(String(m.file), ROOT), ":", m.line)
function argsig(m)
    ts = Base.unwrap_unionall(m.sig).parameters[2:end]
    s = join(map(t -> replace(string(t), "EcoSISTEM." => "", "Diversity." => "",
                              "DimensionalData." => "",
                              "Base.Broadcast." => ""), ts), ", ")
    return length(s) > 110 ? s[1:107] * "…" : s
end

# The prose belongs in the generator, not in the generated file — hand-editing the output would be
# lost the next time this runs, which is the one thing a generated page must not invite.
const NOTES = Dict(:Base => """
               ⭐ **The bulk of this is ONE decision**: `src/collections.jl` gives every collection *and* every
               single member the standard container interface, so a leaf answers as a one-member container. That is
               `keys` `values` `pairs` `iterate` `getindex` `length` `merge` `haskey` `get` `firstindex` `lastindex`
               `propertynames` `getproperty` `hasproperty` `NamedTuple` — fifteen generics, each with one or two
               methods over the two unions, replacing the fifteen bespoke `regimes`/`supply_names`/`named_demands`
               accessors that came before.

               ⚠️ **`eltype` is deliberately NOT in that set.** It is a **leaf** property — the unit frame a layer's
               data is in, and what supplies a nichefit's frame parameter — so a collection has no single one and
               asking is an error rather than a guess. Its fourteen methods are all on leaves.

               ⭐ **`ClimateRaster` opts in to broadcasting** (`size` `axes` `ndims` `length` `broadcastable`
               `BroadcastStyle` `copy` `+` `-`, all in `src/Climate.jl`). That is what lets a `ConstructedRasterSpec`
               combine name no array type on either side — `compress_landcover(lc) .!= landcoverclass(:open_water)`
               goes in and a raster comes out.

               ⭐ **`show` is two methods per type**, plain and `MIME"text/plain"`, for the seven types a user
               actually sees printed.

               ⭐ **`hash` and `==` travel together**, on `ReadKey` (a cache key, so it must) and on
               `AbstractLayerFate`.
               """,
                   :Diversity => """
               🔴 **The `Diversity.API` hooks split exactly along the `[TF-FORWARD]` rule**: `SpeciesList <:
               AbstractTypes` but it *holds* the real types in `sl.types`, so every hook on it delegates there.
               Eleven do. ⚠️ Every breach of that rule was silent, because Diversity's defaults are not neutral —
               `_hassimilarity` defaults to **true**. A test sweeps `Diversity.API` rather than listing the hooks,
               so a hook added upstream cannot be missed.

               ⚠️ **A call graph cannot see these.** Their only in-repo mention is the `import Diversity.API: …`
               line, so each looks dead or single-caller; deleting one would silently break the interface. Check any
               candidate against the names we import before calling it dead.

               ⚠️ **`faith_pd`/`generalisedfaith_pd` live in `Diversity.Ecology` and are NOT re-exported into
               `Diversity`** — `isdefined(Diversity, :faith_pd)` is `false`, so extending them must import from the
               owning submodule or precompilation fails. Both are constrained to a phylogenetic species list through
               `SpeciesList`'s **fourth type parameter**, which is a static constraint; dispatching on Diversity's
               own `Sim <: PhyloBranches` can never match an `Ecosystem`, because `AbstractEcosystem` puts the
               `SpeciesList` in Diversity's hardcoded `Sim` slot.
               """,
                   :EcoBase => """
               🔴 **Two types only, and that is the point.** These methods used to sit on `AbstractRegime` as well,
               where `xmin`/`ymin` were hardcoded `0` and `xcellsize` was `Float64(size / km)` — so a geographic grid
               answered `1.0 ° km⁻¹`, silently. **A layer is not a grid**; `StudyGrid` is, and it is the package's
               only `EcoBase.AbstractGrid`.

               ⚠️ **EcoBase's `indices`/`coordinates` are `(x, y)` columns — the opposite order to this package.**
               Read `EcoBase.convert_to_image`, which uses `indices(grd, 1)` as the matrix *column*. The old
               implementation transposed every grid plotted through EcoBase.
               """,
                   :GeoInterface => """
               ⭐ All four are on `LatLong`, declaring a coordinate pair to be a point so that anything speaking
               GeoInterface can consume one.
               """)

records = collectall()
rows = Any[]
for (f, ms) in records
    o = try
        parentmodule(f)
    catch
        nothing
    end
    (o === nothing || isours(o)) && continue
    f === Core.kwcall && continue      # artefact: every keyword method gets one
    push!(rows,
          (root = Base.moduleroot(o), owner = o, name = string(nameof(f)),
           ms = sort(collect(ms), by = m -> (String(m.file), m.line))))
end
ngen = length(rows)
nmeth = sum(r -> length(r.ms), rows)
pirates = sort(reduce(vcat, (filter(ispirate, r.ms) for r in rows),
                      init = Method[]),
               by = m -> (String(m.file), m.line))

byroot = Dict{Module, Vector{Any}}()
for r in rows
    push!(get!(byroot, r.root, Any[]), r)
end
order = sort(collect(keys(byroot)),
             by = m -> (-sum(r -> length(r.ms), byroot[m]), string(nameof(m))))

open(joinpath(ROOT, "data", "overloads.md"), "w") do io
    println(io,
            """
# Overloads — every generic EcoSISTEM extends but does not own

Generated by reflection: every method whose source file is under `src/` or `ext/`, whose generic's
`parentmodule` is **not** one of ours. Run with the weak dependencies loaded (`MPI`, `Phylo`,
`RasterDataSources`), or the extensions are invisible and the counts come out short.

**$ngen generics, $nmeth methods.**

| owner | generics | methods |
|---|---:|---:|""")
    for root in order
        rs = byroot[root]
        println(io, "| `", nameof(root), "` | ", length(rs), " | ",
                sum(r -> length(r.ms), rs), " |")
    end
    println(io)
    println(io,
            """
🔴 **Type piracy — $(length(pirates)) method$(length(pirates) == 1 ? "" : "s")** — a method on a foreign generic where no argument type is ours either. Listed in full at the foot of this page.

⭐ **How to regenerate.** This page is produced by reflection, from the `data/src` environment, which
carries the weak dependencies so that every extension loads:

```julia
julia --project=data/src data/src/overloads.jl
```

⚠️ **The prose below lives in the generator, not here** — editing this file by hand loses it on the
next run.
""")
    for root in order
        rs = sort(byroot[root], by = x -> (string(nameof(x.owner)), x.name))
        nm = sum(r -> length(r.ms), rs)
        println(io, "\n## ", nameof(root), " — ", length(rs),
                length(rs) == 1 ? " generic, " : " generics, ", nm,
                nm == 1 ? " method\n" : " methods\n")
        haskey(NOTES, nameof(root)) && println(io, NOTES[nameof(root)])
        println(io, "| generic | defined at | on |")
        println(io, "|---|---|---|")
        for r in rs
            for m in r.ms
                pir = ispirate(m) ? " 🔴" : ""
                println(io, "| `", nameof(r.owner), ".", r.name, "`", pir,
                        " | `", loc(m), "` | `", argsig(m), "` |")
            end
        end
    end
    println(io,
            """

## 🔴 Type piracy

⭐ **All of these are the sanctioned exception**, and there is no other. `api.md`'s `@autodocs` cannot
see inside an extension, so a public name whose implementation moves keeps its docstring on a
method-less stub in the parent — but `Base.read`'s dataset methods **cannot take a stub**: the parent
would have to define `read(::Type, …)` and pirate `Base.read` for every type in Julia. So `api.md`
names that one extension module in `Modules` instead, and `test/core_ext.jl` carries it as a named
exception rather than a count.
""")
    for m in pirates
        println(io, "- `", loc(m), "` — `", m.name, "(", argsig(m), ")`")
    end
    return println(io,
                   """

   ## Two loose ends worth knowing

   ⚠️ **`Trapezoid` implements only part of the `Distributions` interface.** ✅ Measured on
   `Trapezoid(0, 1, 2, 3)`: `pdf` and `rand` answer, while `cdf`, `logpdf`, `quantile`, `mean` and
   `var` are each a `MethodError` — Distributions supplies no fallback for any of them. That is enough
   for the suitability code, which only ever evaluates the density, but a `Trapezoid` will not work
   everywhere a `ContinuousUnivariateDistribution` is expected.

   ⚠️ **Several `import` lines extend nothing** — `import Diversity.Gamma`, `import Rasters: Projected`,
   `import ArchGDAL`. 🔴 `import ArchGDAL` is **load-bearing**: it registers the Rasters GDAL backend and
   has no referenced symbol, so a dead-import sweep must not touch it. The others are plain uses, where
   `using` would read more honestly.""")
end
println("wrote data/overloads.md: ", ngen, " generics, ", nmeth, " methods, ",
        length(pirates), " pirate")
