# SPDX-License-Identifier: LGPL-3.0-or-later

module TestEcoSISTEM

using EcoSISTEM
using DimensionalData
using Hwloc
using Test

struct _CacheTestOwner end

@testset "assetdir with owner" begin
    dir = EcoSISTEM.assetdir(owner = _CacheTestOwner)
    @test isdir(dir)
    @test basename(dir) == "_CacheTestOwner"
    @test dirname(dir) == EcoSISTEM.assetdir()
    rm(dir, recursive = true, force = true)
end

@testset "CachedAsset / assetpath" begin
    # A `file://` URL needs no network access, so this test isn't flaky.
    src = tempname()
    write(src, "cached asset test content")
    asset = EcoSISTEM.CachedAsset(_CacheTestOwner, "file://" * src)
    path = EcoSISTEM.assetpath(asset)
    @test isfile(path)
    @test read(path, String) == "cached asset test content"
    @test dirname(path) == EcoSISTEM.assetdir(owner = _CacheTestOwner)
    # a second call reuses the cache without re-downloading — remove the source and confirm the
    # already-cached copy is returned, unchanged
    rm(src)
    @test EcoSISTEM.assetpath(asset) == path
    @test isfile(path)
    rm(dirname(path), recursive = true, force = true)
end

@testset "a failed download leaves no cache residue" begin
    # `isfile(path)` is `assetpath`'s only cache check, so anything left at that path after a failure
    # would be served as valid forever. This test does *not* distinguish the temp-then-rename
    # implementation from writing straight to the destination: `Downloads.download` already deletes
    # its output on any raised error (`ArgTools.arg_write`'s bare `catch`), so both satisfy it. The
    # case the rename actually guards — a process killed outright mid-transfer — cannot be simulated
    # in-process. What is locked in here is the observable invariant: after a failure there is
    # neither a cache entry nor a stray temporary, and the asset still recovers afterwards.
    missing_src = tempname()                       # deliberately never created
    asset = EcoSISTEM.CachedAsset(_CacheTestOwner, "file://" * missing_src)
    dir = EcoSISTEM.assetdir(owner = _CacheTestOwner)
    @test_throws Exception EcoSISTEM.assetpath(asset)
    # No cache entry, and no `tempname` residue either — the destination must be absent, not empty.
    @test !isfile(joinpath(dir, basename(missing_src)))
    @test isempty(readdir(dir))

    # …and the path stays usable afterwards: once the source exists, the same asset downloads
    # cleanly rather than being permanently poisoned by the earlier failure.
    write(missing_src, "recovered")
    path = EcoSISTEM.assetpath(asset)
    @test read(path, String) == "recovered"
    rm(missing_src)
    rm(dir, recursive = true, force = true)
end

@testset "assetpath with a bare relative path" begin
    @test EcoSISTEM.assetpath("some/file.txt") ==
          joinpath(EcoSISTEM.assetdir(), "some/file.txt")
end

@testset "species_blocksize" begin
    bs = EcoSISTEM.species_blocksize()
    # Used as an inner-loop block width in update! (cld(spp, bs), 1:bs:...), so it
    # must be a positive integer.
    @test bs isa Int
    @test bs >= 1
    # It is the detected CPU cache line divided by the abundance element size
    # (Int64), or the 128-byte-line fallback (16) if hwloc detection failed.
    expected = try
        max(1, Hwloc.cachelinesize() ÷ sizeof(Int64))
    catch
        16
    end
    @test bs == expected
end

# Every abstract type this package owns is `public` and not exported — a supertype is something to
# dispatch on or subtype, not everyday vocabulary, so it does not earn a slot in every
# `using EcoSISTEM` namespace. Asserted by *visibility* rather than by membership of the `public`
# block at the foot of `src/EcoSISTEM.jl`, so a hierarchy declared anywhere is still covered.
@testset "abstract types are public, not exported" begin
    # `MPIEcosystem`/`MPIGridLandscape` are the deliberate exceptions: abstract only because their
    # concrete subtypes live in the MPI extension, named directly by users, and exported in v0.4.0.
    exceptions = (:MPIEcosystem, :MPIGridLandscape)

    # Abstract types owned by `m` (not merely visible in it, which would sweep in Base's).
    #
    # A deprecated binding is skipped *before* `getfield` reaches it, and that is the whole reason
    # `Base.isdeprecated` is here: reading one fires its `depwarn`, so a sweep over every name emits
    # one warning per deprecated alias. Measured at **51 warnings from this single line**, which is
    # every deprecation warning `core_test` produced. Permanent log noise trains a reader to ignore
    # warnings, which is what makes it worth a line of code.
    #
    # Skipping them costs nothing here: every one is an alias of a type this sweep sees under its
    # live name anyway, so none could contribute an abstract type of its own.
    function ownabstracts(m)
        return filter(names(m, all = true, imported = false)) do n
            (isdefined(m, n) && !Base.isdeprecated(m, n)) || return false
            v = getfield(m, n)
            return v isa Type && isabstracttype(v) && parentmodule(v) === m
        end
    end

    # Every niche axis is abstract, but the rule above does not apply to them: they are the
    # user-facing vocabulary (`axis = Temperature`), not supertypes users never name, and they are
    # abstract only so that a leaf can be subdivided later without changing the hierarchy. They get
    # their own two-sided rule below instead.
    isaxis(m, n) = getfield(m, n) <: EcoSISTEM.NicheAxis

    for m in (EcoSISTEM, EcoSISTEM.ClimatePref)
        abstracts = ownabstracts(m)
        # A guard on the guard: if this ever finds nothing, the filter has broken rather than the
        # package having no hierarchies.
        @test length(abstracts) > 20 || m === EcoSISTEM.ClimatePref
        for n in abstracts
            (n in exceptions || isaxis(m, n)) && continue
            @test Base.ispublic(m, n)
            @test !Base.isexported(m, n)
        end
    end

    # …and the two exceptions really are still exported, so the carve-out cannot rot into silence.
    for n in exceptions
        @test Base.isexported(EcoSISTEM, n)
    end

    # The axis carve-out gets the same treatment: axes split by whether they are a *grouping* — a
    # `⋯Axis`-suffixed intermediate, public but not exported, since naming one is unusual — or a
    # leaf that a user writes as `axis =`, which must be exported. Nothing may be merely private.
    axes = filter(n -> isaxis(EcoSISTEM, n), ownabstracts(EcoSISTEM))
    @test length(axes) > 40
    for n in axes
        @test Base.ispublic(EcoSISTEM, n)
        if endswith(String(n), "Axis")
            @test !Base.isexported(EcoSISTEM, n)
        else
            @test Base.isexported(EcoSISTEM, n)
        end
    end
end

# **A layer's array type is an implementation detail, and must not cross the public boundary.**
# The package has migrated `AxisArrays` → `DimensionalData` once and may migrate again; every public
# function that hands out a bare array makes that migration *users'* problem too. Asserting it here
# turns the next leak into a failing test rather than something to be noticed in review.
#
# `Union{}` return types are **inference failures on error paths, not leaks** — `GridHabitat`,
# `build_species`, `update!` and friends all appear that way — so the check looks for a genuine
# `AbstractDimArray` in the inferred results and ignores the empty union.
@testset "no public function hands out a naked array" begin
    # `Union{} <: AbstractDimArray` is **true** — the empty union is a subtype of everything — so
    # excluding it is what separates a real leak from an inference failure, and without the guard
    # this reports five. Measured: `GridHabitat`, `build_species`, `update!`,
    # `investigate_study_area` and `emptyMPIgridlandscape` all infer `Union{}` on an error path.
    isdim(T) = T !== Union{} &&
               (T <: DimensionalData.AbstractDimArray ||
                (T isa Union && (isdim(T.a) || isdim(T.b))))

    # **One carve-out, and it is a released API constraint rather than an oversight.** `readfile`
    # returns a bare `DimArray` and has done since v0.4.0, when it was exported by
    # `EcoSISTEM.ClimatePref`; the check only walks `EcoSISTEM` and `EcoSISTEM.Units`, so it never saw
    # it. Dissolving that submodule moved `readfile` into the parent and exported it there, which is
    # what surfaced this. Changing what it returns would break released code, so it is recorded here
    # instead of hidden — a named exception, like the other lists in this file, so it has to be looked
    # at rather than passing silently.
    #
    # `hasdata` was a carve-out too until the combine contract was normalised: every combine now
    # returns a `ClimateRaster`, a mask being one whose element type is `Bool`.
    #
    # VERSION-DEPENDENT, and measured: on Julia 1.11 inference does not surface `readfile`'s
    # `DimArray` return through `Base.return_types`, so it is not detected and the list is empty.
    # The named list is what makes that visible — it FAILED on 1.11 rather than passing quietly,
    # which is the whole point of naming entries instead of counting them. What it also says is that
    # this check catches strictly less on 1.11 than on 1.12: a real leak could go unreported there.
    known_naked = VERSION >= v"1.12" ? [:readfile] : Symbol[]
    leaks = Symbol[]
    for n in names(EcoSISTEM, all = false)
        isdefined(EcoSISTEM, n) || continue
        f = getfield(EcoSISTEM, n)
        f isa Function || continue
        for m in methods(f)
            rts = try
                Base.return_types(f, Base.tuple_type_tail(m.sig))
            catch
                continue
            end
            any(isdim, rts) && (push!(leaks, n); break)
        end
    end
    @test leaks == known_naked

    # And the detector is not vacuous: a function that *does* return one is caught.
    _leaky(r::EcoSISTEM.ClimateRaster) = r.array
    @test any(isdim,
              Base.return_types(_leaky,
                                Tuple{EcoSISTEM.ClimateRaster}))
end

# **A comment between a docstring and its definition silently detaches it.** Julia does not
# error, `@doc` simply returns `nothing`, the package loads and every test passes — only the
# Documenter build notices, and then as an unresolved `@ref` somewhere else entirely. Found the
# hard way: five public docstrings were lost that way in one sitting, including
# `GridHabitat`'s.
#
# A control for the placeholder detector below, and the reason it is out here rather than inside the
# testset: `@doc` has to run at a module's top level to register anything, and `@testset begin` is a
# `let` block. Neither name is exported, so the sweep itself never reaches them.
_placeholder_control_() = nothing
@doc (@doc _no_such_binding_anywhere_) _placeholder_control_
"""    _documented_control_() — a real docstring, so the detector must NOT flag this one. """
_documented_control_() = nothing

# Asserted over the whole public surface rather than the five, because the mistake is easy to
# repeat and costs nothing to check here — `core_test` is a much shorter feedback loop than the
# docs build.
@testset "every public name keeps its docstring" begin
    # Checked through the **binding**, not `@doc(value)`: `@doc` is a macro over an *expression*,
    # so handing it a variable documents the variable and reports every name as undocumented. That
    # mistake made this very test claim the whole package was undocumented.
    hasdoc(m, n) = haskey(Base.Docs.meta(m), Base.Docs.Binding(m, n))
    mods = (EcoSISTEM, EcoSISTEM.ClimatePref, EcoSISTEM.Units)

    undocumented = Symbol[]
    for m in mods, n in names(m, all = false)
        (isdefined(m, n) && n !== nameof(m)) || continue
        v = getfield(m, n)
        # A method-less function is an extension hook or a `@deprecate` shim — not this module's
        # to document. A name re-exported from a submodule is documented *there*, so every
        # module's own table is consulted before calling it a gap.
        (v isa Function && isempty(methods(v))) && continue
        any(mm -> hasdoc(mm, n), mods) || push!(undocumented, n)
    end

    # A *named* list, not a count: a new gap has to be looked at, and an old one disappearing
    # should shrink this rather than pass silently. These are `@deprecate`d shims, whose redirect is
    # their documentation.
    known = [
        # `@deprecate`d names — the redirect the macro generates is their documentation.
        :CHELSA_bioclim, :CHELSA_monthly, :Landcover, :Worldclim_bioclim,
        :Worldclim_monthly, :compressLC, :fitBrownian,
        :bioclimAE, :eraAE, :eraChange, :lcAE, :peakedgradAE, :raingradAE,
        :simplehabitatAE, :simplenicheAE, :tempgradAE, :worldclimAE,
        :worldclimChange, :getbudget, :gethabitat, :gettraitrel, :getsize,
        :getgridsize,
        :resetrate!, :traitpopulate!, :traitrepopulate!,
        # The flatcase conformance pass, 2026-08-25: four names that were exported in v0.4.0 and so
        # kept a redirect. The other twenty renamed names owed no shim, being unexported or new on
        # this branch, so they do not appear here.
        :assign_traits!, :gather_abundance, :gather_diversity, :get_traits,
        # `[C7-VIS]` C, 2026-08-26: `clearcache` gained its `!` (it deletes files) and `searchdir`
        # became private. Both were exported in v0.4.0, so both keep a redirect.
        :clearcache, :searchdir,
        # These two were `@deprecate_binding`s until the two categorical tolerances merged, and
        # a binding inherits its target's docstring. They are now `@deprecate` *functions*, because
        # each has to pin its own `penalty` — which a binding cannot supply — so they land here with
        # the rest of the redirects.
        :DiscreteTrait, :LCtrait,
        # `EcoSISTEM.Units` constants — `Unitful.@unit` and the `Dates` month-number re-exports
        # generate the binding, so there is nowhere to hang a docstring of our own.
        :January, :February, :March, :April, :May, :June, :July, :August,
        :September, :October, :November, :December,
        :january_duration, :february_common_year_duration,
        :february_leap_year_duration, :february_mean_duration,
        :march_duration, :april_duration, :may_duration, :june_duration,
        :july_duration, :august_duration, :september_duration,
        :october_duration, :november_duration, :december_duration,
        :month_mean_duration, :quarter_mean_duration,
        :arcminute, :arcsecond, :day, :week, :year,
        # Genuine gaps, recorded rather than hidden — pre-existing, and outside this pass.
        :ContinuousEvolve, :DiscreteEvolve, :emptyMPIgridlandscape]
    @test isempty(setdiff(undocumented, known))

    # **A name can HAVE a docstring that says it is undocumented**, which is the hole the
    # check above cannot see. `@doc (@doc other) name` republishes whatever `@doc other` returned,
    # and for a binding that does not exist that is not an error but a *placeholder* reading "No
    # documentation found. Binding `X` does not exist." So the name has an entry, `hasdoc` is true,
    # and the gap reaches the docs build. `getlong` shipped that way, its docstring naming a `lat`
    # that had been renamed years earlier.
    #
    # No exception list: publishing a placeholder is always wrong, so a single name here is a
    # failure rather than something to be looked at.
    function placeholding(m, n)
        txt = string(Base.Docs.doc(Base.Docs.Binding(m, n)))
        return occursin("No documentation found", txt) &&
               occursin("does not exist", txt)
    end
    placeholders = [(m, n)
                    for m in mods
                    for n in names(m, all = false)
                    if isdefined(m, n) && n !== nameof(m) && hasdoc(m, n) &&
                           placeholding(m, n)]
    @test isempty(placeholders)

    # **And the detector is not vacuous** — which matters more here than usual, because the
    # failure it guards against is *silent*: a comment between a docstring and its definition
    # detaches the docstring, Julia does not complain, and only the Documenter build notices, as an
    # unresolved `@ref` somewhere else. Five public docstrings were lost that way in one sitting.
    @test !hasdoc(EcoSISTEM, :_no_such_name_exists_)
    for n in (:GridHabitat, :investigate_study_area, :materialise,
        :Varying, :StudyArea)
        @test hasdoc(EcoSISTEM, n)
    end

    # The placeholder detector against its two controls: a name whose docstring was taken from a
    # binding that does not exist, and one with an ordinary docstring. Both `hasdoc`, so only the
    # content tells them apart — which is the whole point.
    @test hasdoc(@__MODULE__, :_placeholder_control_)
    @test placeholding(@__MODULE__, :_placeholder_control_)
    @test !placeholding(@__MODULE__, :_documented_control_)
end

# A docstring that HAS an `# Arguments` block looks exactly like one that COVERS every argument, and
# a partial list is the harder mistake to see: nothing is missing from the page, so nothing looks
# wrong. `force` was undocumented on all five `public` cell-geometry functions at once, and stayed
# that way through a whole documentation pass because each of them did have an `# Arguments` block.
#
# Keywords only, deliberately. A positional argument is named in the signature header a reader
# already has in front of them, while a keyword is invisible unless the prose says it exists — so
# this is where the asymmetry bites, and where the check is worth its exception list.
@testset "every public function documents its keywords" begin
    hasdoc(m, n) = haskey(Base.Docs.meta(m), Base.Docs.Binding(m, n))
    mods = (EcoSISTEM, EcoSISTEM.ClimatePref, EcoSISTEM.Units)

    # Every keyword any of `f`'s **own** methods declares. Methods from elsewhere are skipped: a
    # `@deprecate` shim is generated into `Base`, and its redirect is its documentation.
    ours = ("EcoSISTEM", "ClimatePref", "Units")
    function keywords(f)
        ks = Symbol[]
        for m in methods(f)
            String(nameof(m.module)) in ours || continue
            for k in Base.kwarg_decl(m)
                # A splat is `kw...`, which names nothing a caller could look up.
                endswith(String(k), "...") || push!(ks, k)
            end
        end
        return unique(ks)
    end

    # Backticked, and as a whole token. A bare substring is far too generous: `cut` matches
    # "cutting", `active` matches "inactive", and `force` matches "enforced" — so prose that never
    # mentions the keyword would satisfy it. Every argument this package documents is written in
    # backticks, alone in a bullet or sharing one, so requiring that costs nothing and closes the
    # hole.
    function _documents(txt::AbstractString, k::Symbol)
        return occursin(Regex("`\\s*" * String(k) * "\\s*(`|=|,|\\))"), txt)
    end

    undocumented = Tuple{Symbol, Symbol}[]
    for m in mods, n in names(m, all = false)
        (isdefined(m, n) && n !== nameof(m)) || continue
        f = getfield(m, n)
        (f isa Function && !isempty(methods(f))) || continue
        # Every module's own table, because a re-exported name is documented where it is defined.
        txt = join(string(Base.Docs.doc(Base.Docs.Binding(mm, n)))
                   for mm in mods if hasdoc(mm, n))
        isempty(txt) && continue
        for k in keywords(f)
            _documents(txt, k) || push!(undocumented, (n, k))
        end
    end

    # A named list, on the same terms as the docstring check above: a new gap has to be looked at,
    # and an old one disappearing shrinks this rather than passing silently.
    known = [
        # `src/actions.jl` is held back from the documentation pass, so its two gaps are recorded
        # here rather than fixed. Both are real.
        (:build_species, :verbosity),
        (:simulate!, :intervention),
        # `@deprecate`d readers, whose `cut` forwards to the replacement's. A name being removed is
        # documented by the redirect the macro generates, and describing its arguments is work with
        # a negative lifespan — the same exemption the docstring check above makes.
        (:readCERA, :cut), (:readCRUTS, :cut), (:readERA, :cut)]
    @test isempty(setdiff(unique(undocumented), known))

    # And the detector is not vacuous: it has to find a keyword that genuinely is not written down,
    # or it would pass on a package with no documentation at all.
    docof(n) = string(Base.Docs.doc(Base.Docs.Binding(EcoSISTEM, n)))
    @test !_documents(docof(:getcellareas), :kw_that_is_never_documented)
    for (n, k) in ((:getcellareas, :force), (:build_habitat, :verbosity),
        (:StudyArea, :cellsize), (:readfile, :cut))
        @test _documents(docof(n), k)
    end
    # And it is not a bare substring search: `readfile`'s own prose says "cutting a global layer",
    # which `occursin` would accept as documentation of a keyword called `cutting`.
    #
    # What it does NOT catch, measured: a keyword mentioned only in passing — `` `force = true` ``
    # in prose satisfies it even with the `# Arguments` bullet removed. That is deliberate. The
    # failure this exists for is a keyword written down **nowhere**, which is how `force` stayed
    # invisible on five functions at once; demanding a bullet would mean parsing the block, and
    # would fail on the shared form `` `units`, `x`, `force` ``.
    @test !_documents(docof(:readfile), :cutting)
end

end
