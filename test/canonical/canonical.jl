# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The blessing machinery shared by every canonical test. See README.md in this directory for what a
# canonical test is for and how to re-bless one.

module Canonical

using Test
using TOML

export canonical, canonical_reference, writereference, blessing, heavydata

# Where the blessed numbers live. One file for the whole directory, so a re-blessing produces a single
# reviewable diff rather than a scatter of them.
const REFERENCE = joinpath(@__DIR__, "reference.toml")

# Blessed values are written through **on every call**, not accumulated and flushed at the end.
# That looks wasteful and is deliberate: each canonical test file `include`s this file into its own
# module, so every file gets its *own* copy of any in-memory state — an accumulate-then-flush design
# silently blessed nothing at all, because the runner wrote out its own empty dict rather than the
# values the test files had recorded into theirs. Writing through cannot have that bug.
const RECORDED = Dict{String, Any}()

"""
    heavydata()

May this run download the **large** rasters — the CHELSA `BioClimPlus` layers, tens of gigabytes
between them — as opposed to the modest WorldClim ones every canonical real-data test uses?

`true` when running locally, `false` on a CI runner. **The GitHub runners already download enough
to be close to timing out**, so anything that needs a `BioClimPlus` layer is local-only by default.
Set `ECOSISTEM_HEAVY_DATA=true` to force it on (or `=false` to suppress it locally).

 **Skipping is safe for the reference file, by construction**: `writereference` *merges*, so a
re-blessing run that skipped these tests cannot delete their blessed values — see its own docstring.
That is what lets the heavy tests be blessed once, locally, and simply not exercised on CI.

 It is *not* a substitute for checking the same fact more cheaply. The "when is an accumulated
layer divided?" discriminator needs no download at all — it is a catalogue property, checked in
`test_LayerCatalogue.jl` and so run on every CI build. Only the **read values** are gated here.
"""
function heavydata()
    haskey(ENV, "ECOSISTEM_HEAVY_DATA") &&
        return ENV["ECOSISTEM_HEAVY_DATA"] == "true"
    return !haskey(ENV, "RUNNER_OS")
end

"""
    blessing()

Is this a re-blessing run? True when `ECOSISTEM_BLESS=true`, in which case the canonical tests record
their results instead of checking them.
"""
blessing() = get(ENV, "ECOSISTEM_BLESS", "false") == "true"

# The blessed values, or an empty dict the first time. Read once and cached, since every call needs it.
const _CACHE = Ref{Union{Nothing, Dict{String, Any}}}(nothing)
function canonical_reference()
    isnothing(_CACHE[]) &&
        (_CACHE[] = isfile(REFERENCE) ? TOML.parsefile(REFERENCE) :
                    Dict{String, Any}())
    return _CACHE[]
end

# TOML holds numbers, strings and arrays of them — not `Quantity`s. Refusing a united value here, with
# the fix in the message, is much clearer than letting `TOML.print` fail deep inside a write.
function _plain(name, value)
    value isa Real && return float(value)
    value isa AbstractArray{<:Real} && return float.(collect(value))
    return error("canonical value `$name` is a $(typeof(value)); a blessed value must be a real " *
                 "number or an array of them. If it carries a unit, strip it explicitly — " *
                 "`ustrip(u\"mm/d\", x)` — so the reference file records which unit it is in.")
end

"""
    canonical(name, value; rtol = 1e-8)

Compare `value` against the blessed result for `name`, or record it when re-blessing.

 Strip units before calling: the reference is a plain number, and stripping *explicitly* at the call
site is what pins which unit the blessed figure is in. A blessed `7.5398` means nothing on its own.

`rtol` is deliberately tight by default. A canonical test exists to notice change, so a loose
tolerance defeats it; widen it only where a result is genuinely only reproducible to fewer digits, and
say why at the call site.
"""
function canonical(name::AbstractString, value; rtol = 1e-8)
    key = String(name)
    plain = _plain(key, value)
    RECORDED[key] = plain
    if blessing()
        _writethrough(key, plain)
        return @test true                    # nothing to compare against; this run defines it
    end
    ref = canonical_reference()
    if !haskey(ref, key)
        return @test_broken "no blessed value for `$key` — run the canonical suite with " *
                            "ECOSISTEM_BLESS=true to record one" == ""
    end
    return @test isapprox(plain, ref[key], rtol = rtol)
end

# Merge one blessed value into the reference file. Read-modify-write per call, for the reason given
# at `RECORDED` above; the file is small and blessing is rare.
function _writethrough(key, plain)
    merged = merge(canonical_reference(), Dict(key => plain))
    _CACHE[] = merged
    open(REFERENCE, "w") do io
        println(io,
                "# Blessed canonical results — regenerate with ECOSISTEM_BLESS=true, and read")
        println(io,
                "# test/canonical/README.md before committing a change to this file.")
        return TOML.print(io,
                          Dict(k => merged[k]
                               for k in sort(collect(keys(merged)))))
    end
    return nothing
end

"""
    writereference()

Write everything recorded this run to `reference.toml`. Call once, after all canonical tests.

 **Merges rather than replaces.** A run that executed only some of the canonical files would
otherwise silently delete the blessed values of the rest, turning a partial re-blessing into a
wholesale loss — the sort of damage that only shows up much later, as a test that stopped checking
anything.
"""
function writereference()
    blessing() || return nothing
    merged = merge(canonical_reference(), RECORDED)
    open(REFERENCE, "w") do io
        println(io,
                "# Blessed canonical results — regenerate with ECOSISTEM_BLESS=true, and read")
        println(io,
                "# test/canonical/README.md before committing a change to this file.")
        return TOML.print(io,
                          Dict(k => merged[k]
                               for k in sort(collect(keys(merged)))))
    end
    @info "blessed $(length(RECORDED)) canonical value(s) into $(REFERENCE)"
    return nothing
end

end
