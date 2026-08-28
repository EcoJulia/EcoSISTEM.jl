# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How to declare a niche axis: the interface hooks the macro fills in, and `@nicheaxis` itself.
# This file declares no types at all — the root `NicheAxis` is in `Ecology.jl` with the rest of the
# conceptual framework, and the axes themselves are in `NicheAxis.jl`. What is here is the machinery
# that declares them, which is why it belongs before them rather than with the code that merely
# operates on axes.
#
# The four `function … end` forward declarations are load-bearing: `@nicheaxis` emits its methods
# fully qualified (`EcoSISTEM.canonicalunit(…) = …`), and a qualified definition needs the name to
# exist already, or it is an `UndefVarError` on the name rather than a definition of it.
# ---------------------------------------------------------------------------
# Niche axes
# ---------------------------------------------------------------------------
# A `NicheAxis` names *what* a layer measures — an axis of the (Hutchinsonian) niche — in
# plain English, independent of its physical unit, and orthogonal to `Role` (which says how
# a layer is *used*). A layer and the trait matched to it share the same `NicheAxis`, so
# "same axis" is the matching rule rather than "same unit". Every axis — this package's own and
# anyone else's — is declared with [`@nicheaxis`](@ref), which is the only supported route.

using Unitful
using Unitful.DefaultSymbols

# Forward declarations, and they are load-bearing rather than stylistic. `@nicheaxis` emits its
# methods **fully qualified**, and a qualified definition needs the function to exist already —
# otherwise it is an `UndefVarError` on the name rather than a definition of it. The docstrings and
# the fallback methods stay below, with the rest of the interface.
function canonicalunit end
function bounds end
function supplytype end
function demandtype end

# ---------------------------------------------------------------------------
# `@nicheaxis` — the one way to declare an axis
# ---------------------------------------------------------------------------
# **Why a macro rather than several public functions.** An axis is defined by up to five separate
# methods — `canonicalunit` in two role forms, `supplytype`, `demandtype`, `bounds` — and the three
# resource ones have to agree. Declared by hand they can drift, and a missing method leaves "nobody
# has decided" indistinguishable from "deliberately not a condition". Emitting them together from one
# declaration makes disagreement unrepresentable and turns the four states into syntax.
#
# It also keeps the **public surface smaller**: one macro, rather than `canonicalunit`, `supplytype`,
# `demandtype` and `bounds` all becoming API. `bounds` in particular must not be public, since it is
# exported by both `DimensionalData` and `Rasters`, which this package depends on.
#
# **Compulsory in the sense that nothing else works, rather than by the language.** Julia cannot stop
# someone writing the abstract type by hand, but such a type fails loudly at first use:
# `canonicalunit(::Type{<:NicheAxis})` **errors**, `canonicalunit(::Type{NicheAxis})` answers
# `nothing` for the root alone, and every macro-declared axis emits its own method — a pure grouping
# node by **delegating to its parent** rather than by falling through, which is what would otherwise
# make a group and an undeclared type indistinguishable.

# The ways an axis can answer "what am I", as they appear in a declaration. `condition` and
# `resource` are deliberately parallel: each **names a role and gives that role's canonical unit**,
# so the two sides read alike and an axis states which of the four states it is in. (`unit =` was the
# earlier spelling and was ambiguous about *which* unit — an axis genuinely has two.)
#
# condition = U                            a condition, canonicalised to U
# resource = U, supply = S, demand = D     a resource; all three or none, so they cannot drift
# condition = nothing                      explicitly NOT a condition (a supply-only axis)
# reference                                neither — buildable and materialisable, never simulated
# (nothing at all)                         inherit from the group
#
# `bounds = (lo, hi)` is optional on any of them, and is refused without something to state it in.
# `supply = Supply{ThisAxis}` now restates the axis, and `supply`/`demand` both leave this list
# once **demands** are re-parameterised too — at that point `supplytype(A)` *is* `Supply{A}`,
# `demandtype(A)` is `Demand{A}`, and `resource = U` becomes the whole resource side. Kept for now
# because the demand side still keys on the rate, which is the asymmetry that step removes.
const _AXIS_KEYWORDS = (:condition, :resource, :supply, :demand, :bounds,
                        :densitywidth, :categorical)

# Split `Name <: Parent` into its two halves, erroring helpfully rather than with a `BoundsError`
# on a malformed head — this is the first thing a mistyped declaration hits.
# `ex` is deliberately unannotated although its body reads `.head`/`.args`: the `ex isa Expr`
# check inside is what produces that helpful message, and `::Expr` would replace it with a
# `MethodError` naming an internal function.
function _axishead(ex)
    (ex isa Expr && ex.head === :(<:) && length(ex.args) == 2) ||
        error("@nicheaxis expects `Name <: Parent`, got `$ex`")
    return ex.args[1], ex.args[2]
end

# **These live ABOVE the macro and the declarations because `@nicheaxis` CALLS them.** A
# grouping node's emitted method delegates to its parent, and `_checkaxisunit` reads both
# answers at declaration time, so the root must already answer by then. Moving them down
# again gives `MethodError: no method matching canonicalunit(::Type{NicheAxis})` while the
# first axis is being declared.
# **The root exactly answers `nothing`, and every other undeclared axis errors.** `NicheAxis` is the
# axis of a layer built without one — "I could be anything" — so `nothing` is the honest answer there
# and only there. An exact-`Type{NicheAxis}` method does **not** leak to subtypes, which is what keeps
# those two meanings apart: a single `<:NicheAxis` fallback would conflate them, so a hand-rolled
# `abstract type MyAxis <: NicheAxis` would answer "no axis was named" about an axis that **was**
# named, and then accept any unit for it. Erroring is what makes the extension point safe to
# document. Every macro-declared axis emits its own method, a pure grouping node by delegating to its
# parent, so the only types reaching this are hand-rolled ones.
"""
    canonicalunit(::Type{<:NicheAxis})
    canonicalunit(::Type{<:Role}, ::Type{<:NicheAxis})

The unit a layer on this axis is normalised to when materialised (e.g. temperature → `K`),
or `nothing` to leave values as-is. The 2-argument form is role-aware: a `Condition` (a
niche tolerance — a descriptive climatological normal) and a `Resource` (a literal
consumable pool, replenished over time) can legitimately want a different canonical unit for
the same axis — e.g. `Precipitation` stays a bare depth (`mm`) as a `Condition`, but is a
genuine volumetric flow (`L/day`) as a `Resource`. It defaults to the 1-argument form for
any role/axis combination without a specific override, so existing (implicitly
`Condition`-role) call sites are unaffected.
"""
canonicalunit(::Type{NicheAxis}) = nothing
function canonicalunit(A::Type{<:NicheAxis})
    return error("niche axis $(nameof(A)) declares no canonical unit. Declare it with " *
                 "`@nicheaxis($(nameof(A)) <: SomeParent, condition = <unit>)` — or, if it is a " *
                 "pure grouping node whose leaves each carry their own unit, with " *
                 "`@nicheaxis($(nameof(A)) <: SomeParent)`, which inherits its parent's answer. " *
                 "Declaring the abstract type by hand does not register it: the macro is what emits " *
                 "the unit, bounds and resource methods the rest of the package asks for.")
end
canonicalunit(::Type{<:Role}, A::Type{<:NicheAxis}) = canonicalunit(A)

"""
    densitywidth(::Type{<:NicheAxis})

The **fixed physical width** a continuous suitability density is multiplied by to make it a
dimensionless weight, or `nothing` for an axis that declares none (no scaling is applied).

**This is what stops an axis's `canonicalunit` being a model parameter in disguise.** A
probability density carries `1/x`, so a *stripped* `pdf` is really "probability mass in a window one
canonical unit wide" — re-declaring an axis from `mm/d` to `cm/d` would widen that window tenfold,
multiply every suitability by ten and move every equilibrium. Multiplying by a width that is a
**fixed physical quantity** cancels exactly that, because the width shrinks in the new frame by the
same factor the density grows.

**It must be written as a literal, never derived from [`canonicalunit`](@ref).** A width defined
as "one canonical unit" changes when the canonical unit changes, which is the bug it exists to fix.
Each declaration therefore states its quantity outright, and today's values are exactly one of each
axis's *current* canonical unit — so behaviour is unchanged and `test/canonical/reference.toml` does
not move.

It preserves the required specialist advantage, because the width is the **same for every
species**: a narrower niche still peaks higher. That is what separates it from peak-normalising by
σ, which would leave the specialist strictly dominated.
"""
densitywidth(::Type{NicheAxis}) = nothing
densitywidth(::Type{<:NicheAxis}) = nothing

"""
    iscategorical(axis::Type{<:NicheAxis})

Return whether `axis`' values are **class labels** rather than measurements.

This is a property of the axis, and nothing else needs to declare it. An axis holding Koppen
climate classes or land-cover classes holds codes whose arithmetic mean is meaningless, so anything
resampling such a layer must take the nearest class rather than interpolate. Every other axis holds
numbers that may be averaged.

A day-count is deliberately not modelled as a third case. The catalogue distinguishes one from a
continuous measurement, but nothing in the package acts on the distinction: every consumer asks only
whether a layer is categorical, because that is the only answer that changes what resampling may do.
A count is an ordinary number that may be averaged, exactly as a temperature is.

An axis whose values are class labels says so in its [`@nicheaxis`](@ref) declaration, with
`categorical = true`; every other axis declares nothing and inherits from its parent, ending at
`NicheAxis` itself, which is not categorical. This is the same declared-or-delegated shape as
[`canonicalunit`](@ref) and `densitywidth`, and the macro is the only route — writing a method by
hand covers the exact type named and not its subtypes, so a hand-declared group would silently fail
to pass the property to its own leaves.

The fallback is `false` rather than an error, which is where this differs from
[`canonicalunit`](@ref). Very few axes are categorical, so `false` is a real majority answer that an
axis can be assumed into; nearly every axis has a *different* canonical unit, so there is no answer
to fall back on and not declaring one has to be refused. The same question decides how any future
axis property should behave when undeclared: ask whether a majority answer exists.
"""
iscategorical(::Type{<:NicheAxis}) = false

"""
    @nicheaxis Name <: Parent  [condition = U] [resource = U supply = S demand = D] [reference]
                               [bounds = (lo, hi)] [densitywidth = W] [categorical = true]

Declare a niche axis: its type, and what it means. This is the supported way to add one — the
underlying `canonicalunit`/`supplytype`/`demandtype`/`bounds`/`densitywidth`/`iscategorical` methods
are internal, and are emitted here so that they cannot disagree.

**Every axis is an `abstract type`** — an axis is only dispatched on, never constructed — so there
is no group/leaf distinction to declare. An axis that declares nothing inherits its parent, which is
the common case and quite different from `reference`, which *states* that it has no canonical unit.

`condition` and `resource` are the two roles, and each **names a role and gives that role's
canonical unit** — so an axis is a condition, a resource, both, or neither, and says which.

| what you write | what it means |
| --- | --- |
| `condition = U` | a **condition**: layers are canonicalised to `U`, and a tolerance is built in it |
| `resource = U supply = S demand = D` | a **resource**: species consume it. All three, or none |
| both | a condition *and* a resource, as [`Precipitation`](@ref) is |
| `condition = nothing` | *not* a condition — use with `resource` for a supply-only axis, as `CarbonFlux` is |
| `reference` | **neither** — buildable, materialisable and composable, but never simulated |
| `categorical = true` | values are **class labels**, so a layer on this axis is resampled by nearest class rather than interpolated |
| nothing at all | inherit the group's declarations |

Omitting `condition` and writing `condition = nothing` are **different**: the first inherits
whatever the group declares, the second states that this axis has no condition reading and overrides
the group. Declaring a `resource` implies neither — an axis can be both, and `SolarRadiation` is.

`bounds = (lo, hi)` gives the physical range in the axis's own canonical unit, either end `nothing`.

`categorical` is independent of the roles, and combines with `reference`: a layer of class codes
carried only as ground truth is still put on a grid, and still must not be interpolated between its
classes. Declaring it on one axis covers that axis' descendants, which is why `TypologyAxis` is the
only declaration in this package — `LandCoverTypology` and `ClimateTypology` inherit it.

```julia
@nicheaxis(TemperatureAxis <: NicheAxis, condition = K, bounds = (0.0K, nothing), densitywidth=1.0K)
@nicheaxis(Temperature <: TemperatureAxis)                     # inherits both
@nicheaxis(CumulativeHeat <: TemperatureAxis, condition = K * Unitful.d, densitywidth=1.0K * Unitful.d)
```

Both call forms work — `@nicheaxis Name <: Parent condition = U` parses identically. The
parenthesised, comma-separated one is used throughout this package because JuliaFormatter keeps it
readable, while it reflows the space-separated form into wrapped fragments.
"""
macro nicheaxis(args...)
    isempty(args) && error("@nicheaxis needs at least `Name <: Parent`")
    name, parent = _axishead(args[1])
    opts = Dict{Symbol, Any}()
    isreference = false
    for a in args[2:end]
        if a === :reference
            isreference = true
        elseif a isa Expr && a.head === :(=) && a.args[1] isa Symbol
            a.args[1] in _AXIS_KEYWORDS ||
                error("@nicheaxis $name: unknown option `$(a.args[1])`; expected one of " *
                      "$(join(_AXIS_KEYWORDS, ", ")) or the bare word `reference`")
            haskey(opts, a.args[1]) &&
                error("@nicheaxis $name: `$(a.args[1])` given twice")
            opts[a.args[1]] = a.args[2]
        else
            error("@nicheaxis $name: expected `option = value` or `reference`, got `$a`")
        end
    end

    # The three resource declarations are all-or-nothing **by construction**, which is what keeps
    # `supplytype`, `demandtype` and `canonicalunit(Resource, ·)` in step without anyone having to.
    res = intersect(keys(opts), (:resource, :supply, :demand))
    isempty(res) || length(res) == 3 ||
        error("@nicheaxis $name: a resource needs `resource`, `supply` and `demand` together " *
              "(got $(join(sort(collect(res)), ", "))) — declaring some of them is how they drift apart")
    isresource = length(res) == 3

    # `reference` is a positive statement of "neither", so it cannot be combined with either role.
    isreference && (haskey(opts, :condition) || isresource) &&
        error("@nicheaxis $name: `reference` means the axis is neither a condition nor a resource, " *
              "so it cannot be given a `condition` or a `resource`")
    # A bound is stated *in* the canonical unit, so there must be one to state it in — either here or,
    # for a leaf, inherited from the group, which we cannot see from inside the macro.
    haskey(opts, :bounds) && isreference &&
        error("@nicheaxis $name: `reference` declares no unit, so there is nothing for `bounds` " *
              "to be expressed in")

    # **Every method is emitted fully qualified** (`$M.canonicalunit(…)`), and this is load-bearing
    # rather than tidy: an unqualified definition in a *caller's* module defines that module's own
    # `canonicalunit`, and Julia refuses to extend a `using`-imported function anyway. So without this
    # a user would have to `import EcoSISTEM: canonicalunit, supplytype, demandtype, bounds` before
    # the macro worked — importing the exact internals the macro exists to keep private. Caught by the
    # first test of the macro from outside the package, which is why that test is worth having.
    M = @__MODULE__
    out = Expr(:block)
    # **Always an `abstract type`, never a struct.** An axis is only ever dispatched on, never
    # constructed, so there is nothing for a concrete type to buy — and making every axis abstract
    # means one that is a leaf today can gain children tomorrow without its own declaration changing,
    # which a `struct` would make a breaking change.
    # `Base.@__doc__` is what lets a docstring sit above the macro call. Without it the docsystem
    # sees a `begin` block of several statements and refuses ("cannot document the following
    # expression") — it needs to be told which one the docstring belongs to.
    push!(out.args, :(Base.@__doc__ abstract type $name <: $parent end))
    # `condition = nothing` is a *statement*, not an omission — "this axis has no condition reading" —
    # and is emitted as a method so it beats anything the group declares. That distinction is what
    # makes `nothing` unambiguous later: omitting `condition` inherits, writing `condition = nothing` refuses.
    #
    # **Found by dogfooding, and it would have been silent**: an earlier version emitted
    # `canonicalunit = nothing` for *every* resource. `SolarRadiation` is a resource whose condition
    # unit is inherited from `SolarRadiationAxis` (the group exists precisely to reconcile WorldClim's
    # `kJ·m⁻²` with CHELSA's `MJ·m⁻²`), so that would have overridden the group with `nothing` and
    # unreconciled the two sources — a wrong *number*, not an error. Being a resource says nothing
    # about whether the axis is also a condition; only `condition = nothing` says that.
    if haskey(opts, :condition)
        push!(out.args,
              :($M.canonicalunit(::Type{<:$name}) = $(opts[:condition])))
    elseif isreference
        push!(out.args, :($M.canonicalunit(::Type{<:$name}) = nothing))
    else
        # **Inheritance is emitted explicitly rather than left to the method table**, so that every
        # macro-declared axis owns a `canonicalunit` method. That is what lets the
        # `canonicalunit(::Type{<:NicheAxis})` fallback **error** below: the only types left reaching
        # it are hand-rolled ones that never ran this macro, which is exactly who should be told.
        # Delegation rather than a frozen copy, so that a group and its leaves cannot drift; the
        # chain is at most four deep, on a setup-only path.
        push!(out.args,
              :($M.canonicalunit(::Type{<:$name}) = $M.canonicalunit($parent)))
    end
    # **An axis may not contradict an ancestor's declared unit.** Where a group and its leaves
    # disagree, they are typically not the same quantity in different units at all — `K` against
    # `K·d`, degree-days, against a dimensionless ratio or count — so there is **no meaningful
    # conversion**, and promotion, which rewraps without converting, would silently relabel data onto
    # an axis it is not in. Refusing the *declaration* makes that unrepresentable, rather than
    # something a runtime guard has to catch.
    #
    # A parent that declares **no** unit is the safe case and stays legal: it made no claim to
    # contradict, which is what every pure grouping node does.
    push!(out.args, :($M._checkaxisunit($name, $parent)))
    # The density width follows `canonicalunit`'s shape exactly: declared here, or delegated to the
    # parent so a grouping node's leaves inherit it. It is **not** derived from the canonical unit —
    # see `densitywidth`'s docstring for why that would reintroduce the very bug it fixes.
    if haskey(opts, :densitywidth)
        push!(out.args,
              :($M.densitywidth(::Type{<:$name}) = $(opts[:densitywidth])))
    else
        push!(out.args,
              :($M.densitywidth(::Type{<:$name}) = $M.densitywidth($parent)))
    end
    # Whether the axis' values are class labels rather than measurements, following the same
    # declared-or-delegated shape as `densitywidth`. Only `TypologyAxis` declares `true`; everything
    # else reaches `iscategorical(::Type{<:NicheAxis}) = false` up the chain of parents.
    #
    # The delegation is what makes a single declaration cover a whole branch: `LandCoverTypology`
    # declares nothing and emits `iscategorical(::Type{<:LandCoverTypology}) = iscategorical(TypologyAxis)`,
    # which is `true`. Emitting nothing instead would leave the branch to the method table, which
    # works here but stops working the moment a default other than the parent's is wanted.
    #
    # This is deliberately independent of `reference`. An axis that is neither a condition nor a
    # resource may still hold class codes, and it still needs to say so: a reference layer is
    # materialised onto a grid like any other, and `_resamplemethod` reads this answer to choose
    # `:mode` over `:bilinear`. Interpolating between class codes is meaningless whether or not a
    # species responds to them, so unlike `bounds` — which is stated in a unit a reference axis has
    # shed — there is nothing here for `reference` to take away.
    if haskey(opts, :categorical)
        push!(out.args,
              :($M.iscategorical(::Type{<:$name}) = $(opts[:categorical])))
    else
        push!(out.args,
              :($M.iscategorical(::Type{<:$name}) = $M.iscategorical($parent)))
    end
    if isresource
        push!(out.args,
              :($M.canonicalunit(::Type{$M.Resource},
                                 ::Type{<:$name}) = $(opts[:resource])))
        push!(out.args, :($M.supplytype(::Type{<:$name}) = $(opts[:supply])))
        push!(out.args, :($M.demandtype(::Type{<:$name}) = $(opts[:demand])))
    end
    # A `reference` axis must **shed** its group's bounds, not merely be forbidden its own. It
    # already overrides the group's `canonicalunit` with `nothing`; leaving `bounds` inherited would
    # give it a floor stated in a unit it no longer has — and `_enforcebounds!` runs on every
    # `Condition` layer, so a reference layer under a bounded group would be compared against
    # `0.0K` with nothing to convert. Refusing a written `bounds` (above) does not cover this,
    # because the offending bound is not written here at all.
    if isreference
        push!(out.args, :($M.bounds(::Type{<:$name}) = (nothing, nothing)))
    elseif haskey(opts, :bounds)
        push!(out.args, :($M.bounds(::Type{<:$name}) = $(opts[:bounds])))
    end
    push!(out.args, name)
    return esc(out)
end

# The guard `@nicheaxis` emits after declaring an axis's unit — see the comment at its call site for
# why a contradiction is refused rather than converted. Reads both answers rather than the macro's
# arguments, so it catches an inherited contradiction as well as a written one.
function _checkaxisunit(axis::Type, parent::Type)
    inherited = canonicalunit(parent)
    isnothing(inherited) && return nothing
    own = canonicalunit(axis)
    own === inherited && return nothing
    # **Shedding a claim is safe; making a different one is not**, and that distinction is the whole
    # rule. `reference` — and `condition = nothing` on a supply-only axis — says "I assert no
    # condition unit", so a layer promoted onto it carries its data with nothing false claimed about
    # it. Declaring `K·day` beneath a `K` group asserts something false instead, and no conversion
    # could make it true.
    #
    # A blanket "any difference is a contradiction" rule would break the `reference` case, since such
    # an axis must shed its group's unit *and* its bounds; `test_nicheaxis_macro.jl` guards that with
    # a fixture of its own.
    isnothing(own) && return nothing
    return error("@nicheaxis $(nameof(axis)): its ancestor $(nameof(parent)) declares the canonical " *
                 "unit `$inherited`, and this axis declares `$own` — an axis cannot contradict the " *
                 "unit it inherits, because there is no meaningful conversion between them and " *
                 "promoting a layer from one to the other would relabel its data rather than " *
                 "convert it. Give $(nameof(parent)) no `condition` (making it a pure grouping " *
                 "node, as $(nameof(parent)) may well be) and state the unit on each leaf instead.")
end

# --- Which axis is this on? ---------------------------------------------------
# The typed methods are with their own families (`Layer.jl`, `SpeciesRequirement.jl`, `NicheFit.jl`,
# `Demand.jl`); only the fallback is here, with the rest of the axis interface.

# The niche axis of a threaded trait / layer (`NicheAxis` when none was declared).
# **The fallback is deliberate and load-bearing — do not annotate it.** A collection constructor asks
# `axisof` of every member *before* the guard that rejects a non-member, so this is what lets
# `SpeciesRequirementCollection((tolerance, 1))` reach that guard and fail with a message naming the
# mistake. Annotating it turns that into `MethodError: no method matching axisof(::Int64)`, which
# names an internal function instead of explaining anything, and a test asserts the `ErrorException`.
#
# `NicheAxis` is also the honest answer for a layer built without a declared axis, so this doubles as
# that — but the error path is why it must stay untyped.
axisof(_) = NicheAxis
