# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations — main `EcoSISTEM` module
#
# Every deprecated public trait/nichefit API is collected here (sorted into
# sections by context) and included late in `EcoSISTEM.jl`, after all the types
# it shims. Each shim warns — so downstream code gets a migration message rather
# than a silent `MethodError` — and forwards to the current API. Mirrored by
# `test/test_deprecations.jl`. The `ClimatePref` submodule keeps its own
# deprecations in `src/ClimatePref/deprecations.jl` (a single file cannot span
# two modules).
# ===========================================================================

# ---------------------------------------------------------------------------
# Trait line: `GaussTrait` → `NicheTolerance(A, Normal, …)`
#
# A Gaussian preference is just the `Normal` case of a `NicheTolerance`. The redesign takes
# a trait's unit from its niche axis (not its data), so a *unitful* preference
# must name an axis; only dimensionless (bare) data may stay axis-less.
# ---------------------------------------------------------------------------
"""
    GaussTrait(mean, sd)
    GaussTrait(::Type{A}, mean, sd)

!!! warning "Deprecated"
    `GaussTrait` is deprecated and will be removed. A Gaussian preference is just the `Normal` case of a
    [`NicheTolerance`](@ref), so use one directly — `NicheTolerance(A, Normal, mean, sd)` (one unitful vector per parameter,
    `support` imputed). This shim forwards to that. `GaussTrait(::Type{A}, mean, sd)` names the niche axis
    `A` explicitly; the axis-less `GaussTrait(mean, sd)` form takes dimensionless (bare) data as
    `NicheAxis`, and — for back-compatibility only — *infers* the axis from a unitful vector's unit
    (temperature → `Temperature`, precipitation rate → `Precipitation`) with an extra warning, because a
    [`NicheTolerance`](@ref)'s unit must come from its axis, not its data.
"""
function GaussTrait(::Type{A}, mean::AbstractVector,
                    sd::AbstractVector) where {A <: NicheAxis}
    Base.depwarn("`GaussTrait` is deprecated; use `NicheTolerance($A, Normal, mean, sd)` instead.",
                 :GaussTrait)
    return NicheTolerance(A, Normal, mean, sd)
end

# Axis-less bare (dimensionless) form: a root-axis (`NicheAxis`) NicheTolerance, as before.
function GaussTrait(mean::AbstractVector{<:Real}, sd::AbstractVector)
    return GaussTrait(NicheAxis, mean, sd)
end

# Map a bare unit back to the niche axis it canonically measures. This is a *deprecated* fallback for the
# axis-less unitful `GaussTrait(mean, sd)` form only: inferring a layer's meaning from its unit is exactly
# the anti-pattern the axis redesign removes, so it is intentionally minimal (temperature/precipitation) and
# errors for anything else, telling the caller to name the axis. The two units are read off `canonicalunit`
# rather than written out, so a change there (e.g. `Precipitation` becoming the rate `mm/d`) cannot leave
# this inference silently mapping the old unit.
function _axis_from_unit(u)
    dimension(u) == dimension(canonicalunit(Temperature)) &&
        return Temperature
    dimension(u) == dimension(canonicalunit(Precipitation)) &&
        return Precipitation
    return error("Cannot infer a niche axis from unit `$u`; name the axis explicitly, e.g. " *
                 "`GaussTrait(Temperature, mean, sd)` or `NicheTolerance(axis, Normal, mean, sd)`.")
end

# Axis-less *unitful* form (doubly deprecated): kept working via unit→axis inference so downstream code gets
# a clear migration message instead of a `MethodError`, but the meaning-from-unit inference will be dropped.
function GaussTrait(mean::AbstractVector{<:Unitful.AbstractQuantity},
                    sd::AbstractVector)
    u = unit(eltype(mean))
    axis = _axis_from_unit(u)
    Base.depwarn("`GaussTrait(mean, sd)` with unitful data infers the niche axis from its unit " *
                 "($u → $axis); this meaning-from-unit fallback is deprecated. Name the axis: " *
                 "`GaussTrait($axis, mean, sd)` or `NicheTolerance($axis, Normal, mean, sd)`.",
                 :GaussTrait)
    return NicheTolerance(axis, Normal, mean, sd)
end

# ---------------------------------------------------------------------------
# Trait line: `Gauss` / `Trapeze` / `Unif` relationships → `NicheSuitability`
#
# A Gaussian/trapezoidal/uniform preference is just the `Normal`/`Trapezoid`/`Uniform` case of a `NicheTolerance`, so
# all three relationships are `NicheSuitability` now. Each is kept as a distinct type that warns on construction but
# shares `NicheSuitability`'s 2-arg density functor, so old hand-built ecosystems still build and evaluate.
# ---------------------------------------------------------------------------
"""
    Gauss{A, V} <: AbstractNicheFit{A, V}

!!! warning "Deprecated"
    `Gauss` is deprecated and will be removed; use [`NicheSuitability`](@ref) instead (a Gaussian preference is the
    `Normal` case of a [`NicheTolerance`](@ref)). This shim shares `NicheSuitability`'s 2-argument density functor and also
    retains the legacy 3-argument `(current, opt, sd)` Gaussian call for back-compatibility.
"""
struct Gauss{A, V} <: AbstractNicheFit{A, V}
    function Gauss{A, V}() where {A, V}
        Base.depwarn("`Gauss` is deprecated; use `NicheSuitability{A, V}()` instead.",
                     :Gauss)
        return new{A, V}()
    end
end

"""
    Trapeze{A, V} <: AbstractNicheFit{A, V}

!!! warning "Deprecated"
    `Trapeze` is deprecated and will be removed; use [`NicheSuitability`](@ref) instead (it pairs with a `Trapezoid`
    [`NicheTolerance`](@ref)).
"""
struct Trapeze{A, V} <: AbstractNicheFit{A, V}
    function Trapeze{A, V}() where {A, V}
        Base.depwarn("`Trapeze` is deprecated; use `NicheSuitability{A, V}()` instead.",
                     :Trapeze)
        return new{A, V}()
    end
end

"""
    Unif{A, V} <: AbstractNicheFit{A, V}

!!! warning "Deprecated"
    `Unif` is deprecated and will be removed; use [`NicheSuitability`](@ref) instead (it pairs with a `Uniform`
    [`NicheTolerance`](@ref)).
"""
struct Unif{A, V} <: AbstractNicheFit{A, V}
    function Unif{A, V}() where {A, V}
        Base.depwarn("`Unif` is deprecated; use `NicheSuitability{A, V}()` instead.",
                     :Unif)
        return new{A, V}()
    end
end

# The deprecated shims share `NicheSuitability`'s density functor + trait interface.
for Old in (:Gauss, :Trapeze, :Unif)
    @eval begin
        function (::$Old{A, V})(dist::ContinuousUnivariateDistribution,
                                current::V) where {A, V}
            return pdf(dist, ustrip(current))
        end
        iscontinuous(::$Old) = true
        Base.eltype(::$Old{A, V}) where {A, V} = V
    end
end

# Legacy 3-argument Gaussian functor `Gauss{A, V}()(current, opt, sd)`, restored for back-compatibility: the
# hand-written Gaussian density (a dimensionless `Quantity`) that `Gauss` evaluated before it became a
# `NicheSuitability` shim. New code should build a `Normal` `NicheTolerance` and use `NicheSuitability`'s 2-argument functor instead.
function (::Gauss{A, V})(current::V, opt::V, var::V) where {A, V}
    pref = 1.0 / sqrt(2 * π * var^2) *
           exp(-abs(current - opt)^2 / (2 * var^2))
    return pref * unit(current)
end

# ---------------------------------------------------------------------------
# Arity-numbered collection aliases (v0.5.0): plumbing for the v0.4.0 shims below, not API
#
# Each family now has **one** collection type over a `Tuple` or `NamedTuple` — `LayerCollection`,
# `SpeciesRequirementCollection`, `CombiningFit` — so an arity in a name says nothing the
# backing's own type does not. The numbered names are *not* released: they were themselves new in
# v0.5.0 (v0.4.0 had `HabitatCollection2/3`, `BudgetCollection2`, `TraitCollection2/3`,
# `ReqCollection2`, `multiplicativeTR2/3`, `additiveTR2/3`), so by the usual rule — a name new on the
# branch is changed outright, no shim — they owe nothing and are simply no longer exported.
#
# They survive here only because the v0.4.0 shims need somewhere to point: `Base.@deprecate_binding`
# makes `const old = new`, which cannot retarget `HabitatCollection2{H1, H2}` onto
# `LayerCollection{R, C}` (the first parameter means a different thing), and dropping the parameters
# would break v0.4.0 code that *dispatches* on `::HabitatCollection2{H1, H2}`. So each numbered name
# stays as an alias of the right shape, and goes when its v0.4.0 shim does.
#
# Deliberately **no `depwarn` of their own**: the only ways to reach them are through a v0.4.0 shim,
# which warns already, or by writing `EcoSISTEM.RegimeCollection2` explicitly, which is reaching into
# internals. Adding one would give a caller of the oldest name two warnings naming each other —
# exactly what the `*AE` and gradient sections below take care to avoid.
# ---------------------------------------------------------------------------

"""    RegimeCollection2{M1, M2} — two positional regimes over one grid. Deprecated plumbing: use [`LayerCollection`](@ref). """
const RegimeCollection2{M1, M2} = LayerCollection{Condition, A,
                                                  NamedTuple{N, Tuple{M1, M2}}} where {A,
                                                                                       N}
"""    RegimeCollection3{M1, M2, M3} — three positional regimes over one grid. Deprecated plumbing: use [`LayerCollection`](@ref). """
const RegimeCollection3{M1, M2, M3} = LayerCollection{Condition, A,
                                                      NamedTuple{N,
                                                                 Tuple{M1, M2,
                                                                       M3}}} where {A,
                                                                                    N}
"""    SupplyCollection2{M1, M2} — two positional supplies over one grid. Deprecated plumbing: use [`LayerCollection`](@ref). """
const SupplyCollection2{M1, M2} = LayerCollection{Resource, A,
                                                  NamedTuple{N, Tuple{M1, M2}}} where {A,
                                                                                       N}
"""    ToleranceCollection2{M1, M2} — two positional tolerances. Deprecated plumbing: use [`SpeciesRequirementCollection`](@ref). """
const ToleranceCollection2{M1, M2} = SpeciesRequirementCollection{Condition, A,
                                                                  NamedTuple{N,
                                                                             Tuple{M1,
                                                                                   M2}}} where {A,
                                                                                                N}
"""    ToleranceCollection3{M1, M2, M3} — three positional tolerances. Deprecated plumbing: use [`SpeciesRequirementCollection`](@ref). """
const ToleranceCollection3{M1, M2, M3} = SpeciesRequirementCollection{Condition,
                                                                      A,
                                                                      NamedTuple{N,
                                                                                 Tuple{M1,
                                                                                       M2,
                                                                                       M3}}} where {A,
                                                                                                    N}
"""    DemandCollection2{M1, M2} — two positional demands. Deprecated plumbing: use [`SpeciesRequirementCollection`](@ref). """
const DemandCollection2{M1, M2} = SpeciesRequirementCollection{Resource, A,
                                                               NamedTuple{N,
                                                                          Tuple{M1,
                                                                                M2}}} where {A,
                                                                                             N}
"""    MultiplicativeFit2{M1, M2} — two positional nichefits, multiplied. Deprecated plumbing: use [`MultiplicativeFit`](@ref). """
const MultiplicativeFit2{M1, M2} = MultiplicativeFit{A,
                                                     NamedTuple{N,
                                                                Tuple{M1,
                                                                      M2}}} where {A,
                                                                                   N}
"""    MultiplicativeFit3{M1, M2, M3} — three positional nichefits, multiplied. Deprecated plumbing: use [`MultiplicativeFit`](@ref). """
const MultiplicativeFit3{M1, M2, M3} = MultiplicativeFit{A,
                                                         NamedTuple{N,
                                                                    Tuple{M1,
                                                                          M2,
                                                                          M3}}} where {A,
                                                                                       N}
"""    AdditiveFit2{M1, M2} — two positional nichefits, added. Deprecated plumbing: use [`AdditiveFit`](@ref). """
const AdditiveFit2{M1, M2} = AdditiveFit{A,
                                         NamedTuple{N, Tuple{M1, M2}}} where {A,
                                                                              N}
"""    AdditiveFit3{M1, M2, M3} — three positional nichefits, added. Deprecated plumbing: use [`AdditiveFit`](@ref). """
const AdditiveFit3{M1, M2, M3} = AdditiveFit{A,
                                             NamedTuple{N,
                                                        Tuple{M1, M2,
                                                              M3}}} where {A,
                                                                           N}

# **Each alias pins the arity and the member types, and leaves the axis structure `A` and the names
# `N` free** — which is why every one ends `} where {A, N}`. It reads oddly and the shape is
# load-bearing: a collection's backing is a `NamedTuple` whose names are derived from its members'
# **axes**, so an alias cannot state them, and the axis structure is not something a caller of the
# released API ever knew about.
#
# They discriminate correctly: a two-member fit `isa MultiplicativeFit2{M1, M2}` and is **not**
# `isa MultiplicativeFit2{M2, M1}`, and a three-member one is **not** `isa MultiplicativeFit2`.
#
# **Put a member tuple in the axis slot and an alias silently stops matching what its own constructor
# builds** — which a test comparing only the *bindings*, `additiveTR2 === AdditiveFit2`, cannot catch.
# `test_deprecations.jl` has a round-trip testset that uses them as **types**; keep it.
#
# Each alias fixes the arity but leaves the member parameters free, so Julia generates no
# constructor for it; these are what a v0.4.0 `HabitatCollection2(h1, h2)` call actually reaches.
RegimeCollection2(l1, l2) = LayerCollection((l1, l2))
RegimeCollection3(l1, l2, l3) = LayerCollection((l1, l2, l3))
SupplyCollection2(b1, b2) = LayerCollection((b1, b2))
ToleranceCollection2(t1, t2) = SpeciesRequirementCollection((t1, t2))
ToleranceCollection3(t1, t2, t3) = SpeciesRequirementCollection((t1, t2, t3))
DemandCollection2(r1, r2) = SpeciesRequirementCollection((r1, r2))
MultiplicativeFit2(f1, f2) = CombiningFit(prod, (f1, f2))
MultiplicativeFit3(f1, f2, f3) = CombiningFit(prod, (f1, f2, f3))
AdditiveFit2(f1, f2) = CombiningFit(sum, (f1, f2))
AdditiveFit3(f1, f2, f3) = CombiningFit(sum, (f1, f2, f3))

# ---------------------------------------------------------------------------
# Resource line: `Resource` → `Supply` (v0.4.0 rename; the environment's resource layer)
# ---------------------------------------------------------------------------
# **`SimpleBudget` is REMOVED, not shimmed** — the free/dimensionless supply it named no longer
# exists, so there is nothing to re-point it at. It was exported in v0.4.0, which makes this a
# breaking removal with a NEWS entry rather than a cleanup. Recorded here, where the shim used to
# be, so the absence reads as a decision rather than an oversight.
# These re-point at the **parametric** alias, which was checked before being adopted: a
# `@deprecate_binding` cannot retarget `HabitatCollection2{H1, H2}` onto `LayerCollection{R, C}`
# — the only reason the arity-numbered collection names still exist — but `Supply{A}` is a different
# shape and does work, so no plumbing alias has to survive here.
# --- The v0.4.0 coordinate-less supply -------------------------------------
#
# **These two exist only for the released `SolarBudget(matrix)` path, and nothing else may use them.**
# A released budget *was* a bare `Matrix`, so there are no coordinates for a cell size to be derived
# from and none to attach — which is why such a layer can only report a placeholder size, and why
# nothing built today may take this route.
#
# **Kept here rather than deleted**, because `SolarBudget` was an
# *exported* v0.4.0 name and `@deprecate_binding SolarBudget Supply{SolarRadiation}` resolves a
# released `SolarBudget(mat)` call straight onto the `Matrix` method: removing it would drop released
# public behaviour with no shim. Here it keeps working, and keeps warning.
#
# **A layer built this way still cannot answer `getcellsizes` from its dims** — it has none. That
# is consistent with the rest of this file, which deliberately reproduces pre-projection v0.4.0
# behaviour (the same reason the `*habitat` builders keep `_cellsize`'s metric conversion).
_noLookupYX(M::AbstractMatrix) = DimArray(M, (Y(NoLookup()), X(NoLookup())))

# **Reaches the INNER constructor directly, bypassing `_checkhascoords`** — deliberately, and it is
# the only place that may. The outer `Supply{A}(::DimArray)` now refuses `NoLookup` dims; this path
# has none to give, because a v0.4.0 budget was a bare `Matrix`. So `deprecations.jl` is the sole
# owner of coordinate-less layers, alongside `_noLookupYX` above.
# **The ONE route to a coordinate-less layer, and every deprecated path goes through it.** The
# public constructors refuse `NoLookup` dims (`_checkhascoords`); these reach the *inner* constructor
# directly, which is legitimate here and nowhere else — a v0.4.0 budget was a bare `Array` and had no
# coordinates to give.
# **Written as one helper rather than repeated at each site**, because it was first patched into the
# `Matrix` constructor alone and the **3-D** path (`SolarTimeBudget(stack, time)` → a slice of a
# `NoLookup` stack) still went through the public one and was refused. Caught by
# `test_deprecations.jl` — 16 new errors, all in this file.
# **The REGIME half of this pair is gone.** `_legacyregime`/`_legacycategorical` built a v0.4.0
# layer around a raster's own cells with a *stated* cell size — `_cellsize`'s implicit equal-area
# projection — and their only callers were the four data-backed builders, which have been removed
# because that is precisely what no longer has an equivalent. `_legacysupply` survives because the
# `SolarTimeBudget`/`WaterTimeBudget` family still needs it: a v0.4.0 budget was a bare `Matrix` with
# no coordinates at all, which is a different problem from a raster whose coordinates are geographic.
function _legacysupply(::Type{A}, mat::AbstractArray) where {A <: NicheAxis}
    canon = EcoSISTEM._canonicalresource(mat isa
                                         DimensionalData.AbstractDimArray ?
                                         mat : _noLookupYX(mat), A)
    return ContinuousLayer{Resource, A, eltype(canon), typeof(canon),
                           typeof(_SUPPLY_SIZE)}(canon, _SUPPLY_SIZE,
                                                 NoLayerChange())
end

function Supply{A}(mat::Matrix) where {A <: NicheAxis}
    return _legacysupply(A, mat)
end

Base.@deprecate_binding SolarBudget Supply{SolarRadiation}
Base.@deprecate_binding SolarTimeBudget Supply{SolarRadiation}
Base.@deprecate_binding WaterBudget Supply{Precipitation}
Base.@deprecate_binding WaterTimeBudget Supply{Precipitation}
Base.@deprecate_binding BudgetCollection2 SupplyCollection2
# This forwards to today's `getsupply`, which returns the **layer** rather than v0.4.0's bare
# value matrix. A shim redirects a name, and the name it redirects to is the current one — pinning
# the old return type here would make the shim a second implementation, and would put exactly the
# naked array across the public boundary that `getsupply` was changed to stop handing out.
@deprecate getbudget(eco) getsupply(eco)

# Three names that each described the wrong quantity. `getsize` gave the grid's total *area*;
# `getgridsize` gave one *cell's* side length, not the grid's anything; and the grid's actual size in
# cells was `getdimension`, which was private.
#
# `getgridsize` is redirected rather than reused, which is why the replacement is called
# `getgridshape`. Reusing the name would have left it live with a different meaning, so a v0.4.0 call
# would silently get a `(ny, nx)` tuple where it expected a length -- a wrong answer with no warning,
# and the one thing a shim cannot express. Renaming makes the redirect ordinary.
#
# `first(getcellsizes(eco).y)` is exactly what the old function returned, for every ecosystem that
# can exist: `_checksimulatable` refuses a geographic or non-square grid before one is built, so a
# cell's `y` extent is uniform across the grid and equal to its `x`.
@deprecate getsize(eco) getgridarea(eco)
@deprecate getgridsize(eco) first(getcellsizes(eco).y)

# ---------------------------------------------------------------------------
# Resource line: `Requirement` → `Demand` (v0.4.0 rename; the species' resource need)
# ---------------------------------------------------------------------------
# **`SimpleRequirement` is REMOVED, not shimmed** — the same shape as `SimpleBudget` above, and
# for the same reason: the free/dimensionless demand it named no longer exists, so there is nothing
# to re-point it at. It was exported in v0.4.0, which makes this a breaking removal with a NEWS
# entry. Recorded here, where the shim used to be, so the absence reads as a decision.
# **`NoBoundary` said the opposite of what it meant**: there are boundaries, on all four sides —
# it named the *absence* of wrapping, not the absence of edges. `Island` says it correctly, and the
# grid's topology is now stated one axis at a time (`EdgeTopology{Bounded, Bounded}`).
# Owed a shim: `NoBoundary` was **exported in v0.4.0**. `Torus` and `Cylinder` are *not*
# deprecated — they survive as aliases of the same parametric type, so nothing that used them changes.
Base.@deprecate_binding NoBoundary Island

Base.@deprecate_binding SolarRequirement Demand{SolarRadiation}
Base.@deprecate_binding WaterRequirement Demand{Precipitation}
Base.@deprecate_binding ReqCollection2 DemandCollection2

# ---------------------------------------------------------------------------
# Condition line: `Condition`(-role layer) → `Regime` (v0.4.0 rename; the environment's condition layer)
# ---------------------------------------------------------------------------
Base.@deprecate_binding ContinuousHab ContinuousRegime
Base.@deprecate_binding ContinuousTimeHab ContinuousRegime
# Two renames, one hop. v0.4.0's `DiscreteHab` became `DiscreteRegime`, which v0.5.0 renamed
# again to `CategoricalRegime` — the layer is categorical (class codes that must not be averaged),
# not "discrete" in the shipped `ValueType` sense, where a day count is discrete and interpolates
# perfectly well. The shim points straight at the current name rather than chaining through the
# intermediate one, which never shipped as anything a user could have written.
Base.@deprecate_binding DiscreteHab CategoricalRegime
Base.@deprecate_binding HabitatCollection2 RegimeCollection2
Base.@deprecate_binding HabitatCollection3 RegimeCollection3
@deprecate gethabitat getregime

# ---------------------------------------------------------------------------
# Environment container: `AbioticEnv`/`GridAbioticEnv` → `Condition` (v0.4.0 rename). NB the *condition
# layer* `AbstractHabitat` was renamed to `AbstractRegime`, freeing `AbstractHabitat` for the environment;
# v0.4.0's (unexported) `AbstractHabitat` therefore changes meaning — a NEWS breaking note, not a shim.
# ---------------------------------------------------------------------------
Base.@deprecate_binding GridAbioticEnv GridHabitat

# `clearcache` -> `clearcache!` (v0.5.0): it deletes the JLD2 files a `CachedEcosystem` recorded, so
# the mutation marker is owed. Exported in v0.4.0, hence the shim.
@deprecate clearcache clearcache!

# `searchdir` -> `_searchdir` (v0.5.0): a bare directory listing, private now because an end user has
# no reason to call it. It was exported by `EcoSISTEM.ClimatePref` in v0.4.0, so it is owed a shim
# even though its replacement is internal -- the redirect is what tells a caller where it went.
@deprecate searchdir _searchdir

# ---------------------------------------------------------------------------
# Condition line: `Trait`/`Bin` → `Tolerance`/`NicheTolerance` (v0.4.0 rename; the species' condition response)
# ---------------------------------------------------------------------------
# **These two are no longer bindings, because the type they named is gone and the one that
# replaced it needs an argument they did not supply.** `CategoricalTolerance` and
# `LandCoverTolerance` merged into `SimpleCategoricalTolerance`, whose `penalty` field now carries
# the distinction their *types* used to: `DiscreteTrait` scored `0.5` outside a species' class
# (soft exclusion, via `Match`) and `LCtrait` scored `0.0` (hard, via `LCmatch`). So each shim pins
# its own released value, and neither inherits the new `0.0` default by accident.
# `DiscreteTrait` also took **one class per species** rather than a set; the constructor wraps
# each in a one-element set, which is the same model.
# **Each shim also assigns an axis**, which is now required and which neither released name took:
# `DiscreteTrait` is a generic niche label, so it gets the grouping `TypologyAxis`, while `LCtrait`
# is specifically land cover and gets `LandCoverTypology`.
@deprecate DiscreteTrait(vals) SimpleCategoricalTolerance(vals,
                                                          axis = TypologyAxis,
                                                          penalty = 0.5)
@deprecate LCtrait(vals) SimpleCategoricalTolerance(vals,
                                                    axis = LandCoverTypology,
                                                    penalty = 0.0)
Base.@deprecate_binding TraitCollection2 ToleranceCollection2
Base.@deprecate_binding TraitCollection3 ToleranceCollection3
Base.@deprecate_binding TempBin TempTolerance
Base.@deprecate_binding RainBin RainTolerance
@deprecate traitpopulate!(args...) populate_by_tolerance!(args...)
@deprecate traitrepopulate!(args...) repopulate_by_tolerance!(args...)

# ---------------------------------------------------------------------------
# Condition line: the matcher `TraitRelationship`/`Match`/… → `NicheFit`/`Suitability` (v0.4.0 rename)
# ---------------------------------------------------------------------------
# **Both now name the same type, and that is correct rather than a collision.** `Match` and
# `LCmatch` differed only in the weight they gave a species outside its classes, which is no longer
# a property of the fit at all — it moved to the tolerance's `penalty` (see `DiscreteTrait`/`LCtrait`
# above). What is left of either is "score a categorical match", and there is one of those.
Base.@deprecate_binding Match CategoricalSuitability{NicheAxis}
Base.@deprecate_binding LCmatch CategoricalSuitability{NicheAxis}
Base.@deprecate_binding NoRelContinuous NoFitContinuous
Base.@deprecate_binding NoRelDiscrete NoFitCategorical   # re-pointed; see `DiscreteHab` above
Base.@deprecate_binding multiplicativeTR2 MultiplicativeFit2
Base.@deprecate_binding multiplicativeTR3 MultiplicativeFit3
Base.@deprecate_binding additiveTR2 AdditiveFit2
Base.@deprecate_binding additiveTR3 AdditiveFit3
@deprecate gettraitrel getnichefit

# ---------------------------------------------------------------------------
# Layer dynamics: the v0.4.0 change functions → a declared change
#
# **What a shim is owed for, since this file is where it gets inferred from.** A shim is owed only
# for a released **public** name — one v0.4.0 `export`ed (or declared `public`). Two kinds owe nothing:
# a name **new on this branch**, and a released but **unexported** one, because unexported means
# internal and nothing outside the package could have depended on it.
# This was got backwards once and cost seventeen shims. `HabitatUpdate`, `habitatupdate!`,
# `budgetupdate!`, `AbstractBudget`, `AbstractTimeBudget`, `AbstractRequirement`,
# `Abstract1Requirement`, `Abstract2Requirements`, `AbstractAbiotic`, `AbstractTraits`,
# `ContinuousTrait` and `AbstractTraitRelationship` were all unexported in v0.4.0, and
# `Abstract2Demands`, `ContinuousTimeRegime`, `AbstractTimeSupply`, `SolarTimeSupply` and
# `WaterTimeSupply` were never released at all — all seventeen have been removed (2026-08-07).
# The failure was *method*, not judgement: the shims already here were read as evidence of a rule
# and the rule written down to match them, instead of each being checked against the released tag.
# Check the name itself — `git grep "^ *export" v0.4.0 -- 'src/*.jl'`, remembering that an export list
# wraps across lines and a naive grep undercounts it badly — never the neighbouring lines.
#
# `LayerUpdate` (this branch's rename of v0.4.0's unexported `HabitatUpdate`) is replaced by the
# `AbstractLayerChange` hierarchy (`Layer.jl`), which carries the change's *mode* as a type
# parameter and validates its values against the layer's own unit (`changeunit`, `LayerChange.jl`)
# instead of against a dimension passed in by hand and then discarded. The name survives below as a
# constructor mapping each old change function onto its successor.
# ---------------------------------------------------------------------------

# The v0.4.0 change functions. Each was only ever *named* — stored in a `LayerUpdate` and invoked by
# the update loop — so the deprecated form keeps the old `(eco, layer, timestep)` signature, warns,
# and defers to whatever change the layer now holds. `NoChange` is the one whose name is reused as a
# type, so its shim is a three-argument constructor.

"""
    TempChange(eco::AbstractEcosystem, layer::AbstractLayer, timestep::Unitful.Time)

Apply one timestep of `layer`'s own change. Deprecated: a temperature increase is now declared as
`IncrementBy(rate)`, in the layer's unit per unit time, rather than named as a change function
that hard-codes `K`.
"""
function TempChange(eco::AbstractEcosystem, layer::AbstractLayer,
                    timestep::Unitful.Time)
    Base.depwarn("`TempChange` is deprecated; declare `IncrementBy(rate)` on the layer instead.",
                 :TempChange)
    return _applychange!(layer.change, layer, simulationtime(eco), timestep)
end

"""
    RainfallChange(eco::AbstractEcosystem, layer::AbstractLayer, timestep::Unitful.Time)

Apply one timestep of `layer`'s own change. Deprecated: a rainfall change is now declared as
`IncrementBy(rate)`, in the layer's unit per unit time, rather than named as a change function
that hard-codes `mm` — which stopped being precipitation's unit when it became the rate `mm/d`.
"""
function RainfallChange(eco::AbstractEcosystem, layer::AbstractLayer,
                        timestep::Unitful.Time)
    Base.depwarn("`RainfallChange` is deprecated; declare `IncrementBy(rate)` on the layer " *
                 "instead.", :RainfallChange)
    return _applychange!(layer.change, layer, simulationtime(eco), timestep)
end

"""
    TempFluct(eco::AbstractEcosystem, layer::AbstractLayer, timestep::Unitful.Time)

Apply one timestep of `layer`'s own change. Deprecated: an oscillation is now declared as
`IncrementBy(PatternedChange(amplitude, timescale))`.
"""
function TempFluct(eco::AbstractEcosystem, layer::AbstractLayer,
                   timestep::Unitful.Time)
    Base.depwarn("`TempFluct` is deprecated; declare " *
                 "`IncrementBy(PatternedChange(amplitude, timescale))` on the layer instead.",
                 :TempFluct)
    return _applychange!(layer.change, layer, simulationtime(eco), timestep)
end

"""
    cyclic_change(eco::AbstractEcosystem, layer::AbstractLayer, timestep::Unitful.Time)

Apply one timestep of `layer`'s own change. Deprecated: walking a stored time series is now the
layer's declared change rather than a named change function. `eraChange`/`worldclimChange` are
aliases of this — the logic never differed by source or role.
"""
function cyclic_change(eco::AbstractEcosystem, layer::AbstractLayer,
                       timestep::Unitful.Time)
    Base.depwarn("`cyclic_change` is deprecated; a time-varying layer now carries its own " *
                 "series-valued change.", :cyclic_change)
    return _applychange!(layer.change, layer, simulationtime(eco), timestep)
end

const eraChange = cyclic_change
const worldclimChange = cyclic_change

# `NoChange` is now an `AbstractChangeMode`, so this three-argument form is a *constructor* on
# the mode type rather than the old change function. A static layer is `NoLayerChange()`.
function NoChange(eco::AbstractEcosystem, layer::AbstractLayer,
                  timestep::Unitful.Time)
    Base.depwarn("`NoChange(eco, layer, timestep)` is deprecated; `NoChange` is now a change " *
                 "*mode*, and a layer that never changes holds `NoLayerChange()`.",
                 :NoChange)
    return nothing
end

"""
    HabitatLoss(eco::AbstractEcosystem, regime::AbstractRegime, timestep::Unitful.Time)

Destroy regime for one timestep of the ecosystem. The regime's change carries a loss rate (per
unit time); over `timestep` it gives the per-cell probability that an active cell is lost. That
many active cells are drawn at random and have their supply and abundances zeroed.

Deprecated: habitat loss changes the *ecosystem*, not a layer — it never marks the lost cells
inactive, so they stay dispersal targets and can be lost again, and it draws from the global RNG,
so it neither replays reproducibly nor agrees across MPI ranks. It is superseded by an explicit
cell-deactivating intervention.
"""
function HabitatLoss(eco::AbstractEcosystem, regime::AbstractRegime,
                     timestep::Unitful.Time)
    Base.depwarn("`HabitatLoss` is deprecated: it mutates the ecosystem from a layer's change " *
                 "slot, leaves the cells it destroys active, and draws from the global RNG. Use " *
                 "a cell-deactivating intervention instead.", :HabitatLoss)
    # Loss rate × timestep is the (dimensionless) probability of losing a cell.
    prob = uconvert(NoUnits, regime.change.rate * timestep)
    pos = findall(vec(eco.habitat.active))
    smp = sample(pos, rand(Binomial(length(pos), prob)), replace = false)
    eco.habitat.supply.matrix[smp] .= zero(eltype(eco.habitat.supply.matrix))
    eco.abundances.matrix[:, smp] .= 0
    return eco
end

# Map a v0.4.0 change function to the change type that replaces it. The rate is stored as given
# rather than converted to the layer's `changeunit`: with no layer in hand this constructor cannot
# know the unit the rate *should* be in, which is precisely the check the `_attachchange` path adds
# and this one structurally cannot. Applying it still works — `matrix .+= rate .* timestep`
# reconciles any per-time unit — it simply is not validated.
_legacychange(::Type{NoChange}, rate) = NoLayerChange()
function _legacychange(::typeof(TempChange), rate)
    return SteadyLayerChange{typeof(rate)}(rate)
end
function _legacychange(::typeof(RainfallChange), rate)
    return SteadyLayerChange{typeof(rate)}(rate)
end
function _legacychange(::typeof(TempFluct), rate)
    @warn "`TempFluct`'s semantics have changed. It now oscillates as a function of elapsed " *
          "simulation time, with a one-year period, rather than feeding the layer's own values " *
          "back into `sin` (a path-dependent walk whose cycle came from the values themselves). " *
          "Runs using it will NOT reproduce previous results."
    return PatternedLayerChange{RateChange, typeof(sinusoidal), typeof(rate),
                                Nothing}(sinusoidal, rate, uconvert(s, 1.0year),
                                         nothing)
end
_legacychange(::typeof(HabitatLoss), rate) = LegacyLoss(rate)
# `cyclic_change` walked a 3-D layer's own stored stack. There is no layer here to take that stack
# from — `LayerUpdate` is built before the layer it will drive — so this cannot be materialised into
# a `SeriesLayerChange` and says so. The stack is installed by whichever constructor reads it.
function _legacychange(::typeof(cyclic_change), rate)
    return error("`cyclic_change` is no longer a change a layer can be given: a stored series now " *
                 "lives in the layer's own change, installed from the stack when the layer is " *
                 "built, and there is no stack to take here. Build the layer from its stack (e.g. " *
                 "`SourceSpec(WorldClim{Climate}, :wind, month = 1:12)`), or pass " *
                 "`ReplaceWith(SeriesChange(stack))` to `setchange!`.")
end
function _legacychange(changefun, rate)
    return error("`$changefun` is not a recognised layer change function; build an " *
                 "`AbstractLayerChange` (e.g. `IncrementBy($rate)`) instead")
end

"""
    LayerUpdate(changefun, rate)
    LayerUpdate(changefun, rate, ::Type{<:Unitful.Dimensions})

Construct the [`AbstractLayerChange`](@ref) that replaces the v0.4.0 change function `changefun`.
The third argument, a dimension to check `rate` against, is accepted and ignored — a change is now
checked against the layer it is attached to instead.
"""
function LayerUpdate(changefun, rate)
    Base.depwarn("`LayerUpdate(changefun, rate)` is deprecated; layer change is now an " *
                 "`AbstractLayerChange` (`NoLayerChange()`, `IncrementBy(rate)`, …) attached to a " *
                 "layer, which validates it against that layer's own unit.",
                 :LayerUpdate)
    return _legacychange(changefun, rate)
end
function LayerUpdate(changefun, rate, ::Type{<:Unitful.Dimensions})
    return LayerUpdate(changefun, rate)
end

# `resetrate!` → `setchange!`. The old name could only ever install a constant rate on the *whole*
# regime; `setchange!` takes any change and addresses one layer, so it also reaches a sub-layer of a
# collection, which `resetrate!` could not. It stays exported, as it was in v0.4.0; `setchange!` is
# `public` rather than exported, per the convention for new supported API.
@deprecate resetrate!(eco::AbstractEcosystem,
                      rate::Unitful.Quantity) setchange!(eco.habitat.regime,
                                                         IncrementBy(rate))

# ---------------------------------------------------------------------------
# Environment constructors: the `*habitat` family → `GridHabitat`
#
# `simplehabitat` / `simplenichehabitat` / `erahabitat` / `worldclimhabitat` / `bioclimhabitat` /
# `landcoverhabitat` (and the older `*AE` spellings of each) all do the same thing in six shapes:
# take a *value or a raster*, a grid `dimension`, an `area` and a `maxsupply`, and hand back a
# finished `GridHabitat`. `GridHabitat(; regime, supply, area)` takes layer **recipes** on a
# grid that a `StudyArea` has already decided, which is what lets the same layer be inspected,
# masked, reprojected and compared before anything is built on it — none of which a builder that
# invents its own grid from a `dimension` can take part in.
#
# The whole family and its private helpers live here as one unit, out of `src/GridHabitat.jl`, which
# is left holding `GridHabitat` itself. `_maxsupply_env` in particular is now used only by
# deprecated code (these builders and the gradient family above), which is the clearest sign the
# `maxsupply`-and-`area` shape has been superseded: `GridHabitat` takes a supply *spec* per
# cell, not one total to spread.
#
# As with the gradient family, these forward to the private bodies rather than to
# `GridHabitat`. That is now a choice rather than a blocker: the reason recorded above —
# that a spec had no way to declare a change — has since gone (`Varying(spec, IncrementBy(rate))`).
# What remains is that `GridHabitat` would put these on a *different grid*: a `StudyArea`
# derives cells from an extent and a cell size, where these derive a cell size from `dimension` and
# `area` and place the layer at no real-world position at all. Forwarding would therefore change
# results silently, which is worse for a shim than not forwarding.
#
# The `*AE` names resolve straight onto the shared bodies rather than chaining through their
# `*habitat` renames, so a caller of the oldest name gets one warning naming the real target, not
# two naming each other — the same handling as the gradient section.
# ---------------------------------------------------------------------------

# Shared tail of the maximum-supply `*AE` constructors: total the per-area maximum supply
# `maxsupply` (an areal rate) over `area`, spread it uniformly across the grid of `regime`,
# pick the supply type from its axis (`supplytype`),
# and assemble the `GridHabitat` with regime `regime`. `active` is passed through as given (a
# plain `Matrix{Bool}` or an already-`DimArray` mask) — `GridHabitat`'s own constructor wraps/
# validates it against `regime`'s grid.
# The v0.4.0 `BUDGETDICT` in the one place it still belongs: inside the deprecation shims.
#
# **The unit table is gone from the live API, not from history.** These builders' *documented*
# v0.4.0 contract was to pick the budget from `maxsupply`'s unit (`kJ => SolarBudget`,
# `mm => WaterBudget`, `NoUnits => SimpleBudget`), so reproducing that behaviour is what a shim is
# for. Keeping it here rather than as a general method confines unit-based inference to the code
# whose job is to be backwards-compatible, and leaves the live supply path dispatching on the axis.
#
# `NoUnits => SimpleBudget` has **no replacement**: the free supply family was removed, so that
# case now errors instead of silently building something else.
# The element type each `cancel` method converts into, derived from `canonicalunit` rather than
# spelling the literal unit here — so a change to what a resource is measured in cannot leave these
# converting to the old one. Moved here from `Layer.jl` (2026-08-20): the three `cancel` methods
# below are their only callers, and a v0.4.0 conversion belongs with the rest of the shims.
const _SolarRate = typeof(1.0 * canonicalunit(Resource, SolarRadiation))
const _WaterRate = typeof(1.0 * canonicalunit(Resource, Precipitation))
const _CarbonRate = typeof(1.0 * canonicalunit(Resource, CarbonFlux))

# The two-argument, dimension-dispatched `cancel` — v0.4.0's areal-rate × cell-area conversion,
# kept here for the same reason as `_v040supplytype` below and moved beside it (2026-08-09): its
# only caller is `_maxsupply_env`, and choosing a unit from a unit is exactly what the live path
# stopped doing. The three-argument axis form in `GridHabitat.jl` is what everything else uses.
# Verified against `dimension(...)` before wiring: solar's `kJ/m²/day` (𝐌𝐓⁻³) × m² (𝐋²) →
# `kJ/day` (𝐋²𝐌𝐓⁻³, `Unitful.Power`); water's `L/m²/day` (𝐋𝐓⁻¹) × m² → `L/day` (𝐋³𝐓⁻¹,
# `VolumeFlow`); carbon's `g/m²/day` (𝐌𝐋⁻²𝐓⁻¹) × m² → `g/day` (𝐌𝐓⁻¹, `Unitful.MassFlow`).
# `test_rasters.jl`'s wind-speed regression test also asks these directly — deliberately, as
# proof of what the deleted unit table would have said — so they are not callable only from here.
function cancel(a::Quantity{<:Real, 𝐌 * 𝐓^-3}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(unit(_SolarRate), a * b)
end
function cancel(a::Quantity{<:Real, 𝐋 * 𝐓^-1}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(unit(_WaterRate), a * b)
end
function cancel(a::Quantity{<:Real, 𝐌 * 𝐋^-2 * 𝐓^-1}, b::Quantity{<:Real, 𝐋^2})
    return uconvert(unit(_CarbonRate), a * b)
end

# The niche axis a v0.4.0 `maxsupply` implies, chosen from its unit exactly as `_v040supplytype`
# always has. Split out because the reroute below names an *axis* in a spec where the old code named
# a supply *type* directly.
function _v040supplyaxis(B)
    d = dimension(B)
    d == dimension(1.0u"kJ/d") && return SolarRadiation
    d == dimension(1.0u"L/d") && return Precipitation
    d == dimension(1.0u"g/d") && return CarbonFlux
    return error("this deprecated builder chose its supply from `maxsupply`'s unit, and $(unit(B)) " *
                 "no longer maps to one. The dimensionless/free supply (v0.4.0's `SimpleBudget`) " *
                 "has been removed; use `GridHabitat` with a spec naming an axis instead.")
end

# **The bridge from a v0.4.0 grid description to a `StudyArea`, and the reason these shims still
# work at all.** A v0.4.0 builder states its grid as `(dimension, total area)`; a `StudyArea` states
# it as `(extent, cellsize)`. Translating here is exact rather than approximate: `_gridgeometry`
# recovers `sqrt(a / (ny·nx))` from the extent and cell size, which is the very quantity computed
# below — measured over all five surviving families, bit-identical for three and equal to within
# 1e-15 relative for the other two (pure float round-off, the same interpolation reached by a
# different order of operations).
#
# `extent` and `dimension` are both `(y, x)`, so they pass through in order. They used to be crossed
# here, because `extent` alone was `(x, y)`.
# This is still `deprecations.jl` stating a cell size rather than deriving one — the exemption
# `CLAUDE.md` grants this file — but it now states it *through* `StudyArea`, so the grid is described
# once and the equal-area arithmetic v0.4.0 depended on is preserved.
function _v040area(dimension::Tuple{Int64, Int64}, area::Unitful.Area, active)
    cs = sqrt(uconvert(km^2, area) / (dimension[1] * dimension[2]))
    return StudyArea(extent = (dimension[1] * cs, dimension[2] * cs),
                     cellsize = cs, within = active, verbosity = :silent)
end

# Build a v0.4.0-shaped environment from a *spec* plus the old `(dimension, maxsupply, area)` triple.
# Replaces the former `_maxsupply_env`, which took an already-built regime and assembled a
# `GridHabitat` from parts — a construction route that no longer exists, `GridHabitat` now being
# built from specs and a `StudyArea` alone.
function _v040env(spec, dimension::Tuple{Int64, Int64}, active,
                  maxsupply::Unitful.Quantity{Float64},
                  area::Unitful.Area)
    # The supply is stated per unit area and multiplied by the cell area, which is precisely what
    # v0.4.0 did by hand as `cancel(maxsupply, area) / countsubcommunities`.
    axis = _v040supplyaxis(cancel(maxsupply, area))
    return GridHabitat(regime = spec,
                       supply = UniformSpec(maxsupply, axis = axis),
                       area = _v040area(dimension, area, active))
end

# **The three data-regime helpers that used to be here are gone with the builders that used
# them** (`_timeclimate_regime`, `_continuousregime`, `_categoricalregime`). They existed to build a
# regime from a raster's own cells, preserving `_cellsize`'s area-preserving metric conversion — the
# thing that made the data-backed builders unreproducible through a `StudyArea`, and so the reason
# those builders were removed rather than rerouted.

function _simplenichehabitat(numniches::Int64,
                             dimension::Tuple,
                             maxsupply::Unitful.Quantity{Float64},
                             area::Unitful.Area{Float64},
                             active::Matrix{Bool};
                             axis::Type{<:NicheAxis} = TypologyAxis)
    return _v040env(NicheSpec(numniches, axis = axis), dimension, active,
                    maxsupply, area)
end
function _simplenichehabitat(numniches::Int64,
                             dimension::Tuple,
                             maxsupply::Unitful.Quantity{Float64},
                             area::Unitful.Area{Float64};
                             axis::Type{<:NicheAxis} = TypologyAxis)
    active = Matrix{Bool}(undef, dimension)
    fill!(active, true)
    return _simplenichehabitat(numniches, dimension, maxsupply, area, active,
                               axis = axis)
end

# **The four data-backed environment builders are REMOVED, and only their names remain** — see
# `_removedbuilder` below for the message and the reasoning. Their implementations built a regime
# **directly from an in-memory raster with no resampling**, which is not something any construction
# route still does: a `StudyArea` decides the grid and the data is sampled onto it.
#
# Nothing dispatches here any more, so there is no `_erahabitat`/`_worldclimhabitat` body at all —
# the generated shims error before reaching one. The declarations that used to be here (three
# method-less stubs whose methods lived in `EcoSISTEMRasterDataSourcesExt`) are gone with them, and so
# are those extension methods.

function _simplehabitat(val::Union{Float64, Unitful.Quantity{Float64}},
                        dimension::Tuple{Int64, Int64},
                        maxsupply::Unitful.Quantity{Float64},
                        area::Unitful.Area{Float64},
                        active::Matrix{Bool};
                        axis::Type{<:NicheAxis} = NicheAxis)
    return _v040env(UniformSpec(val, axis = axis), dimension, active,
                    maxsupply, area)
end

function _simplehabitat(val::Union{Float64, Unitful.Quantity{Float64}},
                        dimension::Tuple{Int64, Int64},
                        maxsupply::Unitful.Quantity{Float64},
                        area::Unitful.Area{Float64};
                        axis::Type{<:NicheAxis} = NicheAxis)
    active = fill(true, dimension)
    return _simplehabitat(val, dimension, maxsupply, area, active, axis = axis)
end

"""
    simplehabitat(val, dimension, maxsupply, area[, active]; axis = NicheAxis)
    simplenichehabitat(numniches, dimension, maxsupply, area[, active]; axis = TypologyAxis)

!!! warning "Deprecated"
    Both are deprecated and will be removed; use
    `GridHabitat(regime = UniformSpec(val; axis), supply = …, area = …)` — or
    [`NicheSpec`](@ref) in place of [`UniformSpec`](@ref) for the random-niche form — with a
    [`StudyArea`](@ref) deciding the grid. `maxsupply` has no direct equivalent: a supply spec gives
    the resource *per cell*, rather than a total spread across the grid.
"""
simplehabitat

@doc (@doc simplehabitat) simplenichehabitat

"""
    erahabitat(era, maxsupply, area[, active]; axis = NicheAxis)
    erahabitat(era, supply[, active]; axis = NicheAxis)
    worldclimhabitat(worldclim, maxsupply, area[, active]; axis = NicheAxis)
    worldclimhabitat(worldclim, supply[, active]; axis = NicheAxis)
    bioclimhabitat(bioclim, maxsupply, area[, active]; axis = NicheAxis)
    bioclimhabitat(bioclim, supply[, active]; axis = NicheAxis)
    landcoverhabitat(landcover, maxsupply, area[, active]; axis = NicheAxis)
    landcoverhabitat(landcover, supply[, active]; axis = NicheAxis)

!!! warning "Deprecated"
    All four are deprecated and will be removed; use
    `GridHabitat(regime = SourceSpec(source, code), supply = …, area = …)`, with a
    [`StudyArea`](@ref) deciding the grid from the data's own CRS, extent and resolution. A
    [`ConstructedSpec`](@ref) wraps an already-read [`ClimateRaster`](@ref) where one is in hand.
    Unlike these builders, that route reports what the grid costs each layer, masks and reprojects,
    and can be inspected with [`investigate_study_area`](@ref) before anything is built.
"""
erahabitat

@doc (@doc erahabitat) worldclimhabitat
@doc (@doc erahabitat) bioclimhabitat
@doc (@doc erahabitat) landcoverhabitat

# The `*AE` names are re-exported here because they no longer go through the symbol-form
# `@deprecate` that used to do it for them.
export simplehabitatAE, simplenicheAE, eraAE, worldclimAE, bioclimAE, lcAE

# The builders that survive: still deprecated, still working, and now reaching `GridHabitat`
# through a synthetic `StudyArea` rather than assembling a habitat from parts.
for (old, since, body, target) in ((:simplehabitat, "v0.5.0", :_simplehabitat,
    "UniformSpec(val, axis = …)"),
    (:simplehabitatAE, "v0.4.0", :_simplehabitat,
    "UniformSpec(val, axis = …)"),
    (:simplenichehabitat, "v0.5.0",
    :_simplenichehabitat,
    "NicheSpec(numniches, axis = …)"),
    (:simplenicheAE, "v0.4.0",
    :_simplenichehabitat,
    "NicheSpec(numniches, axis = …)"))
    @eval function $old(args...; kwargs...)
        Base.depwarn("`$($(QuoteNode(old)))` is deprecated (since $($since)); use " *
                     "`GridHabitat(regime = $($target), supply = …, area = …)` instead, " *
                     "with a `StudyArea` deciding the grid.", $(QuoteNode(old)))
        return $body(args...; kwargs...)
    end
end

# **The four data-backed builders are REMOVED — the names remain only to explain themselves.**
#
# **Why an error rather than a deletion.** Four of these eight names (`eraAE`, `worldclimAE`,
# `bioclimAE`, `lcAE`) were **exported in v0.4.0**, so a released script calls them. Deleting the
# binding gives `UndefVarError: eraAE not defined`, which says nothing about what happened or what to
# do; the name kept, erroring, is the difference between a five-minute fix and an afternoon.
#
# **No `Base.depwarn` first, deliberately.** A deprecation warning followed by an error is two
# messages for one problem, and the warning is the less useful one — these are not deprecated (still
# working, to be removed later), they are *gone*.
#
# **Each message names `in_memory_raster` as well as the `SourceSpec` form.** Every one of these
# builders took an *already-read* `ERA`/`ClimateRaster`, so anyone reaching this error is by
# definition holding a raster; pointing only at `SourceSpec` would answer a question they did not ask.
function _removedbuilder(name, since, what, sourcespec)
    return error("`$name` has been removed (it was deprecated in $since). It built an environment " *
                 "directly from an in-memory $what, sampling **no** grid: the raster's own cells " *
                 "became the simulation's cells. Nothing does that now — a `StudyArea` decides the " *
                 "grid and the data is sampled onto it — so there is no shim that would give the " *
                 "same answer, and results from a replacement will differ from v0.4.0's.\n" *
                 "  • Reading from the catalogued source:\n" *
                 "        area = StudyArea(regime = $sourcespec)\n" *
                 "        env  = GridHabitat(regime = $sourcespec, supply = …, area = area)\n" *
                 "  • A raster you already hold (what this builder took):\n" *
                 "        spec = in_memory_raster(raster, axis = …)\n" *
                 "        env  = GridHabitat(regime = spec, supply = …, area = StudyArea(regime = spec))")
end

for (old, since, what, sourcespec) in ((:erahabitat, "v0.5.0", "`ERA` raster",
    "SourceSpec(…)"),
    (:eraAE, "v0.4.0", "`ERA` raster", "SourceSpec(…)"),
    (:worldclimhabitat, "v0.5.0", "WorldClim climate raster",
    "SourceSpec(WorldClim{Climate}, code)"),
    (:worldclimAE, "v0.4.0", "WorldClim climate raster",
    "SourceSpec(WorldClim{Climate}, code)"),
    (:bioclimhabitat, "v0.5.0", "WorldClim bioclim raster",
    "SourceSpec(WorldClim{BioClim}, code)"),
    (:bioclimAE, "v0.4.0", "WorldClim bioclim raster",
    "SourceSpec(WorldClim{BioClim}, code)"),
    (:landcoverhabitat, "v0.5.0", "land-cover raster",
    "SourceSpec(EarthEnv{LandCover})"),
    (:lcAE, "v0.4.0", "land-cover raster",
    "SourceSpec(EarthEnv{LandCover})"))
    @eval function $old(args...; kwargs...)
        return _removedbuilder($(QuoteNode(old)), $since, $what, $sourcespec)
    end
end

# ---------------------------------------------------------------------------
# The rainfall-gradient builders: `raingrad` / `raingradhabitat` / `raingradAE` → `GridHabitat`
#
# The whole rainfall-gradient family is deprecated together, and its bodies move here out of
# `src/GridHabitat.jl` and `src/Layer.jl`: `GridHabitat(regime = GradientSpec(low, high;
# axis = Precipitation), …)` expresses a gradient on any axis, so a rain-specific builder — and
# the rain-specific generator behind it — earn nothing. `raingrad`'s only callers were
# `raingradhabitat`'s two real methods, so the pair moves as a unit.
#
# These forward to a *private* body rather than to `GridHabitat`, unlike their `*AE`
# siblings above — because they would land on a **different grid**. A `StudyArea` derives its cells
# from an extent and a cell size; these builders derive a cell size from `dimension` and `area`, so
# a literal forward would change results silently. Behaviour is preserved here instead and the
# warning names the migration target.
#
# **The `rate` argument is no longer the blocker, and was until 2026-08-03.** It used to be:
# every layer `GridHabitat` built was armed at rate 0, so forwarding turned a wetting climate
# static. `Varying(spec, IncrementBy(rate))` now expresses exactly this, so only the grid remains.
#
# `raingradAE` calls that same body directly instead of chaining through `raingradhabitat`, so a
# caller of the oldest name gets one warning naming the real target, not two naming each other.
#
# Being wholly inside the deprecated unit is also what lets these be typed for `Precipitation`'s
# canonical unit — the rate `mm/d` — without touching live code: as `Unitful.Length{Float64}` they
# could not be called on their own default axis at all.
# ---------------------------------------------------------------------------

# The rainfall-gradient regime itself. `minR`/`maxR` are Condition
# values on the `Precipitation` axis, so they carry its canonical *rate* unit; `rate` is the change
# of that per unit time, so both are general `Quantity`s and `LayerUpdate`'s own `rate * 1month_mean_duration`
# check does the validating.
function _raingrad(minR::Unitful.Quantity{Float64},
                   maxR::Unitful.Quantity{Float64},
                   size::Unitful.Length{Float64},
                   dim::Tuple{Int64, Int64},
                   rate::Unitful.Quantity{Float64})
    dim[1] > 1 ||
        error("First dimension should be greater than 1 for rainfall gradient")
    M = Array{typeof(minR)}(undef, dim)
    total = dim[1]
    rain_range = collect(range(minR, stop = maxR, length = total))
    map(1:total) do seq
        return M[seq, :] .= rain_range[seq]
    end
    # Built directly rather than through the deprecated `LayerUpdate`, which would add a second
    # warning on top of the builder's own. The rate is unvalidated here for the same reason it is
    # there: the layer it will drive does not exist yet.
    regimeupdate = SteadyLayerChange{typeof(rate)}(rate)
    return ContinuousRegime(M, size, regimeupdate)
end

# The `maxsupply` form: a rainfall gradient regime plus a supply
# filled from a separate maximum supply value. Which supply that is comes from `maxsupply`'s
# *unit*, via `_v040supplytype` — the v0.4.0 contract these shims exist to reproduce, and the one
# place unit-based inference legitimately survives.
function _raingradhabitat(minR::Unitful.Quantity{Float64},
                          maxR::Unitful.Quantity{Float64},
                          dimension::Tuple{Int64, Int64},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area{Float64},
                          rate::Unitful.Quantity{Float64},
                          active::Matrix{Bool};
                          axis::Type{<:NicheAxis} = Precipitation)
    return _v040env(Varying(GradientSpec(minR, maxR, axis = axis),
                            IncrementBy(rate)), dimension, active, maxsupply,
                    area)
end

# The water-supply form: the gradient's own rainfall values are the supply, so there is no
# separate `maxsupply`.
#
# **This body was broken before this change, and reachable** — `raingradhabitat(minR, maxR, dim,
# area, rate)` is the exported 5-argument form that lands here. It passed the regime's matrix
# straight to `WaterSupply`, but the regime holds `Precipitation`'s *condition* unit (`mm/day`, an
# areal rate) while a supply is `L/day` **per cell**, so the call had no method at all. v0.4.0 could
# do it because both sides were bare `mm`; the v0.5.0 unit change moved the regime to `mm/day` and
# the supply to an absolute `L/day`, and this line was not moved with them. No test reaches it,
# which is why it survived — `raingrad`'s coverage all goes through the `maxsupply` form.
# The fix is what `cancel` exists for and what the other supply paths already do: the areal rate
# against this grid's own cell area. That is the faithful reading of "the rainfall itself is the
# water budget" in units where a supply is per cell.
function _raingradhabitat(minR::Unitful.Quantity{Float64},
                          maxR::Unitful.Quantity{Float64},
                          dimension::Tuple{Int64, Int64},
                          area::Unitful.Area{Float64},
                          rate::Unitful.Quantity{Float64},
                          active::Matrix{Bool};
                          axis::Type{<:NicheAxis} = Precipitation)
    # The regime spec is used **twice**: once as the Condition and once as the Resource. That is
    # exactly what "the rainfall itself is the water budget" means, and it is what the hand-built
    # `cancel.(regime.matrix, cellarea, Precipitation)` computed — a supply spec is stated per unit
    # area and multiplied by the cell area, which is the same arithmetic in one step.
    # Only the regime carries the rate: a declared change belongs to one layer, and it is the
    # *condition* that v0.4.0 drifted.
    grad = GradientSpec(minR, maxR, axis = axis)
    return GridHabitat(regime = Varying(grad, IncrementBy(rate)),
                       supply = GradientSpec(minR, maxR,
                                             axis = Precipitation),
                       area = _v040area(dimension, area, active))
end

# `active`-defaulting convenience forms, matching the two above.
function _raingradhabitat(minR::Unitful.Quantity{Float64},
                          maxR::Unitful.Quantity{Float64},
                          dimension::Tuple{Int64, Int64},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area{Float64},
                          rate::Unitful.Quantity{Float64};
                          axis::Type{<:NicheAxis} = Precipitation)
    return _raingradhabitat(minR, maxR, dimension, maxsupply, area, rate,
                            fill(true, dimension), axis = axis)
end
function _raingradhabitat(minR::Unitful.Quantity{Float64},
                          maxR::Unitful.Quantity{Float64},
                          dimension::Tuple{Int64, Int64},
                          area::Unitful.Area{Float64},
                          rate::Unitful.Quantity{Float64};
                          axis::Type{<:NicheAxis} = Precipitation)
    return _raingradhabitat(minR, maxR, dimension, area, rate,
                            fill(true, dimension), axis = axis)
end

"""
    raingrad(minR, maxR, size, dim, rate)

!!! warning "Deprecated"
    `raingrad` is deprecated and will be removed. Build the gradient as a layer recipe instead —
    `GridHabitat(regime = GradientSpec(minR, maxR, axis = Precipitation), …)` — which does
    the same job on any niche axis, with `Varying(spec, IncrementBy(rate))` for the rate. This shim
    keeps the old behaviour rather than forwarding, because a `StudyArea` would put the layer on a
    different grid.
"""
function raingrad(minR::Unitful.Quantity{Float64},
                  maxR::Unitful.Quantity{Float64},
                  size::Unitful.Length{Float64},
                  dim::Tuple{Int64, Int64},
                  rate::Unitful.Quantity{Float64})
    Base.depwarn("`raingrad` is deprecated; use `GridHabitat(regime = " *
                 "GradientSpec(minR, maxR, axis = Precipitation), …)` instead. " *
                 "(This shim preserves the old behaviour: a `StudyArea` grid differs.)",
                 :raingrad)
    return _raingrad(minR, maxR, size, dim, rate)
end

"""
    raingradhabitat(minR, maxR, dimension[, maxsupply], area, rate[, active]; axis = Precipitation)

!!! warning "Deprecated"
    `raingradhabitat` is deprecated and will be removed; use
    `GridHabitat(regime = GradientSpec(minR, maxR, axis = Precipitation), supply = …,
    area = …)` instead. `minR`/`maxR` are Condition values on the [`Precipitation`](@ref) axis,
    whose canonical unit is the **rate** `mm/d` — an accumulated depth (`mm`) is rejected. Note
    declare the rate with `Varying(spec, IncrementBy(rate))`. This shim preserves the old behaviour
    rather than forwarding.
"""
raingradhabitat

# `raingradAE` (the pre-v0.4.0 name) resolves straight onto the shared body — see the
# one-warning-not-two note above. Exported explicitly because it no longer goes through the
# symbol-form `@deprecate` that used to re-export it for us; `raingrad`/`raingradhabitat` keep the
# exports they already have beside their original `include`s.
export raingradAE
for (old, since) in ((:raingradhabitat, "v0.5.0"), (:raingradAE, "v0.4.0"))
    @eval function $old(args...; kwargs...)
        Base.depwarn("`$($(QuoteNode(old)))` is deprecated (since $($since)); use " *
                     "`GridHabitat(regime = GradientSpec(minR, maxR, axis = Precipitation), " *
                     "supply = …, area = …)` instead. (This shim preserves the old behaviour: a " *
                     "shim preserves the old behaviour.)", $(QuoteNode(old)))
        return _raingradhabitat(args...; kwargs...)
    end
end

# ---------------------------------------------------------------------------
# The temperature-gradient builders: `tempgrad` / `tempgradhabitat` / `peakedgradhabitat` (and the
# older `tempgradAE` / `peakedgradAE`) → `GridHabitat`
#
# The same treatment as the rainfall-gradient family above, for the same reasons — see that
# section's comment for the rationale, the private-body forwarding and the one-warning-not-two
# handling of the `*AE` names. `GradientSpec`/`PeakedSpec` on any niche axis replace them, and
# `_tempgrad` backs both builders, `peakedgradhabitat` mirroring its
# lower half to make the peak.
#
# Unlike the rain family these are *not* retyped: `Temperature`'s canonical unit is still `K`, so
# they were never broken — this is a pure move, and their `Unitful.Temperature`/`𝚯𝐓⁻¹` signatures
# stay exactly as released.
# ---------------------------------------------------------------------------

# The temperature-gradient regime itself; `_tempgradhabitat` and
# `_peakedgradhabitat` are its only callers.
function _tempgrad(minT::Unitful.Temperature{Float64},
                   maxT::Unitful.Temperature{Float64},
                   size::Unitful.Length{Float64},
                   dim::Tuple{Int64, Int64},
                   rate::Quantity{Float64, 𝚯 * 𝐓^-1})
    dim[1] > 1 ||
        error("First dimension should be greater than 1 for temperature gradient")
    M = Array{typeof(minT)}(undef, dim)
    total = dim[1]
    temp_range = collect(range(minT, stop = maxT, length = total))
    map(1:total) do seq
        return M[seq, :] .= temp_range[seq]
    end
    # Built directly rather than through the deprecated `LayerUpdate`, which would add a second
    # warning on top of the builder's own. The rate is unvalidated here for the same reason it is
    # there: the layer it will drive does not exist yet.
    regimeupdate = SteadyLayerChange{typeof(rate)}(rate)
    return ContinuousRegime(M, size, regimeupdate)
end

# A linear temperature gradient plus a supply filled from `maxsupply`, whose type comes from that
# value's unit — see `_v040supplytype`.
function _tempgradhabitat(minT::Unitful.Temperature{Float64},
                          maxT::Unitful.Temperature{Float64},
                          dimension::Tuple{Int64, Int64},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area{Float64},
                          rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                          active::Matrix{Bool};
                          axis::Type{<:NicheAxis} = Temperature)
    return _v040env(Varying(GradientSpec(minT, maxT, axis = axis),
                            IncrementBy(rate)), dimension, active, maxsupply,
                    area)
end
function _tempgradhabitat(minT::Unitful.Temperature{Float64},
                          maxT::Unitful.Temperature{Float64},
                          dimension::Tuple{Int64, Int64},
                          maxsupply::Unitful.Quantity{Float64},
                          area::Unitful.Area{Float64},
                          rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                          axis::Type{<:NicheAxis} = Temperature)
    return _tempgradhabitat(minT, maxT, dimension, maxsupply, area, rate,
                            fill(true, dimension), axis = axis)
end

# As `_tempgradhabitat`, but peaking at `maxT` in the middle: build a gradient over twice the
# range, then mirror its lower half onto the upper.
function _peakedgradhabitat(minT::Unitful.Temperature{Float64},
                            maxT::Unitful.Temperature{Float64},
                            dimension::Tuple{Int64, Int64},
                            maxsupply::Unitful.Quantity{Float64},
                            area::Unitful.Area{Float64},
                            rate::Quantity{Float64, 𝚯 * 𝐓^-1},
                            active::Matrix{Bool};
                            axis::Type{<:NicheAxis} = Temperature)
    return _v040env(Varying(PeakedSpec(minT, maxT, axis = axis),
                            IncrementBy(rate)), dimension, active, maxsupply,
                    area)
end
function _peakedgradhabitat(minT::Unitful.Temperature{Float64},
                            maxT::Unitful.Temperature{Float64},
                            dimension::Tuple{Int64, Int64},
                            maxsupply::Unitful.Quantity{Float64},
                            area::Unitful.Area{Float64},
                            rate::Quantity{Float64, 𝚯 * 𝐓^-1};
                            axis::Type{<:NicheAxis} = Temperature)
    return _peakedgradhabitat(minT, maxT, dimension, maxsupply, area, rate,
                              fill(true, dimension), axis = axis)
end

"""
    tempgrad(minT, maxT, size, dim, rate)

!!! warning "Deprecated"
    `tempgrad` is deprecated and will be removed. Build the gradient as a layer recipe instead —
    `GridHabitat(regime = GradientSpec(minT, maxT, axis = Temperature), …)` — which does the
    same job on any niche axis, with `Varying(spec, IncrementBy(rate))` for the rate. This shim
    keeps the old behaviour rather than forwarding, because a `StudyArea` would put the layer on a
    different grid.
"""
function tempgrad(minT::Unitful.Temperature{Float64},
                  maxT::Unitful.Temperature{Float64},
                  size::Unitful.Length{Float64},
                  dim::Tuple{Int64, Int64},
                  rate::Quantity{Float64, 𝚯 * 𝐓^-1})
    Base.depwarn("`tempgrad` is deprecated; use `GridHabitat(regime = " *
                 "GradientSpec(minT, maxT, axis = Temperature), …)` instead. " *
                 "(This shim preserves the old behaviour: a `StudyArea` grid differs.)",
                 :tempgrad)
    return _tempgrad(minT, maxT, size, dim, rate)
end

"""
    tempgradhabitat(minT, maxT, dimension, maxsupply, area, rate[, active]; axis = Temperature)
    peakedgradhabitat(minT, maxT, dimension, maxsupply, area, rate[, active]; axis = Temperature)

!!! warning "Deprecated"
    Both are deprecated and will be removed; use
    `GridHabitat(regime = GradientSpec(minT, maxT, axis = Temperature), supply = …,
    area = …)` — or [`PeakedSpec`](@ref) in place of [`GradientSpec`](@ref) for the peaked form.
    Declare the rate with `Varying(spec, IncrementBy(rate))`. These shims preserve the old behaviour
    rather than forwarding.
"""
tempgradhabitat

@doc (@doc tempgradhabitat) peakedgradhabitat

# The `*AE` names resolve straight onto the shared bodies, so each warns once naming the real
# target rather than chaining through its `*habitat` rename. Exported here because they no longer
# go through the symbol-form `@deprecate` that used to re-export them.
export tempgradAE, peakedgradAE

for (old, since, body, spec) in ((:tempgradhabitat, "v0.5.0", :_tempgradhabitat,
    "GradientSpec"),
    (:tempgradAE, "v0.4.0", :_tempgradhabitat,
    "GradientSpec"),
    (:peakedgradhabitat, "v0.5.0",
    :_peakedgradhabitat, "PeakedSpec"),
    (:peakedgradAE, "v0.4.0", :_peakedgradhabitat,
    "PeakedSpec"))
    @eval function $old(args...; kwargs...)
        Base.depwarn("`$($(QuoteNode(old)))` is deprecated (since $($since)); use " *
                     "`GridHabitat(regime = $($spec)(minT, maxT; axis = Temperature), " *
                     "supply = …, area = …)` instead. (This shim preserves the old behaviour: a " *
                     "shim preserves the old behaviour.)", $(QuoteNode(old)))
        return $body(args...; kwargs...)
    end
end

# ---------------------------------------------------------------------------
# Member-by-name accessors (v0.5.0): `getregime(r, ::Symbol)` / `getpref(t, ::Symbol)` → `getproperty`
#
# These existed because a fixed-arity collection's members were reachable no other way. A collection
# now forwards property access to its backing, so `regime.rainfall` (or `regime.two` on a positional
# one) reaches a member directly, `regimes`/`tolerances` give them all, and `named_regimes`/
# `named_tolerances` give them with their names. Both shims were `getfield(x, name)` verbatim, so
# `getproperty` is exactly what they meant.
#
# Only the `::Symbol` methods are deprecated. Their numeric siblings mean something entirely
# different and stay: `getregime(regime, pos::Int64)` is the regime *value* in grid cell `pos`, and
# `getpref(tolerance, sp::Int64)`/`getdist(tolerance, sp::Int64)` are keyed by *species* — all three
# are hot-loop calls in `_suitability`. Positional member access is `values(lc)[i]`, on the plain
# tuple, which is where indexing a collection is meant to be spelled out.
#
# `getrelationship(nichefit, ::Symbol)` was the third of the set. It was never exported, so it owes
# no shim and has simply gone; use `nichefit.two` or `NamedTuple(nichefit)`.
# ---------------------------------------------------------------------------

"""
    getregime(regime::AbstractRegime, name::Symbol)

!!! warning "Deprecated"
    Deprecated and will be removed; use `getproperty` — `regime.rainfall`, or `regime.two` on a
    positional collection — or `values`/`NamedTuple` for all of them at once.
    Note `getregime(regime, pos::Int64)` is **not** deprecated and means something different: the
    regime value in grid cell `pos`.
"""
function getregime(regime::H, name::Symbol) where {H <: AbstractRegime}
    Base.depwarn("`getregime(regime, ::Symbol)` is deprecated; use `regime.$name`, or " *
                 "`values(regime)`/`NamedTuple(regime)` for every layer at once.",
                 :getregime)
    return getproperty(regime, name)
end

"""
    getpref(tolerance::AbstractTolerance, name::Symbol)

!!! warning "Deprecated"
    Deprecated and will be removed; use `getproperty` — `tolerance.rainfall`, or `tolerance.two` on
    a positional collection — or `values`/`NamedTuple`. Note
    `getpref(tolerance, sp::Int64)` is **not** deprecated and means something different: species
    `sp`'s own preference.
"""
function getpref(tolerance::T, name::Symbol) where {T <: AbstractTolerance}
    Base.depwarn("`getpref(tolerance, ::Symbol)` is deprecated; use `tolerance.$name`, or " *
                 "`values(tolerance)`/`NamedTuple(tolerance)` for every one at once.",
                 :getpref)
    return getproperty(tolerance, name)
end

# ---------------------------------------------------------------------------
# Synthetic layer generators (v0.5.0): `randomniches` → `NicheSpec`
#
# The same treatment the gradient generators got: a layer is described by a *spec* and materialised
# onto a decided grid, so a bare generator returning a loose regime has no way to take part. Unlike
# `tempgrad`/`raingrad`, though, the body does not move here — `NicheSpec` materialises through it
# (`_syntheticregime`, `rasters.jl`; `_syntheticniche`, `StudyArea.jl`), so it stays live in
# `Layer.jl` as `_randomniches` and this is a shim over it.
# ---------------------------------------------------------------------------

"""
    randomniches(dimension, types, clumpiness, weights, gridsquaresize)

!!! warning "Deprecated"
    `randomniches` is deprecated and will be removed. Describe the layer as a recipe instead —
    `GridHabitat(regime = NicheSpec(numniches, axis = …), …)` — which places it on the study area's
    decided grid rather than at a dimension passed in by hand.
"""
function randomniches(dimension::Tuple,
                      types::Vector{Int64},
                      clumpiness::Float64,
                      weights::Vector,
                      gridsquaresize::Unitful.Length)
    Base.depwarn("`randomniches` is deprecated; use `GridHabitat(regime = " *
                 "NicheSpec(numniches, axis = …), …)` instead.", :randomniches)
    return _randomniches(dimension, types, clumpiness, weights, gridsquaresize)
end

# ---------------------------------------------------------------------------
# Naming standardisation (v0.5.0): function names that were `CamelCase` (reading as types) or ran
# an acronym into lowercase are moved to `snake_case`. The symbol-form `@deprecate` forwards
# args/kwargs and (re-)exports the old name, so callers keep working with a deprecation warning.
# ---------------------------------------------------------------------------
@deprecate ContinuousEvolve continuous_evolve
@deprecate DiscreteEvolve discrete_evolve
@deprecate emptyMPIgridlandscape empty_mpi_gridlandscape

# ---------------------------------------------------------------------------
# Layer time series (v0.5.0): a layer no longer holds a stack and a cursor into it
# ---------------------------------------------------------------------------
# A layer holds one grid of values; which stored slice that is comes from elapsed time, via a
# `SeriesLayerChange`. The `*Time*` layer types are therefore gone — a time-varying layer is the same
# type as a static one — and the constructors that took `(stack, time)` build the pair instead.

# Build the 2-D layer + series a v0.4.0 3-D layer becomes. The old `time` cursor is honoured rather
# than dropped: the layer holds the slice it pointed at, and the series is anchored so that elapsed
# time zero selects that same slice, which is what the cursor meant.
# Declares `UndatedSeries` explicitly rather than relying on the default. A v0.4.0 stack is a bare
# cursor into slices with no calendar identity whatever — that is exactly why `origin` is the right
# knob for it, and `origin` is accepted only for an undated series. Saying so here also means a later
# change to what an unmarked source infers cannot silently break the shim.
function _legacyseries!(layer::AbstractLayer, stack, time::Integer)
    origin = _seriescoords(stack, nothing, UndatedSeries()).times[time]
    setchange!(layer,
               ReplaceWith(SeriesChange(stack, origin = origin,
                                        atend = RepeatAtEnd(),
                                        calendar = UndatedSeries())))
    return layer
end

function ContinuousRegime(matrix::AbstractArray{C, 3}, time::Integer,
                          size::Unitful.Length,
                          change::AbstractLayerChange) where {C}
    Base.depwarn("`ContinuousRegime(stack, time, size, change)` is deprecated (since v0.5.0): a " *
                 "regime no longer holds a stack of slices and a cursor into it, but the current " *
                 "slice plus a `SeriesLayerChange` over the stack, indexed by elapsed time. This " *
                 "builds that pair; `change` is ignored, since the series is the change.",
                 :ContinuousRegime)
    stack = _asdimstack(matrix)
    # **The slice is given real coordinates, unlike the deprecated *supply* path.** A `NoLookup`
    # slice would be refused by `_checkhascoords`, and here it need not be: this constructor is handed
    # a `size`, so the coordinates can be derived exactly as `[CELL-DO]` 4a does for every other
    # `Matrix`-built regime. Only the deprecated supply stays coordinate-less, because a v0.4.0
    # budget carried no size to derive from.
    layer = ContinuousRegime(_sizedyx(parent(stack[:, :, time]), size),
                             NoLayerChange())
    return _legacyseries!(layer, stack, time)
end

# Keyed on the **axis**, not on a supply name, since there are no per-resource supply names left:
# `SolarTimeBudget`/`WaterTimeBudget` are bindings onto `Supply{SolarRadiation}`/`Supply{Precipitation}`,
# so these are the methods those bindings resolve to. The old name is still what the warning says,
# because that is what the caller wrote.
for (axis, old) in ((:SolarRadiation, :SolarTimeBudget),
    (:Precipitation, :WaterTimeBudget))
    @eval function Supply{$axis}(mat::AbstractArray{<:Number, 3},
                                 time::Integer)
        Base.depwarn("`$($(QuoteNode(old)))(stack, time)` is deprecated (since v0.5.0): a " *
                     "supply no longer holds a stack of slices and a cursor into it, but the " *
                     "current slice plus a `SeriesLayerChange` over the stack, indexed by elapsed " *
                     "time. This builds that pair.", $(QuoteNode(old)))
        stack = _asdimstack(mat)
        return _legacyseries!(_legacysupply($axis, stack[:, :, time]), stack,
                              time)
    end
end

# A bare 3-D `Array` carries no grid, so it is wrapped on a fresh `NoLookup` `(Y, X, Ti)` — the
# synthetic case, exactly as the 2-D constructors do. A `DimArray` is already on its own grid.
_asdimstack(stack::DimensionalData.AbstractDimArray) = stack
function _asdimstack(stack::AbstractArray{<:Any, 3})
    return DimArray(stack, (Y(NoLookup()), X(NoLookup()), Ti(NoLookup())))
end

# ---------------------------------------------------------------------------
# Ecosystem-level scenarios (v0.5.0): callbacks become declarations
# ---------------------------------------------------------------------------
# `AbstractScenario` carried a *function reference*, so nothing about a scenario could be
# dispatched on, validated, reported or reproduced — and `RandHabitatLoss!`/`ClustHabitatLoss!` drew
# from the **global** RNG, which meant no run using one was reproducible and MPI ranks silently
# diverged. `Intervention(schedule, region, operation)` states the same three things as types.
#
# **The old types are GONE** (v0.5.0), not deprecated: `AbstractScenario`, `SimpleScenario`,
# `FluctScenario`, `MultiScenario` and `runscenario!` are removed outright, with no shim. There was
# nothing to point one at — a callback cannot be mechanically turned into a declaration, since what
# it does is invisible until it runs, and `Intervention` is a plain struct whose parts are the
# abstract types (`AbstractSchedule`, `AbstractRegion`, `AbstractOperation`), so the new mechanism
# has no single supertype standing where the old one stood.
#
# The `scenario` argument is gone from `simulate!`, `simulate_record!`, `simulate_record_diversity!`
# and `simulate_action!` too. Use the `intervention` keyword.

# ---------------------------------------------------------------------------
# The accessor interface (v0.5.0): five ambiguously-indexed names become a 2×2
# ---------------------------------------------------------------------------
# Each of the old names took a second argument whose *meaning was invisible in the name*, and they
# did not all mean the same thing: `getregime(r, 3)` meant cell 3, `getpref(t, 3)` meant species 3.
# That became actively dangerous once collections gained members, since `getregime(r, 3)` acquired a
# second plausible reading ("layer 3"). The replacements say what they are indexed by.

"""
    getdispersaldist(eco, sp)

Deprecated. Use [`EcoSISTEM.speciesdispersal`](@ref), which returns the dispersal **kernel** itself
rather than one of its parameters — more honest (the distance distribution is Rayleigh, not Normal)
and directly reusable as [`AddSpecies`](@ref)' `dispersal` keyword.
"""
function getdispersaldist(eco::AbstractEcosystem, sp)
    Base.depwarn("`getdispersaldist` is deprecated; use `EcoSISTEM.speciesdispersal(eco, sp)`, " *
                 "which returns the kernel itself.", :getdispersaldist)
    return getkernels(eco.spplist.movement)[_speciesindex(eco, sp)].dist
end

"""
    getdispersalvar(eco, sp)

Deprecated. Use [`EcoSISTEM.speciesdispersal`](@ref).

Its formula did not match the kernel it claimed to describe — see the deferred item on dispersal
parameterisation in the master plan — which is part of why the replacement returns the kernel and
lets the caller read whatever it actually holds.
"""
function getdispersalvar(eco::AbstractEcosystem, sp)
    Base.depwarn("`getdispersalvar` is deprecated; use `EcoSISTEM.speciesdispersal(eco, sp)`, " *
                 "which returns the kernel itself.", :getdispersalvar)
    k = getkernels(eco.spplist.movement)[_speciesindex(eco, sp)]
    return (k.dist)^2 * pi / 4
end

# **Hand-written rather than `@deprecate`, and that is the whole point**: the old `getregime` takes
# a *linear* position and returns a **bare value**, while `cellregime` takes `(y, x)` and returns a
# `NamedTuple`. A mechanical redirect would silently change both the indexing and the return type, so
# the shim owns the conversion and keeps the old return.
"""
    getregime(regime::AbstractRegime, pos::Int64)

Deprecated. Use [`EcoSISTEM.cellregime`](@ref)`(eco, y, x)`, which is indexed by cell coordinates
rather than by a linear position and returns a `NamedTuple` naming each regime layer.
"""
function getregime(regime::AbstractLayer{Condition}, pos::Int64)
    Base.depwarn("`getregime(regime, pos)` is deprecated; use " *
                 "`EcoSISTEM.cellregime(eco, y, x)`, which is indexed by cell coordinates and " *
                 "names each layer. This shim keeps the old bare return.",
                 :getregime)
    height = _getdimension(regime)[1]
    (y, x) = convert_coords(pos, height)
    return _cellvalue(regime, y, x)
end

# ---------------------------------------------------------------------------
# Flatcase conformance, 2026-08-25
# ---------------------------------------------------------------------------
# Public function names are flatcase unless flattening hurts readability, matching the packages
# EcoSISTEM sits under: Diversity, Phylo and EcoBase are together about 92% flatcase, and Julia's own
# style guide lumps words together without underscores, reserving them for names that would otherwise
# be hard to read. Twenty names moved; these four are the ones that were exported in v0.4.0 and so owe a
# shim. The rest were unexported or new on this branch, which owes nothing.
#
# Each is a pure rename with an unchanged signature, so a mechanical redirect is exactly right here --
# contrast `getregime` above, where the shim has to own a conversion.
@deprecate assign_traits! assigntraits!
@deprecate gather_abundance gatherabundance
@deprecate gather_diversity gatherdiversity
@deprecate get_traits gettraits

# ---------------------------------------------------------------------------
# The deprecated names' own export/public declarations
# ---------------------------------------------------------------------------
# Declared here rather than in `src/EcoSISTEM.jl` (v0.5.0), so that file's declarations read as a
# map of where each LIVE name is defined. Every name below is a shim defined in this file.

export Gauss, Trapeze, Unif

export tempgrad, raingrad

export simplenichehabitat, tempgradhabitat, raingradhabitat, peakedgradhabitat,
       simplehabitat, erahabitat, worldclimhabitat, bioclimhabitat,
       landcoverhabitat

export getdispersaldist, getdispersalvar, resetrate!

# The v0.4.0 change functions these replace were exported, so their shims stay exported too; they
# are defined in `deprecations.jl`, included below.
export TempChange, RainfallChange, TempFluct, cyclic_change, eraChange,
       worldclimChange

export randomniches

export GaussTrait

# ---------------------------------------------------------------------------
# The netCDF climate sources became SOURCES rather than CONTAINERS (v0.5.0)
# ---------------------------------------------------------------------------
# `ERA`, `CERA` and `CRUTS` used to be struct wrappers holding an array. They are now fieldless
# `EcoSISTEMSource` markers, and a reader returns a `ClimateRaster{ERA}` — the same shape every other
# data source already had. That is what makes them usable in a `SourceSpec`.
#
# The constructor call is redirected, so `ERA(array)` keeps working with a warning.
# 🔴 What CANNOT be shimmed is a **type test or a dispatch** — `x isa ERA` was true of a container
# and is now true of nothing, and `f(::ERA)` no longer matches a read result. Both are recorded in
# `NEWS.md` as breaking; there is no shim that can make a type mean two things at once.
@deprecate ERA(array::DimensionalData.AbstractDimArray) ClimateRaster(ERA,
                                                                      array)
@deprecate CERA(array::DimensionalData.AbstractDimArray) ClimateRaster(CERA,
                                                                       array)
@deprecate CRUTS(array::DimensionalData.AbstractDimArray) ClimateRaster(CRUTS,
                                                                        array)
