# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The abiotic environment as the hot loop sees it: a regime, a supply and an `active` mask, all on
# one `(Y, X)` grid, plus the study area that says where that grid is.

using Diversity

using Diversity.API

using EcoBase

using Unitful

using Unitful.DefaultSymbols

using DimensionalData

using EcoSISTEM.Units

# **The unions are in the signature so that `methods(GridHabitat)` shows what it takes** — see
# [`LayerInput`](@ref), including the cost that buys (a `TypeError` rather than a message with a
# remedy, since Julia does not dispatch on keyword types).
# **Both roles admit an already-built layer**, and for the same reason: a `ContinuousLayer{Condition,
# A}` and a `Supply{A}` each carry their axis in their own type, so each has said what it means —
# which is exactly what a bare raster cannot do. Only on a *data-driven* area: a built layer has to
# have been built on the grid it is passed to, and a synthetic area's builder generates its layers
# rather than accepting them.
"""
    GridHabitat(; regime::Union{LayerInput, AbstractRegime},
        supply::Union{LayerInput, AbstractSupply}, area::StudyArea,
        topology::AbstractTopology = Island())

Build an abiotic environment on an already-decided grid.

All three inputs are **required**. `regime` (the environmental **Condition** layer) and `supply`
(the **Resource** layer) are each an [`AbstractSpec`](@ref) — a synthetic layer
([`UniformSpec`](@ref)/[`GradientSpec`](@ref)/[`PeakedSpec`](@ref)/[`NicheSpec`](@ref)), a
[`SourceSpec`](@ref), a [`ConstructedSpec`](@ref), or a tuple of 2–3 of these for a multi-variable
environment. `area` is the [`StudyArea`](@ref) that decided the grid — its CRS, extent, cell size and
active mask — which is where every geometric choice lives:

```julia
area = StudyArea(regime = cultivated, within = coastline, crs = EPSG(27700), cellsize = 5km)
env  = GridHabitat(regime = cultivated, supply = sunlight, area = area)
```

Use [`investigate_study_area`](@ref) first to see the grid, what it costs each layer and what it
warns about, before committing to it. Omitting a required input errors — use
[`build_habitat`](@ref) to fill omissions with announced defaults.

**The area fixes the frame; the layers can only shrink activity.** Nothing here chooses a grid: this
samples `regime` and `supply` onto the area's, and where a layer has no data the cell is marked
inactive. A layer that was not named when the area was decided can therefore remove cells but never
move or resize the grid — so the area you inspected is the area you get — and doing so warns. Name a
layer in `StudyArea(...)` if you want it to shape the grid; omit it and it can only subtract.

**This is the type's only constructor.** A habitat is described by the layers it should have and
the grid they sit on; there is no form taking already-assembled parts, so `new` is reachable from
exactly one place.

# Fields

  - `regime`: the **Condition** layer(s) of type `H` — what each cell is *like*.
  - `supply`: the **Resource** layer(s) of type `B` — what each cell *provides*.
  - `active`: a boolean `(Y, X)` `DimArray` marking which cells are simulated.
  - `topology`: how the grid's edges join — see [`EdgeTopology`](@ref).
  - `area`: the [`StudyArea`](@ref) this habitat sits on, carrying both the
    [`StudyAreaReport`](@ref) it was built from — stamped [`AsBuilt`](@ref) and refined by whatever
    the build changed — and the [`StudyGrid`](@ref) it was actually built on. The report's
    `active` is the mask the habitat *started* with, so comparing it against the live `active` above
    says what the simulation has done since.
"""
struct GridHabitat{H, B, L} <: AbstractHabitat{H, B, L}
    regime::H
    supply::B
    active::DimensionalData.AbstractDimArray{Bool, 2}
    # **The grid's topology lives here, beside `active`** — the other whole-grid fact the dispersal
    # path consults. It belongs to the map rather than to a species: two species on one grid cannot
    # disagree about whether that grid wraps.
    # Defaults to `Island` on every grid: a study area is a window on a larger world, so an
    # individual dispersing past its margin has left the simulation rather than reappearing opposite.
    # `StudyArea` is where a caller chooses otherwise, and where a wrapping choice on a real-world
    # grid earns its warning; this constructor is reached directly only by synthetic and deprecated
    # paths, which have no geography to misstate.
    topology::AbstractTopology
    # **The grid this habitat sits on, and how it was built — one field, because they are one
    # fact.** The `StudyArea` the caller passed belongs to them and may be gone; and the build itself
    # *changes* things that area could not know — a layer passed here but never named when the area
    # was decided can cost cells it listed as active. So this is a **new** area, carrying that report
    # with the changes written over it (`AsBuilt`) beside the [`StudyGrid`](@ref) actually built on.
    # **`area.report.active` is not the same array as `active` above, and that is the point.**
    # `active` is **live**: an intervention deactivates cells as the run proceeds. The report's is the
    # mask this habitat *started* with, so the difference between the two is the record of what the
    # simulation has done — which no other object can report.
    # It costs ~0.125 bytes/cell: one `BitMatrix` over the grid's own coordinates, read cache
    # dropped. The grid beside it adds nothing measurable — it holds the *same* dimension objects
    # `active` is indexed by, not a copy of them.
    # **`L` is what makes this a typed field rather than an `Any` one**, and it is why `StudyArea`
    # is parameterised at all: it flows straight into `AbstractPartition{L}`, so `getcoords(habitat)`
    # infers concretely instead of inheriting the report's deliberately abstract fields.
    area::StudyArea{L}
    # **The SOLE constructor, and it builds from specs and a `StudyArea`.** There is no
    # parts-taking form: a habitat is described by what its layers should *be* and the grid they sit
    # on, and `new` is reachable from exactly one place. The auto-generated positional constructor is
    # suppressed by this one existing, so `GridHabitat(regime, active, supply)` does not exist either.
    #
    # `area` can be annotated because a type annotation is resolved when the method is defined, and
    # `StudyArea`'s **struct** is declared before this file. Everything that *decides* an area still
    # lives below.
    #
    # Julia does not dispatch on keyword types, so a wrong `regime`/`supply` gives a `TypeError`
    # naming the keyword rather than a message with a remedy — an accepted trade for the signature
    # documenting itself. A bad *element* inside a tuple is still caught later, with the fuller
    # message.
    function GridHabitat(;
                         regime::Union{LayerInput,
                                       AbstractRegime} = _require(:regime),
                         supply::Union{LayerInput,
                                       AbstractSupply} = _require(:supply),
                         area::StudyArea = _requirearea(),
                         topology::AbstractTopology = Island())
        _checktopology(topology, area)
        # `_desugarsupply` here as well as in `StudyArea`, because a supply need not have been
        # named when the area was decided — this accepts one the area never saw. Doing it in both
        # places is idempotent, and the difference is only whether the divisor got to shape the grid.
        parts = _buildonarea(regime, _desugarsupply(supply), area, topology)
        reg, sup = parts.regime, parts.supply
        yx = _yx(reg)
        act = parts.active isa DimensionalData.AbstractDimArray ? parts.active :
              DimArray(parts.active, yx)
        # These four were the old parts-constructor's guards on what a *caller* handed in. They are
        # kept as internal consistency checks: they are cheap, they run once per habitat, and they are
        # the only thing between a builder bug and a silently wrong grid.
        _yxcompatible(dims(act, (Y, X)), yx) ||
            error("`active`'s grid does not match the regime's grid.")
        _yxcompatible(_yx(sup), yx) ||
            error("`supply`'s grid does not match the regime's grid.")
        countsubcommunities(reg) == countsubcommunities(sup) ||
            error("Condition and supply must have same dimensions")
        # Bounds before gaps: `_zerogaps` writes zeros, which would mask nothing here (a `NaN`
        # compares false against any bound) but would make the error message describe a matrix the
        # caller never wrote. A negative supply is the caller's data or spec being wrong.
        cleaned = _zerogaps(_checksupplybounds(sup))
        # **After `_zerogaps`, matching what the builder did before this move** — see
        # `_buildonarea` for why the ordering is load-bearing rather than incidental.
        _applydeclared!(cleaned, parts.supplychange)
        # **The grid is taken from `act`, not from the area's report, and the two can differ.**
        # `act` is the mask this habitat is actually indexed by; the report's is what the area
        # proposed. They share their dimensions today, but reading the habitat's own array is what
        # makes that a property rather than a promise — and it is the same rule as everywhere else
        # here, that a grid is read off the coordinates rather than stated beside them.
        grid = StudyGrid(area.report.crs, act)
        # A **new** area, not the caller's: theirs was not built on, and this one may be narrower.
        built = StudyArea(_refinedreport(area.report, parts.active,
                                         parts.problems), grid)
        return new{typeof(reg), typeof(cleaned), typeof(grid)}(reg, cleaned,
                                                               act, topology,
                                                               built)
    end
end

# Length of a degree of latitude (≈ constant); a degree of longitude is this
# times cos(latitude).
const LONGITUDE_DEGREE_LENGTH = 111.32km / °

# The keyword spelling of a topology, whose types are declared in `src/Topology.jl` — a topology is a
# property of the map rather than of a species, so it is a `GridHabitat` field, and the types must be
# declared before this file for that to be possible.

# **This does not remove the positional spelling, and cannot.** `Torus` and friends are aliases of
# `EdgeTopology{…, …}`, so `Torus()` *is* the parametric constructor and `EdgeTopology{Bounded,
# Periodic}()` stays writable. This is the safer default spelling, not an enforcement — an accepted
# departure from "one spelling per construction", forced by the aliases rather than chosen.
# Takes the boundary conditions as **types**, matching how axes are already passed
# (`axis = Temperature`), rather than as instances.
function EdgeTopology(; y::Type{<:AbstractBoundaryCondition},
                      x::Type{<:AbstractBoundaryCondition})
    return EdgeTopology{y, x}()
end

# ══ GridHabitat ════════════════════════════════════════════════════════════════════════════════════
# **A habitat prints as a summary, not as its type.** Without this it falls back to Julia's
# struct dump, which for a type parameterised over `DimArray`s means thousands of characters of
# nested parameters before any content — nearly eight thousand characters for a 4 by 6 grid, which is
# unusable at exactly the size where one would want to look.
#
# **The line no other object can print is the `active` one.** `active` is live and `report.active`
# is the mask this habitat was built with, so their difference is what the simulation has done —
# cells an intervention has deactivated since. That is the whole argument for keeping the report.
#
# **Two methods, deliberately.** The compact `show(io, x)` is what runs when a habitat appears
# *inside* something else (an `Ecosystem`, a container), and it must stay to one line; the
# `MIME"text/plain"` one is the REPL display. `StudyAreaReport` aliases the two, which is wrong for
# anything nestable — and a habitat is nestable.
function Base.show(io::IO, h::GridHabitat)
    ny, nx = Base.size(h.active)
    return print(io,
                 "GridHabitat($(ny) × $(nx), $(count(h.active)) active)")
end

function Base.show(io::IO, ::MIME"text/plain", h::GridHabitat)
    ny, nx = Base.size(h.active)
    live, built = count(h.active), count(h.area.report.active)
    println(io, "GridHabitat")
    println(io, "  grid      $(ny) × $(nx) cells of ",
            _cellsizetext(h.area.report.cellsize), ", ",
            _crsname(h.area.report.crs))
    lost = built - live
    println(io, "  active    $(live) of $(length(h.active)) cells",
            iszero(lost) ? " (as built)" :
            " — as built: $(built), $(lost) lost since")
    println(io, "  regime    ", _layersummary(h.regime))
    println(io, "  supply    ", _layersummary(h.supply))
    print(io, "  topology  ", _topologyname(h.topology))
    return nothing
end

"""
    totalsupply(habitat::GridHabitat)

Return the **total** resource of a [`GridHabitat`](@ref) — one number per resource, not a per-cell
grid — as a `NamedTuple`, one entry per
resource, keyed by the supply's own **axis names** — `(SolarRadiation = …, Precipitation = …)` for a
multi-resource habitat, and `(SolarRadiation = …,)` for a single one.
Each total is over the habitat's **active** cells only: an inactive cell's resource cannot be used by
anything, so counting it would overstate what the landscape can support.
"""
function totalsupply(habitat::GridHabitat)
    # The mask is applied here rather than being baked into the supply. A habitat keeps the values
    # of its inactive cells so that a cell can later be reactivated — see `_zerogaps` for why a zeroed
    # supply is worse than useless — which makes this the one place that has to do the masking.
    #
    # `NamedTuple(x)` on the container gives every member by name — a single layer answers as a
    # one-member container, so one `map` covers
    # (single, collection)
    # to a `NamedTuple`, so one `map` covers them and the result is keyed the same way every other
    # supply accessor is. It is also concretely typed, which a `Vector` of mixed units could not be:
    # solar `kJ/day` beside water `L/day` can only meet at an abstract `Quantity{Float64}`.
    #
    # Indexes the `parent`: a `DimArray` indexed by a plain `Matrix{Bool}` has no dimensions to
    # match the mask against.
    active = parent(habitat.active)
    return map(NamedTuple(habitat.supply)) do supply
        return sum(parent(supply.matrix)[active])
    end
end

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

# Turn a per-area *rate* (an areal flux — energy/water/carbon per unit area per unit time) into
# an absolute Resource quantity (per cell, not per area) by multiplying by that cell's `area`, and
# state the result in the **axis's** canonical resource unit. This is the last step of the supply
# path, and asking the axis is the same rule as everywhere else — so a new resource axis needs no new
# method here at all. A two-argument, dimension-dispatched form also exists, for the deprecated
# callers that cannot name an axis; it lives in `deprecations.jl`, with them.
# `a::Number` rather than `a::Quantity`: an areal **space** value is a fraction, `m²/m²` and so
# `NoDims`, which Julia represents as a bare `Float64`. Every other areal rate carries a unit, which
# is why the wider annotation is needed.
function cancel(a::Number, b::Quantity{<:Real, 𝐋^2}, axis::Type{<:NicheAxis})
    return _canonicalresource(a * b, axis)
end

# The grid's own cell size, in whatever unit its coordinates are in — `sqrt(Δlat · Δlong)`, and
# nothing else. A projected grid gives a length; a **geographic** grid gives an **angle**, which is
# what its cells genuinely are.
#
# **Distinct from `_cellsize` above, deliberately.** That one converts degrees to kilometres with a
# hardcoded degree length and a cosine at the mean latitude — an *implicit equal-area projection*,
# chosen on the caller's behalf. It predates the package being able to project at all, and it is
# still right for the two places that need a real physical area (`_cellareas`, and the deprecated
# v0.4.0 constructors that reproduce pre-projection behaviour). It is **wrong for a layer's
# `size`**, which should say what the grid is, not what a simulator would like it to be — and a
# geographic grid cannot be simulated anyway (`build_ecosystem` refuses it).
function _gridcellsize(lats, longs)
    return sqrt(abs(lats[2] - lats[1]) * abs(longs[2] - longs[1]))
end

function _gridcellsize(A)
    return _gridcellsize(parent(DimensionalData.lookup(A, Y)),
                         parent(DimensionalData.lookup(A, X)))
end

# The ratio of a cell's *true* area to the nominal one `_cellsize` computes, as a column indexed by
# latitude. On a projected grid it is exactly 1 — a scalar, so nothing downstream grows an
# array — and on a geographic grid it is
#
# f(φ) = [sin(φ + Δφ/2) − sin(φ − Δφ/2)] / (Δφ_rad · cos φ_centre)
#
# **Only the area becomes per-cell, never the cell *side*.** `_cellsize` also feeds dispersal
# (`_uniformcellside` → `genlookups`), where a per-cell value would mean a lookup table per cell rather
# than per species. So the supplies become exact and dispersal keeps the scalar approximation,
# which is what `_cellsize`'s own `@info` now says.
#
# Both `R` and the `111.32 km/°` constant cancel out of `f`, so it inherits none of `_cellsize`'s
# approximations and is exact for a sphere (sphere-vs-ellipsoid is ≲0.3% against the 8–20% being
# corrected). It is applied as `nominal × f` rather than rebuilding `R²·Δλ·[sin…]` from scratch, so
# that area and dispersal keep deriving from *one* Earth radius — and they do: `111.32 km/°` *is*
# a radius, 6378.2 km, because `°` is dimensionless in Unitful.
#
# **`f` is NOT exactly 1 on the centre row**, and an earlier version of this comment claimed it
# was. At `φ = φ_centre` the factor is `sin(h)/h` for a half-cell `h`, i.e. `1 − h_rad²/6` — a
# 1.1e-5 shortfall for a 0.95° cell (measured against `sin(h)/h` to 1e-15). That is the *nominal*
# `Δφ·L` overstating a meridional arc, and correcting it is right; it is simply not the property to
# assert.
# **The property that IS exact, and the one to test**: the per-cell areas sum to the exact
# spherical band area — verified to 1e-12 over an 11 × 7 grid.
# And the reason this defect stayed invisible: over that same grid the old scalar was **0.139%
# out in total** while individual cells were **−10.3% to +13.9%** out. Anything checking a total
# would have passed.
#
# The column is returned as `(ny, 1)` so it broadcasts across X — and across the `Ti` of a
# monthly stack — without a full matrix ever being built. It follows the grid's own row order, so a
# north-up and a south-up grid are both correct with nothing to reverse.
#
# **It takes each cell's own `(lo, hi)` rather than a centre.** Rebuilding the edges as
# `centre ± half` would assume the lookup labels its cells by their centres — true under
# `Intervals(Center)` and false under `Intervals(Start)` — so a labelling change would silently shift
# every correction by half a cell. Asking for the interval is blind to that, and **exact**: there is
# no half-step reconstruction from a differenced step.
_areafactor(::AbstractVector{<:Tuple{Unitful.Length, Unitful.Length}}) = 1.0

function _areafactor(intervals)
    # Each cell's own span, and the grid's centre latitude — the latter is what `_cellsize`'s nominal
    # area is stated at, so the factor is relative to the same place.
    lo, hi = minimum(first, intervals), maximum(last, intervals)
    dlat = abs(hi - lo) / length(intervals)
    φcentre = (lo + hi) / 2
    band = [sin(maximum(iv)) - sin(minimum(iv)) for iv in intervals]
    return reshape(band ./ (ustrip(u"rad", dlat) * cos(φcentre)), :, 1)
end

# Real area of every cell of the `(Y, X)`-dimensioned array `A`: a scalar on a
# projected grid, and an `(ny, 1)` column on a geographic one, where a cell's area
# falls towards the poles.
#
# This is what an areal rate is multiplied by to give a per-cell supply, so it is
# the quantity `cancel` needs; [`_cellsize`](@ref) remains the scalar cell *side*
# that dispersal is expressed against.
#
# It reads the array's own dims rather than taking coordinate vectors, so it asks each axis for the
# thing it means — cell **midpoints** for the nominal size, cell **intervals** for the latitude
# correction, and is blind to how the lookup labels its cells.
function _cellareas(A)
    return _cellareasyx(_gridyx(A))
end

# --- The cell side from its coordinates ---------------------------------------------------
# The implicit equal-area projection, and `src/deprecations.jl` depends on this BEHAVIOUR
# surviving: a v0.4.0 budget was a bare `Matrix` with no coordinates, so `_legacysupply` is the one
# place that may STATE a cell size rather than derive one. Everything else reads the grid.

# Area-preserving physical side of a grid cell whose cell-centre latitudes and
# longitudes are `lats` and `longs`. The behaviour depends on the units of the
# coordinates (multiple dispatch):
#
# - **Geographic (angle) coordinates** (degrees): the north–south side is
# `Δlat × 111.32 km`; the east–west side uses the true length of a degree of
# longitude at the grid's **centre** latitude, `Δlong × 111.32 km × cos φ`.
# Returns the geometric mean of the two, and reports (`@info`) how much the
# east–west length varies from the top to the bottom of the grid.
#
# !!! warning "One side for the whole grid, so it is exact only at the centre"
# This is a single scalar, so on a geographic grid it is right at the
# centre latitude and progressively wrong away from it. That is why it is
# never used on its own: both callers square it and then correct it per
# row with `_areafactor`, giving each cell its true area. It is the
# *nominal* side the correction is applied to, not an answer in itself.
#
# **Dispersal does not use it.** A kernel is expressed against
# `_uniformcellside`, which reads the grid's own step — and a geographic
# grid cannot be simulated at all, `_checksimulatable` refusing one before
# an ecosystem exists.
# - **Projected (length) coordinates** (`m`, `km`, …): the grid is already metric,
# so there is no spherical adjustment — the cell side is simply
# `sqrt(Δlat × Δlong)`.
function _cellsize(lats, longs)
    dlat = abs((lats[2] - lats[1]))
    dlong = abs((longs[2] - longs[1]))
    φtop, φbot = maximum(lats), minimum(lats)
    ns = dlat * LONGITUDE_DEGREE_LENGTH
    ew = dlong * LONGITUDE_DEGREE_LENGTH * cos((φtop + φbot) / 2)
    ewtop = dlong * LONGITUDE_DEGREE_LENGTH * cos(φtop)
    ewbot = dlong * LONGITUDE_DEGREE_LENGTH * cos(φbot)
    isapprox(ewtop, ewbot, rtol = 1.0e-2) ||
        @info "East–west cell length varies with latitude across this grid: " *
              "$(round(typeof(1.0km), ewtop, digits = 2)) at $(round(typeof(1.0°), φtop, digits = 1)) (top), " *
              "$(round(typeof(1.0km), ew, digits = 2)) at the centre, " *
              "$(round(typeof(1.0km), ewbot, digits = 2)) at $(round(typeof(1.0°), φbot, digits = 1)) (bottom); " *
              "using the area-preserving cell size $(round(typeof(1.0km), sqrt(ns * ew), digits = 2)) " *
              "as the nominal, which is then corrected per row. Supplies are unaffected — each " *
              "cell is scaled by its own true area — and dispersal does not use this figure at " *
              "all. It is reported because it is the scale the grid is described at."
    return sqrt(ns * ew)
end

# Projected (length) coordinates: already metric, so no spherical adjustment. The
# result keeps its native length unit; `genlookups` makes the cell-size/dispersal
# ratio dimensionless explicitly, so no unit normalisation is needed here.
function _cellsize(lats::AbstractVector{<:Unitful.Length},
                   longs::AbstractVector{<:Unitful.Length})
    return sqrt(abs(lats[2] - lats[1]) * abs(longs[2] - longs[1]))
end

function _cellsize(A)
    return _cellsize(parent(DimensionalData.lookup(A, Y)),
                     parent(DimensionalData.lookup(A, X)))
end

# The same computation reached from the `(Y, X)` dims alone, for callers that hold a grid rather than
# an array. `_cellareas` delegates here so there is one route to a cell's true area, not two: the
# supply constructor and `getcellareas` must not be able to disagree about how much ground a cell
# covers.
function _cellareasyx(yx)
    lats = _dimcentres(yx[1])
    longs = _dimcentres(yx[2])
    return _cellsize(lats, longs)^2 .*
           _areafactor(DimensionalData.Lookups.intervalbounds(DimensionalData.lookup(yx[1])))
end

# Each cell's midpoint along one dim, whatever the lookup labels its cells by -- the locus-blind
# reading `_cellcentres` gives for an array.
function _dimcentres(d)
    return [(lo + hi) / 2
            for (lo, hi) in DimensionalData.Lookups.intervalbounds(DimensionalData.lookup(d))]
end

# Turn a supply's data gaps into zero resource. This is where "no data" finally becomes "no
# resource": supply constructors preserve their `NaN`s so the grid builder can *see* the gaps (they
# mark cells inactive, warn, and trim the grid), and only once that has happened is the information
# safe to discard. It keeps the defensive property that no `NaN` ever reaches the resource arithmetic
# or `populate!`'s supply-weighted placement, where it would throw or silently spread every
# individual onto nothing.
#
# **It deliberately does NOT zero merely-inactive cells** — that is the `active` mask's job, and
# `active` is read live. Zeroing them here destroyed the supply of any cell that might later be
# *reactivated*, and a zero supply is not inert: `death_resource` is `E/K`, so `K = 0` in a cell that
# holds individuals gives `Inf` and a death probability of exactly 1. A reactivated cell would not
# merely fail to support life, it would kill everything that dispersed into it, silently. Every
# consumer that needs the mask already applies it for itself — the hot loop gates on
# `active && totalE > 0` before reading any supply, `populate!`/`repopulate!` (and their MPI
# counterparts) re-derive `activity` and zero their own working copy, and `totalsupply` passes
# the mask in. (R3/R6, decided with the user 2026-08-05.)
#
# Rebuilds rather than mutating: `ContinuousLayer` is mutable, so zeroing in place would reach back
# into the caller's own supply object.
function _zerogaps(supply::AbstractLayer{Resource})
    (any(isnan, supply.matrix) || _hasgaps(supply.change)) || return supply
    # Applied to the layer's own matrix *and* to any slices its change is holding to write later — a
    # series supply keeps the values it will install on future steps, and they carry the same gaps.
    function zeronan!(values)
        values[isnan.(values)] .= zero(eltype(values))
        return values
    end
    cleaned = copy(supply.matrix)
    # Work on the parent array throughout: the values are all this needs to touch.
    zeronan!(parent(cleaned))
    # `typeof(supply)` already carries every type parameter, and `copy` preserves the array type, so
    # the layer is rebuilt without having to name (or risk mis-stating) its role/axis/eltype.
    return typeof(supply)(cleaned, supply.size,
                          _cleanstored(supply.change, zeronan!))
end

# A collection is cleaned sub-layer by sub-layer. `map` over the backing preserves its names, so a
# named supply collection stays named.
function _zerogaps(supply::LayerCollection{Resource})
    return LayerCollection(map(_zerogaps, getfield(supply, :nt)))
end

# **A habitat owns its layers: give it a copy of anything it did not build itself.**
#
# **Why this is correctness and not hygiene.** `_applychange!` writes `layer.matrix .= …` on **every
# timestep**, so a layer stored by reference makes the caller's own object live simulation state —
# and two habitats handed the same pre-built layer share it. `_zerogaps` above already names this
# hazard, but rebuilds only when there are *gaps* to clean, so a gap-free layer passed straight
# through. The spec path was never affected: `materialise` builds fresh arrays on every call.
#
# **A shallow copy, and the change is deliberately shared.** Nothing mutates a change in place:
# `_repointseries!` rebuilds it, `_cleanstored` copies the slices before cleaning, and `setchange!`
# and `_applydeclared!` replace the field rather than writing through it. So the mutable state is
# exactly `matrix` and the `change` **binding**, both per layer, while a deep copy would clone a
# `SeriesLayerChange`'s whole slice stack — twelve months of a full grid — for nothing. Should a
# change ever write `change.slices` in place, this must become a deep copy.
#
# `typeof(layer)` carries every parameter, so the role, axis, eltype and array type are preserved
# without being named — the same trick `_zerogaps` uses.
function _ownlayer(layer::AbstractLayer)
    return typeof(layer)(copy(layer.matrix), layer.size,
                         layer.change)
end

function _ownlayer(coll::LayerCollection)
    return LayerCollection(map(_ownlayer,
                               getfield(coll, :nt)))
end

# `nameof(typeof(Island()))` is `EdgeTopology` — the useful names are `const` aliases for
# parameterisations, so they have to be matched rather than read off the type.
function _topologyname(::EdgeTopology{Bounded, Bounded})
    return "Island — hard edges on all four sides"
end

function _topologyname(::EdgeTopology{Periodic, Periodic})
    return "Torus — both pairs of edges join"
end

function _topologyname(::EdgeTopology{Bounded, Periodic})
    return "Cylinder — east–west edges join"
end

_topologyname(t::AbstractTopology) = string(nameof(typeof(t)))

# One line naming what a layer (or a collection of them) is on: its axis and unit, which is what
# distinguishes two layers of the same shape.
function _layersummary(l::AbstractLayer)
    return "$(nameof(_layeraxis(l))) ($(unit(eltype(l))))"
end

function _layersummary(c::LayerCollection)
    return join(map(_layersummary, values(c)), " + ")
end

# --- Getting the report out of whatever carries it ---------------------------
# A report, a `StudyArea` or a `GridHabitat` — the habitat is the last of the three, so the family
# lives here. `_basis` returns the *inputs* an analysis is re-run from; `_basereport` the report
# itself, and conflating the two is why the copy branch once never fired.

# Resolve a `base` (a `StudyArea`, a `StudyAreaReport`, or nothing) plus keywords into the layers,
# constraints and cache an analysis needs. Both public entry points go through this, so refining an
# area and re-investigating one behave identically.
function _basis(::Nothing)
    return (specs = NamedTuple(), constraints = NamedTuple(),
            cache = nothing)
end

function _basis(base::StudyAreaReport)
    return (specs = base.specs,
            constraints = base.constraints,
            cache = base.cache)
end

# Anything else carrying a report — a `StudyArea` — delegates to it.
# **`base` is deliberately UNANNOTATED, and a signature audit will keep suggesting `::StudyArea`.**
# It is duck-typed on purpose: the three methods cover a report, a habitat and anything else that
# carries a `.report`, so pinning the fallback to `StudyArea` would refuse the very cases it exists
# to catch. The same applies to `_basereport` below.
_basis(base) = _basis(base.report)

# A habitat carries its report one level deeper, inside the `StudyArea` it was built into, so it
# needs a method of its own rather than falling through to the `.report` duck-type above.
_basis(base::GridHabitat) = _basis(base.area.report)

# The report a base carries, whatever face it wears — a report itself, a `StudyArea` or a
# `GridHabitat`. Distinct from `_basis` above, which returns the *inputs* an analysis is re-run
# from rather than the report itself; conflating the two is why the copy branch below first never
# fired.
_basereport(base::StudyAreaReport) = base

_basereport(base) = _basereport(base.report)

_basereport(base::GridHabitat) = base.area.report

# --- What a habitat refuses to guess -----------------------------------------
# `GridHabitat` has one construction route and no defaults, so these are what turn an omitted or
# unusable input into a message naming the remedy. `build_habitat` (`actions.jl`) is the forgiving
# face that supplies a value instead.

# Sentinel default for a required keyword: errors (naming the missing input) when the keyword is
# omitted from a strict `GridHabitat(…)` / `build_species(…)` call. The `DefaultEcosystem`
# builders fill omissions instead.
function _require(field)
    return error("the required keyword `$field` was not passed — pass it explicitly, or use " *
                 "the `DefaultEcosystem()` builder (e.g. `build_habitat(DefaultEcosystem(); …)`) " *
                 "to fill omitted inputs with announced defaults.")
end

# `area`'s own version of `_require`. Which grid to build on is the one required input whose answer is
# not obvious from the keyword's name, so the error says how to produce one rather than only naming it.
function _requirearea()
    return error("the required keyword `area` was not passed: `GridHabitat` builds on a grid " *
                 "that has already been decided, and never chooses one itself. Decide it first — " *
                 "`StudyArea(; regime, supply)` takes the grid from the layers, and " *
                 "`investigate_study_area(...)` shows what it would be before you commit — then pass " *
                 "`area = …`. Or use `build_habitat()` for a toy grid.")
end

# Warn when a **wrapping** axis is chosen on a grid that has a real-world position. The default is
# `Island` on every grid — one answer everywhere, deliberately, because a default that varied with the
# CRS would mean the same call meant different things in different files. What the CRS decides is only
# whether a *deviation* from that default is announced.
#
# A synthetic grid is silent whatever it chooses: with no real-world position there is no geography
# to misstate, and a torus is simply true of it — there is no outside to leave to.
#
# **Stated per axis, which is what makes it exact rather than a special case.** The rule says what it
# means and every exemption falls out of it:
#
#   * a periodic **X** is right only where the grid spans the whole longitude sweep — east-west
#     really does join, there;
#   * a periodic **Y** is *never* right, however wide the grid: latitude does not wrap, and a step
#     past the pole is not a step to the far south.
#
# Stating it per axis also covers the Y-wrapping topology, which no check keyed on the named
# combinations could see.
_checktopology(::EdgeTopology{Bounded, Bounded}, area::StudyArea) = nothing

function _checktopology(::EdgeTopology{BCY, BCX},
                        area::StudyArea) where {BCY, BCX}
    isnothing(area.report.crs) && return nothing        # synthetic: nothing to misstate
    BCY === Periodic &&
        @warn "`topology` wraps the **y** axis on a study area with a real-world position: an " *
              "individual dispersing past the northern edge reappears in the south, which no " *
              "geography does — latitude does not wrap, whatever the grid's extent. The run is " *
              "still built; `Bounded` on y (as `Island()` and `Cylinder()` both use) is the " *
              "geography as it is."
    BCX === Periodic && !_spansglobe(area) &&
        @warn "`topology` wraps the **x** axis on a study area with a real-world position that " *
              "does not span the whole longitude sweep: an individual dispersing past the eastern " *
              "edge reappears in the west, which this grid's geography does not do. The run is " *
              "still built; use `Island()` for the geography as it is. (On a grid that *does* span " *
              "all longitudes, `Cylinder()` is exactly right and is not warned about.)"
    return nothing
end

# Does the grid's longitude span the whole globe? Judged in **degrees**, so a projected grid — whose
# X is a length — answers `false`: a projected CRS cannot wrap, whatever its extent.
#
# The unit test here means what it says: a synthetic grid's coordinates are plainly `km`, so it is
# refused for the same reason a projected one is — its X is a length, and a length cannot wrap.
function _spansglobe(area::StudyArea)
    lo,
    hi = DimensionalData.Lookups.bounds(DimensionalData.lookup(area.report.active,
                                                               X))
    unit(typeof(lo)) == u"°" || return false
    return (hi - lo) ≥ 359.9u"°"
end
