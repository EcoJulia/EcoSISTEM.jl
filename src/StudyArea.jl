# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Where the cells actually are - the CRS and the very `Y`/`X` dimensions a habitat's `active` mask
# is indexed by - and the area that carries both that grid and the report saying how it was decided.

using DimensionalData

using EcoBase

using Unitful

using Rasters

import Rasters: Projected

"""
    StudyGrid{C, YD, XD} <: EcoBase.AbstractGrid

The geometry of the grid a [`GridHabitat`](@ref) was built on: where its cells are, how big they are
and what coordinate reference system they are stated in.

This is the package's [`EcoBase`](https://github.com/EcoJulia/EcoBase.jl) location data, so anything
that speaks that interface - `EcoBase`, `SpatialEcology`, the plot recipes - can ask a habitat where
its cells are. Every answer is **unitful**, in the grid's own units: a projected grid answers in `km`
(or whatever it was built in) and a geographic one in degrees, neither converted into the other.

**It is never constructed by hand and never carries its own copy of the coordinates**: it holds the
very `Y` and `X` dimensions the habitat's `active` mask is indexed by, so the two cannot drift apart.
Reach it through the habitat's [`StudyArea`](@ref).

# Fields

  - `crs`: the coordinate reference system the coordinates are stated in, or `nothing` for a synthetic
    grid, which has none.
  - `y`: the grid's `Y` dimension - the **first** array dimension throughout this package.
  - `x`: the grid's `X` dimension.
"""
struct StudyGrid{C, YD, XD} <: EcoBase.AbstractGrid
    crs::C
    y::YD
    x::XD

    # **The sole constructor, and it validates what the grid must be able to say.** A `NoLookup`
    # dimension is bare cell indices and an unitless one is numbers dressed as coordinates; neither
    # can state where a cell is or how big it is, which is the whole job of this type. The same rule
    # `_checkhascoords` applies to a layer, applied here for the same reason.
    # Cheap enough to be unconditional: it runs once per habitat, never per cell.
    # It cannot fire on anything a `GridHabitat` can currently be built from - the grid is read
    # off the **regime**, and no regime constructor accepts `NoLookup` (`_checkhascoords` refuses it);
    # only the deprecated `Supply{A}(::Matrix)` makes such dims, and a supply is not what is read
    # here. So this guards the type's own invariant rather than a live route, and is exercised
    # directly rather than incidentally - see the `StudyGrid` testsets in `test_Layer.jl`.
    function StudyGrid(crs, y, x)
        for (d, name) in ((y, "Y"), (x, "X"))
            _isreallookup(d) ||
                error("a `StudyGrid`'s `$name` coordinates are `NoLookup` - bare cell indices, " *
                      "which cannot say where the cells are. This grid came from a habitat built " *
                      "without real coordinates; build it on a `StudyArea`.")
            eltype(parent(DimensionalData.lookup(d))) <: Unitful.Quantity ||
                error("a `StudyGrid`'s `$name` coordinates must carry a unit - these are bare " *
                      "numbers, so how far apart the cells are is unstated.")
        end
        return new{typeof(crs), typeof(y), typeof(x)}(crs, y, x)
    end
end

# == Functions ==================================================================================

# The grid an already-built `(Y, X)` array sits on. This is the form every caller in the package
# uses: the array is the habitat's own `active` mask, so the grid is guaranteed to be the one the
# habitat is indexed by rather than a second description of it.
function StudyGrid(crs, active::DimensionalData.AbstractDimArray)
    yx = dims(active, (Y, X))
    return StudyGrid(crs, yx[1], yx[2])
end

"""
    CellNames{G} <: AbstractVector{String}

The names of a grid's cells - each the cell's own extent - computed on demand rather than stored.

**Lazy because the eager form is what a grid this size cannot afford.** Measured on a
1.2 million-cell grid: storing the names costs ~33 MB, against **8 bytes** for this. That matters
here specifically - a `GridHabitat`'s `active` mask was deliberately squeezed from 10.3 MB to
0.14 MB for the same reason, and [`getspeciesstorage`](@ref) exists because these runs are bound by
memory. Nothing in the simulation reads a cell's name; only output does.

# Fields

  - `grid`: the [`StudyGrid`](@ref) the names describe.
"""
struct CellNames{G <: StudyGrid} <: AbstractVector{String}
    grid::G
end

"""
    StudyArea(base = missing; regime::Union{LayerInput, Missing, Nothing} = missing,
        supply::Union{LayerInput, Missing, Nothing} = missing, within = missing,
              crs = missing, cellsize = missing, extent = missing, align = missing,
              simulate_safely = missing, verbosity = :normal)

Decide the grid a simulation will run on, before anything is built on it.

Nothing is compulsory: the grid is decided from whatever is given, so naming only `regime` means only
the regime shapes it. **A layer that is not named here can never move or resize the grid** - when
[`GridHabitat`](@ref) later samples it, it can only mark cells inactive (and says so). Name a
layer to let it shape the grid; omit it and it can only subtract.

The optional positional `base` refines an existing area - `StudyArea(area, cellsize = 1km)` - keeping
everything not overridden, including its cache of reads, so trying several grids does not re-read the
data. A [`StudyAreaReport`](@ref) or a [`GridHabitat`](@ref) may be given in its place.

**A habitat's grid is copied rather than re-derived.** A habitat can end up *narrower* than the
area it was built on, because a layer named only to [`GridHabitat`](@ref) can cost cells the area
listed as active - and nothing but the habitat's own report records that. So a base whose report is
[`AsBuilt`](@ref), with no other keyword given, is taken **verbatim**: `StudyArea(habitat)` is the
grid the habitat is actually on, not the one originally investigated. Naming any grid keyword
alongside re-derives instead, which loses the narrowing; that is deliberate, and not warned about,
because overriding is an explicit act.

**What is copied is the report, not the built grid.** `StudyArea(habitat)` describes the habitat's
grid but is itself a *fresh* area with nothing built on it, so its `builtgrid` is `nothing` while its
report stays [`AsBuilt`](@ref). The two say different things - the report's stage describes what the
*report* is, `builtgrid` what has been done to *this area* - which is why neither can be derived from
the other.

Keywords are tri-state: `missing` (the default) means "not specified - derive it, or take it from
`base`", `nothing` means "explicitly cleared, ignore any inherited value", and any other value is used
as given.

# Arguments

  - `regime`, `supply`: the layers the simulation will use. Naming one here lets it **shape the
    grid** - its extent, its CRS and its resolution all become candidates. A layer not named here can
    still be used later, but can only mark cells inactive, never move or resize the grid.
  - `within`: **what positions the area** - a [`ShapeSpec`](@ref), a [`CircleMaskSpec`](@ref), a
    [`LatLong`](@ref) box, or `EcoSISTEM.boundingbox("Scotland")`. It both restricts which cells are
    active and, where it can state an extent, sets the grid's. This is the argument to reach for
    when a global dataset would otherwise give you the globe.
  - `crs`: the coordinate reference system to work in. **A projected CRS is required to simulate**:
    dispersal is expressed against one cell size, which only a projected grid has. Omitted, it is
    adopted from the layers, and a geographic result is warned about.
  - `cellsize`: the length of a cell's side. Omitted, it is taken from the `align` layer's own step.
  - `extent`: **a size, not a bounding box** - how big the area is, never where it is, given as a
    `(y, x)` tuple of lengths: north-south first, then east-west, the dimension order used
    throughout the package. It is therefore only meaningful for a *synthetic* area, and combining
    it with data layers is an **error**: those already carry their own extent, and a second
    unanchored size cannot be reconciled with it. Use `within` to position and limit an area built
    from data.
  - `align`: names the layer whose own grid is preserved exactly. By default whichever layer is
    already in the target CRS (finest first), since that one needs no reprojection.
  - `simulate_safely`: whether a cell must be **wholly** inside every layer's data to be simulated.
    `true` (the default) is the safe reading - a cell the data only partly describes is dropped, and
    the grid crops inwards rather than carrying it - since nothing can put data where a file has
    none. `false` restores the older rule, under which a cell survives if its centre has a value
    (i.e. if it is more than half covered) and is then given a whole cell's worth of supply over
    ground that is partly missing; that is announced, with a count. The same rule applies to
    layers named only later, at [`GridHabitat`](@ref).
  - `verbosity`: see below.

With no layers that can shape a grid, `extent` and `cellsize` build a synthetic one with no CRS -
which is also what a *synthetic* layer ([`UniformSpec`](@ref), [`GradientSpec`](@ref), ...) needs, since
it is generated at whatever shape it is handed and so has no CRS, extent or resolution to contribute.
Naming one here is harmless but decides nothing.

`verbosity` is `:silent` (errors only), `:normal` (announce every guessed value and every lossy step,
and warn about grids that will work but are probably not what you want) or `:verbose` (aliases
`:full`/`:debug`; show the whole [`investigate_study_area`](@ref) report).

# Fields

  - `report`: the [`StudyAreaReport`](@ref) recording how this grid was decided - the CRS, cell
    size and extent it settled on, which layer decided each, and every [`Problem`](@ref) raised on
    the way.
  - `builtgrid`: the [`StudyGrid`](@ref) a [`GridHabitat`](@ref) was built on here, or `nothing`
    until one has been. `nothing` does **not** mean the grid is undecided - an investigated area
    has decided one; that is what investigation is. It means nothing has been *built* on it yet, so
    the cells it lists are a prediction: a layer named only to [`GridHabitat`](@ref) can still
    narrow them.
"""
struct StudyArea{G <: Union{StudyGrid, Nothing}}
    report::StudyAreaReport
    builtgrid::G

    # The same unions as `investigate_study_area`, which shares this signature - the two are one
    # analysis behind two faces, so a difference between them could only be a mistake.
    function StudyArea(base = missing;
                       regime::Union{LayerInput, Missing, Nothing} = missing,
                       supply::Union{LayerInput, Missing,
                                     Nothing} = missing,
                       within = missing, crs = missing, cellsize = missing,
                       extent = missing, align = missing,
                       simulate_safely::Union{Bool, Missing,
                                              Nothing} = missing,
                       verbosity::Symbol = :normal)
        # **A COPY CONSTRUCTOR, and the stage is what makes it possible.** An `AsBuilt` report
        # describes a grid a `GridHabitat` was actually built on - its `active` has been narrowed by
        # layers the area never saw, and that narrowing is recorded nowhere else. Re-deriving from
        # `specs`/`constraints` would answer the *original* question and hand back the wider,
        # as-investigated mask: measured at **48 cells against the habitat's 46**, silently
        # discarding the narrowing *and* the `Problem` that recorded it.
        # So an `AsBuilt` base with **no** keyword is taken **verbatim** - that is the whole point
        # of copying rather than rebuilding, and it is why `[AREA-PROV]` Part C added the stage.
        # **The copy keeps `AsBuilt`** (user's decision): the flag says what a report *describes*,
        # and resetting it to protect the connotation that "a `StudyArea` is a proposal" would be a
        # lie that also breaks the copy on its second hop.
        # **Any keyword sends it down the ordinary path**, with no warning: overriding is an
        # explicit act, and someone writing `cellsize =` beside a built base is presumed to mean it.
        if _iscopy(base) &&
           !any(!ismissing,
                (regime, supply, within, crs, cellsize, extent,
                 align, simulate_safely))
            return new{Nothing}(_copyablereport(_basereport(base)),
                                nothing)
        end
        inputs = _resolveinputs(base; regime, supply, within, crs, cellsize,
                                extent, align, simulate_safely)
        # "No compulsory arguments" is not "any combination works": with neither layers nor a
        # synthetic geometry there is nothing whatever to decide a grid from.
        (all(isnothing, inputs.layers) &&
         (isnothing(inputs.constraints.extent) ||
          isnothing(inputs.constraints.cellsize))) &&
            error("a `StudyArea` needs something to work from: either at least one layer " *
                  "(`regime`/`supply`), or both `extent` and `cellsize` for a synthetic grid.")
        report = _analyse(inputs.layers; inputs.constraints...,
                          cache = inputs.cache)
        _emit(report, verbosity)
        return new{Nothing}(report, nothing)
    end

    # **The one route by which an area gains its grid, and it is internal.** `GridHabitat` calls
    # it with the report it refined and the grid it actually built on. It builds a *new* area
    # rather than mutating the caller's: theirs was not built on, and the habitat's grid may be
    # narrower than the one they inspected, so overwriting it would misreport both.
    function StudyArea(report::StudyAreaReport, builtgrid::StudyGrid)
        return new{typeof(builtgrid)}(report, builtgrid)
    end
end

# How much slack to leave around the window, as a fraction of its own span. The grid is snapped
# outwards onto cell boundaries, widened to at least two cells and resampled with a method that reads
# neighbours, so a window cut exactly to the mask would starve those edges. Generous on purpose: the
# saving is already three orders of magnitude, and a window too tight is a wrong answer where a
# window too loose is merely a slower one.
const _WINDOW_PAD = 0.1

# --- The decision core -----------------------------------------------------
# Deliberately written over plain numbers rather than rasters: what a grid choice costs is
# arithmetic on cell sizes and origins, and keeping it separable makes it directly testable without
# constructing rasters. `_analyse` supplies the numbers; `StudyArea` and `investigate_study_area`
# share the answer, so the report can never describe something different from what is built.

# Cell-boundary alignment is judged as a fraction of a source cell, so the tolerance is absolute.
const _ORIGIN_ATOL = 1e-6

# Cell counts beyond which a grid is more likely a mistake than an intention: `_BIG_GRID` cells is
# already ~80 MB per layer and grows per species, and an aggregation factor this large means each
# target cell swallows ~10^6 source cells, where a majority class stops meaning anything.
const _BIG_GRID = 10_000_000

const _HUGE_FACTOR = 1000

# Below this fraction of active cells the grid is mostly dead space - usually a wrong mask or extent.
const _SPARSE_ACTIVE = 0.05

# Warn when the grid sits substantially outside the ground its CRS was designed for - British
# National Grid stretched across continental Europe, or one UTM zone spanning several. Such a
# projection still *computes*, silently, while distorting distance and area badly enough to change
# what a simulation means, so it is exactly the sort of thing that should not pass unremarked.
#
# Judged on the fraction of the study extent falling outside the declared area of use, with a wide
# tolerance: a real study area routinely overhangs its CRS's box a little (BNG's stops at the UK
# coastline), and warning on that would be noise.
const _CRS_OUTSIDE_FRAC = 0.2

# == StudyGrid and CellNames ========================================================================
function Base.show(io::IO, grid::StudyGrid)
    return print(io, "StudyGrid(", length(grid.y), " × ", length(grid.x),
                 " cells of ", EcoBase.ycellsize(grid), " × ",
                 EcoBase.xcellsize(grid),
                 isnothing(grid.crs) ? "" : ", $(_crsname(grid.crs))", ")")
end

# == StudyArea ======================================================================================
function Base.show(io::IO, a::StudyArea)
    ny, nx = Base.size(a.report.active)
    nwarn = count(p -> p.severity isa ProblemWarning, a.report.problems)
    return print(io, "StudyArea($(ny) × $(nx) cells of $(a.report.cellsize)",
                 isnothing(a.report.crs) ? "" : ", $(_crsname(a.report.crs))",
                 ", $(count(a.report.active)) active",
                 nwarn == 0 ? "" : ", $nwarn warning$(nwarn == 1 ? "" : "s")",
                 ")")
end

function Base.show(io::IO, ::MIME"text/plain", a::StudyArea)
    return show(io, MIME"text/plain"(), a.report)
end

"""
    getspeciesstorage(x)

Return the bytes **one** species' abundances occupy on `x`'s grid, or `nothing` where `x` has no grid
to answer for.

- `x` - anything that knows the grid: a [`StudyArea`](@ref), its [`StudyAreaReport`](@ref), a
  [`GridHabitat`](@ref), an [`Ecosystem`](@ref), a layer, a [`ClimateRaster`](@ref) or a raster.
  A spec that is a rule rather than data has no grid and answers `nothing`.

Multiply by a species count to size a run before building it - which is the point, since it can be
asked of a `StudyAreaReport` from [`investigate_study_area`](@ref) before anything is constructed.

It counts the **whole** grid, not the active cells: [`GridLandscape`](@ref)'s matrix is allocated
over all of it, so an inactive cell costs exactly as much as an active one. And it is one array -
a run holds several (abundances and net migration serially, more when distributed), so this is a
per-array figure to be multiplied, not a total.

See also [`getcellsizes`](@ref) and [`getcellareas`](@ref), which ask the same grid other questions.
"""
function getspeciesstorage(x)
    yx = _gridyx(x)
    isnothing(yx) && return nothing
    # Routed through `_gridfootprint` so the report's `footprint.perspecies` and this can never
    # disagree about what a cell costs.
    return _gridfootprint(length(yx[1]), length(yx[2])).perspecies
end

# --- The `EcoBase.AbstractGrid` interface ----------------------------------
#
# Seven methods, all reading the dimensions rather than a stored number, so there is nothing to
# disagree with the coordinates. `cellsize`, `cells`, `xmax`/`ymax` and `xrange`/`yrange` are
# EcoBase's own derivations from these and need nothing here - measured to work unchanged on
# unitful values, angular ones included, which is why widening EcoBase was not needed.

# The rank of array position `i` along `d`, counting in **increasing coordinate** order.
#
# A raster's `Y` commonly runs north to south, so its array rows descend, while EcoBase's `yrange`
# ascends by construction (`ymin:ycellsize:ymax`). Ranking rather than handing the array row straight
# out is what keeps such a grid the right way up once something plots it - the row and the rank agree
# on every ascending grid, which is every grid this package builds, and differ on one that arrives
# from elsewhere.
function _ascendingrank(d, i::Int)
    lk = DimensionalData.lookup(d)
    isreverse = DimensionalData.Lookups.order(lk) isa
                DimensionalData.Lookups.ReverseOrdered
    return isreverse ? length(lk) - i + 1 : i
end

# --- Naming a cell by where it is ------------------------------------------

# Is this coordinate an angle - i.e. is the grid geographic? Unitful makes angles dimensionless, so
# a coordinate with no dimension but a real unit is one; a metre or a kilometre has a dimension, and
# a genuinely unitless coordinate is refused by the constructor above.
_isangular(v) = Unitful.dimension(v) === Unitful.NoDims

# Enough decimal places to tell neighbouring cells apart, taken from the step rather than fixed: a
# 1 km grid needs one, a 30-arcsecond one (0.008333°) needs three. At least one either way, so a
# coordinate always reads as a coordinate rather than as a count.
function _extentdigits(step)
    return max(1, ceil(Int, -log10(abs(Unitful.ustrip(step)))))
end

# One coordinate as it appears in a cell-extent label: stripped of its unit and rounded, because the
# unit is stated once for the whole label rather than on every number in it.
_extentnumber(v, digits) = string(round(Unitful.ustrip(v), digits = digits))

# One cell's half-open extent along one axis, as `[lo, hi)`. `intervalbounds` is locus-blind and
# order-aware, so this is right whether the grid labels its cells by their corner or their centre and
# whether the axis ascends or descends - the reason `_cellintervals` uses it too.
function _extenttext(d, i::Int, digits::Int)
    bounds = DimensionalData.Lookups.intervalbounds(DimensionalData.lookup(d),
                                                    i)
    lo, hi = minmax(bounds...)
    return "[" * _extentnumber(lo, digits) * ", " * _extentnumber(hi, digits) *
           ")"
end

# The name of cell `i` - its own extent, `Y` first, so a cell can be identified from its name alone
# rather than only by counting.
#
# **Two shapes, because the two kinds of grid read differently.** A geographic grid names its axes
# `N` and `E`, which says which is which and which way is positive (a negative longitude is still
# `°E`, just west of the meridian); a projected one shares one unit between both axes and states it
# once, `Y` before `X` as everywhere else in this package.
function _cellname(grid::StudyGrid, i::Int)
    ny = length(grid.y)
    r = ((i - 1) % ny) + 1
    c = div(i - 1, ny) + 1
    dy, dx = _extentdigits(_axisstep(grid.y)), _extentdigits(_axisstep(grid.x))
    ytext, xtext = _extenttext(grid.y, r, dy), _extenttext(grid.x, c, dx)
    yunit = Unitful.unit(first(DimensionalData.lookup(grid.y)))
    xunit = Unitful.unit(first(DimensionalData.lookup(grid.x)))
    _isangular(yunit) && return "$ytext$(yunit)N × $xtext$(xunit)E"
    return "$ytext × $xtext $yunit"
end

# --- The study area --------------------------------------------------------

# == Deciding the grid ==========================================================================
#
# Everything that turns a set of layers into a decided grid: a CRS, an extent, a cell size and
# an active mask, with a report saying how each was chosen. Sampling layers onto the result is
# `materialise.jl`, included next.

# --- Windowing and the grid decision ------------------------------------------

# Refuse a combine over layers that are not on one grid, before the combine sees them.
#
# Combining *reads* is only meaningful where the reads share a lattice. `gsp` and `gsl` do - one
# dataset, one resolution - but a WorldClim layer against an EarthEnv one does not, and handing
# those to a combine either throws something obscure from inside the broadcast or, where the shapes
# happen to agree, lines up cells describing different ground. The second is the reason this is a
# guard rather than a comment.
#
# Checked here rather than only on the `CombineOnSourceGrid` path, because this route is also how
# `_analyse` asks a `ConstructedSpec` where it is - so the hazard predates early collapse, and
# fixing it in one place covers both. `_yxcompatible` is the same agreement test the habitat builder
# applies to a regime and its supply.
function _checksourcegrid(spec::ConstructedSpec, layers,
                          idx = eachindex(layers))
    length(layers) < 2 && return nothing
    yx = dims(first(layers).array, (Y, X))
    for i in 2:length(layers)
        theirs = dims(layers[i].array, (Y, X))
        _yxcompatible(theirs, yx) && continue
        return error("a `ConstructedSpec`'s layers must be on one grid to be combined, but " *
                     "$(_speclabel(spec.layers[idx[i]])) is $(join(length.(theirs), "×")) cells where " *
                     "$(_speclabel(spec.layers[first(idx)])) is $(join(length.(yx), "×")) - they are " *
                     "different ground, so combining them cell by cell would be meaningless. Read " *
                     "them as separate layers of the environment (a tuple `regime`/`supply`), " *
                     "which puts each on the study grid first.")
    end
    return nothing
end

# --- Windowing: read only the ground the study area can possibly use ---------
#
# The point of the exercise: a Scotland-sized area has no use for a global land-cover layer, and
# reading one costs 1752 × 4320 cells where 76 × 96 would do. The obstacle is an ordering one - the
# window must be known *before* the read, but it depends on the target CRS, which is decided from the
# layers. It is broken by asking each source only for its **CRS**, which a lazy open answers from the
# file header without touching a pixel (`EcoSISTEM.sourcecrs`), and which the aggregation `scale`
# cannot affect. Everything else still comes from the real, now much smaller, read.

# A spec's CRS from metadata alone, or `nothing` when it cannot be had without reading - in which
# case the caller simply does not window. A `ConstructedSpec`'s `combine` is an opaque closure whose
# output CRS is not predictable from its inputs, so it declines.
_probecrs(raster::ClimateRaster) = _rastercrs(raster)

function _probecrs(spec::SourceSpec)
    return sourcecrs(spec.source, spec.code; spec.readkw...)
end

_probecrs(::Any) = nothing

# The target CRS decided without reading anything: an explicit `crs` needs no probe at all, and
# otherwise the same staged rule `_targetcrs` applies, fed from headers. `nothing` if any layer
# declines to say, since adopting a CRS from an incomplete picture could pick the wrong one.
function _probetargetcrs(layers::NamedTuple, crs)
    isnothing(crs) || return crs
    specs = [s
             for (_, ss) in ((n, _expandspecs(s))
                             for (n, s) in pairs(layers) if !isnothing(s))
             for s in ss if _shapesgrid(s)]
    isempty(specs) && return nothing
    crss = map(_probecrs, specs)
    any(isnothing, crss) && return nothing
    projected = filter(_isprojectedcrs, crss)
    unique_projected = [projected[i]
                        for i in eachindex(projected)
                        if all(j -> !_samecrs(projected[i], projected[j]),
                               1:(i - 1))]
    return length(unique_projected) == 1 ? only(unique_projected) : first(crss)
end

# Do the target's cell boundaries fall on the source's? Only when they do can the source be copied or
# aggregated exactly; otherwise every target cell straddles source cells and must be interpolated.
#
# Exact on an angular lattice, where "is this offset a whole number of cells?" is integer
# divisibility and needs no tolerance - the same reason `_snap` works on arcseconds. The `_ORIGIN_ATOL`
# path remains for projected grids, where there is no integer lattice to appeal to.
function _originaligned(sourceorigin, targetorigin, sourcestep)
    so, to, s = _arcsecs(sourceorigin), _arcsecs(targetorigin),
                _arcsecs(sourcestep)
    if !isnothing(so) && !isnothing(to) && !isnothing(s) && !iszero(s)
        return (to - so) % s == 0
    end
    offset = (targetorigin - sourceorigin) / sourcestep
    return isapprox(offset, round(offset), atol = _ORIGIN_ATOL)
end

# An angle as a whole number of arcseconds, or `nothing` if it is not one - either because it is a
# projected coordinate (a length, where arcseconds mean nothing) or because it genuinely does not sit
# on the arcsecond lattice.
#
# Used for *origins* as well as steps, so `0` is a legitimate answer - the equator and the prime
# meridian are both exactly zero. A sub-arcsecond **step** is still rejected, but by the ULP bound
# rather than by a zero check: `eps(0.0)` is 5e-324, so only an exactly-zero value passes, and a 0.4
# arcsec step fails as it must. Callers that divide by the result guard against zero themselves.
#
# The comparison cannot be `==`, and this is a property of `Float64` rather than a fudge: a whole
# number of arcseconds has **no exact `Float64` degree representation**, so a genuine value arrives up
# to **1 ULP** off. Worse, *which* way it lands depends on how it was built - `uconvert` scales by a
# precomputed factor and so rounds twice, while a hand-written `n / 3600` rounds once, and the two
# disagree in the last bit for 29 of the first 3600 counts. Both occur here (`_snaparcsec` produces
# the first, a caller writing `cellsize = (30 / 3600)°` the second), so an equality test against
# either one silently rejects values built the other way.
#
# The bound below is measured, not chosen: over all 3600 counts and both construction routes the worst
# deviation is 1.0 ULP, while the nearest *wrong* answer is a whole arcsecond away - relatively
# 2.8e-4, twelve orders of magnitude further. There is no value this can misclassify.
#
# Rounding to an `Int`-valued arcsecond quantity keeps units through the comparison and hands back the
# plain `Int` that `_stepratio`'s `%`/`÷` arithmetic needs, with no second rounding to get it.
#
# `float` is not decoration: `30arcsecond` is an `Int`-valued quantity, and `eps` has no `Integer`
# method, so an caller writing `cellsize = 30arcsecond` - exactly the idiom these units exist to
# support - would otherwise hit a bare `MethodError`.
function _arcsecs(angle)
    _isangle(angle) || return nothing
    n = round(typeof(1arcsecond), angle)
    return abs(angle - n) <= 4eps(float(uconvert(arcsecond, angle))) ?
           ustrip(n) : nothing
end

# How two cell sizes relate: equal, finer, and - the one that matters - whether the target is a whole
# multiple of the source.
#
# Computed on **integer arcseconds** whenever both grids are angular, because the float ratio is not
# reliable for this question even with exactly-snapped steps: 21 arcsec into 7 is exactly 3, but
# `(21/3600)/(7/3600)` is 3.0000000000000004, and a sweep of every pair up to 200 arcsec finds 107
# such cases. Testing those with `==` would report a genuinely exact aggregation as a resample; testing
# them with a tolerance is the arbitrary fudge this replaced. Integers simply answer it. Projected
# grids fall back to the float ratio, where there is no integer lattice to appeal to.
function _stepratio(source, target)
    s, t = _arcsecs(source), _arcsecs(target)
    # `s == 0` would be a sub-arcsecond source grid; there is no integer lattice to compare on, and
    # the `%` below would divide by zero, so fall through to the float ratio.
    if !isnothing(s) && !isnothing(t) && !iszero(s)
        return (equal = t == s, finer = t < s,
                factor = (t % s == 0) ? t ÷ s : nothing)
    end
    r = target / source
    n = round(Int, r)
    return (equal = r == 1, finer = r < 1,
            factor = (n > 0 && r == n) ? n : nothing)
end

# What building `target`-sided cells out of a `source`-sided layer costs - see
# [`AbstractLayerFate`](@ref) for the three answers.
function _resamplecost(samecrs::Bool, source, target, aligned::Bool)
    samecrs || return LayerResampled("it is in a different CRS")
    rel = _stepratio(source, target)
    if rel.equal
        aligned && return LayerKeptExactly()
        return LayerResampled("its cell boundaries are offset from the target grid's")
    end
    rel.finer &&
        return LayerResampled("the target grid is finer than it is, so its values would be interpolated up")
    !isnothing(rel.factor) && aligned && return LayerAggregated(rel.factor)
    return LayerResampled(aligned ?
                          "the target cell size is not a whole multiple of its own" :
                          "its cell boundaries are offset from the target grid's")
end

# Pick the layer whose grid is preserved exactly. The rule (user's choice) is "whichever is already in
# the target CRS", since that is the one needing no reprojection, which maximises how many layers
# survive untouched. Ties are broken by finest resolution, then by the order the layers were given
# (regime before supply) so the choice is deterministic. Returns `nothing` when no layer matches the
# target CRS - nothing can then be preserved exactly, and the caller reports that.
function _choosealign(facts, tcrs)
    candidates = filter(f -> _samecrs(f.crs, tcrs), facts)
    isempty(candidates) && return nothing
    finest = minimum(f -> f.step, candidates)
    for f in candidates                       # first (given-order) among the finest
        f.step == finest && return f
    end
    return first(candidates)
end

# What a grid of `ny × nx` cells costs to hold. Reported for information only (at `:verbose`, and in
# the investigator's report) and used by the "this grid will not survive the simulation" warning.
#
# Two honest figures, since the third is not knowable here:
# * one environment layer is `Float64` per cell - a regime, a supply, one time slice of either;
# * abundances are `Int64` per species per cell (`GridLandscape.matrix`; its `grid` field is a
# *reshape* of the same memory, so it does not double), hence a per-species figure the caller can
# multiply once a species count exists.
# Dispersal lookups are deliberately *not* estimated: they depend on each species' kernel radius, not
# just the cell count, so any figure here would be invented.
function _gridfootprint(ny::Integer, nx::Integer)
    cells = ny * nx
    return (cells = cells, layer = cells * sizeof(Float64),
            perspecies = cells * sizeof(Int64))
end

# The footprint as a line of prose. `Base.format_bytes` gives the familiar KiB/MiB/GiB rendering.
function _footprintmessage(fp::NamedTuple)
    return "$(fp.cells) cells: $(Base.format_bytes(fp.layer)) per environment layer, " *
           "$(Base.format_bytes(fp.perspecies)) per species of abundances"
end

# The integer counterpart of a rounding direction: `fld`/`cld` are floor/ceil division, exact on
# integers where `floor(a / b)` on floats is not.
_intdiv(::typeof(floor)) = fld

_intdiv(::typeof(ceil)) = cld

# Snap a bound outwards onto the alignment layer's own cell boundaries, so that every target cell is a
# whole number of source cells and the layer can be aggregated (or copied) rather than interpolated.
# `dir` is `floor` for a lower bound and `ceil` for an upper one, so snapping only ever grows the box -
# shrinking it could drop data the mask asked for.
#
# Done in **integer arcseconds** whenever the lattice is angular, because the float form is wrong
# for the case that matters most: a bound *already on* a cell boundary. `(bound - origin) / step`
# should then be a whole number, but lands either side of one, so `floor` and `ceil` disagree and the
# box grows by a spurious cell. Swept over every whole-cell bound of a global 30 arcsec lattice,
# **7664 of 21601** did that; on integers it is 0. This is not fixed by the origin and step being
# exact - they already are - because `range` elements in degrees cannot themselves be exact.
#
# Whether a *bound* is on the lattice decides which route is exact. When it is, `fld`/`cld` on
# integers need no tolerance at all. When it is not - a shapefile envelope, a reprojected box - no
# whole-cell answer exists to be lost, so plain `floor`/`ceil` of the float quotient is right; only
# the origin and step need be integral to keep the reconstruction exact.
function _snap(bound, origin, step, dir)
    o, s = _arcsecs(origin), _arcsecs(step)
    (isnothing(o) || isnothing(s) || iszero(s)) &&
        return origin + dir((bound - origin) / step) * step
    b = _arcsecs(bound)
    k = isnothing(b) ? Int(dir((ustrip(arcsecond, bound) - o) / s)) :
        _intdiv(dir)(b - o, s)
    # Given back in the origin's own unit, so an arcsecond-denominated lattice stays in arcseconds.
    return uconvert(unit(origin), float(o + k * s) * arcsecond)
end

# `bounds` snapped onto the grid implied by `origin`/`step`. With no alignment layer there is nothing
# to snap to, so the bounds pass through unchanged.
function _snapbounds(bounds::Extents.Extent, origin, step)
    isnothing(origin) && return bounds
    yo, xo = origin
    return _extentof(_snap(bounds.Y[1], yo, step, floor),
                     _snap(bounds.Y[2], yo, step, ceil),
                     _snap(bounds.X[1], xo, step, floor),
                     _snap(bounds.X[2], xo, step, ceil))
end

# The report's one grid-and-mask array: `mask`'s values on `template`'s coordinates.
#
# **`BitMatrix`, and it is free here**: the report is never mutated, so the one thing a `BitArray`
# is bad at - element writes, which are a read-modify-write of a whole `UInt64` word - never happens
# to it. **`GridHabitat.active` is deliberately NOT backed this way**: interventions write it
# (`Deactivate`/`Reactivate`), so it stays `Matrix{Bool}`, where a write is a plain byte store. The
# shared `AbstractDimArray{Bool, 2}` annotation covers both backings, which is why it is abstract.
#
# **`rebuild`, not `DimArray(bits, dims(template))` - the data path's template is a
# `Rasters.Raster`, and constructing a plain `DimArray` from its dims SILENTLY DROPS THE CRS.**
# `_fullycovered` opens with `tcrs = Rasters.crs(target)` and returns "covered everywhere" when that
# is `nothing`, so the loss did not error - it quietly disabled `simulate_safely`, taking a coverage
# count from 4 cells to 9. Caught by `test_StudyArea.jl`'s own assertion; it would have been
# invisible in a smaller test. `rebuild` keeps the wrapper type, its dims and its metadata, swapping
# only the payload.
#
# The synthetic path passes a bare dims tuple instead - it has no CRS to keep
# (`NoRealWorldPosition`), so a plain `DimArray` is right there.
_activegrid(template::Tuple, mask) = DimArray(BitMatrix(mask), template)

function _activegrid(template, mask)
    return DimensionalData.rebuild(template, data = BitMatrix(mask))
end

# Refine a study area's report into the one a habitat keeps: the same decisions, with what the build
# actually changed written over them.
#
# **Only three fields can differ**, and this is the whole list: `active` (the data branch ANDs the
# area's mask with each layer's coverage), `problems` (a layer the area never saw can cost cells, and
# that is learned here), and `layers` (a layer passed to the builder but not named in the area has no
# `LayerPlan`). Everything else - CRS, cell size, footprint, specs, constraints, every `*source` -
# was settled before the habitat existed.
#
# **The cache is dropped, not carried.** See the field's own comment.
function _refinedreport(r::StudyAreaReport, active, problems)
    return StudyAreaReport(r.crs, r.crssource, r.cellsize, r.cellsizesource,
                           r.align, _activegrid(r.active, active),
                           r.simulate_safely, r.layers, r.footprint,
                           vcat(r.problems, problems), r.specs, r.constraints,
                           nothing, AsBuilt())
end

# Gather what the grid decision needs to know about one materialised layer: its extent both in the
# target CRS (so layers in different CRSs can be positioned against each other at all) and in its
# **own** (where it is always exact - see `_planbounds` for why that matters).
function _layerfacts(name::Symbol, raster::ClimateRaster, tcrs)
    return (name = name, raster = raster, crs = _rastercrs(raster),
            step = _rastercellstep(raster),
            bounds = _dimsextent(raster.array, tcrs),
            native = _dimsextent(raster.array, nothing))
end

# A synthetic layer is generated at whatever shape it is handed, so it has no CRS, extent or
# resolution of its own and can never shape a study area. Skipping it (rather than erroring) means
# `StudyArea(; regime, supply)` reads the same whichever kinds of spec those happen to be, and a
# study area with *only* synthetic layers falls through to the synthetic plan, which asks for the
# `extent`/`cellsize` it genuinely needs.
_shapesgrid(::AbstractSyntheticSpec) = false

_shapesgrid(::Any) = true

# **A zero period must give `NaN`, and plain division does not.** `x/0` is `Inf` in floating
# point - only `0/0` is `NaN` - and `Inf` is far worse than `NaN` here: `_coverage` marks
# a cell inactive on `NaN`, so an `Inf` sails through as a *valid infinite supply* in a cell that
# has no growing season at all. Physically absurd, and silent. The design note that this case
# "needs no policy of its own - the division gives `NaN` and the cell is marked inactive, for free"
# was simply wrong about the arithmetic; the policy is one line, and this is it.
function _perperiod(amount, period)
    quotient = amount / period
    return iszero(period) ? NaN * Unitful.unit(quotient) : quotient
end

_expandspecs(spec::Tuple) = map(_unwrapspec, spec)

# A *named* tuple of specs means the same thing, with the caller's names attached - those are what
# `_expandednames` reports and aligns by. `NamedTuple <: Tuple` is `false`, so this needs its own
# method or a named regime would be taken for a single spec.
_expandspecs(spec::NamedTuple) = map(_unwrapspec, spec)

_expandspecs(spec) = (_unwrapspec(spec),)

# The name each expanded layer is reported and aligned by: the plain keyword for the usual
# single-layer case, numbered (`:regime1`, `:regime2`, ...) when it is a positional multi-variable
# tuple, which has no names of its own.
function _layername(name::Symbol, i::Integer, n::Integer)
    return n == 1 ? name :
           Symbol(name, i)
end

# The names of every layer a `regime`/`supply` argument expands to, in `_expandspecs`'s order. A
# *named* tuple carries the caller's own names - which is what makes `align = :temperature` and
# "layer `:rainfall` is resampled" say something recognisable - and anything else falls back to
# the numbering above.
_expandednames(::Symbol, spec::NamedTuple) = keys(spec)

function _expandednames(name::Symbol, spec::Tuple)
    return ntuple(i -> _layername(name, i, length(spec)), length(spec))
end

_expandednames(name::Symbol, spec) = (name,)

# `align` given by name - the layer must be one that was actually supplied.
function _namedalign(facts, align::Symbol)
    i = findfirst(f -> f.name === align, facts)
    isnothing(i) &&
        error("`align = :$align` names a layer that was not given; the layers here are " *
              join(("`:$(f.name)`" for f in facts), ", ") * ".")
    return facts[i]
end

# The extent before snapping: the mask's if it states one (the mask leads), else the layers' shared
# footprint. Either way it is intersected with every layer, so the grid never covers ground no layer
# reaches - and where that intersection actually cuts into what the mask asked for, it is reported,
# since nothing later can see it (there is no grid beyond the layers' edge to mark as missing).
function _planbounds(facts, maskextent, tcrs, cellsize, problems)
    known = [f for f in facts if !isnothing(f.bounds)]
    isempty(known) &&
        error("no layer carries a CRS, so there is no real-world extent to build a study area on.")
    bounds = isnothing(maskextent) ? first(known).bounds : maskextent
    for f in known
        bounds = _clipto(bounds, f, tcrs)
    end
    isnothing(maskextent) || _maskclamped!(problems, maskextent, bounds,
                  cellsize)
    return bounds
end

# `bounds` (in the target CRS) clipped to the ground layer `f` actually covers.
#
# The clip happens in the **layer's own** CRS, not the target's. Re-expressing a small local box in
# a big layer's CRS is well conditioned, whereas re-expressing a *global* layer's four corners in a
# local projected CRS is meaningless - the corners of the world have no honest British National Grid
# easting, and the envelope of whatever comes back would carve the study area down to nothing in
# particular. Done this way round, a layer that simply contains the area cuts nothing at all: the
# round trip is skipped outright in that case, so the common "the data covers this" path stays exactly
# lossless instead of being nudged outwards by two envelope transforms.
function _clipto(bounds, f::NamedTuple, tcrs)
    inlayer = _bboxin(tcrs, f.crs, bounds)
    clipped = _intersectbounds(inlayer, f.native,
                               "the study area and layer `:$(f.name)`")
    clipped == inlayer && return bounds
    return _bboxin(f.crs, tcrs, clipped)
end

# Report a mask reaching past the ground the layers cover, so the grid had to be clamped to the data.
# Half a cell of slack, because a mask envelope that covers the outermost cells exactly still sits up
# to half a cell outside their centres - that is coverage, not a cut.
function _maskclamped!(problems, requested::Extents.Extent,
                       achieved::Extents.Extent, cellsize)
    over = maximum(abs,
                   _extentvalues(achieved) .- _extentvalues(requested))
    over > cellsize / 2 || return problems
    push!(problems,
          Problem(ProblemWarning(), :mask_clamped,
                  "the mask extends up to " *
                  "$(round(typeof(1.0 * unit(over)), over, digits = 1)) beyond the ground the " *
                  "layers cover, so the grid has been clamped to the data. Widen the layers' " *
                  "coverage, or narrow the mask, if that is unintended."))
    return problems
end

# Report cells the mask marks active where some layer has no data at all. Measured against the mask's
# own **active cells** rather than its bounding box: a mask is usually irregular, so its box
# legitimately contains cells it marks inactive and measuring against it would warn spuriously.
function _masklost!(problems, mask, active, cellsize)
    lost = count(mask) - count(active)
    lost > 0 || return problems
    push!(problems,
          Problem(ProblemWarning(), :mask_lost,
                  "$lost of $(count(mask)) cells inside the mask have no data in at least one " *
                  "layer, so they are inactive$(_areaphrase(lost, cellsize)). Check the layers' " *
                  "real extent against the mask if that is unexpected."))
    return problems
end

# How much ground a cell count covers, as a parenthesised aside. Omitted on a geographic grid, whose
# ° "cell size" is not a physical length at all - a degree cell's real area varies with latitude, so
# any figure derived from it here would be invented. In km^2 once it reaches that scale, since a real
# grid's raw m^2 figure runs to ten digits and more.
function _areaphrase(cells, cellsize)
    dimension(cellsize) === 𝐋 || return ""
    raw = cells * cellsize^2
    area = raw >= 1u"km^2" ? uconvert(u"km^2", raw) : raw
    return " ($(round(typeof(1.0 * unit(area)), area, digits = 1)) of the masked area)"
end

# Sample every layer onto `grid`, intersect the mask with their real coverage, and recut to the
# result - the same "mask AND coverage, then tighten" rule the environment builder uses. The mask is
# returned alongside the outcome so `_masklost!` can say what the coverage took away from it, and the
# full-coverage test with it so `_partialcover!` can cost out what `simulate_safely = false` let
# through. Both are recomputed on each recursion, so all four describe the grid finally returned.
#
# **`simulate_safely` decides one `.&`, and that is the whole of it here.** With it a cell has to be
# wholly inside every layer's footprint to be simulated, so the recut crops inwards to ground the data
# actually describes; without it the rule is what it has always been - a cell survives if the
# resampler gave its centre a value, i.e. if it is more than half covered - and such a cell is then
# handed a *whole* cell's worth of supply over ground that is partly not there.
function _coveredgrid(rasters, payload, grid; simulate_safely::Bool)
    # **`categorical = false` here, and it is only ever asked where the values fall, not what they
    # are.** These are the study area's *data layers*, and all this does with them is find which cells
    # came back non-`NaN` - so the interpolation is never read, and the choice cannot move a class
    # code. Sampling them as `:mode` instead genuinely changes the answer: nearest-class fills cells
    # that bilinear leaves `NaN`, which moves the coverage, the recut and so the study area's own
    # size. Measured - three `test_StudyArea` assertions went from equal to `NaN`-bearing.
    sampled = map(r -> _sampledata(r, grid, name = "layer",
                                   categorical = false), rasters)
    mask = _rastermaskonly(payload, grid, first(sampled))
    full = Matrix{Bool}(reduce(.&, map(r -> _fullycovered(r, grid), rasters)))
    usable = mask .& reduce(.&, map(_nanfree, sampled))
    active = Matrix{Bool}(simulate_safely ? usable .& full : usable)
    rows, cols = _activerange(active, simulate_safely)
    rows = _atleast2(rows, Base.size(active, 1))
    cols = _atleast2(cols, Base.size(active, 2))
    (length(rows) == Base.size(active, 1) &&
     length(cols) == Base.size(active, 2)) &&
        return (grid = grid, active = active, mask = Matrix{Bool}(mask),
                full = full)
    return _coveredgrid(rasters, payload, grid[Y(rows), X(cols)],
                        simulate_safely = simulate_safely)
end

# Say how many cells `simulate_safely = false` let through - cells the simulation will run on although
# the data does not describe all of them, and which are handed a whole cell's worth of supply anyway.
#
# **Nothing is said on the `true` path** (user, 2026-08-13): a cell is dropped at the edge of almost
# every re-gridded area, so a notice there would fire on nearly all of them and mean nothing.
function _partialcover!(problems, covered::NamedTuple, cellsize,
                        simulate_safely::Bool)
    simulate_safely && return problems
    partial = count(covered.active .& .!covered.full)
    partial > 0 || return problems
    push!(problems,
          Problem(ProblemNotice(), :partly_covered,
                  "$partial of $(count(covered.active)) active cells are only partly backed by " *
                  "data and are being simulated anyway$(_areaphrase(partial, cellsize)), because " *
                  "`simulate_safely = false`. Each is given a whole cell's worth of supply. Drop " *
                  "the keyword to simulate only cells the data covers completely."))
    return problems
end

# The non-`NaN` cells of a sampled layer - its real coverage, before any mask is applied. A 3-D
# (monthly) read is judged on its first time slice, the convention used throughout.
function _nanfree(sampled)
    slice = ndims(sampled) == 2 ? sampled : view(sampled, :, :, 1)
    return Matrix{Bool}(.!isnan.(slice))
end

# The mask alone, rasterised onto `grid`. `_rastermask` takes a regime so it can fall back to that
# regime's coverage when there is no mask; here coverage is handled separately, so the no-mask case is
# simply "everything", and every other payload delegates to the existing rasteriser.
_rastermaskonly(::Nothing, grid, reference) = trues(Base.size(reference)[1:2])

_rastermaskonly(payload, grid, reference) = _rastermask(payload, nothing, grid)

# The target grid itself. Where the chosen cell size *is* the alignment layer's own step, that layer's
# grid is **adopted and merely cropped** rather than rebuilt - which is what makes `_resamplecost`'s
# `:exact` verdict true, and what stops a mask narrowing the grid from also shifting every cell by
# part of a step. Synthesising instead would force *square* cells of the layer's north-south step,
# silently re-gridding any source whose east-west step differs. Any other cell size is a deliberate
# re-gridding, so it does synthesise.
function _targetgrid(chosen::Union{NamedTuple, Nothing}, tcrs, bounds,
                     cellsize)
    (isnothing(chosen) || isnothing(chosen.step) ||
     chosen.step != cellsize) &&
        return _crstemplate(tcrs, bounds, cellsize)
    return _cropgrid(_owngrid(chosen.raster), bounds)
end

# The contiguous index range of `vals` lying within `[lo, hi]`. A *range*, never `findall`'s vector: a
# `Y`/`X` lookup is monotone so the hits are contiguous either way, and indexing with a range keeps a
# `Regular` span where a vector of indices would drop it to `Irregular` - which `Rasters.resample`
# rejects as a `to =` target.
#
# Compared in **integer arcseconds** on an angular lattice, for the same reason `_snap` is: a
# `lo <= v <= hi` test on floats is decided by the last bit. The bounds arrive exactly on the lattice
# (they come from `_snap`), but a layer's own coordinates need not - aggregating a *cropped* raster
# yields values ~1 ULP from those of the same cells aggregated whole. That is enough for the first
# and last cell to fail `<=` and be dropped, so a windowed read of a coarsened layer produced a grid
# two rows and one column smaller than a whole read of it. Integers make the endpoints exact.
function _indexrange(vals, lo, hi)
    a, b = _arcsecs(lo), _arcsecs(hi)
    ns = (isnothing(a) || isnothing(b)) ? nothing : map(_arcsecs, vals)
    idx = (isnothing(ns) || any(isnothing, ns)) ?
          findall(v -> lo <= v <= hi, vals) : findall(n -> a <= n <= b, ns)
    return isempty(idx) ? (1:0) : first(idx):last(idx)
end

# `grid` cropped to the unitful extent by plain integer indexing, which preserves
# its cell registration exactly (and needs no span rebuild, so it survives cropping to a single row).
function _cropgrid(grid, bounds::Extents.Extent)
    ylo, yhi, xlo, xhi = _extentvalues(bounds)
    # Cell **centres**, asked for by name - this selects the cells whose centre falls in the box,
    # which is what the message below already describes and what it has always done. Locus-blind:
    # the raw lookup values were only the centres by convention, and stopped being once the grid moved
    # to `Intervals(Start)`.
    yx = _cellcentres(grid)
    rows = _indexrange(yx.lat, ylo, yhi)
    cols = _indexrange(yx.long, xlo, xhi)
    (isempty(rows) || isempty(cols)) &&
        error("the study area's extent selects no cell of the grid: it falls between cell centres, " *
              "or outside the data. Bounds were $(bounds).")
    return grid[Y(_atleast2(rows, Base.size(grid, Y))),
                X(_atleast2(cols, Base.size(grid, X)))]
end

# One layer's entry in the plan: where it is, and what the chosen grid costs it.
#
# Alignment is only asked about a layer already in the target CRS. Its `bounds` are stated in the
# *target's* units while its `step` is in its own, so across a reprojection the two are not even
# dimensionally comparable (degrees against metres) - and the answer would be ignored regardless,
# since a different CRS always costs a resample.
function _planfor(f::NamedTuple, tcrs, cellsize, origin)
    same = _samecrs(f.crs, tcrs)
    aligned = same && !isnothing(origin) && !isnothing(f.step) &&
              _originaligned(f.bounds.Y[1], first(origin), f.step)
    cost = isnothing(f.step) ?
           LayerResampled("it has no resolvable cell size") :
           _resamplecost(same, f.step, cellsize, aligned)
    return LayerPlan(f.name, f.crs, f.step, f.native, cost)
end

# The warnings, all of which describe grids that work but are probably not what was wanted. They are
# collected rather than emitted so that `investigate_study_area` can list them and `StudyArea` can
# choose how loudly to say them.
# The WGS84 degree box a CRS was actually defined for, straight from PROJ's own database (via GDAL's
# OSR, so no dependency beyond ArchGDAL). `nothing` when the database declares none.
function _areaofuse(crs)
    return ArchGDAL.importCRS(crs, order = :trad) do sr
        w, s, e, n = Ref{Cdouble}(), Ref{Cdouble}(), Ref{Cdouble}(),
                     Ref{Cdouble}()
        name = Ref{Cstring}()
        ArchGDAL.GDAL.osrgetareaofuse(sr.ptr, w, s, e, n, name) == 1 ||
            return nothing
        return (west = w[], south = s[], east = e[], north = n[],
                name = unsafe_string(name[]))
    end
end

# Record a problem when the grid reaches outside the CRS's declared area of use. A projection is only
# accurate within it, and exceeding it distorts silently rather than failing - which is why this is a
# reported problem rather than an error.
function _crsareaofuse!(problems, tcrs, bounds)
    use = _areaofuse(tcrs)
    isnothing(use) && return problems
    # The extent is in the target's own coordinates; the area of use is always WGS84 degrees.
    wgs = _bboxin(tcrs, Rasters.EPSG(4326), bounds)
    la1, la2 = ustrip(°, wgs.Y[1]), ustrip(°, wgs.Y[2])
    lo1, lo2 = ustrip(°, wgs.X[1]), ustrip(°, wgs.X[2])
    span = max((la2 - la1) * (lo2 - lo1), eps())
    inside = max(0.0, min(la2, use.north) - max(la1, use.south)) *
             max(0.0, min(lo2, use.east) - max(lo1, use.west))
    outside = 1 - inside / span
    outside > _CRS_OUTSIDE_FRAC || return problems
    push!(problems,
          Problem(ProblemWarning(), :crs_area_of_use,
                  "$(round(Int, 100 * outside))% of this grid lies outside the area " *
                  "$(_crsname(tcrs)) was defined for ($(use.name); longitude " *
                  "$(use.west)...$(use.east), latitude $(use.south)...$(use.north)). A projection " *
                  "used far from its intended ground distorts distance and area badly - " *
                  "$(_crsadvice(wgs))."))
    return problems
end

# Everything wrong with a decided grid, gathered in one pass so a report names all of it at once
# rather than one fault per re-run. Mutates `problems`, which the report then owns.
function _collectproblems!(problems, plans, tcrs, active, fp::NamedTuple,
                           chosen::Union{NamedTuple, Nothing}, bounds)
    # On a geographic target the bounds are already degrees, so they name the UTM zone directly -
    # the same concrete suggestion the metric-`cellsize` error makes, so the fix is one paste either
    # way.
    if !_isprojectedcrs(tcrs)
        push!(problems,
              Problem(ProblemWarning(), :geographic,
                      "this is a geographic (° coordinate) grid, where a cell's physical size " *
                      "varies with latitude, so it cannot be simulated - `build_ecosystem` " *
                      "requires a projected grid, because dispersal assumes one uniform cell " *
                      "size. To build a metric grid, $(_crsadvice(bounds))."))
    else
        _crsareaofuse!(problems, tcrs, bounds)
    end
    for p in plans
        p.kind isa LayerResampled && occursin("finer", p.kind.reason) &&
            push!(problems,
                  Problem(ProblemWarning(), :upsampling,
                          "layer `:$(p.name)` is being interpolated up onto a finer grid than it " *
                          "has, which invents detail it never measured."))
        p.kind isa LayerAggregated && p.kind.factor >= _HUGE_FACTOR &&
            push!(problems,
                  Problem(ProblemWarning(), :extreme_aggregation,
                          "layer `:$(p.name)` is aggregated $(p.kind.factor)× per cell, which is enough " *
                          "to stop the result representing the source."))
    end
    # Degrading the layer you probably care about: something finer than the aligned layer is
    # resampled. Only layers in the target CRS can be compared - elsewhere a step is stated in
    # different units, and "finer" has no meaning across degrees and metres.
    if !isnothing(chosen)
        for p in plans
            p.kind isa LayerResampled && !isnothing(p.step) &&
                _samecrs(p.crs, tcrs) && p.step < chosen.step &&
                push!(problems,
                      Problem(ProblemWarning(), :degraded_finest,
                              "layer `:$(p.name)` is finer than the layer being aligned to " *
                              "(`:$(chosen.name)`), so the more detailed layer is the one being " *
                              "resampled. Consider `align = :$(p.name)`."))
        end
    end
    fp.cells > _BIG_GRID &&
        push!(problems,
              Problem(ProblemWarning(), :grid_too_big,
                      "this grid may not survive the simulation - $(_footprintmessage(fp))."))
    Base.size(active, 1) >= 2 && Base.size(active, 2) >= 2 ||
        push!(problems,
              Problem(ProblemWarning(), :grid_too_small,
                      "this grid is $(join(Base.size(active), "×")) cells, too few for dispersal to " *
                      "mean anything."))
    frac = count(active) / length(active)
    frac < _SPARSE_ACTIVE &&
        push!(problems,
              Problem(ProblemWarning(), :sparse_active,
                      "only $(round(100 * frac, digits = 1))% of this grid is active, so most of it " *
                      "is dead space - usually a sign the mask or extent is not what was meant."))
    return problems
end

# A study area with no data layers at all: the synthetic grid, positioned nowhere in particular, whose
# geometry comes entirely from `extent` (a `(y, x)` tuple of lengths) and `cellsize`.
function _syntheticplan(extent, cellsize, within, problems, layers, cons,
                        simulate_safely, cache)
    (isnothing(extent) || isnothing(cellsize)) &&
        error("a study area with no layers needs both `extent` (a `(y, x)` tuple of lengths) and " *
              "`cellsize`; got extent = $(repr(extent)), cellsize = $(repr(cellsize)).")
    geometry = _gridgeometry(extent, cellsize)
    active = _resolvesyntheticactive(within, geometry.dim, geometry.cellsize)
    fp = _gridfootprint(geometry.dim...)
    # Recorded, not applied: a synthetic area has no data footprint for a cell to fall outside of,
    # so every cell is covered by construction and the flag is a no-op here. It is still stored, so
    # `StudyArea(synthetic, regime = <data>)` refines with the value that was asked for.
    return StudyAreaReport(nothing, NoRealWorldPosition(), geometry.cellsize,
                           GivenByUser(), nothing,
                           _activegrid(_syntheticyx(geometry), active),
                           simulate_safely,
                           LayerPlan[], fp,
                           problems, layers, cons, cache, AsInvestigated())
end

# --- Reporting -------------------------------------------------------------

# A CRS as a human would name it. `EPSG` objects print as `GeoFormatTypes.EPSG{1}((27700,))`, and a
# WKT string runs to hundreds of characters, neither of which belongs in a message meant to be read.
_crsname(::Nothing) = "none - an index grid"

function _crsname(crs::Rasters.EPSG)
    return "EPSG:" * join(Rasters.GeoFormatTypes.val(crs), ", ")
end

function _crsname(crs)
    text = string(Rasters.GeoFormatTypes.val(crs))
    name = match(r"\"([^\"]+)\"", text)          # the first quoted name in a WKT string
    return isnothing(name) ? first(text, 60) : name.captures[1]
end

# A length floored to one significant figure - 14218.48 m becomes 10 km, 927.66 m becomes 900 m.
#
# Floored, never rounded to nearest, so the chosen grid is never *coarser* than the size measured
# off the projection and no resolution is lost that the source data carried. The cost is cells: the
# worst case is a value just under a power of ten (19999 m -> 10000 m), which nearly halves the cell
# edge and so gives ~4× the cells. That is deliberate and visible - the report prints the grid
# dimensions and the memory both.
function _floor1sf(size)
    m = ustrip(u"m", uconvert(u"m", size))
    m > 0 || return size
    decade = 10.0^floor(Int, log10(m))
    floored = floor(m / decade) * decade
    # Carried in the unit it reads best in, not just shown in it, so a caller reaching for
    # `report.cellsize` gets `10.0 km` rather than `10000.0 m`.
    return isinteger(floored / 1000) ? (floored / 1000) * u"km" : floored * u"m"
end

# The cell size to actually adopt, given what `_inferredcellsize` measured.
#
# Only a `MeasuredAcrossProjection` size is rounded. That one was measured across the projection for a grid that
# is *synthesised* - no layer is in the target CRS, so every layer is resampled whatever size we pick,
# and rounding forfeits nothing that was not already forfeit. A `GivenByUser` size is the user's, and an
# `TakenFromAlignedLayer`/`AgreedByAllLayers` one is a real layer's own step, where changing it would destroy the exact
# alignment that made it worth adopting.
#
# Lengths only. The angular analogue would be flooring to a whole number of arcseconds, which
# belongs with the arcsecond lattice work rather than here; a geographic target is refused by
# `build_ecosystem` anyway.
function _roundedcellsize(size, source, problems)
    (source isa MeasuredAcrossProjection && dimension(size) === Unitful.𝐋) ||
        return (size, source)
    rounded = _floor1sf(size)
    rounded == size && return (size, source)
    push!(problems,
          Problem(ProblemNotice(), :cellsize_rounded,
                  "the cell size measured across the projection was " *
                  "$(round(typeof(1.0u"m"), size, digits = 2)); it has been rounded down to " *
                  "$(_cellsizetext(rounded)) so the grid is described in whole units. Rounding " *
                  "down never coarsens the grid, but it does mean more cells than the measured " *
                  "size implies - pass `cellsize` to choose your own."))
    return (rounded, RoundedFromMeasurement())
end

# A cell size as a human would write it. Geographic grids are laid out in whole arcseconds, so a
# WorldClim cell is 10 arcmin and an EarthEnv one 30 arcsec - but in degrees those read
# `0.16666666666666666°` and `0.008333333333333333°`, which say nothing to a reader and invite the
# suspicion that the grid is ragged. Shown in the coarsest subdivision that keeps the number whole.
# Anything off the arcsecond lattice, and any projected size, is left exactly as it is rather than
# dressed up.
function _cellsizetext(size)
    n = _arcsecs(size)
    if !isnothing(n)
        n % 3600 == 0 && return string((n ÷ 3600) * °)
        n % 60 == 0 && return string((n ÷ 60) * arcminute)
        return string(n * arcsecond)
    end
    dimension(size) === Unitful.𝐋 || return string(size)
    # A whole number of kilometres reads better as one, and a whole number of metres should not carry
    # a pointless `.0`. Anything else is left exactly as it is rather than given spurious precision.
    m = ustrip(u"m", uconvert(u"m", size))
    isinteger(m / 1000) && return string(Int(m ÷ 1000) * u"km")
    isinteger(m) && return string(Int(m) * u"m")
    return string(size)
end

# One line per layer: what it is, and what the grid costs it.
function _layerline(p::LayerPlan)
    step = isnothing(p.step) ? "no resolvable cell size" :
           "$(_cellsizetext(p.step)) cells"
    return "  :$(p.name) - $step, $(_fatephrase(p.kind))"
end

# Public keyword -> what `_analyse` expects. `missing` ("not specified") and `nothing` ("explicitly
# cleared") both mean "no constraint" once there is no base area left to inherit from; they differ
# only during refinement, which `_inherit` resolves first.
_constraint(x) = ismissing(x) ? nothing : x

# Refinement: `missing` takes the base area's value, anything else (including `nothing`, meaning
# "clear it") overrides.
_inherit(given, base) = ismissing(given) ? base : given

# Is this base a report of a grid that was actually *built*, rather than merely investigated? Only
# such a base can be copied verbatim - everything else has nothing to preserve that re-deriving would
# not reproduce.
function _iscopy(base)
    (ismissing(base) || isnothing(base)) && return false
    return _basereport(base).stage isa AsBuilt
end

# **A copied report needs a WORKING cache, and an as-built one has none.** A `GridHabitat`
# deliberately discards its report's reads (`cache === nothing`, see `StudyAreaReport`) - they were
# consumed inputs, and keeping them would pin every raster the build touched. But a `StudyArea` is a
# thing you *build on*, and `_materialisefield` hands `area.report.cache` straight to
# `_asraster`/`_combineon`, which have methods for a `LayerCache` and none for `nothing`.
# So copying the report verbatim made an as-built area unable to **read data at all**: rebuilding
# with a synthetic spec worked and with a data-backed one died on a `MethodError` naming a private
# function. Measured, and it is the whole of the defect.
# A *fresh, empty* cache rather than the original: the reads are genuinely gone, so this says
# "nothing cached yet" - which is true, costs one re-read, and is exactly what `_resolveinputs` does
# for the re-derive path (`isnothing(old.cache) ? LayerCache() : old.cache`). The copy branch
# returns early and so bypasses `_resolveinputs` entirely, which is the one place that rule was
# missed.
# Nothing the copy exists to preserve is touched - `active`, `problems` and `stage` all carry
# across - because a read cache is not part of what a report *describes*.
function _copyablereport(r::StudyAreaReport)
    isnothing(r.cache) || return r
    return StudyAreaReport(r.crs, r.crssource, r.cellsize, r.cellsizesource,
                           r.align, r.active, r.simulate_safely, r.layers,
                           r.footprint, r.problems, r.specs, r.constraints,
                           LayerCache(), r.stage)
end

# Fold a `base` area and the caller's keywords into one set of inputs, with `missing` meaning
# "inherit from base" and `nothing` meaning "clear what would have been inherited". Keeping that
# resolution here is what lets refining an area and building one fresh behave identically.
function _resolveinputs(base; regime, supply, within, crs, cellsize, extent,
                        align, simulate_safely)
    old = _basis(ismissing(base) ? nothing : base)
    prev(field) = get(old.specs, field, nothing)
    prevcon(field) = get(old.constraints, field, nothing)
    # The supply is desugared **here**, at the one point where both roles are named, so everything
    # downstream sees the rewritten spec: the divisor layer joins the CRS probe and the grid decision,
    # and receives the same read window and `cut` as its numerator. Desugaring later would leave `gsl`
    # out of the grid it is meant to be divided on.
    layers = map(_constraint,
                 (regime = _inherit(regime, prev(:regime)),
                  supply = _desugarsupply(_inherit(supply, prev(:supply)))))
    cons = map(_constraint,
               (within = _inherit(within, prevcon(:within)),
                crs = _inherit(crs, prevcon(:crs)),
                cellsize = _inherit(cellsize, prevcon(:cellsize)),
                extent = _inherit(extent, prevcon(:extent)),
                align = _inherit(align, prevcon(:align)),
                simulate_safely = _inherit(simulate_safely,
                                           prevcon(:simulate_safely))))
    return (layers = layers, constraints = cons,
            cache = isnothing(old.cache) ? LayerCache() : old.cache)
end

# Say what was decided, at the requested volume. Pure presentation over one already-computed report,
# so a level can never report something the build does not do.
function _emit(report::StudyAreaReport, verbosity::Symbol)
    verbosity === :silent && return nothing
    if verbosity in (:verbose, :full, :debug)
        show(stdout, MIME"text/plain"(), report)
        println(stdout)
        return nothing
    end
    verbosity === :normal ||
        error("`verbosity` must be `:silent`, `:normal` or `:verbose` (aliases `:full`, `:debug`); " *
              "got `:$verbosity`.")
    # Announce every value that was guessed rather than given - the antidote to hidden automation.
    # `NoRealWorldPosition` is not a guess: a synthetic area has no real-world position for a CRS to
    # describe.
    report.crssource isa Union{GivenByUser, NoRealWorldPosition} ||
        @info "study area CRS: $(_crsname(report.crs)) ($(_sourcephrase(report.crssource)))."
    report.cellsizesource isa GivenByUser ||
        @info "study area cell size: $(_cellsizetext(report.cellsize)) ($(_sourcephrase(report.cellsizesource)))."
    isnothing(report.align) ||
        @info "study area aligned to layer `:$(report.align)`, which is kept exactly."
    # ...and every layer the grid costs something.
    for p in report.layers
        p.kind isa LayerAggregated &&
            @info "layer `:$(p.name)` is aggregated $(p.kind.factor)× onto the study grid (exact)."
        p.kind isa LayerResampled &&
            @info "layer `:$(p.name)` is resampled onto the study grid, because $(p.kind.reason)."
    end
    for p in report.problems
        _report(p.severity, p.message)
    end
    return nothing
end

# --- The study area itself -------------------------------------------------
#
# **The `StudyArea` type itself is declared at the head of this file, well before this point.**
# `GridHabitat` names an `area::StudyArea{L}` in a **field annotation**, which Julia resolves when the
# struct is defined, so the type must exist before `GridHabitat.jl` is included - the same constraint
# that put `StudyAreaReport` in a file of its own. Everything that *decides* an area is below.

# --- Inspecting a layer on a decided grid ----------------------------------

# **A member that is itself a tuple is the removed `(source, code)` pair form, and dispatch cannot
# see it**: `materialise(::Tuple, ...)` is the *container* method, so a nested tuple recurses into it
# and reports something about `WorldClim{BioClim}` and `1` as though they were two layers. Checking
# the member here is what makes the refusal say what is actually wrong.
_checklayerspec(member::Tuple) = _sourcepairnotaspec(member)

_checklayerspec(member) = member

# Whether `NaN` can be stored in values of this type at all - the mechanical form of the categorical
# exception above, which also covers a `Bool`- or `Int`-valued combine that never called itself
# categorical. A layer that cannot hold `NaN` cannot express "no data" either, so there is nothing to
# express it with rather than a rule being waived.
_cannan(::Type{<:AbstractFloat}) = true

_cannan(::Type{<:Unitful.AbstractQuantity{T}}) where {T} = _cannan(T)

_cannan(::Type) = false

# `NaN` (in whatever unit the values carry) wherever `full` is false. Two methods rather than a
# branch on `ndims`, because a monthly read is `(Y, X, time)` and the whole time series of an
# uncovered cell goes together.
function _blankuncovered!(values::AbstractArray{T, 2}, full) where {T}
    for i in Base.axes(values, 1), j in Base.axes(values, 2)
        full[i, j] || (values[i, j] = NaN * oneunit(T))
    end
    return values
end

function _blankuncovered!(values::AbstractArray{T, 3}, full) where {T}
    for i in Base.axes(values, 1), j in Base.axes(values, 2)
        full[i, j] || (view(values, i, j, :) .= NaN * oneunit(T))
    end
    return values
end

# --- Building an environment on a decided area -----------------------------

# --- CRS advice: what to suggest when the grid is probably not what was wanted ---

# Whether a CRS measures its own coordinates as a projected length rather than geographic degrees.
# Decided by *dimension* of the CRS's declared coordinate unit (via `_crsunit`, reused from
# `rasters.jl` so there is a single authority on "what unit is this CRS in") - 𝐋 => projected,
# dimensionless => degrees. Dimension rather than an exact unit match, so a CRS in feet (or any other
# length `_crsunit` learns) classifies correctly without touching this.
_isprojectedcrs(crs) = dimension(_crsunit(crs)) === 𝐋

function _crscandidates()
    isnothing(_CRS_CANDIDATES[]) || return _CRS_CANDIDATES[]
    out = NamedTuple[]
    count = Ref{Cint}(0)
    list = ArchGDAL.GDAL.osrgetcrsinfolistfromdatabase("EPSG", C_NULL, count)
    try
        for i in 1:count[]
            info = unsafe_load(unsafe_load(list, i))
            (info.bBboxValid == 1 && info.bDeprecated == 0 &&
             info.eType == ArchGDAL.GDAL.OSR_CRS_TYPE_PROJECTED) || continue
            info.dfWestLongitudeDeg < info.dfEastLongitudeDeg || continue
            method = info.pszProjectionMethod == C_NULL ? "" :
                     unsafe_string(info.pszProjectionMethod)
            push!(out,
                  (code = parse(Int, unsafe_string(info.pszCode)),
                   name = unsafe_string(info.pszName), method = method,
                   west = info.dfWestLongitudeDeg,
                   south = info.dfSouthLatitudeDeg,
                   east = info.dfEastLongitudeDeg,
                   north = info.dfNorthLatitudeDeg))
        end
    finally
        ArchGDAL.GDAL.osrdestroycrsinfolist(list)
    end
    _CRS_CANDIDATES[] = out
    return out
end

# Cylindrical Mercator's area distortion grows without bound towards the poles, so it is never a
# defensible grid for an ecological simulation, where a cell's area *is* the habitat available in it.
# Excluded from the suggestions rather than merely ranked down, since the whole point is a paste-able
# answer.
#
# Only the *cylindrical* variants - "Mercator (variant A/B)" and "Popular Visualisation Pseudo
# Mercator". **Transverse** and **Oblique** Mercator are different projections that merely share the
# name: they are near-area-true within the narrow zone they are defined for, and excluding them would
# throw away almost every good local answer, British National Grid and every UTM zone included.
function _areadistorting(method)
    return occursin("Mercator", method) &&
           !occursin("Transverse", method) && !occursin("Oblique", method)
end

# A projected CRS actually suited to a WGS84 degree extent: the one with the
# smallest declared area of use that still wholly contains it. "Smallest containing" is what makes
# the answer specific - British National Grid for Scotland rather than a continental projection - and
# an extent too wide for any local CRS (Kenya spans two UTM zones) falls back to a continental or
# global one by the same rule, with no special case. `nothing` if the database offers nothing, which
# should not happen while global equal-area CRSs exist.
function _suggestcrs(e::Extents.Extent)
    la1, la2 = ustrip(°, e.Y[1]), ustrip(°, e.Y[2])
    lo1, lo2 = ustrip(°, e.X[1]), ustrip(°, e.X[2])
    fits = [c
            for c in _crscandidates()
            if c.west <= lo1 && c.east >= lo2 && c.south <= la1 &&
                   c.north >= la2 && !_areadistorting(c.method)]
    isempty(fits) && return nothing
    return argmin(c -> (c.east - c.west) * (c.north - c.south), fits)
end

# The suggestion as a paste-able clause, naming the CRS so the reader can judge it rather than
# pasting a bare number. Falls back to bare prose when the database offers nothing.
function _crsadvice(e::Extents.Extent)
    best = _suggestcrs(e)
    isnothing(best) && return "pass a projected `crs`"
    return "for this extent `crs = Rasters.EPSG($(best.code))` ($(best.name)) is a reasonable choice"
end
