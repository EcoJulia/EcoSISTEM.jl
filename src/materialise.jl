# SPDX-License-Identifier: LGPL-3.0-or-later
#
# PUTTING LAYERS ON A DECIDED GRID. The study area decides where the cells are; this samples layers
# onto them. Three things that must stay together:
#
# the read cache      raw source reads, kept so that refining an area re-reads nothing
# inspection          `materialise` — what one layer looks like on the grid, without building
# building            what `GridHabitat`'s constructor calls to assemble the real thing
#
# Why one file. `_materialiseon` (building) and `_materialisefield` (inspection) implement the same
# thing method for method — synthetic, data-backed, `ConstructedSpec`, `ClimateRaster` — differing
# only in their return wrapper, and they have drifted apart three separate times: united `(Y, X)`
# dims on one path and not the other, a hardcoded sampling mode, and a `NicheSpec` method one had
# and the other did not. Each was found by accident. `materialise`'s docstring promises that its
# whole point is to see what the builder will get, so a divergence does not merely duplicate code —
# it makes the inspection function lie, and nothing can catch it: there is no shared thing to
# compare against and both paths pass their own tests. Keeping them in one file is the cheapest
# guard there is, and it is what makes uniting them later a local change.
#
# `_unitedyx` is here for the same reason: it is the shared fix from the first of those three
# divergences, and both paths call it.

# Both of these name `SourceSpec`, which is not in scope where `ReadKey`/`LayerCache` are now
# defined — so the *constructor* and the lookup stay here with the code that uses them. A method can
# live anywhere; a type cannot.
using Rasters

import Rasters: Projected

using DimensionalData

using Unitful

ReadKey(spec::SourceSpec) = ReadKey(spec.source, spec.code, spec.readkw)

# `::AbstractSpec` rather than the full [`LayerInput`](@ref): the tuple/named-tuple forms are the
# *separate* method below, so admitting them here would make the two ambiguous.
"""
    materialise(spec::AbstractSpec, area::StudyArea; role = missing)
    materialise(specs::Union{Tuple, NamedTuple}, area::StudyArea; role = missing)

Put `spec` — or a tuple/named tuple of specs — onto `area`'s grid and hand back the **layer** it
becomes, so it can be plotted or checked before a simulation is built from it.

Works the same way for every kind of spec — a data source, a `ConstructedSpec`, or a synthetic
gradient — so layers can be compared cell for cell on the grid they will actually share. Data-backed
specs are *sampled onto* the grid and synthetic ones *generated at* its shape, and the members of a
tuple may be **mixed**. The result is an [`AbstractLayer`](@ref) whose `matrix` is a `(Y, X)` array
carrying the area's real coordinates; several specs give a [`LayerCollection`](@ref) keeping the
caller's names, exactly as the builder produces.

This runs the same code the builder does, so what you see is what the simulation gets. It is also the
way to check a *synthetic* layer's layout, which is otherwise invisible: a gradient's direction
follows real-world north (not array row order), and its endpoints span the grid's bounding rectangle,
so on a non-rectangular study area the extremes may fall on inactive cells.

# Arguments

  - `spec`: **any** layer spec — a [`SourceSpec`](@ref), a [`ConstructedSpec`](@ref), or a synthetic
    one ([`UniformSpec`](@ref), [`GradientSpec`](@ref), [`PeakedSpec`](@ref), [`NicheSpec`](@ref)) —
    or a **tuple/named tuple of them** for a multi-layer regime, which is materialised member by
    member and keeps its names. The members may be **mixed**: data-backed specs are *sampled onto*
    the grid and synthetic ones *generated at* its shape, so both can share one grid.
  - `area`: the [`StudyArea`](@ref) whose grid to put it on.
  - `role`: `missing` (the default) gives a regime layer holding the values as they are — its spatial
    pattern. **A spec does not know its own role**: `UniformSpec(298K)` is a regime and
    `UniformSpec(1kJ/km^2/day)` a supply only by virtue of which keyword it is passed to. Pass
    `role = Condition` to canonicalise to the axis's unit, or `role = Resource` for a
    [`Supply`](@ref) whose per-area rate has become an absolute per-cell one — giving exactly the
    layer the simulation uses.
"""
function materialise(spec::AbstractSpec, area::StudyArea; role = missing)
    # A declared change is ignored here: `materialise` shows what a spec becomes on this grid, and a
    # change is a separate declaration that `GridHabitat` hangs on the layer afterwards — the
    # same reason `StudyArea` ignores it.
    spec = _unwrapspec(spec)
    return _applyrole(_materialisefield(spec, area), role, _specaxis(spec),
                      area)
end

# **A multi-layer regime, materialised member by member.** Without this a tuple fell through to the
# data-backed path and failed trying to coerce the first member into a `SourceSpec` — even though
# `GridHabitat` has always built exactly such a regime, including a **mix** of real and
# synthetic layers. So the capability existed and only inspection could not reach it, which is the
# wrong way round: the whole point of `materialise` is to see what the builder will get.
#
# `map` over a named tuple keeps its names, so `materialise((temp = …, rain = …), area)` gives a
# `LayerCollection` carrying them — the same convention `GridHabitat` uses, and what makes the
# members identifiable rather than positional. It is a *collection*, not a bare tuple, because that
# is what the builder produces: reach the members with `regimes`/`named_regimes` (or
# `supplies`/`named_supplies`), which work identically on a single layer.
function materialise(specs::Union{Tuple, NamedTuple}, area::StudyArea;
                     role = missing)
    return LayerCollection(map(s -> materialise(_checklayerspec(s), area,
                                                role = role), specs))
end

# **Anything that is not a spec, refused by naming what is wrong with it.** `_specaxis` already
# says the useful thing — *"a supply given as a bare value … has no niche axis, and its meaning
# cannot be taken from its unit"* — so this reaches it rather than restating it.
#
# **Ordering, and it is load-bearing.** The axis must be resolved *before* any conversion, or the
# refusal falls to dispatch and a tailored message becomes a bare `MethodError`.
# It matters most **inside a container**, where a signature cannot help: a multi-layer supply
# legitimately *is* a tuple, so `supply = (3.0m/s, SUP)` reaches the element one at a time and only a
# per-element check can say which element is wrong.
function materialise(x, area::StudyArea; role = missing)
    _specaxis(x)                       # raises, naming the value and why it has no axis
    return error("`materialise` needs a spec (an `AbstractSpec`), or a tuple of them; got a " *
                 "$(nameof(typeof(x))).")
end

# **An already-built layer passes through, with its grid checked.** The builder routes through
# `materialise`, so this method is what serves `GridHabitat(regime = <a built layer>)`.
#
# **The grid check is the whole of the work.** A built layer is used as it stands — never
# resampled, because its values are already canonical — so it has to have been built on this area.
# `role` is ignored, deliberately: the layer already *is* a regime or a supply, and re-applying a
# role would canonicalise values that are canonical already, i.e. convert them twice.
function materialise(layer::AbstractLayer, area::StudyArea; role = missing)
    _yxcompatible(_yx(layer), dims(area.report.active, (Y, X))) ||
        error("an already-built layer is $(size(layer.matrix)) on its own grid, but this study " *
              "area's grid is $(size(area.report.active)). A built layer is used as it stands — it " *
              "is never resampled, because its values are already canonical — so it has to have " *
              "been built on this area. Pass a spec instead if it should be sampled onto the grid.")
    # **Copied, not passed through — this is the ONE place a layer the caller owns enters the
    # package.** Returning it as-is made a habitat alias it, and `_applychange!` then wrote
    # `layer.matrix .= …` into the caller's own object on every timestep; two habitats given one
    # layer shared it. See `_ownlayer` for why a *shallow* copy is enough.
    # **`materialise`'s promise gets stronger, not weaker**: it says it shows what the builder will
    # get, and the builder now gets exactly this — an independent layer. Every other `materialise`
    # method already returns something freshly generated, so this is the only one that had to change.
    return _ownlayer(layer)
end

# **A raster reaches the shared refusal, not the generic one.** `_specaxis(::ClimateRaster)`
# *succeeds* — it resolves an axis from the layer code, defaulting to `NicheAxis` — so the fallback
# above would fall through to its own message and say nothing about the remedy. This is the fourth
# entry point onto `_rasternotaspec`, and it is here so that all four say the same thing.
function materialise(raster::ClimateRaster, area::StudyArea; role = missing)
    return _rasternotaspec(raster)
end

# `_asraster` through a cache. The uncached methods stay as they are, so every existing caller is
# unaffected; these add the memoised path used while resolving a study area.
#
# A `ConstructedSpec` is deliberately *not* cached at this level — its `combine` is an opaque
# closure that two separately-built specs can never share, so a top-level entry could never hit.
# Its children are `SourceSpec`s and *are* cached, which is where the cost actually is.
#
# `cut`, when given, restricts the read to that box — the windowing described under `_readwindow`.
# It is folded into the **cache key**, not applied afterwards: a windowed read and a whole-world
# read of the same layer are genuinely different results, and returning one for the other would hand
# back the wrong ground. A `cut` the spec already carries always wins, so an explicit user window is
# never silently widened or narrowed.
function _asraster(raster::ClimateRaster, ::LayerCache; cut = nothing)
    return _rasternotaspec(raster)
end

function _asraster(spec::SourceSpec, cache::LayerCache; cut = nothing)
    kw = (isnothing(cut) || haskey(spec.readkw, :cut)) ? spec.readkw :
         merge(spec.readkw, (cut = cut,))
    spec.code isa AbstractVector && return _stackcached(spec, cache, kw)
    return get!(cache.reads, ReadKey(spec.source, spec.code, kw)) do
        return _read(spec; kw...)
    end
end

# The cached twin of `_asraster(::Tuple)` — a `(source, code)` pair, refused with the spelling
# that replaces it. Both entry points must refuse it, or deciding a grid and building a layer would
# disagree about what is accepted.
_asraster(spec::Tuple, ::LayerCache; cut = nothing) = _sourcepairnotaspec(spec)

# **Two passes, because a synthetic layer has no grid of its own to be read on.** The data layers
# are read first and must agree (`_checksourcegrid`); the synthetic ones are then generated **at the
# grid they agree on**. That is the whole of the rule: a generated layer needs shape and orientation,
# never coordinates, so it can adopt any grid it is given — it just cannot supply one.
#
# **All-synthetic is refused, and it is the one case that cannot work.** There is no grid to adopt
# and none can be invented; the combine has nowhere to happen. Caught here rather than at
# construction because a `ConstructedSpec` is also reached this way by `_analyse`, so this is the
# choke point both routes pass through.
#
# **Order is preserved throughout.** The layers are positional arguments to a combine that is *user*
# code, so the array handed to it must be in the order they were written, not data-first.
function _asraster(spec::ConstructedSpec, cache::LayerCache; cut = nothing)
    issynth = map(_issyntheticspec, spec.layers)
    # `!isempty` first, and it is load-bearing: `all` of an **empty** collection is `true`, so a
    # nullary thunk — a spec with no layers at all, which is what `in_memory_raster` builds — would
    # otherwise be refused for having no grid when it supplies its own raster outright.
    !isempty(issynth) && all(issynth) &&
        error("a `ConstructedSpec` whose layers are all synthetic has no grid to be combined on: " *
              "a generated layer takes its shape from a grid rather than supplying one. Give it at " *
              "least one data layer, or leave the spec on the default `CombineOnTargetGrid`, where " *
              "the study area's grid decides the shape.")
    layers = Vector{Any}(undef, length(spec.layers))
    for (i, l) in pairs(spec.layers)
        issynth[i] || (layers[i] = _asraster(l, cache, cut = cut))
    end
    readidx = [i for i in eachindex(layers) if !issynth[i]]
    read = [layers[i] for i in readidx]
    # The indices travel with the rasters: the error names a layer, and after filtering out the
    # synthetic ones a position in `read` is no longer a position in `spec.layers`.
    _checksourcegrid(spec, read, readidx)
    if any(issynth)
        grid = _owngrid(first(read))
        for (i, l) in pairs(spec.layers)
            issynth[i] && (layers[i] = _materialiseon(l, grid, cache))
        end
    end
    return _combined(spec.combine(layers...), spec)
end

# A multi-layer read, assembled from **per-layer** cache entries rather than cached whole.
#
# The point is that the entries are shared. Keyed on the whole request, a study area whose regime
# is `ConstructedSpec(EarthEnv{LandCover}) do lc … end` and whose supply is
# `SourceSpec(EarthEnv{LandCover}, 7)` would read file 7 twice: the keys differ, so neither hits.
# Keyed per layer, the second request finds what the first already read. Narrowing the *key*
# instead would be wrong — it must describe what was asked for, or one request is served another's
# data — so the fix is to decompose the read, which is only possible because a spec's code list is
# now explicit rather than a `nothing` sentinel.
#
# Only reached for a spec whose layers can share one array; `_read` refuses the rest, and this
# path goes through `_read` per layer, so each keeps its own correct unit before stacking.
function _stackcached(spec::SourceSpec, cache::LayerCache, kw::NamedTuple)
    codes = spec.code
    layers = map(codes) do c
        one = SourceSpec(spec.source, c; kw...)
        return get!(cache.reads, ReadKey(spec.source, c, kw)) do
            return _read(one)
        end
    end
    length(layers) == 1 && return only(layers)
    stacked = _stacklayers([l.array for l in layers],
                           Dim{:layer}(map(c -> _axisname(spec.source,
                                                          c),
                                           codes)))
    return ClimateRaster(spec.source, stacked, collect(codes))
end

# The name a layer appears under on the stacked `:layer` axis. Not the code itself: assembling the
# stack here rather than in `_readmultilayer` must not relabel it, and `getraster` keys its layers by
# canonical *name*, so `EarthEnv{LandCover}` code `7` has always shown as `:cultivated_and_managed`.
# The shipped table lists both spellings as aliases, numeric first, so the name is the non-numeric one
# — falling back to the code for a layer that has only the one alias.
function _axisname(source, code)
    named = filter(a -> !all(isdigit, a),
                   layerinfo(source, code).aliases)
    return isempty(named) ? code : Symbol(first(named))
end

# Put a **child** spec's data on `target`, combining a nested `ConstructedSpec` **after** its own
# children are sampled rather than before.
#
# **This is the `ConstructedSpec` sub-layer path, and only that** — the three kinds a spec's *layers*
# can be. The builder itself routes through `materialise`, so nothing here duplicates it. It takes a
# bare `target` grid and a cache rather than a `StudyArea`, which is why it cannot simply *be*
# `materialise`: on the `CombineOnSourceGrid` path the grid is the layers' own, and `_analyse` reaches
# it before any study area exists.
#
# **The ordering is the whole point.** Sampling first and combining second interpolates the
# percentages, which is what they are for, and takes the argmax on the target grid. Combining first
# would have `compress_landcover` argmax EarthEnv's per-class percentage bands into class codes and
# then leave the resampler interpolating *between codes*, which is meaningless. This is not a
# workaround but the better numerics: it keeps the sub-cell mix that interpolating codes destroys and
# that nearest-class discards.
#
# The children are handed over as `ClimateRaster`s; only their **grid** changes.
# That is deliberate: `compress_landcover(lc)`, `landcoverclass` comparisons and every user combine
# keep working unchanged, where passing bare arrays would have broken all of them.
#
# Valid only for a combine that **commutes with regridding**, which is what the default assumes.
# Cell-wise is necessary but not sufficient: a combine that looks at neighbours, crops or aggregates
# fails it outright, and so does a cell-wise but *nonlinear* one, since the ratio of two
# interpolations is not the interpolation of the ratio. Either says so with
# `combinestage = CombineOnSourceGrid()`, and `_combineon` below runs it on the layers' own grid
# instead.
function _materialiseon(spec::SourceSpec, target, cache::LayerCache)
    read = _asraster(spec, cache)
    return ClimateRaster(spec.source,
                         _sampledata(read, target, name = "layer",
                                     categorical = iscategorical(read,
                                                                 _specaxis(spec))),
                         spec.code)
end

# **A synthetic spec on a positioned grid — generated at the target's shape, not sampled onto it.**
# A synthetic layer needs **shape and orientation**, never coordinates, so there is nothing about a
# real-world grid it cannot be built on: `examples/ScottishCultivatedLand.jl` puts a synthetic solar
# gradient over a real land-cover grid.
#
# **The masking question that made it look risky resolves itself.** A regime normally narrows
# `active` by where its data actually exist (`_coverage`), and a generated layer has no gaps to
# narrow by — so it contributes nothing, and the area's own `active` (plus any data-backed supply's
# coverage) still decides. That is the correct semantics rather than a concession: a uniform layer
# genuinely does cover every cell.
#
# Wrapped as a code-less `ClimateRaster` on the target's own dims, because that is what a combine
# takes; the spec's declared axis is what gives it meaning, exactly as for a `ConstructedSpec`.
function _materialiseon(spec::AbstractSyntheticLayerSpec, target,
                        ::LayerCache)
    yx = dims(target, (Y, X))
    # `_specfield`, not `_syntheticsupplyfield`, so the **whole** synthetic family is reachable —
    # including `NicheSpec`, which the latter answers only through its `fill(float(spec), …)`
    # fallback and so dies on `float(::NicheSpec)`. Sharing the field generator covers a sub-layer
    # and a top-level one alike.
    field = _specfield(spec, length.(yx), _rowsincreasenorth(yx))
    # **United, CRS-bearing dims, not the target template's own.** `target`'s lookups are
    # deliberately bare numbers, but `_cropto` reads a raster's coordinates with `ustrip(crsunit, …)`
    # — so a raster carrying bare dims fails on the unit rather than being recognised as already on
    # the grid. Rebuilt in the same form a reprojected raster has, which makes the subsequent crop a
    # pure index.
    # **Shared with `_materialisefield` deliberately.** As two near-identical blocks they would
    # diverge silently over whether to unite their dims — see `_unitedyx`. One function, so it cannot
    # happen
    # The sampling comes from the template rather than being stated here, or a synthetic raster
    # would disagree with the grid it was built on.
    united = _unitedyx(yx, Rasters.crs(target))
    # **`SyntheticData`, because that is what this is**: the field is generated by `_specfield` from
    # a spec, with no dataset, no layer and nothing read. Naming a real dataset here would be a claim
    # the values cannot support.
    return ClimateRaster(SyntheticData, DimArray(field, united))
end

# **Sampled afterwards, and that is not redundant** — the same reason `_materialisefield` does it.
# A spec with **no** layers is a nullary thunk (what `in_memory_raster` builds): `_combineon` has
# nothing to put on the target, so it returns the combine's own raster at its own resolution. As a
# *layer of another spec* that matters, because the enclosing combine broadcasts its layers together
# — measured, a `_reg(raster)` layer at 4×4 against a synthetic one generated at the target's 2×2
# is a `DimensionMismatch`. Where the layers were already put on the target, this is a no-op:
# `_cropto` recognises a raster on the target grid and crops instead of resampling.
function _materialiseon(spec::ConstructedSpec, target, cache::LayerCache)
    out = _combineon(spec.combinestage, spec, target, cache)
    return ClimateRaster(_sourceof(out),
                         _sampledata(out, target, name = "layer",
                                     categorical = iscategorical(out,
                                                                 _specaxis(spec))),
                         out.code)
end

# The `(Y, X)` dims of a purely synthetic grid: real coordinates from an origin of `(0, 0)`, in the
# cell size's own length unit, `Intervals(Start)` and a `Regular` span.
#
# **A synthetic grid has no *position*, but it does have real *spacing*.** Cells 1 and 3 really are
# two cell-widths apart, and that follows from `cellsize`, which the area knows. Only the origin is a
# convention — and `(0, 0)` is a stated one, not an invented fact.
# **`Sampled`, not `Projected`, and deliberately so**: there is genuinely no CRS
# (`NoRealWorldPosition`), and claiming one would be the lie this avoids. It costs nothing — a
# synthetic grid is never resampled, which is the only thing a CRS is needed for.
# **`Start`, so the grid runs `[0, extent)`** — with `Center` the same labels would give
# `(-cell/2, extent - cell/2)`, i.e. the grid would start half a cell before its own origin.
# Two independently built areas of the same extent and cell size produce **identical** dims, so their
# layers compare equal; and a `km`-labelled grid compares equal to the same grid in `m`, so the unit
# costs nothing.
function _syntheticyx(geometry::NamedTuple)
    s = geometry.cellsize
    return map((geometry.dim, (Y, X))...) do n, D
        return D(DimensionalData.Lookups.Sampled((0:(n - 1)) .* s,
                                                 order = DimensionalData.Lookups.ForwardOrdered(),
                                                 span = DimensionalData.Lookups.Regular(s),
                                                 sampling = DimensionalData.Lookups.Intervals(DimensionalData.Lookups.Start())))
    end
end

# The area's own `(Y, X)` dims, as a materialised layer carries them.
#
# **Just the grid's own dims.** A unitless template would need its unit *put back* on the way out,
# in one place per path — and two copies of one re-attachment, disagreeing about whether to run, is
# exactly the kind of drift that costs cells silently. The grid carries its unit instead, and the
# strip happens once, in
# `_reproject`, at the only boundary that cannot take a unit.
# Kept as a named function rather than inlined at its three call sites: it says *why* those dims
# are the right ones, and it is the natural place to reinstate a check if the grid ever stops being
# the single source of truth.
_unitedyx(yx, tcrs) = yx

# Late: every layer is put on the target grid, and the combine runs there. The default, and today's
# behaviour exactly — see `_materialiseon` above for why it is the right ordering whenever the
# combine commutes with regridding.
function _combineon(::CombineOnTargetGrid, spec::ConstructedSpec, target,
                    cache::LayerCache)
    # Nothing to stamp on the result: what it *is* comes from `spec.axis`, and the callers that
    # need to know ask `iscategorical(raster, axis)`. On this path the layers were sampled before
    # the combine ran, so no resampling decision is left to make here at all.
    return _combined(spec.combine(map(l -> _materialiseon(l, target, cache),
                                      spec.layers)...), spec)
end

# Early: the combine runs on the layers' own shared grid (which `_asraster` has just checked they
# share) and its *result* is what gets sampled.
#
# **The spec's axis has to travel with it**, and this is the path that needs it: the combine's
# result is sampled *here*, so whether it holds class codes decides `:mode` against `:bilinear`. A
# derived layer has no catalogue row, so the axis is the only thing that can say.
function _combineon(::CombineOnSourceGrid, spec::ConstructedSpec, target,
                    cache::LayerCache)
    return _sampledeclared(_asraster(spec, cache), target, spec.axis)
end

# The WGS84 box worth reading, or `nothing` to read everything. `nothing` whenever the answer is not
# certain — no mask, a mask that cannot state an extent, or a target CRS that could not be probed —
# because reading too little silently truncates the study area, while reading too much only costs
# time.
# Returns the *prepared mask* alongside the window, so it can be reused rather than prepared twice.
# Preparing a `ShapeSpec` downloads and reads a shapefile and a `ConstructedSpec` reads its source
# layers, so doing it once is the whole reason `_preparemask` exists as a separate step — computing
# the window must not quietly undo that.
function _readwindow(within, layers::NamedTuple, crs)
    probed = _probetargetcrs(layers, crs)
    (isnothing(within) || isnothing(probed)) &&
        return (cut = nothing, prepared = nothing, tcrs = probed)
    prepared = try
        _preparemask(within, probed)
    catch                    # a mask that cannot be prepared yet is simply not a window
        nothing
    end
    (isnothing(prepared) || isnothing(prepared.extent)) &&
        return (cut = nothing, prepared = prepared, tcrs = probed)
    return (cut = _padded(prepared.extent, probed), prepared = prepared,
            tcrs = probed)
end

# `extent` (in `tcrs`) as a padded WGS84 box for `cut`, or `nothing` when there is nothing to gain.
#
# Longitude is clamped to *just inside* -180°: `LatLong` accepts (-180°, 180°], half-open, so a
# box padded out to exactly -180° is rejected outright. A globe-spanning mask hits this immediately —
# and is also precisely the case where a window saves nothing — so it declines instead, leaving the
# ordinary unwindowed read to do what it already does well.
function _padded(extent, tcrs)
    wgs = _bboxin(tcrs, Rasters.EPSG(4326), extent)
    ylo, yhi, xlo, xhi = _extentvalues(wgs)
    la1, la2 = ustrip(°, ylo), ustrip(°, yhi)
    lo1, lo2 = ustrip(°, xlo), ustrip(°, xhi)
    dla, dlo = _WINDOW_PAD * (la2 - la1), _WINDOW_PAD * (lo2 - lo1)
    la1, la2 = max(-90.0, la1 - dla), min(90.0, la2 + dla)
    lo1, lo2 = max(nextfloat(-180.0), lo1 - dlo), min(180.0, lo2 + dlo)
    (la2 - la1 >= 179.0 && lo2 - lo1 >= 359.0) && return nothing
    return Extents.Extent(Y = (la1 * °, la2 * °), X = (lo1 * °, lo2 * °))
end

# Materialise every named layer through the cache, dropping the ones not given and the ones that
# cannot shape a grid. Order matters: it is the tie-break for `_choosealign`, so `regime`
# deliberately comes before `supply`.
function _materialiselayers(layers::NamedTuple, cache::LayerCache;
                            cut = nothing)
    return [(name = nm, raster = _asraster(s, cache, cut = cut))
            for (n, spec) in pairs(layers) if !isnothing(spec)
            for (nm, s) in zip(_expandednames(n, spec), _expandspecs(spec))
            if _shapesgrid(s)]
end

# Decide the grid for a study area and report what it costs.
#
# This is the single engine behind both [`StudyArea`](@ref) and [`investigate_study_area`](@ref): the
# former builds from its result, the latter formats it, so a report can never describe a grid other
# than the one that would be built. Every argument is optional — the grid is decided from whatever is
# given — and `nothing` means "not specified" throughout (the `missing`/`nothing` distinction belongs
# to the public constructors, and is resolved before reaching here).
function _analyse(layers::NamedTuple; within = nothing, crs = nothing,
                  cellsize = nothing, extent = nothing, align = nothing,
                  simulate_safely = nothing,
                  cache::LayerCache = LayerCache())
    problems = Problem[]
    cons = (within = within, crs = crs, cellsize = cellsize, extent = extent,
            align = align, simulate_safely = simulate_safely)
    # `nothing` here is "not specified" (`_constraint` maps both `missing` and a cleared value onto
    # it), and the default this resolves to is `true` — never simulate ground the data does not
    # describe unless asked to. `cons` keeps it *as given*, like every other constraint.
    safely = isnothing(simulate_safely) ? true : simulate_safely
    # Window the reads before making them: only the ground this area could possibly use is read.
    window = _readwindow(within, layers, crs)
    materialised = _materialiselayers(layers, cache, cut = window.cut)

    # A synthetic area — no layers at all — is decided entirely by `extent` and `cellsize`.
    isempty(materialised) &&
        return _syntheticplan(extent, cellsize, within, problems, layers, cons,
                              safely, cache)

    # **`extent` is meaningless once data positions the grid, and must be refused rather than
    # ignored.** It says how *big* an area is and never *where*, so with layers present the extent
    # already comes from the data, optionally narrowed by `within`, and a second unanchored size
    # cannot be reconciled with it. Silently dropping it is worse than it sounds:
    # `StudyArea(supply = <a global dataset>, crs = EPSG(27700), cellsize = 25km,
    # extent = (350km, 500km))` then gives the **identical** grid to the same call without `extent` —
    # 953 by 6 cells, northings from −19 315 km, far outside the projection's valid domain, and
    # noticed. Refusing names the argument that *does* position an area, which is the thing the
    # caller actually wanted.
    isnothing(extent) ||
        error("`extent` sizes a **synthetic** area and cannot be combined with data layers: those " *
              "already carry their own extent, and `extent` says how big an area is, never where " *
              "it is. Use `within` to position and limit it — a bounding box " *
              "(`EcoSISTEM.boundingbox(\"Scotland\")`), a shape file, or a mask spec — or drop " *
              "the layers to decide a synthetic area from `extent` and `cellsize` alone.")

    rasters = Tuple(m.raster for m in materialised)
    # Stage 1, unchanged from the existing engine: an explicit `crs` wins, else a single projected
    # input CRS is adopted, else the first layer's own.
    tcrs = _targetcrs(rasters, crs, cellsize)
    crssource = isnothing(crs) ? AdoptedFromLayers() : GivenByUser()
    facts = [_layerfacts(m.name, m.raster, tcrs) for m in materialised]

    # Stage 2: which layer is preserved exactly. An explicit `align` names a layer; otherwise the rule
    # is "whichever is already in the target CRS" (finest first, then given order).
    chosen = isnothing(align) ? _choosealign(facts, tcrs) :
             _namedalign(facts, align)
    if isnothing(chosen)
        push!(problems,
              Problem(ProblemWarning(), :no_alignment,
                      "no layer is in the target grid's CRS, so no layer can be kept exactly — " *
                      "every one of them will be resampled. Choosing a `crs` that matches one of " *
                      "them would avoid that."))
    end

    # Stage 3: the cell size. An explicit `cellsize` wins; else the alignment layer's own step; else
    # the existing unanimity rule (which fails closed on disagreement).
    cs, cssource = if !isnothing(cellsize)
        (cellsize, GivenByUser())
    elseif !isnothing(chosen) && !isnothing(chosen.step)
        (chosen.step, TakenFromAlignedLayer())
    else
        _roundedcellsize(_inferredcellsize(rasters, tcrs)..., problems)
    end

    # Stage 4: the extent. The mask leads where it can state one, intersected with every layer's
    # footprint, then snapped outwards onto the alignment layer's cell boundaries so the alignment is
    # real rather than nominal.
    # Reuse the mask `_readwindow` already prepared, when it prepared one against the CRS we
    # actually settled on. Preparing a `ShapeSpec` reads a shapefile and a `ConstructedSpec` reads its
    # source layers, so doing it twice would undo a good part of what the windowing just saved.
    prepared = (!isnothing(window.prepared) && _samecrs(window.tcrs, tcrs)) ?
               window.prepared : _preparemask(within, tcrs)
    bounds = _planbounds(facts, prepared.extent, tcrs, cs, problems)
    # The alignment layer's own lower-left corner, which `_snapbounds` snaps onto. Reading it by
    # name rather than as elements 1 and 3 of a tuple is the point of the extent type.
    origin = isnothing(chosen) ? nothing :
             (chosen.bounds.Y[1], chosen.bounds.X[1])
    snapped = _snapbounds(bounds, origin, cs)
    grid = _targetgrid(chosen, tcrs, snapped, cs)

    # Stage 5: sample the layers on, take the coverage, apply the mask, and recut to what is genuinely
    # usable — the same rules the builder applies, run here so the reported grid is the final one.
    covered = _coveredgrid(rasters, prepared.payload, grid,
                           simulate_safely = safely)

    plans = [_planfor(f, tcrs, cs, origin) for f in facts]
    fp = _gridfootprint(Base.size(covered.active)...)
    _collectproblems!(problems, plans, tcrs, covered.active, fp, chosen, bounds)
    _masklost!(problems, covered.mask, covered.active, cs)
    _partialcover!(problems, covered, cs, safely)
    return StudyAreaReport(tcrs, crssource, cs, cssource,
                           isnothing(chosen) ? nothing : chosen.name,
                           _activegrid(covered.grid, covered.active), safely,
                           plans, fp, problems, layers, cons, cache,
                           AsInvestigated())
end

# A data-driven spec is read (through the area's cache) and sampled onto the grid; a synthetic one is
# generated at the grid's shape. Both end up on the area's own dims.
#
# Whether the values are **class codes** travels with them, because that is what decides the kind
# of layer they become and neither the values nor the spec can be asked afterwards: a synthetic
# `NicheSpec` and a categorical raster arrive by different routes at the same answer.
function _materialisefield(spec::AbstractSyntheticLayerSpec, area::StudyArea)
    grid = area.report.active
    yx = dims(grid, (Y, X))
    field = _specfield(spec, length.(yx), _rowsincreasenorth(yx))
    # **United dims, so a synthetic layer and a data layer on one area agree** — see `_unitedyx`.
    return (values = DimArray(field, _unitedyx(yx, area.report.crs)),
            categorical = spec isa NicheSpec)
end

function _materialisefield(spec, area::StudyArea)
    raster = _asraster(spec, area.report.cache)
    categorical = iscategorical(raster, _specaxis(spec))
    values = _sampledata(raster, area.report.active, name = "layer",
                         categorical = categorical)
    return (values = _restricttocovered(values, raster, area, categorical),
            categorical = categorical)
end

# **A `ConstructedSpec` goes through `_combineon`, so inspection obeys `combinestage`.**
#
# The generic method above reads every layer on its **own** grid, combines there, and samples the
# result — which is `CombineOnSourceGrid`'s behaviour. Applied to a spec that says
# `CombineOnTargetGrid`, that is the opposite ordering to the one the builder uses, and the two give
# different numbers for any combine that does not commute with regridding. `materialise` would then
# show a layer the model would never use, which is the one thing an inspection function must not do.
#
# **`_combineon`'s result still needs sampling, and that is not belt-and-braces.** A spec with no
# layers — a nullary thunk, which is what `in_memory_raster` builds — has nothing for `_combineon` to
# sample, so it hands back the combine's own raster on whatever grid that raster has. Measured: a
# `_reg(raster)` layer came out 4×4 on a 2×2 study area when this step was omitted. The builder has
# always sampled afterwards for the same reason; this mirrors it.
# Where the layers *were* sampled first, the second pass is a no-op: `_cropto` recognises a raster
# already on the target and crops rather than resampling.
function _materialisefield(spec::ConstructedSpec, area::StudyArea)
    out = _combineon(spec.combinestage, spec, area.report.active,
                     area.report.cache)
    categorical = iscategorical(out, _specaxis(spec))
    values = _sampledata(out, area.report.active, name = "layer",
                         categorical = categorical)
    return (values = _restricttocovered(values, out, area, categorical),
            categorical = categorical)
end

# **This is what carries `simulate_safely` into `GridHabitat`** (user, 2026-08-13: *"it
# should apply equally if possible"*). A layer named on the `StudyArea` already shaped the grid, so
# every remaining cell is covered and this does nothing; a layer introduced only at build time never
# faced that test, and blanking the cells it does not reach makes `_coverage` — which already reads
# `NaN` — mark them inactive. No new plumbing, and `materialise` still shows exactly what is built.
#
# **A categorical layer is the documented exception**: its class codes are integers and cannot
# carry `NaN`, which `_coverage` already notes (such a layer always reports full coverage). Harmless
# for a layer named in the area — the grid was cropped to it — so only one introduced at build time
# escapes the rule. Recorded as a limitation rather than answered with a parallel mask.
function _restricttocovered(values, raster, area::StudyArea, categorical::Bool)
    (categorical || !_cannan(eltype(values)) ||
     !area.report.simulate_safely) && return values
    return _blankuncovered!(values, _fullycovered(raster, area.report.active))
end

# The role-specific final step, and only that. `missing` and `Condition` both give a **regime** —
# `Condition` canonicalises to the axis's unit first, `missing` leaves the values as written — and
# `Resource` gives a `Supply`.
#
# **This IS the builder's final step, not merely the same shape as it.** `GridHabitat` calls
# `materialise`, so there is one chain and nothing left for inspection and building to agree about.
function _applyrole(f::NamedTuple, ::Missing, axis, area::StudyArea)
    return _asregime(f.values, f.categorical, axis)
end

function _applyrole(f::NamedTuple, ::Type{Condition}, axis, area::StudyArea)
    f.categorical &&
        return _asregime(f.values, true, axis)
    return _asregime(_canonical.(f.values, Ref(axis)), false, axis)
end

function _applyrole(f::NamedTuple, ::Type{Resource}, axis, area::StudyArea)
    return _wrapsupply(f.values, _inspectioncellareas(area), axis)
end

# The area of each of `area`'s cells, as `cancel` needs it to turn a per-area rate into a per-cell
# one: the nominal `cellsize^2`, scaled by the latitude factor where there is a latitude.
#
# The nominal area stays the **report's** own `cellsize^2`, rather than being recomputed from the
# grid as `_cellareas` would: the report's is the cell size this area was actually decided on, and
# keeping it means the only thing the latitude correction moves is the factor itself.
#
# **No guard is needed for a synthetic area, because `_areafactor` dispatches on the coordinates.**
# `_areafactor(::AbstractVector{<:Unitful.Length}) = 1.0`, and a synthetic grid's coordinates are
# lengths (`km`) exactly as a projected grid's are (`m`) — so both take that method and only a
# genuinely angular grid reaches the latitude correction.
#
# **No `isnothing(area.report.crs)` guard is needed, and adding one would be weaker than what is
# here.** Both this and `_buildhabitat` read the **coordinates themselves** to decide whether a grid
# is angular, so inspection and building cannot disagree about the kind of area; a guard on the CRS
# would be a second rule that has to be remembered and kept in step.
# **This one stays metric, and that is not the same question as the layer's `size`.** A layer's
# `size` is the grid's own cell size — an angle on a geographic grid, and reported as such. A cell's
# **area** here is a real *physical* area, because that is what `cancel` needs to turn a per-area rate
# into a per-cell one.
# **Squaring the report's `cellsize` on a geographic grid would give `degree²`, and the reason that
# is wrong is NOT that a `degree²` is meaningless** — square degrees are a perfectly good unit of
# angular area. It is that `Δlat · Δlong` is a **coordinate extent**, not a physical area and not even
# the cell's solid angle: the true solid angle is `Δλ · (sin φ₂ − sin φ₁)`, which *varies with
# latitude* while the coordinate product does not. `_areafactor` is precisely the ratio between the
# two — its numerator is the solid-angle term — so the correction this function applies is exactly
# what stands between the coordinate product and a real area. See `getcellareas`, which offers both
# readings and refuses to conflate them.
function _inspectioncellareas(area::StudyArea)
    yx = _cellcentres(area.report.active)
    nominal = (area.report.cellsize isa Unitful.Length ?
               area.report.cellsize :
               _cellsize(yx.lat, yx.long))^2
    # The latitude correction is asked for the cells' **own top and bottom**, so
    # it no longer reconstructs them as `centre ± half` from a differenced step — exact, and blind to
    # how the lookup happens to label its cells.
    return nominal .* _areafactor(_cellintervals(area.report.active, Y))
end

# The environment on an already-decided grid. Nothing here chooses a grid: the area did that, and this
# only samples layers onto it — which is the whole point of separating the two.
#
# **The frame is fixed; only activity can shrink.** A layer that was not named when the area was
# decided may still have gaps, and those gaps make cells inactive — but they must never move or resize
# the grid, or the area that was inspected would not be the area that gets built. Such a layer is
# reported rather than silently narrowing the result.
# Strip the declared changes, build the habitat exactly as before, then apply them to the layers
# that were built. Splitting here — before the synthetic/data fork on the next line — is what keeps
# every path downstream free of wrappers: `materialise` can never see one, on either fork, so
# nothing below here needs a guard.
# **Returns assembled PARTS, not a habitat**, because `GridHabitat`'s sole constructor is now the
# specs-and-area one and `new` has to be the last thing that happens. The declared changes are
# applied here, to the layers, before anything is constructed from them.
#
# **The supply's declared change is applied by the constructor, AFTER `_zerogaps`, not here** —
# and the asymmetry is deliberate. `_zerogaps` *rebuilds* a supply layer that has gaps, passing its
# change through `_cleanstored`; `_hasgaps` recurses into a `SumOfLayerChanges`, so attaching the
# declared change first would let `_zerogaps` fire on a gap-free data supply merely because a
# *declared* series carried `NaN`. That may well be the better behaviour, but it is a behaviour
# change in the supply path and belongs to whoever decides it, not to a constructor move.
function _buildonarea(regime, supply, area::StudyArea,
                      topology::AbstractTopology = Island())
    reg = _splitvarying(regime)
    sup = _splitvarying(supply)
    parts = _buildhabitat(reg.spec, sup.spec, area, topology)
    _applydeclared!(parts.regime, reg.change)
    return (regime = parts.regime, active = parts.active,
            supply = parts.supply, supplychange = sup.change,
            problems = parts.problems)
end

# **Both branches build their layers with `materialise`**, the very function a user calls to see what
# the builder will get, so its documented promise holds by construction rather than by agreement. A
# near-copy of that chain here would drift — over uniting dims, over the sampling to declare, over
# which spec kinds it handles — and each drift would be found by accident.
#
# What is left here is what is genuinely the *builder's* and not inspection's: the multi-layer
# arity check, the coverage AND that turns data gaps into inactive cells, and the assembly into a
# `GridHabitat`.
function _buildhabitat(regime, supply, area::StudyArea,
                       topology::AbstractTopology = Island())
    isnothing(area.report.crs) &&
        return _syntheticonarea(regime, supply, area, topology)
    (regime isa Union{Tuple, NamedTuple} && length(regime) < 2) &&
        error("a tuple `regime` describes a multi-layer environment, so it needs at least 2 " *
              "elements; got $(length(regime)). Pass the spec on its own for a single layer.")
    # `role`, because the builder canonicalises: a regime written in °C must reach the layer in K,
    # and a supply's per-area rate must become the absolute per-cell one the simulation divides a
    # demand by. The role-free `materialise(spec, area)` deliberately does neither.
    reg = materialise(regime, area, role = Condition)
    sup = materialise(supply, area, role = Resource)
    # The supply's coverage counts as much as the regime's: a data-driven supply with a hole leaves
    # cells that cannot be simulated just as surely as a regime with one. A synthetic layer of either
    # role has no gaps, so it never restricts.
    active = Matrix{Bool}(parent(area.report.active) .& _coverage(reg) .&
                          _coverage(sup))
    return (regime = reg, active = active, supply = sup,
            problems = _warnextralayers(area, active))
end

# A synthetic area has no CRS and no real-world position, so its layers are generated at its shape
# rather than sampled onto it — which only a synthetic spec, or a layer already built on this very
# grid, can do.
function _syntheticonarea(regime, supply, area::StudyArea,
                          topology::AbstractTopology = Island())
    specs = _expandspecs(regime)
    (regime isa Union{Tuple, NamedTuple} && length(specs) < 2) &&
        error("a tuple `regime` describes a multi-layer environment, so it needs at least 2 " *
              "elements; got $(length(specs)). Pass the spec on its own for a single layer.")
    # **An already-BUILT layer is now admitted here, and this is a CORRECTION of an expired
    # premise rather than a reversal of the guard 6d tightened.** The old `_refusebuilt` turned both
    # roles away outright, on the ground that this builder *generates* its layers at the area's shape
    # and a synthetic area had nothing for a layer to have been built on. That was true when a
    # synthetic grid carried `NoLookup` dims; `[S4]` gave it real `Sampled`/`Intervals(Start)`
    # coordinates, so a built layer now either sits on this exact grid — **provably**, by
    # `_yxcompatible` — or it does not, and `materialise`'s `AbstractLayer` method decides precisely
    # that, refusing a mismatch by size and copying the rest in. So there is nothing left for a
    # separate guard to say that is not said better one layer down.
    # **The general claim below is untouched**: arbitrary data still cannot be sampled onto a
    # synthetic grid. What changed is only the case where the values are *already on it*, which is
    # what lets `build_habitat` reseed a habitat from its own as-built area.
    #
    # The check runs over the *expanded* specs, never the `regime` as passed. `_shapesgrid` is
    # written for a single spec — its `::Any` fallback says "data-backed" — so handing it a
    # multi-variable tuple classified an entirely synthetic regime as data-backed and refused it with
    # an error blaming the caller for something they had not done. The data-driven branch of
    # `_buildhabitat` already expands first for exactly this reason.
    # A built layer reaches that same `::Any` fallback, so it has to be excluded here explicitly:
    # without this, relaxing the rule above would achieve nothing, because a built layer would simply
    # be refused by the next line instead.
    any(s -> !(s isa AbstractLayer) && _shapesgrid(s), specs) &&
        error("a data-backed `regime` cannot be built on a synthetic study area: the area has no " *
              "CRS, so there is nothing to say where the data's values fall on it. Name the layer " *
              "when deciding the area (`StudyArea(regime = …)`) so it takes its grid from the " *
              "data. A layer that is *already built* is accepted, provided it was built on this " *
              "very grid.")
    # **The same `materialise` a user calls**, with `role = Condition` because the builder
    # canonicalises: a spec written in °C must reach the layer in K, which the role-free form
    # deliberately does not do.
    reg = materialise(regime, area, role = Condition)
    # `copy`, and it is load-bearing: a `StudyArea` is reusable by design, so passing its own
    # `active` straight through gave **every environment built from one area the same array object**.
    # An intervention that deactivates cells in one ecosystem then silently deactivated them in every
    # other — and in the area itself, so the next `GridHabitat` inherited the damage. Harmless
    # until something mutated `active`, which is why it survived until interventions arrived; the
    # data-driven branch above never had it, because `Matrix{Bool}(… .& …)` allocates.
    # **The supply goes through the same `materialise` the regime does**, so that a synthetic
    # environment's regime and supply describe their shared grid the same way. Two routes could agree
    # on every *value* and still disagree on the **dims**, and that stays invisible while either side
    # carries `NoLookup`, because `_yxcompatible` then falls back to a size-only check.
    # `materialise` handles a `Tuple`/`NamedTuple` supply too, returning a `LayerCollection` — the
    # same shape `_synthetic_env`'s multi-supply method built by hand.
    sup = materialise(supply, area, role = Resource)
    # A synthetic area has no data footprint, so no layer can cost it a cell — hence no problems.
    # **`Matrix{Bool}`, not the report's `BitMatrix`.** The report's mask is never written, so a
    # `BitMatrix` is free there; the habitat's is **live** — `Deactivate`/`Reactivate` write it — and a
    # `BitArray` element write is a read-modify-write of a whole `UInt64` word. Letting the bit
    # backing through here also changed simulation results, which `extras_canonical` caught.
    return (regime = reg, active = Matrix{Bool}(parent(area.report.active)),
            supply = sup,
            problems = Problem[])
end

# Say so when the layers handed to the builder cost activity the study area did not know about — the
# safe half of the frame-vs-activity rule, made visible.
function _warnextralayers(area::StudyArea, active)
    lost = count(area.report.active) - count(active)
    lost > 0 || return Problem[]
    message = "$lost cell$(lost == 1 ? "" : "s") active in the study area have no data in the layers " *
              "passed to `GridHabitat`, so they are inactive in the simulation" *
              (area.report.simulate_safely ?
               " (`simulate_safely` counts a cell the layer only partly covers as having no data)" :
               "") *
              ". The grid itself " *
              "is unchanged — a layer not named when the area was decided can only remove cells, never " *
              "move or resize it. Name it in `StudyArea(...)` if it should shape the grid."
    @warn message
    # **Returned as well as warned.** This is the one thing the builder learns that the study area
    # could not know, and a line on stderr is seen once and then gone. Returning it lets the habitat's
    # own report carry it, so a habitat can say *why* it has fewer active cells than the area it was
    # built on.
    return [Problem(ProblemWarning(), :layers_lost_cells, message)]
end
