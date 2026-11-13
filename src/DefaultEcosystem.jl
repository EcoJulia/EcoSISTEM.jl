# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Where a builder's unnamed inputs come from: `DefaultEcosystem()`, every value the builders fill in
# for a first run - a small synthetic toy grid and its species - or an existing habitat, whose own
# layers, grid and topology a rebuild carries across. Each value is announced as it is chosen. The
# builders in `actions.jl` ask for them; the values live here, and nowhere else.

"""
    DefaultEcosystem

Marker selecting opt-in defaults for the builders. `build_habitat(DefaultEcosystem(); ...)`,
`build_species(DefaultEcosystem(); ...)` and `build_ecosystem(DefaultEcosystem(); ...)` fill any omitted
**required** input from a defaults table - announcing each choice with an `@info` - then delegate to
the strict builder. Use it to spin up a fully-parameterised runnable toy ecosystem, or to fill in
only the required inputs you do not care about while passing the rest explicitly.
"""
struct DefaultEcosystem end

# As a placement, the one a builder uses when given none: `BySupply()`'s default numbers, which the
# builder's announcement has already stated, so without `BySupply()`'s own notice of them.
function (p::DefaultEcosystem)(species::SpeciesList, habitat::GridHabitat,
                               context::PlacementContext)
    counted = _drawingcounts(p, species, habitat, context)
    return _placeweighted(counted,
                          _drawingweights(p, species, habitat, context),
                          species, habitat, context)
end

# == Functions ==================================================================================

# The default placement's counts: the default number of individuals for each active species, split
# across them.
function _drawingcounts(::DefaultEcosystem, species, habitat, context)
    k = count(_activespecies(species, context))
    counts = _splittotal(_default(Val(:individuals)).value * k, species,
                         context)
    return (counts = counts, problems = Problem[])
end

# The default placement's weights: `BySupply()`'s.
function _drawingweights(::DefaultEcosystem, species, habitat, context)
    return _drawingweights(BySupply(), species, habitat, context)
end

# The default placement, drawn straight into the abundances as `BySupply()` is.
function _placeinto!(eco::AbstractEcosystem, placement::DefaultEcosystem,
                     context::PlacementContext, verbosity::Symbol; clear::Bool,
                     keep::Bool = true, undated::Bool = false)
    return _placedirectly!(eco, placement, context, verbosity, clear = clear,
                           keep = keep)
end

_placementphrase(::DefaultEcosystem) = "DefaultEcosystem()"

# Announce a default that was filled in, unless silenced. Takes `StudyArea`'s vocabulary: `:silent`
# says nothing, anything else announces. Announcing is the point of these builders - a default you
# cannot see is the hidden automation this API exists to avoid - so `:silent` is for callers who
# have read them once already, such as an example running under a harness that checks nothing is
# emitted.
function _announce(verbosity::Symbol, msg::AbstractString)
    verbosity in (:silent, :normal, :verbose, :full, :debug) ||
        error("`verbosity` must be `:silent`, `:normal` or `:verbose` (aliases `:full`, `:debug`); " *
              "got `:$verbosity`.")
    verbosity === :silent || @info msg
    return nothing
end

# One keyword a builder was not given, from where it is **coming from**, announced as `builder` fills
# it in: one function for every builder and keyword, so a new source or a new keyword is one added
# method. `what` is the keyword, as a `Val` written literally at the call, so dispatch is resolved
# when the builder compiles. A builder asks for it in its keyword defaults, which are evaluated only
# for a keyword not given, so only what was filled in is announced.
function _getdefaultvalue(what::Val{K}, ::DefaultEcosystem, verbosity::Symbol,
                          builder::AbstractString) where {K}
    default = _default(what)
    phrase = haskey(default, :phrase) ? default.phrase :
             sprint(show, default.value)
    _announce(verbosity, "$builder: `$K` defaulted to $phrase.")
    return default.value
end

# A source nothing else handles. A fallback rather than a `Union` annotation because the open source
# is the whole point - a new kind of source must be one added method - but a bare `MethodError` on a
# private helper tells a caller nothing useful.
function _getdefaultvalue(::Val, source, ::Symbol, builder::AbstractString)
    return error("`$(first(split(builder, '(')))` cannot take defaults from a " *
                 "$(nameof(typeof(source))). Pass " *
                 "`DefaultEcosystem()` for the announced defaults, or an existing `GridHabitat` to " *
                 "rebuild from one - or name every input explicitly.")
end

# How each builder taking defaults names itself when it announces one.
const _HABITAT = "build_habitat(DefaultEcosystem())"

const _DEFAULTSPECIES = "build_species(DefaultEcosystem())"

const _POOLSPECIES = "build_species(pool)"

# --- Every value a builder fills in --------------------------------------------------------------
# Each its value, by keyword, and how an announcement describes it where its own `show` would not.

# Stated in °C because that is what the announcement says and what a reader recognises as a
# temperate default. `Temperature`'s canonical unit is K, so it converts affinely on the way in -
# 25.0 °C is 298.15 K, not 25 K of anything.
function _default(::Val{:regime})
    return (value = UniformSpec(25.0°C, axis = Temperature),
            phrase = "a uniform 25.0 °C temperature")
end

function _default(::Val{:supply})
    return (value = UniformSpec(1.0 * canonicalunit(Resource, SolarRadiation) /
                                km^2,
                                axis = SolarRadiation),
            phrase = "a uniform canonical solar rate")
end

# The synthetic toy grid is a sensible default only because the default `regime` and `supply` are
# synthetic too. An explicitly passed **data** regime still needs its own `area`: no toy grid could
# suit unseen data, and the builder will not guess one.
function _default(::Val{:area})
    return (value = StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                              verbosity = :silent),
            phrase = "a synthetic 10.0 km × 10.0 km grid of 1.0 km cells")
end

# Announced even though `GridHabitat`'s own default is this same `Island()`: every value this builder
# chose has to be visible, and one silent choice would be an exception a reader must know about.
function _default(::Val{:topology})
    return (value = Island(),
            phrase = "`Island()` - the grid's edges do not wrap")
end

_default(::Val{:numspecies}) = (value = 10,)

function _default(::Val{:tolerance})
    return (value = ToleranceSpec(Temperature, μ = 298.0K, σ = 2.0K),
            phrase = "a (298.0 K, 2.0 K) temperature niche")
end

function _default(::Val{:demand})
    return (value = DemandSpec(SolarRadiation,
                               10.0 * canonicalunit(Resource, SolarRadiation)),
            phrase = "10× the canonical solar rate")
end

function _default(::Val{:dispersal})
    return (value = 10.0km, phrase = "a mean distance of 10.0 km")
end

function _default(::Val{:movement})
    return (value = MovementSpec(InsectDispersal(),
                                 meandist = _default(Val(:dispersal)).value),
            phrase = "`InsectDispersal()` with a mean distance of 10.0 km")
end

_default(::Val{:birth}) = (value = 0.6 / year,)

_default(::Val{:death}) = (value = 0.6 / year,)

_default(::Val{:longevity}) = (value = 1.0,)

_default(::Val{:survival}) = (value = 0.2,)

# How many individuals of each active species a weighted placement given no numbers places.
_default(::Val{:individuals}) = (value = 1000,)

function _default(::Val{:placement})
    return (value = DefaultEcosystem(),
            phrase = "`DefaultEcosystem()`, $(_default(Val(:individuals)).value) individuals " *
                     "per active species, split among them at random and placed by supply")
end

# --- An existing habitat ----------------------------------------------------------------------
# The same four questions asked of an existing habitat, which answers with what it actually has.
# That is what makes `build_habitat(h, supply = ...)` a rebuild with one thing changed rather than a
# fresh start: everything not named is carried across, the grid included. The layers come back
# **built** and are handed to a `GridHabitat` on the habitat's own as-built area, which `materialise`
# copies rather than aliases. Nothing is resampled, because nothing moved.
function _getdefaultvalue(::Val{:regime}, source::GridHabitat,
                          verbosity::Symbol,
                          ::AbstractString)
    _announce(verbosity,
              "build_habitat(habitat): `regime` taken from the habitat as built.")
    return source.regime
end

function _getdefaultvalue(::Val{:supply}, source::GridHabitat,
                          verbosity::Symbol,
                          ::AbstractString)
    _announce(verbosity,
              "build_habitat(habitat): `supply` taken from the habitat as built.")
    return source.supply
end

function _getdefaultvalue(::Val{:topology}, source::GridHabitat,
                          verbosity::Symbol, ::AbstractString)
    _announce(verbosity,
              # Interpolated directly rather than through `nameof(typeof(...))`: `Island` and
              # `Torus` are aliases for a parametric `EdgeTopology`, so the type's name is
              # `EdgeTopology` while its `show` says `Island()`, which is what a reader wants.
              "build_habitat(habitat): `topology` taken from the habitat ($(source.topology)).")
    return source.topology
end

# A habitat hands back the grid it was **actually built on**, narrowing included, which is why the
# `AsBuilt` stage exists. Re-deriving from the report's `specs` would answer the original question
# and hand back the wider as-investigated mask, silently discarding whatever cells the build cost.
function _getdefaultvalue(::Val{:area}, source::GridHabitat, verbosity::Symbol,
                          ::AbstractString)
    _announce(verbosity,
              "build_habitat(habitat): `area` taken from the habitat's own as-built grid.")
    return StudyArea(source, verbosity = :silent)
end
