# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The marker that asks `build_habitat`/`build_species`/`build_ecosystem` for their own defaults —
# a small synthetic toy grid, announced value by value as it is chosen.

"""
    DefaultEcosystem

Marker selecting opt-in defaults for the builders. `build_habitat(DefaultEcosystem(); …)`,
`build_species(DefaultEcosystem(); …)` and `build_ecosystem(DefaultEcosystem(); …)` fill any omitted
**required** input from a defaults table — announcing each choice with an `@info` — then delegate to
the strict builder. Use it to spin up a fully-parameterised runnable toy ecosystem, or to fill in
only the required inputs you do not care about while passing the rest explicitly.
"""
struct DefaultEcosystem end

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════
#
# Filling in what a caller did not name. Every method dispatches on `DefaultEcosystem` or on the
# `GridHabitat` a rebuild takes its values from: the builders in `actions.jl` state the workflow, and
# the values behind it live here, with the type that selects them.

# Announce a default that was filled in, unless silenced. Takes `StudyArea`'s vocabulary: `:silent`
# says nothing, anything else announces. Announcing is the point of these builders — a default you
# cannot see is the hidden automation this API exists to avoid — so `:silent` is for callers who
# have read them once already, such as an example running under a harness that checks nothing is
# emitted.
function _announce(verbosity::Symbol, msg::AbstractString)
    verbosity in (:silent, :normal, :verbose, :full, :debug) ||
        error("`verbosity` must be `:silent`, `:normal` or `:verbose` (aliases `:full`, `:debug`); " *
              "got `:$verbosity`.")
    verbosity === :silent || @info msg
    return nothing
end

# One value, chosen by what is **wanted** against where it is **coming from** — one function rather
# than a `_getregime`/`_getsupply` family, so a new source or a new keyword is one added method
# rather than a new name, and the call site reads as what it asks
# (`_getdefaultvalue(AbstractSupply, source, verbosity)`).
#
# This first method is the fallback for a source nothing else handles. It is a fallback rather than a
# `Union` annotation because the open first argument is the whole point — a new kind of source must
# be one added method — but a bare `MethodError` on a private helper tells a caller nothing useful.
function _getdefaultvalue(::Type, source, ::Symbol)
    return error("`build_habitat` cannot take defaults from a $(nameof(typeof(source))). Pass " *
                 "`DefaultEcosystem()` for the announced defaults, or an existing `GridHabitat` to " *
                 "rebuild from one — or name every input explicitly and call `GridHabitat` itself.")
end

function _getdefaultvalue(::Type{AbstractRegime}, ::DefaultEcosystem,
                          verbosity::Symbol)
    _announce(verbosity,
              "build_habitat(DefaultEcosystem()): `regime` defaulted to a uniform 25.0 °C temperature.")
    # Stated in °C because that is what the announcement says and what a reader recognises as a
    # temperate default. `Temperature`'s canonical unit is K, so it converts affinely on the way in
    # — 25.0 °C is 298.15 K, not 25 K of anything.
    return UniformSpec(25.0°C, axis = Temperature)
end

function _getdefaultvalue(::Type{AbstractSupply}, ::DefaultEcosystem,
                          verbosity::Symbol)
    _announce(verbosity,
              "build_habitat(DefaultEcosystem()): `supply` defaulted to a uniform canonical solar rate.")
    return UniformSpec(1.0 * canonicalunit(Resource, SolarRadiation) / km^2,
                       axis = SolarRadiation)
end

# Announced even though `GridHabitat`'s own default is this same `Island()`: every value this builder
# chose has to be visible, and one silent choice would be an exception a reader must know about.
function _getdefaultvalue(::Type{AbstractTopology}, ::DefaultEcosystem,
                          verbosity::Symbol)
    _announce(verbosity,
              "build_habitat(DefaultEcosystem()): `topology` defaulted to `Island()` — the grid's edges do not wrap.")
    return Island()
end

# The synthetic toy grid is a sensible default only because `DefaultEcosystem`'s `regime` and
# `supply` are synthetic too. An explicitly passed **data** regime still needs its own `area`: no toy
# grid could suit unseen data, and the builder will not guess one.
function _getdefaultvalue(::Type{StudyArea}, ::DefaultEcosystem,
                          verbosity::Symbol)
    _announce(verbosity,
              "build_habitat(DefaultEcosystem()): `area` defaulted to a synthetic 10.0 km × 10.0 km grid of 1.0 km cells.")
    return StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                     verbosity = :silent)
end

# The same four questions asked of an existing habitat, which answers with what it actually has.
# That is what makes `build_habitat(h, supply = …)` a rebuild with one thing changed rather than a
# fresh start: everything not named is carried across, the grid included. The layers come back
# **built** and are handed to a `GridHabitat` on the habitat's own as-built area, which `materialise`
# copies rather than aliases. Nothing is resampled, because nothing moved.
function _getdefaultvalue(::Type{AbstractRegime}, source::GridHabitat,
                          verbosity::Symbol)
    _announce(verbosity,
              "build_habitat(habitat): `regime` taken from the habitat as built.")
    return source.regime
end

function _getdefaultvalue(::Type{AbstractSupply}, source::GridHabitat,
                          verbosity::Symbol)
    _announce(verbosity,
              "build_habitat(habitat): `supply` taken from the habitat as built.")
    return source.supply
end

function _getdefaultvalue(::Type{AbstractTopology}, source::GridHabitat,
                          verbosity::Symbol)
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
function _getdefaultvalue(::Type{StudyArea}, source::GridHabitat,
                          verbosity::Symbol)
    _announce(verbosity,
              "build_habitat(habitat): `area` taken from the habitat's own as-built grid.")
    return StudyArea(source, verbosity = :silent)
end
