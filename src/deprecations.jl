# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations - main `EcoSISTEM` module
#
# Every deprecated public API is collected here, in sections, and included late in
# `EcoSISTEM.jl`, after all the types it shims. Each shim warns - so downstream code gets a
# migration message rather than a silent `MethodError` - and forwards to the current API. A released
# form that cannot be forwarded without changing what it does is instead removed, keeping a method
# that raises an error saying how to rewrite the call. Mirrored by `test/test_deprecations.jl`. The
# shims whose signatures name a `RasterDataSources` type are in
# `ext/EcoSISTEMRasterDataSourcesExt/deprecations.jl`, since they cannot be defined without it; it
# holds none at present.
#
# ## Every section says which release deprecated or removed it
#
# Each section header ends with a `Deprecated in vX.Y.Z` or `Removed in vX.Y.Z` line saying which
# release the shim belongs to, and **the sections are ordered newest release first**, so the oldest
# deprecations sit at the foot of the file and dropping a release's worth is deleting from the bottom
# up to the first section carrying a newer label. `clean_Deprecations.jl` asserts every section has a
# label naming a release it knows and that the labels never get newer going down the file.
#
# Two things a whole-section delete does not reach, and both fail quietly:
#
#   - an `export` or `public` declaration for a deprecated name that sits anywhere but beside its
#     shim. A stale `export` of a deleted binding is a precompile *warning*, not an error.
#   - `test/test_deprecations.jl`, which mirrors this file and must lose the same sections.
# ===========================================================================

# ---------------------------------------------------------------------------
# A species list no longer holds where its individuals start: the placement `build_ecosystem` is
# given says so. A form that took starting abundances, or read them from a species list, could go on
# working only by carrying them on the list, where they would change what other code does, so each
# raises an error saying how to rewrite the call.
#
# Removed in v0.11.0.
# ---------------------------------------------------------------------------
function build_species(numspecies::Integer; kw...)
    return error(_speciesrewrite(numspecies, values(kw)))
end

function SpeciesList(::Int64, ::AbstractTolerance, ::AbstractVector{<:Integer},
                     ::AbstractDemand, ::AbstractVector, ::AbstractParams,
                     ::AbstractVector{Bool}; kw...)
    return error(_REMOVEDLISTFORM)
end

function SpeciesList(::Int64, ::Int64, ::AbstractVector{<:Integer},
                     ::AbstractDemand, ::AbstractVector, ::AbstractParams,
                     ::AbstractVector{Bool}; kw...)
    return error(_REMOVEDLISTFORM)
end

function SpeciesList(::Int64, ::Int64, ::AbstractVector{<:Integer},
                     ::AbstractDemand, ::AbstractVector, ::AbstractParams,
                     ::AbstractVector{Bool}, ::Vector{Float64}; kw...)
    return error(_REMOVEDLISTFORM)
end

function SpeciesList(::Int64, ::Int64, ::AbstractVector{<:Integer},
                     ::AbstractDemand, ::AbstractVector, ::AbstractTypes,
                     ::AbstractParams, ::AbstractVector{Bool}; kw...)
    return error(_REMOVEDLISTFORM)
end

function SpeciesList{TL, DM, MO, T, P}(::Vector{String}, ::TL,
                                       ::AbstractVector{<:Integer}, ::DM, ::T,
                                       ::MO, ::P,
                                       ::AbstractVector{Bool};
                                       kw...) where {TL, DM, MO, T, P}
    return error(_REMOVEDLISTFORM)
end

function SpeciesList{TL, DM, MO, T, P}(::TL, ::AbstractVector{<:Integer}, ::DM,
                                       ::T, ::MO, ::P,
                                       ::AbstractVector{Bool}) where {TL, DM,
                                                                      MO, T,
                                                                      P}
    return error(_REMOVEDLISTFORM)
end

function Ecosystem(::Function, ::SpeciesList, ::GridHabitat,
                   ::AbstractNicheFit; kw...)
    return error(_removedbuild("Ecosystem(popfun, spplist, habitat, nichefit)"))
end

function Ecosystem(::SpeciesList, ::GridHabitat, ::AbstractNicheFit; kw...)
    return error(_removedbuild("Ecosystem(spplist, habitat, nichefit)"))
end

function MPIEcosystem(::Function, ::SpeciesList, ::GridHabitat, ::Any; kw...)
    return error(_removedbuild("MPIEcosystem(popfun, spplist, habitat, nichefit)"))
end

function MPIEcosystem(::SpeciesList, ::GridHabitat, ::Any; kw...)
    return error(_removedbuild("MPIEcosystem(spplist, habitat, nichefit)"))
end

"""
    populate!(ml, spplist, habitat, nichefit, rngs)

Removed in v0.11.0: an ecosystem's individuals are placed by the placement [`build_ecosystem`](@ref)
is given, and placed again by [`repopulate!`](@ref); [`BySupply`](@ref) places by supply.
"""
function populate!(::Any, ::SpeciesList, ::AbstractHabitat, ::Any, ::Any)
    return error(_removedpopulate("populate!", "BySupply()"))
end

"""
    populate_by_tolerance!(ml, spplist, habitat, nichefit, rngs)

Removed in v0.11.0: an ecosystem's individuals are placed by the placement [`build_ecosystem`](@ref)
is given, and placed again by [`repopulate!`](@ref); [`BySuitability`](@ref) places by suitability.
"""
function populate_by_tolerance!(::Any, ::SpeciesList, ::AbstractHabitat, ::Any,
                                ::Any)
    return error(_removedpopulate("populate_by_tolerance!", "BySuitability()"))
end

"""
    emptypopulate!(ml, spplist, habitat, nichefit, rngs)

Removed in v0.11.0: an ecosystem's individuals are placed by the placement [`build_ecosystem`](@ref)
is given; [`NoAbundance`](@ref) places none.
"""
function emptypopulate!(::Any, ::SpeciesList, ::AbstractHabitat, ::Any, ::Any)
    return error(_removedpopulate("emptypopulate!", "NoAbundance()"))
end

export populate!, populate_by_tolerance!, emptypopulate!

const _REMOVEDLISTFORM = "`SpeciesList`'s `abun` and `native` arguments were removed in " *
                         "v0.11.0: a species list no longer holds where its individuals " *
                         "start. Drop the two - `SpeciesList(numspecies, tolerance, demand, " *
                         "movement, params)`, or `SpeciesList(numspecies, numtraits, demand, " *
                         "movement, params)` for a list evolved along a tree - and place the " *
                         "individuals with the placement `build_ecosystem` is given: " *
                         "`BySupply(each = row -> abun[row.key])` places these counts by " *
                         "supply. A species that was not native is now one deactivated in " *
                         "its pool."

# The error a removed way of building an ecosystem raises, naming the form it was given.
function _removedbuild(form::AbstractString)
    return "`$form` was removed in v0.11.0, with the population functions it placed by: build " *
           "the ecosystem with `build_ecosystem(placement, species, habitat; nichefit, seed, " *
           "distributed)`, which places its individuals by the placement - `DefaultEcosystem()` " *
           "($(_default(Val(:individuals)).value) per species, by supply), " *
           "`BySupply(each = row -> abun[row.key])` for counts of your own, `NoAbundance()` for " *
           "none."
end

# The error a removed population function raises, naming the placement that places as it did.
function _removedpopulate(name::AbstractString, placement::AbstractString)
    return "`$name` was removed in v0.11.0: a species list no longer holds starting abundances " *
           "to place. Place an ecosystem's individuals with the placement `build_ecosystem` is " *
           "given, or again with `repopulate!(eco, placement)`: `$placement` places as " *
           "`$name` did."
end

# Keywords only the removed `build_species(numspecies; ...)` took.
const _REMOVEDKEYWORDS = (:toleranceaxis, :demandaxis, :dispersal, :pthresh,
                          :disperse_safely, :abundance, :native, :names)

# Refuse a keyword of `build_species(DefaultEcosystem(); ...)` that only the removed
# `build_species(numspecies; ...)` took, or a requirement in its shape, with the call's rewrite.
function _refuseremovedkeywords(numspecies, tolerance, demand, kw)
    movement = get(kw, :movement, nothing)
    (isempty(intersect(keys(kw), _REMOVEDKEYWORDS)) &&
     _specshaped(tolerance) && _specshaped(demand) &&
     (isnothing(movement) || movement isa Union{MovementSpec, ChosenBy})) &&
        return nothing
    return error(_speciesrewrite(numspecies,
                                 merge((tolerance = tolerance, demand = demand),
                                       values(kw)), shorthand = true))
end

# Whether a requirement is given as a pool's build takes one: a spec, several, or one built already.
_specshaped(::Union{AbstractTraitSpec, AbstractSpeciesRequirement}) = true

_specshaped(x::Union{Tuple, NamedTuple}) = all(_specshaped, values(x))

_specshaped(::Any) = false

# The rewrite of a call to the removed `build_species(numspecies; ...)`, from the keywords it gave,
# or, for the `shorthand` `build_species(DefaultEcosystem(); ...)`, of the keywords only that took.
function _speciesrewrite(numspecies, kw::NamedTuple; shorthand::Bool = false)
    given(k) = get(kw, k, nothing)
    notes = String[]
    args = String[]
    for (k, rewrite, axis) in ((:tolerance, _tolerancerewrite, :toleranceaxis),
                               (:demand, _demandrewrite, :demandaxis))
        shorthand && _specshaped(given(k)) && continue
        push!(args, "$k = " * rewrite(given(k), given(axis), notes))
    end
    if !isnothing(given(:dispersal)) || given(:movement) isa AbstractMovement
        mode = something(given(:movement), _default(Val(:movement)).value.mode)
        distance = something(given(:dispersal), _default(Val(:dispersal)).value)
        push!(args,
              "movement = MovementSpec($(_modephrase(mode)), meandist = " *
              "$(_valuephrase(distance, :meandist, notes)))")
    end
    isnothing(given(:disperse_safely)) ||
        push!(notes,
              "`disperse_safely` is the dispersal mode's own: `InsectDispersal(disperse_safely = " *
              "false)`.")
    for k in (:birth, :death, :longevity, :survival, :seed, :provenance)
        isnothing(given(k)) ||
            push!(args, "$k = $(_valuephrase(given(k), k, notes))")
    end
    isnothing(given(:native)) ||
        push!(notes,
              "`native` has no replacement: a species to leave out is deactivated in its pool, " *
              "`deactivate_species!(pool, name)`.")
    pool = !isnothing(given(:names)) ? "SpeciesPool($(repr(given(:names))))" :
           shorthand ? "DefaultEcosystem(), numspecies = $numspecies" :
           "SpeciesPool($numspecies)"
    abundance = given(:abundance)
    placement = isnothing(abundance) ? "DefaultEcosystem()" :
                abundance isa AbstractVector ?
                "BySupply(each = row -> abundance[row.key])" :
                "BySupply(total = $abundance)"
    abundance isa AbstractVector &&
        push!(notes, "`abundance` is the vector of counts this call gave.")
    lead = shorthand ?
           "`build_species(DefaultEcosystem(); ...)` builds species from a `SpeciesPool` since " *
           "v0.11.0, taking its specs, and where their individuals start is the placement " *
           "`build_ecosystem` is given." :
           "`build_species($numspecies; ...)` was removed in v0.11.0: species are built from a " *
           "`SpeciesPool`, and where their individuals start is the placement `build_ecosystem` " *
           "is given."
    separator = startswith(pool, "DefaultEcosystem") ? ", " : "; "
    return lead * " Rewritten, this call is\n\n    species = build_species(" *
           pool *
           (isempty(args) ? "" : separator * join(args, ", ")) * ")\n    " *
           "ecosystem = build_ecosystem($placement, species, habitat)\n" *
           (isempty(notes) ? "" : "\n" * join(notes, "\n"))
end

# The removed form's tolerance as a pool's spec: a spec or a built one as it is, a `(mean, width)`
# pair on its axis, or a tuple of pairs, one per regime, on theirs.
function _tolerancerewrite(value, axis, notes)
    _specshaped(value) && return "the same tolerance"
    value isa Union{Tuple, NamedTuple} && first(value) isa Tuple &&
        return _membersrewrite(value, axis, _tolerancespec, notes)
    return _tolerancespec(axis, value, notes)
end

# The removed form's demand as a pool's spec: a spec or a built one as it is, a value on its axis, or
# a tuple of values, one per supply, on theirs.
function _demandrewrite(value, axis, notes)
    _specshaped(value) && return "the same demand"
    value isa Union{Tuple, NamedTuple} &&
        return _membersrewrite(value, axis, _demandspec, notes)
    return _demandspec(axis, value, notes)
end

# Several members' specs, each on its own axis, keeping their names.
function _membersrewrite(value, axis, spec, notes)
    axes = axis isa Union{Tuple, NamedTuple} ? values(axis) :
           ntuple(_ -> axis, length(value))
    parts = [spec(a, v, notes) for (a, v) in zip(axes, values(value))]
    value isa NamedTuple || return "(" * join(parts, ", ") * ")"
    return "(" * join(("$k = $p" for (k, p) in zip(keys(value), parts)), ", ") *
           ")"
end

function _tolerancespec(axis, value, notes)
    (value isa Tuple && length(value) == 2) ||
        return "ToleranceSpec($(_axisname(axis)), μ = ..., σ = ...)"
    return "ToleranceSpec($(_axisname(axis)), μ = " *
           "$(_valuephrase(value[1], :μ, notes)), σ = $(_valuephrase(value[2], :σ, notes)))"
end

function _demandspec(axis, value, notes)
    return "DemandSpec($(_axisname(axis)), $(_valuephrase(value, :demand, notes)))"
end

# A dispersal mode as written in the rewrite: `InsectDispersal()` as it is usually written, rather
# than as it prints.
function _modephrase(mode)
    return mode == InsectDispersal() ? "InsectDispersal()" : sprint(show, mode)
end

# An axis as written in the rewrite.
_axisname(axis::Type) = string(nameof(axis))

_axisname(::Nothing) = "axis"

# A value as written in the rewrite: as it is, or for one per species a column of the pool's table.
function _valuephrase(value::AbstractVector, name::Symbol, notes)
    push!(notes,
          "One value per species becomes a column of the pool's table, `SpeciesPool(DataFrame(" *
          "species = ..., $name = ...))`, read by `TraitColumn(:$name)`.")
    return "TraitColumn(:$name)"
end

_valuephrase(value, ::Symbol, notes) = sprint(show, value)

# ---------------------------------------------------------------------------
# Placing again by suitability is `repopulate!` given the suitability placement, which places over the
# active cells only, places none of a species no cell suits rather than spreading it evenly, and takes
# its numbers from the placement rather than re-splitting the species' totals at random.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
@deprecate(repopulate_by_tolerance!(eco::AbstractEcosystem),
           repopulate!(eco, BySuitability()), false)

export repopulate_by_tolerance!

# ---------------------------------------------------------------------------
# The kernels take their distances by keyword, each named for what it is, since the positional
# distance is a Gaussian kernel's mean but neither the mean nor the scale of a long-tailed one. Each
# shim must build the kernel its positional form denotes to the bit, or every lookup and canonical
# reference built from one moves: the long-tailed distance is the scale `(2 * dist) / sqrt(pi)`.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
@deprecate(GaussianKernel(dist::Unitful.Length{Float64}, thresh::Float64),
           GaussianKernel(meandist = dist, thresh = thresh), false)
@deprecate(LongTailKernel(dist::Unitful.Length{Float64}, shape::Float64,
                          thresh::Float64),
           LongTailKernel(scale = (2 * uconvert(km, dist)) / sqrt(pi),
                          shape = shape, thresh = thresh), false)

# ---------------------------------------------------------------------------
# A species list's movement is one `SpeciesDispersal` per species, a dispersal mode with the species'
# kernel, where it was one object for the whole list. Each list-wide movement type is a function
# returning that vector, in the mode it named, each kernel's type its family and each species' flag
# its dead-cell fate.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
function BirthOnlyMovement(kernels::AbstractVector{<:AbstractKernel},
                           disperse_safely::AbstractVector{Bool} = fill(true,
                                                                        length(kernels)))
    Base.depwarn("`BirthOnlyMovement(kernels, disperse_safely)` is deprecated: give each " *
                 "species its movement, as `SeedDispersal(family; disperse_safely).(kernels)`.",
                 :BirthOnlyMovement)
    return _listmovement(SeedDispersal, kernels, disperse_safely)
end

function AlwaysMovement(kernels::AbstractVector{<:AbstractKernel},
                        disperse_safely::AbstractVector{Bool} = fill(true,
                                                                     length(kernels)))
    Base.depwarn("`AlwaysMovement(kernels, disperse_safely)` is deprecated: give each " *
                 "species its movement, as `AlwaysDispersal(family; disperse_safely).(kernels)`.",
                 :AlwaysMovement)
    return _listmovement(AlwaysDispersal, kernels, disperse_safely)
end

function NoMovement(kernels::AbstractVector{<:AbstractKernel},
                    ::AbstractVector{Bool} = fill(true, length(kernels)))
    Base.depwarn("`NoMovement(kernels)` is deprecated: give each species its movement, " *
                 "`NoDispersal()()`.", :NoMovement)
    return fill(NoDispersal()(), length(kernels))
end

# A list-wide movement's species, each in `mode` with its own kernel and dead-cell fate.
function _listmovement(mode, kernels, disperse_safely)
    length(disperse_safely) == length(kernels) ||
        error("`disperse_safely` has $(length(disperse_safely)) entries but there are " *
              "$(length(kernels)) dispersal kernels: it is a per-species setting, so it needs one " *
              "entry per species.")
    return [SpeciesDispersal(mode(typeof(k), disperse_safely = s), k)
            for (k, s) in zip(kernels, disperse_safely)]
end

export BirthOnlyMovement, AlwaysMovement, NoMovement

# ---------------------------------------------------------------------------
# `InputRecord` is `ProvenanceRecord`: a record describes a file however it enters a run, an
# output of the run that made it included, so its name no longer calls it an input. The binding
# stays `public`, as the old name was.
#
# Deprecated in v0.10.0.
# ---------------------------------------------------------------------------
Base.@deprecate_binding InputRecord ProvenanceRecord false
public InputRecord
