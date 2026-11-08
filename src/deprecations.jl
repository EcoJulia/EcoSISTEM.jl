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

function AddSpecies(; kw...)
    return error(_addspeciesrewrite(values(kw)))
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

# The error the removed `AddSpecies(; tolerance, ..., abundance, name)` raises, with the call that does
# its job built from the keywords it was given.
function _addspeciesrewrite(kw::NamedTuple)
    name = repr(String(something(get(kw, :name, nothing), "arrival")))
    total = something(get(kw, :abundance, nothing), 1000)
    return "`AddSpecies(; tolerance, demand, dispersal, birth, death, abundance, name)` was removed " *
           "in v0.11.0: a run's species are declared in its pool before it is built, so declare " *
           "the arrival there - its traits in the pool's table, or `copy_species!` for a copy of " *
           "one - and leave it out of the run with `deactivate_species!`. Then " *
           "`AddSpecies($name, BySupply(total = $total))` brings it in, placed by supply over the " *
           "intervention's region."
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
# A run's species are all declared in its pool before it is built, so a recording has a row for each
# from the start, an arrival's included, and no room needs leaving for one.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
function _maxspeciesdeprecated(caller::Symbol)
    return Base.depwarn("`$caller`'s `maxspecies` is deprecated, and ignored: every species a run " *
                        "may hold is in its pool from the start, so the storage has a row for " *
                        "each already.", caller)
end

# ---------------------------------------------------------------------------
# Individuals of a species arriving are placed by the placement its `AddSpecies` is given, and are
# added to a species already in the run by `AddAbundance`. Adding a number of the last species at
# random cells with any supply keeps its body.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
function repopulate!(eco::Ecosystem, abun::Int64)
    Base.depwarn("`repopulate!(eco, abun)` is deprecated: an arriving species' individuals are " *
                 "placed by its `AddSpecies` placement, and individuals are added to a species " *
                 "already in the run by `AddAbundance`.", :repopulate!)
    grid, activity = _gridactivity(eco.habitat)
    # Set up copy of supply
    b = reshape(parent(copy(_getsupply(eco.habitat.supply))), length(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    # Draw locations from the last species' own RNG stream
    pos = sample(getrng(eco, lastindex(eco.rngs)), grid[b .> (0 * units)], abun)
    # Add individual to this location
    map(pos) do p
        return eco.abundances.matrix[end, p] += 1
    end
end

# The linear indices of every cell of `habitat`'s grid, and a flattened copy of its active mask, which
# `repopulate!(eco, abun)` zeroes the supply of inactive cells by, so that none is placed outside it.
function _gridactivity(habitat::AbstractHabitat)
    dim = getgridshape(habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # `reshape` on a `DimArray` isn't a clean zero-copy operation the way a plain `Array`
    # reshape is - unwrap via `parent()` first.
    activity = reshape(copy(parent(habitat.active)), len)
    return grid, activity
end

# ---------------------------------------------------------------------------
# A table of dated additions is a placement: `AddAbundance(FromTable(...))` serves its rows as they
# fall due, naming species and places as every placement does. The table operation, which also takes a
# species by its position and a cell by its linear index, forwards there, its rows converted when the
# intervention resolves them; its fractional counts of one species, cell and time are summed and
# rounded, ties to even. `build_abundance_table` writes that table from names and `(y, x)` positions,
# which `FromTable` reads as they are given (`CartesianIndex(y, x)` cells), and keeps its body.
#
# Deprecated in v0.11.0.
# ---------------------------------------------------------------------------
function AddAbundanceTable(table)
    Base.depwarn("`AddAbundanceTable(table)` is deprecated: give `AddAbundance(FromTable(table))`, " *
                 "naming each row's species by name and its place by a `cell` grid key or a " *
                 "location, in an `Intervention` of its own with `EveryStep()`.",
                 :AddAbundanceTable)
    return AddAbundance(_TableRows(table))
end

# The rows of an `AddAbundanceTable` as a dated placement: a species given by its position in the
# species list as its name, a cell given by its linear index as its grid key, and the counts of one
# species, cell and time summed and rounded to whole individuals, ties to even.
struct _TableRows{T}
    table::T
end

_isdated(::_TableRows) = true

_placementphrase(::_TableRows) = "the rows of an `AddAbundanceTable`"

function (p::_TableRows)(species::SpeciesList, habitat::GridHabitat,
                         ::PlacementContext)
    shape = Tuple(getgridshape(habitat))
    grid = CartesianIndices(shape)
    totals = Dict{Tuple{String, CartesianIndex{2}, Any}, Float64}()
    for (n, row) in enumerate(Tables.rows(p.table))
        columns = Tables.columnnames(row)
        for name in (:species, :cell, :time)
            name in columns ||
                error("row $n of the table has no `$name` column; an `AddAbundanceTable` needs " *
                      "`species`, `cell` and `time`, and optionally `count`.")
        end
        name = _tablespeciesname(Tables.getcolumn(row, :species), species.names,
                                 n)
        cell = grid[_tablecell(Tables.getcolumn(row, :cell), shape, n)]
        count = :count in columns ? Tables.getcolumn(row, :count) : 1
        time = Tables.getcolumn(row, :time)
        _tableelapsedtype(time, n)
        key = (name, cell, _normaltime(time))
        totals[key] = get(totals, key, 0.0) + _tablecount(count, n)
    end
    rows = [(k, round(Int, v)) for (k, v) in totals if round(Int, v) > 0]
    return (species = String[k[1] for (k, _) in rows],
            cell = CartesianIndex{2}[k[2] for (k, _) in rows],
            count = Int[n for (_, n) in rows],
            time = Any[k[3] for (k, _) in rows])
end

# A table row's species by name: its own, or the name of the species at its position in the list.
function _tablespeciesname(species::Union{AbstractString, Symbol}, ::Any, ::Int)
    return String(species)
end

function _tablespeciesname(species::Integer, names::AbstractVector, row::Int)
    1 <= species <= length(names) ||
        error("row $row of the abundance table names species $species, but the species list " *
              "has $(length(names)).")
    return names[species]
end

function _tablespeciesname(species, ::Any, row::Int)
    return error("row $row of the abundance table gives its species as `$(repr(species))`: give " *
                 "a name or an index.")
end

function build_abundance_table(eco::AbstractEcosystem, table)
    Base.depwarn("`build_abundance_table(eco, table)` is deprecated: give the table to `FromTable` " *
                 "as it is, naming each row's species by name and its place by a `cell` - " *
                 "`CartesianIndex(y, x)` - or by `latitude` and `longitude`, and add it with " *
                 "`AddAbundance(FromTable(table))`.", :build_abundance_table)
    lookup = _emptynamelookup()
    shape = Base.size(parent(eco.habitat.active))
    rows = Tables.rows(table)
    first = iterate(rows)
    isnothing(first) &&
        error("the table has no rows to build an abundance table from.")
    row, state = first
    return _collectabundances(rows, state,
                              _positionedrow(row, 1, eco, lookup, shape), eco,
                              lookup, shape)
end

# An empty lookup from species names to their indices, for `_speciesnamelookup!` to fill.
_emptynamelookup() = Ref{Union{Nothing, Dict{String, Int}}}(nothing)

# Each species' name and its index in the species list, built the first time a row names a species
# and kept in `lookup` for the rest of the table, so a table of indices never builds it.
function _speciesnamelookup!(lookup::Base.RefValue{Union{Nothing,
                                                         Dict{String, Int}}},
                             eco::AbstractEcosystem)
    if isnothing(lookup[])
        lookup[] = Dict(name => i for (i, name) in enumerate(eco.spplist.names))
    end
    return lookup[]
end

# A table row's species as its index in the species list, refused where there is no such species.
function _tablespecies(species::Integer, eco::AbstractEcosystem, _, row::Int)
    count = length(eco.spplist.names)
    1 <= species <= count ||
        error("row $row of the abundance table names species $species, but the species list " *
              "has $count.")
    return Int(species)
end

function _tablespecies(species::Union{AbstractString, Symbol},
                       eco::AbstractEcosystem, lookup, row::Int)
    index = get(_speciesnamelookup!(lookup, eco), String(species), nothing)
    isnothing(index) &&
        error("row $row of the abundance table names species `$species`, which is not in the " *
              "species list.")
    return index
end

function _tablespecies(species, ::AbstractEcosystem, _, row::Int)
    return error("row $row of the abundance table gives its species as `$(repr(species))`: give " *
                 "a name or an index.")
end

# A table row's cell as a linear index into a grid of `shape`, refused off the grid.
function _tablecell(cell::Integer, shape::Tuple{Int, Int}, row::Int)
    1 <= cell <= prod(shape) ||
        error("row $row of the abundance table names cell $cell, but the grid has " *
              "$(prod(shape)) cells.")
    return Int(cell)
end

function _tablecell(cell::CartesianIndex{2}, shape::Tuple{Int, Int}, row::Int)
    checkbounds(Bool, CartesianIndices(shape), cell) ||
        error("row $row of the abundance table names the cell at $(Tuple(cell)), which is off " *
              "the $(shape[1]) by $(shape[2]) grid.")
    return LinearIndices(shape)[cell]
end

function _tablecell(cell, ::Tuple{Int, Int}, row::Int)
    return error("row $row of the abundance table gives its cell as `$(repr(cell))`: give a " *
                 "linear index or a `CartesianIndex` `(y, x)`.")
end

# A table row's count, refused unless it is a finite number, zero or more.
function _tablecount(count::Real, row::Int)
    isfinite(count) && count >= 0 ||
        error("row $row of the abundance table adds $count individuals: a count must be a finite " *
              "number, zero or more.")
    return Float64(count)
end

function _tablecount(count, row::Int)
    return error("row $row of the abundance table gives its count as `$(repr(count))`: give a " *
                 "number.")
end

# One row of a table given to `build_abundance_table`: its species as an index, its `(y, x)` as a
# linear cell index, its count as given and its time, an elapsed one in seconds.
function _positionedrow(row, n::Int, eco::AbstractEcosystem, lookup,
                        shape::Tuple{Int, Int})
    columns = Tables.columnnames(row)
    for name in (:species, :y, :x, :time)
        name in columns ||
            error("row $n of the table has no `$name` column; `build_abundance_table` needs " *
                  "`species`, `y`, `x` and `time`, and optionally `count`.")
    end
    y = Tables.getcolumn(row, :y)
    x = Tables.getcolumn(row, :x)
    y isa Integer && x isa Integer ||
        error("row $n of the table gives its position as ($(repr(y)), $(repr(x))): `y` and `x` " *
              "are integer indices into the grid.")
    count = :count in columns ? Tables.getcolumn(row, :count) : 1
    _tablecount(count, n)
    time = Tables.getcolumn(row, :time)
    _tableelapsedtype(time, n)
    return (species = _tablespecies(Tables.getcolumn(row, :species), eco,
                                    lookup, n),
            cell = _tablecell(CartesianIndex(y, x), shape, n), count = count,
            time = _normaltime(time))
end

# A row's time as `build_abundance_table` keeps it: a date as it is, an elapsed time in seconds so
# that the same instant in two units is one key.
_normaltime(time::Unitful.Time) = uconvert(s, float(time))

_normaltime(date::Dates.TimeType) = date

# Refuse a time that is neither an elapsed time nor a date.
_tableelapsedtype(::Union{Unitful.Time, Dates.TimeType}, ::Int) = nothing

function _tableelapsedtype(time, row::Int)
    return error("row $row of the table gives its time as `$(repr(time))`: give an elapsed time " *
                 "or a date.")
end

# Sum the rows of a table, the first already read, by time, species and cell, and return them in that
# order as the columns of an abundance table.
function _collectabundances(rows, state, record::NamedTuple,
                            eco::AbstractEcosystem, lookup,
                            shape::Tuple{Int, Int})
    T = typeof(record.time)
    totals = Dict{Tuple{T, Int, Int}, Float64}()
    whole = true
    n = 1
    while true
        key = (record.time, record.species, record.cell)
        totals[key] = get(totals, key, 0.0) + record.count
        whole &= record.count isa Integer
        next = iterate(rows, state)
        isnothing(next) && break
        row, state = next
        n += 1
        record = _positionedrow(row, n, eco, lookup, shape)
        record.time isa T ||
            error("row $n of the table gives its time as a $(typeof(record.time)), where the " *
                  "rows before it gave a $T.")
    end
    ordered = sort!(collect(keys(totals)))
    counts = [totals[key] for key in ordered]
    return (species = [key[2] for key in ordered],
            cell = [key[3] for key in ordered],
            count = whole ? round.(Int, counts) : counts,
            time = T[key[1] for key in ordered])
end

export AddAbundanceTable, build_abundance_table

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
