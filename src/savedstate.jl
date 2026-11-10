# SPDX-License-Identifier: LGPL-3.0-or-later
#
# A run's state saved to a file, and an ecosystem started from one. The file holds plain data only -
# arrays, numbers and strings, never a package's type - so it is read without the types that wrote it.

using JLD2: jldopen
using Random

# What a state file says it is, and the version of its layout.
const _STATEFORMAT = "EcoSISTEM state"

const _STATEVERSION = 1

"""
    write_state(path::AbstractString, eco::AbstractEcosystem)

Save the state `eco` has reached to `path`, a JLD2 file a later run can start from with
`build_ecosystem(species, habitat; from = path)`, with the run's account beside it, as `path` with
`.provenance.toml` in place of `.jld2`; and return `path`.

The state is what a run carries forward that its inputs do not give again: every species' abundance
in every cell, which species are in the run, each species' random stream, the clock, the run's seed,
epoch and calendar, which cells are active, and every layer's current values. It names its grid, its
species and its layers, so a run starting from it can be checked against them. A layer's change rule
is code rather than data, so it is not saved; its description is, and a run starting from the state
warns where its own differs.

Under MPI this is a collective call: every rank sends its own species' abundances and streams, and the
first rank writes the file.

# Arguments

  - `path`: the file to write, ending `.jld2`.
  - `eco`: the ecosystem whose state is saved.
"""
function write_state(path::AbstractString, eco::AbstractEcosystem)
    endswith(path, ".jld2") ||
        error("a state is saved as a JLD2 file, so its `path` ends `.jld2`; got `$path`.")
    abundances = _wholeabundances(eco)
    streams = _wholestreams(eco)
    _isrootrank(eco) || return path
    mkpath(dirname(abspath(path)))
    _writestate(path, abundances, streams, eco.spplist, eco.habitat,
                eco.elapsed, eco.seed, eco.epoch, eco.calendar)
    write_provenance(_stateaccount(path), provenance(eco))
    return path
end

# Where a state's account is written: beside it, `.provenance.toml` in place of `.jld2`.
_stateaccount(path::AbstractString) = chop(path, tail = 5) * ".provenance.toml"

# The whole abundance matrix and every species' stream, serially as they are; under MPI the
# extension's methods gather them on the first rank.
_wholeabundances(eco::Ecosystem) = eco.abundances.matrix

_wholestreams(eco::Ecosystem) = _streamwords(eco.rngs)

# Each stream's words, one column a stream.
function _streamwords(rngs::AbstractVector{Random.Xoshiro})
    return UInt64[getfield(rng, f)
                  for f in fieldnames(Random.Xoshiro), rng in rngs]
end

# Set each stream whose saved column `columns` names (zero for none) to that column's words.
function _setstreams!(rngs::AbstractVector{Random.Xoshiro},
                      words::AbstractMatrix{UInt64}, columns)
    for (rng, column) in zip(rngs, columns)
        iszero(column) && continue
        for (i, f) in enumerate(fieldnames(Random.Xoshiro))
            setfield!(rng, f, words[i, column])
        end
    end
    return rngs
end

# Write a state's plain data to `path`: the abundances (species by cells), the streams' words (one
# column a species), the species, the clock and the run's dates, the grid and its active cells, and
# every layer's values, unit and change, by role and name.
function _writestate(path::AbstractString, abundances::AbstractMatrix,
                     streams::AbstractMatrix{UInt64}, spplist::SpeciesList,
                     habitat::GridHabitat, elapsed::Unitful.Time, seed::Integer,
                     epoch, calendar::AbstractRunCalendar)
    yx = _gridyx(habitat)
    jldopen(path, "w") do f
        f["format"] = _STATEFORMAT
        f["version"] = _STATEVERSION
        f["abundances"] = Matrix{Int64}(abundances)
        f["species/names"] = Vector{String}(spplist.names)
        f["species/keys"] = Vector{Int}(spplist.traits.key)
        f["species/active"] = Vector{Bool}(spplist.traits.active)
        f["streams/words"] = Matrix{UInt64}(streams)
        f["streams/fields"] = String[string(n)
                                     for n in fieldnames(Random.Xoshiro)]
        f["seed"] = UInt64(seed)
        f["elapsed"] = Float64(ustrip(s, elapsed))
        f["epoch"] = isnothing(epoch) ? "" : string(Dates.DateTime(epoch))
        f["calendar"] = string(nameof(typeof(calendar)))
        f["grid/shape"] = collect(size(parent(habitat.active)))
        f["grid/crs"] = _crsidentity(_gridcrs(habitat))
        for (axis, d) in (("y", yx[1]), ("x", yx[2]))
            coords = collect(DimensionalData.lookup(d))
            f["grid/$axis"] = Float64.(ustrip.(coords))
            f["grid/$(axis)unit"] = _unitname(eltype(coords))
        end
        f["grid/active"] = Matrix{Bool}(parent(habitat.active))
        for (role, layers) in (("regime", habitat.regime),
            ("supply", habitat.supply))
            names = String[string(n) for n in keys(layers)]
            f["$role/names"] = names
            for (name, layer) in zip(names, values(layers))
                held = parent(layer.matrix)
                f["$role/$name/values"] = _plainvalues(held)
                f["$role/$name/unit"] = _unitname(eltype(held))
                f["$role/$name/change"] = _changeidentity(layer.change)
                f["$role/$name/changephrase"] = _changephrase(layer.change)
            end
        end
    end
    return path
end

# A grid's CRS in full, to tell it from any other, empty for a synthetic grid with none.
_crsidentity(::Nothing) = ""

_crsidentity(crs) = sprint(show, crs)

# A unit's name, empty for plain numbers; and an array's values with their unit stripped.
_unitname(::Type{<:Unitful.Quantity{T, D, U}}) where {T, D, U} = string(U())

_unitname(::Type) = ""

_plainvalues(A::AbstractArray{<:Unitful.Quantity}) = Array(ustrip.(A))

_plainvalues(A::AbstractArray) = Array(A)

# A state file's contents as plain data, refusing a file that is not one, with how to carry an older
# file across.
function _readstate(path::AbstractString)
    isfile(path) || error("there is no saved state at `$path`.")
    return jldopen(path, "r") do f
        (haskey(f, "format") && f["format"] == _STATEFORMAT) ||
            error(_notastate(path))
        f["version"] == _STATEVERSION ||
            error("`$path` is a state of layout $(f["version"]); this version of EcoSISTEM " *
                  "reads layout $_STATEVERSION.")
        function layers(role)
            return Dict(name => (values = f["$role/$name/values"],
                         unit = f["$role/$name/unit"],
                         change = f["$role/$name/change"],
                         changephrase = f["$role/$name/changephrase"])
                        for name in f["$role/names"])
        end
        epoch = f["epoch"]
        return (abundances = f["abundances"], names = f["species/names"],
                keys = f["species/keys"], active = f["species/active"],
                streams = f["streams/words"], fields = f["streams/fields"],
                seed = f["seed"], elapsed = f["elapsed"] * s,
                epoch = isempty(epoch) ? nothing : Dates.DateTime(epoch),
                calendar = f["calendar"], shape = f["grid/shape"],
                crs = f["grid/crs"], y = f["grid/y"], yunit = f["grid/yunit"],
                x = f["grid/x"], xunit = f["grid/xunit"],
                cells = f["grid/active"], regime = layers("regime"),
                supply = layers("supply"))
    end
end

# What a file not a state is, and how one saved before states were becomes one.
function _notastate(path::AbstractString)
    return "`$path` is not a saved state. Files `SaveAbundance` and a `CachedEcosystem` wrote before " *
           "v0.11.0 hold abundances alone; to start from one, build its ecosystem with " *
           "`NoAbundance()`, copy the saved matrix into `eco.abundances.matrix` and save it with " *
           "`write_state`, which gives it fresh streams and a clock at zero."
end

# The file a state is read from, and the record of it the run keeps: a path, or the record of one.
function _statefile(path::AbstractString)
    return (path = String(path),
            record = ProvenanceRecord(role = :state, source = "saved state",
                                      path = basename(path)))
end

function _statefile(record::ProvenanceRecord)
    _checktabledata(record, "`from`")
    return (path = record.location, record = record)
end

# An ecosystem of `species` on `habitat` starting from the state saved at `from`, a path or a record
# of one: checked against the state's grid, species, epoch and calendar, then given its abundances,
# which species are in the run, which cells are active and every layer's values, and on the state's
# own seed each saved species' stream; its clock at `elapsed`, the state's own unless given, the
# series following it. Under MPI the first rank reads the file and every rank takes its own rows.
function _startfromstate(species::SpeciesList, habitat::GridHabitat, from;
                         nichefit = nothing, seed = nothing,
                         distributed = :auto,
                         epoch::Union{Nothing, Dates.TimeType} = nothing,
                         calendar::Union{Nothing, AbstractRunCalendar} = nothing,
                         elapsed::Union{Nothing, Unitful.Time} = nothing,
                         changing = (),
                         provenance = ProvenanceRecord[],
                         verbosity::Symbol = :normal)
    file = _statefile(from)
    changed = _changingnames(changing, habitat)
    saved = Ref{Any}(nothing)
    state = _madeonroot() do
        read = _readstate(file.path)
        saved[] = read.abundances
        return Base.structdiff(read, NamedTuple{(:abundances,)})
    end
    _checkstategrid(state, habitat, file.path)
    columns = _statecolumns(state, species, file.path)
    exact = isnothing(seed) || seed % UInt64 == state.seed
    exact && _checkstreamfields(state, file.path)
    eco = build_ecosystem(NoAbundance(), species, habitat, nichefit = nichefit,
                          seed = exact ? state.seed : seed,
                          distributed = distributed,
                          epoch = _stateepoch(state, epoch, file.path),
                          calendar = _statecalendar(state, calendar,
                                                    file.path),
                          provenance = vcat(_recordvector(provenance),
                                            file.record),
                          verbosity = verbosity)
    _restoreground!(species, habitat, state, columns, file.path)
    eco.elapsed = isnothing(elapsed) ? state.elapsed :
                  float(uconvert(s, elapsed))
    _primeseries!(habitat.regime, eco.elapsed)
    _primeseries!(habitat.supply, eco.elapsed)
    exact && _setstreams!(eco.rngs, state.streams, columns)
    _placesaved!(eco, saved[], columns)
    eco.placement = nothing
    eco.savedchanges = _savedchanges(state, file.path, changed)
    invalidatecaches!(eco)
    return eco
end

# Give the saved species whether each is in the run, and the habitat its active cells and every layer's
# values, as the state has them.
function _restoreground!(species::SpeciesList, habitat::GridHabitat, state,
                         columns, path::AbstractString)
    for (row, column) in enumerate(columns)
        iszero(column) || (species.traits.active[row] = state.active[column])
    end
    parent(habitat.active) .= state.cells
    _restorelayers!(habitat.regime, state.regime, "regime", path)
    _restorelayers!(habitat.supply, state.supply, "supply", path)
    return habitat
end

# Refuse streams saved from a generator with other fields than this Julia's.
function _checkstreamfields(state, path::AbstractString)
    Tuple(state.fields) == Tuple(string.(fieldnames(Random.Xoshiro))) &&
        return nothing
    return error("`$path` saved its streams from a generator with the fields " *
                 "$(join(state.fields, ", ")), and this Julia's has others; give `seed` for " *
                 "fresh streams.")
end

# Refuse a habitat on another grid than the state's: its shape, CRS and cell coordinates.
function _checkstategrid(state, habitat::GridHabitat, path::AbstractString)
    yx = _gridyx(habitat)
    differs = String[]
    collect(size(parent(habitat.active))) == state.shape ||
        push!(differs,
              "shape $(join(state.shape, " by ")) against " *
              join(size(parent(habitat.active)), " by "))
    _crsidentity(_gridcrs(habitat)) == state.crs ||
        push!(differs,
              "another CRS, this habitat's $(_crsname(_gridcrs(habitat)))")
    for (axis, d, saved, unitname) in (("Y", yx[1], state.y, state.yunit),
                                       ("X", yx[2], state.x, state.xunit))
        coords = collect(DimensionalData.lookup(d))
        (_unitname(eltype(coords)) == unitname &&
         length(coords) == length(saved) &&
         isapprox(Float64.(ustrip.(coords)), saved, rtol = 1.0e-9)) ||
            push!(differs, "its $axis coordinates")
    end
    isempty(differs) && return nothing
    return error("`$path` was saved on another grid than this habitat's: " *
                 join(differs, "; ") *
                 ". A run starts from a state on the grid it was saved on.")
end

# Each of the run's species' column in the state, matched by key and name, zero for a species the
# state does not hold; refusing a species of the state the run does not hold.
function _statecolumns(state, species::SpeciesList, path::AbstractString)
    saved = Dict((k, n) => c
                 for (c, (k, n)) in enumerate(zip(state.keys, state.names)))
    columns = [get(saved, (k, n), 0)
               for (k, n) in zip(species.traits.key, species.names)]
    held = Set(c for c in columns if !iszero(c))
    absent = [c for c in eachindex(state.names) if !(c in held)]
    isempty(absent) && return columns
    return error("`$path` holds species this run does not: " *
                 _namelist(["`$(state.names[c])` (key $(state.keys[c]))"
                            for c in absent],
                           :normal) *
                 ". A run starting from a state holds every species it saved, under the same " *
                 "name and key in its pool, and may hold more.")
end

# The run's epoch and calendar: the state's, refusing others given.
function _stateepoch(state, epoch, path::AbstractString)
    isnothing(epoch) && return state.epoch
    given = Dates.DateTime(epoch)
    given == state.epoch && return state.epoch
    return error("`$path` was saved by a run whose epoch was " *
                 "$(something(state.epoch, "none")), and this one was given $given; a run " *
                 "starting from a state keeps its dates.")
end

function _statecalendar(state, calendar, path::AbstractString)
    saved = _CALENDARS[state.calendar]
    isnothing(calendar) && return saved
    calendar == saved && return saved
    return error("`$path` was saved by a run on `$(state.calendar)()`, and this one was given " *
                 "`$(nameof(typeof(calendar)))()`; a run starting from a state keeps its dates.")
end

# The run calendars by name, as a state names them.
const _CALENDARS = Dict("ExactDates" => ExactDates(),
                        "MeanMonths" => MeanMonths())

# Write each saved layer's values into the habitat's layer of the same role and name, refusing a
# habitat whose layers are not the state's, in name, shape or unit.
function _restorelayers!(layers, saved::AbstractDict, role::AbstractString,
                         path::AbstractString)
    names = Set(string(n) for n in keys(layers))
    names == Set(keys(saved)) ||
        error("`$path` was saved with the $role layers $(join(sort(collect(keys(saved))), ", ")) " *
              "and this habitat has $(join(sort(collect(names)), ", ")); a run starts from a " *
              "state on the habitat it was saved on.")
    for (name, layer) in zip(keys(layers), values(layers))
        kept = saved[string(name)]
        held = parent(layer.matrix)
        (size(held) == size(kept.values) &&
         _unitname(eltype(held)) == kept.unit) ||
            error("`$path` saved its $role layer `$name` in $(kept.unit) by " *
                  "$(join(size(kept.values), " by ")), and this habitat's is in " *
                  "$(_unitname(eltype(held))) by $(join(size(held), " by ")).")
        held .= _withunit(kept.values, eltype(held))
    end
    return layers
end

# Saved values given back their unit.
function _withunit(values, ::Type{<:Unitful.Quantity{T, D, U}}) where {T, D, U}
    return values .* U()
end

_withunit(values, ::Type) = values

# The saved state's record of each layer's change rule, kept for `_resumerun!` to compare with the
# rules in force after the run's first step, together with the layer names given in `changing`.
function _savedchanges(state, path::AbstractString,
                       changing::Vector{Symbol})
    recorded(saved) = Dict(name => _recordedchange(was)
                           for (name, was) in saved)
    return _SavedChanges(path, recorded(state.regime), recorded(state.supply),
                         changing)
end

# The `changing` keyword as a list of layer names without repeats, each of which must name a regime
# or a supply layer of the habitat.
function _changingnames(name::Symbol, habitat::GridHabitat)
    return _changingnames((name,), habitat)
end

function _changingnames(names::Union{Tuple, AbstractVector},
                        habitat::GridHabitat)
    layers = vcat([(n, "regime") for n in keys(habitat.regime)],
                  [(n, "supply") for n in keys(habitat.supply)])
    given = Symbol[]
    for name in names
        name isa Symbol ||
            error("`changing` names layers, as symbols such as `:Temperature`; got $(repr(name)).")
        any(l -> first(l) === name, layers) ||
            error("`changing` names `$name`, which is no layer of the habitat; its layers are " *
                  join(("`$n` ($role)" for (n, role) in layers), ", ") * ".")
        name in given || push!(given, name)
    end
    return given
end

function _changingnames(names, ::GridHabitat)
    return error("`changing` names the layers whose change rule the run means to change, as a " *
                 "symbol or a tuple of them, `(:Temperature,)`; got $(repr(names)).")
end

# One layer's change rule as a state recorded it, without its values.
_recordedchange(was) = (change = was.change, phrase = was.changephrase)

# For each layer, its change rule after the run's first step beside the rule it was saved with: its
# role and name, both rules as phrases, whether they differ, and whether `changing` named it.
function _comparechanges(saved::_SavedChanges, regime, supply)
    compared = @NamedTuple{role::String, layer::Symbol, saved::String,
                           now::String, changed::Bool, declared::Bool}[]
    for (role, layers, recorded) in (("regime", regime, saved.regime),
                                     ("supply", supply, saved.supply))
        for (name, layer) in zip(keys(layers), values(layers))
            was = recorded[string(name)]
            push!(compared,
                  (role = role, layer = name, saved = was.phrase,
                   now = _changephrase(layer.change),
                   changed = _changeidentity(layer.change) != was.change,
                   declared = name in saved.changing))
        end
    end
    return compared
end

# A changed layer as the run's provenance records it.
function _rulechange(c::NamedTuple)
    return (role = c.role, layer = string(c.layer), saved = c.saved,
            now = c.now, declared = c.declared)
end

# Warn of each layer whose change rule changed although `changing` did not name it, and of each layer
# `changing` named whose rule did not change.
function _warnchanges(saved::_SavedChanges, compared::AbstractVector)
    for c in compared
        if c.changed && !c.declared
            @warn "`$(saved.path)` was saved while its $(c.role) layer `$(c.layer)` changed by " *
                  "$(c.saved), and this run starts with it changing by $(c.now). A state keeps " *
                  "a layer's values but not its change rule, which a run takes from its habitat " *
                  "and from the interventions it is given, those due by its start installed at " *
                  "its first step. Where this is meant, name the layer in `build_ecosystem`'s " *
                  "`changing`."
        elseif c.declared && !c.changed
            @warn "`changing` names `$(c.layer)`, but this run starts with its $(c.role) layer " *
                  "changing by $(c.now), as it did when `$(saved.path)` was saved: give the " *
                  "habitat the rule it is meant to change by, or drop `$(c.layer)` from " *
                  "`changing`."
        end
    end
    return nothing
end

# A change rule as text to tell it from another by, its arrays summarised so that a series of any
# length gives a short one.
_changeidentity(change) = sprint(show, change, context = :limit => true)

# A change rule as a reader would name it: a steady rate by its value a year, any other by its kind.
_changephrase(::NoLayerChange) = "nothing"

function _changephrase(change::SteadyLayerChange)
    peryear = unit(change.value) * s / Unitful.yr
    return "a steady $(round(typeof(1.0 * peryear), change.value, sigdigits = 4))"
end

_changephrase(change) = "a $(nameof(typeof(change)))"

# Write the state's abundances into the run's rows, each row from its column of the state, a species
# the state does not hold left empty. Only the first rank read the state, so a serial ecosystem built
# under MPI, which every rank builds, is sent it whole; a distributed one's rows are the extension's.
function _placesaved!(eco::Ecosystem, saved, columns)
    whole = _madeonroot(() -> saved)
    for (row, column) in enumerate(columns)
        iszero(column) && continue
        eco.abundances.matrix[row, :] .= view(whole, column, :)
    end
    return eco
end

# Whether an ecosystem started from a saved state.
_fromstate(eco::AbstractEcosystem) = any(r -> r.role === :state, eco.inputs)
