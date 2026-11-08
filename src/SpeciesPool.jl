# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The species a run may hold, declared before it is built: a pool made from a table, a list of names
# or a count, each species once, standardised through an optional name authority, with the values
# known of each, the names its data used and the records of what it was read from. Changed only
# before the build, by the verbs below.

using DataFrames: DataFrames, DataFrame, nrow, eachrow
using Tables: Tables

"""
    AbstractMergeRule

What a [`SpeciesPool`](@ref) does where several rows of a table reach one species - names that a
name authority standardises to one, or one name written twice: [`RefuseMerge`](@ref), the default,
or [`MergeWith`](@ref) a function choosing the values to keep. A type rather than a symbol, so an
unknown rule is refused where it is written.
"""
abstract type AbstractMergeRule end

"""
    RefuseMerge()

Stop where several rows of a table reach one species, listing each species and the names that
reached it, so the choice between their values is made by whoever knows which to keep.
"""
struct RefuseMerge <: AbstractMergeRule end

"""
    MergeWith(f)

Keep one set of values where several rows of a table reach one species, chosen by `f`, which is
given those rows as a table and returns a row of values for the species - one of the rows it was
given (`MergeWith(first)` keeps the first), or a row it builds from them. Each merge is reported as
a notice.

# Fields

  - `f`: the function choosing the values.
"""
struct MergeWith{F} <: AbstractMergeRule
    f::F
end

Base.show(io::IO, m::MergeWith) = print(io, "MergeWith(", m.f, ")")

"""
    SpeciesPool(numspecies::Integer; verbosity = :normal)
    SpeciesPool(names::AbstractVector{<:AbstractString}; authority = nothing, verbosity = :normal)
    SpeciesPool(table; column = :species, authorship = nothing, authority = nothing,
                merge = RefuseMerge(), provenance = ProvenanceRecord[], verbosity = :normal)

Every species a run may hold, each once, declared before the run is built. A species may be
deactivated - kept, but marked to hold no individuals until it is activated - or removed, leaving it
out altogether.

With a name authority, every name the pool is given is standardised through it
([`EcoSISTEM.resolve_name`](@ref)): a synonym becomes its accepted species, and a name the
authority does not recognise is discarded with a warning. Several of a table's rows reaching one
species are merged by `merge`. Without one, names are matched exactly and nothing is discarded.

A pool is a table, one row per species: its `name` (the accepted name where an authority recognised
it), its `id` (the authority's identifier, `missing` without one), whether it is `active`, what it
is a copy of where it is one ([`copy_species!`](@ref): `copyof`, `copyofid`), then the values known
of it - the table's other columns and those [`add_trait!`](@ref) adds, `missing` where unknown.
[`name_table`](@ref) says how each name of its data was matched, and [`provenance`](@ref) what it
was read from. [`deactivate_species!`](@ref), [`activate_species!`](@ref) and
[`remove_species!`](@ref) change it.

# Arguments

  - `numspecies`: a pool of species named `"1"` to `"numspecies"`.
  - `names`: a pool of the species these name, a name given twice taken once.
  - `table`: a table of the species, one row each: a path to a CSV file, a
    [`ProvenanceRecord`](@ref) of one (its record kept), or any Tables.jl table.
  - `column`: the table's column of names.
  - `authorship`: the table's column of each name's authors, which pick among a homonym's rows,
    or `nothing`, the default, for none.
  - `authority`: the [name authority](@ref EcoSISTEM.AbstractNameAuthority) names are standardised
    through, such as a [`WFOPlantList`](@ref), or `nothing`, the default, for none.
  - `merge`: what to do where several rows reach one species, an
    [`EcoSISTEM.AbstractMergeRule`](@ref); [`RefuseMerge`](@ref) by default.
  - `provenance`: a [`ProvenanceRecord`](@ref), or a vector of them, for what the table was made
    from, where nothing else records it - a table in memory.
  - `verbosity`: `:silent` reports no notices, `:normal` reports them with the first few names,
    `:verbose` with every name; warnings and errors are always reported.

# Fields

  - `authority`: the name authority, or `nothing`.
  - `species`: the species, one row each: the identity columns above and the values, with each
    species' fixed `key` - its position in the table the pool was made from, a copy's after the
    last - which removing another species never changes, and for a copy its original's
    (`copykey`).
  - `names`: every name the pool's data used - its own table's and each [`add_trait!`](@ref)
    table's - with what it was matched to, how, and the species' key.
  - `inputs`: the [`ProvenanceRecord`](@ref)s of what the pool was read from.
  - `nextkey`: the key the next copy takes.
"""
mutable struct SpeciesPool{A <: Union{Nothing, AbstractNameAuthority}}
    authority::A
    species::DataFrame
    names::DataFrame
    inputs::Vector{ProvenanceRecord}
    nextkey::Int
end

function SpeciesPool(numspecies::Integer; verbosity::Symbol = :normal)
    numspecies >= 0 ||
        error("a pool holds a whole number of species; got `$numspecies`.")
    return SpeciesPool(string.(1:numspecies), verbosity = verbosity)
end

function SpeciesPool(names::AbstractVector{<:AbstractString};
                     authority::Union{Nothing, AbstractNameAuthority} = nothing,
                     verbosity::Symbol = :normal)
    return _madeonroot() do
        table = DataFrame(species = String.(names))
        return _newpool(authority, table, :species, nothing, MergeWith(first),
                        "names", ProvenanceRecord[], verbosity,
                        notifymerge = false)
    end
end

function SpeciesPool(table; column::Symbol = :species,
                     authorship::Union{Nothing, Symbol} = nothing,
                     authority::Union{Nothing, AbstractNameAuthority} = nothing,
                     merge::AbstractMergeRule = RefuseMerge(),
                     provenance = ProvenanceRecord[],
                     verbosity::Symbol = :normal)
    return _madeonroot() do
        source = _readtable(table, :species, column = column)
        return _newpool(authority, source.table, column, authorship, merge,
                        source.label,
                        vcat(source.records, _recordvector(provenance)),
                        verbosity)
    end
end

function Base.show(io::IO, pool::SpeciesPool)
    print(io, "SpeciesPool(", length(pool), " species, ",
          count(pool.species.active), " active")
    isnothing(pool.authority) || print(io, ", authority = ", pool.authority)
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", pool::SpeciesPool)
    println(io, "SpeciesPool")
    println(io, "  species   ", length(pool), " (", count(pool.species.active),
            " active)")
    copies = count(!ismissing, pool.species.copyof)
    copies > 0 && println(io, "  copies    ", copies)
    println(io, "  authority ",
            isnothing(pool.authority) ? "none" : sprint(show, pool.authority))
    traits = _traitnames(pool)
    println(io, "  traits    ",
            isempty(traits) ? "none" : join(string.(traits), ", "))
    print(io, "  inputs    ", length(pool.inputs), " record",
          length(pool.inputs) == 1 ? "" : "s")
    return nothing
end

Base.length(pool::SpeciesPool) = nrow(pool.species)

Tables.istable(::Type{<:SpeciesPool}) = true

Tables.columnaccess(::Type{<:SpeciesPool}) = true

Tables.columns(pool::SpeciesPool) = Tables.columns(_pooltable(pool))

# == Functions ==================================================================================

"""
    deactivate_species!(pool::SpeciesPool, names...; rank = "species", verbosity = :normal)
    deactivate_species!(f, pool::SpeciesPool; verbosity = :normal)

Deactivate species of `pool`, and return the pool: those `names` name, every species where no name
is given, or those `f` chooses. A deactivated species stays in the pool, marked to hold no
individuals, until [`activate_species!`](@ref) brings it back. How many species changed is noted.

# Arguments

  - `pool`: the pool to change.
  - `names`: names, or vectors of them, matched as the pool matches names - a copy by its label, any
    other name through the pool's authority, or exactly without one; a name that matches nothing is
    reported.
  - `rank`: the rank `names` are at, as the pool's authority writes it: `"species"`, the default, or
    a group such as `"genus"` or `"family"`, which selects every species the group holds by its
    classification ([`EcoSISTEM.lineage`](@ref)), a copy by its original's. A name found only at
    another rank is warned of. Any rank but `"species"` needs the pool's authority.
  - `f`: a function of a species' row of the pool, its identity and values, returning `true` for
    each species to change; with a do-block, `deactivate_species!(pool) do row ... end`.
  - `verbosity`: `:silent` notes nothing, `:normal` how many species changed, `:verbose` names every
    name that matched nothing.
"""
function deactivate_species!(pool::SpeciesPool, names...;
                             rank::AbstractString = "species",
                             verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    return _poolshared!(pool) do pool
        rows = isempty(names) ? collect(1:length(pool)) :
               _selectrows(pool, names, rank, "deactivate", verbosity)
        return _setactive!(pool, rows, false, verbosity)
    end
end

function deactivate_species!(f::Function, pool::SpeciesPool;
                             verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    return _poolshared!(pool) do pool
        return _setactive!(pool, _chosenrows(f, pool), false, verbosity)
    end
end

"""
    activate_species!(pool::SpeciesPool, names...; rank = "species", verbosity = :normal)
    activate_species!(f, pool::SpeciesPool; verbosity = :normal)

Activate species of `pool`, and return the pool: those `names` name, every species where no name is
given, or those `f` chooses - species [`deactivate_species!`](@ref) deactivated. How many species
changed is noted.

# Arguments

  - `pool`: the pool to change.
  - `names`, `rank`, `f`, `verbosity`: as for [`deactivate_species!`](@ref); with a do-block,
    `activate_species!(pool) do row ... end`.
"""
function activate_species!(pool::SpeciesPool, names...;
                           rank::AbstractString = "species",
                           verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    return _poolshared!(pool) do pool
        rows = isempty(names) ? collect(1:length(pool)) :
               _selectrows(pool, names, rank, "activate", verbosity)
        return _setactive!(pool, rows, true, verbosity)
    end
end

function activate_species!(f::Function, pool::SpeciesPool;
                           verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    return _poolshared!(pool) do pool
        return _setactive!(pool, _chosenrows(f, pool), true, verbosity)
    end
end

"""
    remove_species!(pool::SpeciesPool, names...; keep_active = false, rank = "species",
                    verbosity = :normal)
    remove_species!(f, pool::SpeciesPool; keep_active = false, verbosity = :normal)

Remove species from `pool`, and return the pool: those `names` name, or those `f` chooses; with
`keep_active`, only the inactive ones among them, or with no name given every inactive species - so
that `deactivate_species!(pool)`, `activate_species!(pool, "Acacia", rank = "genus")` and
`remove_species!(pool, keep_active = true)` leave the genus alone. A pool is never left empty.
Every other species keeps its key, its position in the table the pool was made from. The names a
removed species' data used stay in the [`name_table`](@ref), since they were used. How many species
were removed is noted, and the active species `keep_active` kept are named.

# Arguments

  - `pool`: the pool to change.
  - `names`, `rank`, `f`, `verbosity`: as for [`deactivate_species!`](@ref); with a do-block,
    `remove_species!(pool) do row ... end`.
  - `keep_active`: `true` removes only inactive species; `false`, the default, removes every species
    named or chosen, and is refused with no name, which would remove them all.
"""
function remove_species!(pool::SpeciesPool, names...; keep_active::Bool = false,
                         rank::AbstractString = "species",
                         verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    return _poolshared!(pool) do pool
        if isempty(names)
            keep_active ||
                error("`remove_species!` with no names would remove every species; name them, or " *
                      "pass `keep_active = true` to remove only the inactive ones.")
            rows = findall(!, pool.species.active)
        else
            rows = _selectrows(pool, names, rank, "remove", verbosity)
        end
        return _remove!(pool, rows, keep_active, verbosity)
    end
end

function remove_species!(f::Function, pool::SpeciesPool;
                         keep_active::Bool = false, verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    return _poolshared!(pool) do pool
        return _remove!(pool, _chosenrows(f, pool), keep_active, verbosity)
    end
end

"""
    copy_species!(pool::SpeciesPool, name::AbstractString; as)

Add a copy of the species `name` names to `pool`, labelled `as`, and return the pool: the same
species as a separate population, with its original's values, a key of its own after the last, and
its original's key beside it. The copy has no identifier of the authority's: it records its
original's name and identifier instead (`copyof`, `copyofid`), and is matched by its label alone,
never through the authority and never merged with its original. A copy of a copy records the first
original.

# Arguments

  - `pool`: the pool to add to.
  - `name`: the species to copy, matched as the pool matches names.
  - `as`: the copy's label, required: unique in the pool, and not a name the pool's authority
    recognises.
"""
function copy_species!(pool::SpeciesPool, name::AbstractString;
                       as::AbstractString = _requireas())
    return _poolshared!(pool) do pool
        label = String(as)
        label in pool.species.name &&
            _report(ProblemError(),
                    "`$label` already names a species of the pool, so it cannot label a copy.")
        if !isnothing(pool.authority)
            resolved = resolve_name(pool.authority, label)
            isnothing(resolved.name) ||
                _report(ProblemError(),
                        "`$label` is a name $(_authorityphrase(pool.authority)) recognises (as " *
                        "`$(resolved.name)`), so it cannot label a copy.")
        end
        rows = _poolrows(pool, (name,), "copy", :normal)
        isempty(rows) &&
            error("there is no species `$name` in the pool to copy.")
        original = copy(pool.species[only(rows), :])
        # A copy of a copy records the first original, as its own copy would.
        iscopy = !ismissing(original.copykey)
        push!(pool.species,
              merge(original,
                    (key = pool.nextkey, name = label, id = missing,
                     active = true,
                     copyof = iscopy ? original.copyof : original.name,
                     copyofid = iscopy ? original.copyofid : original.id,
                     copykey = iscopy ? original.copykey : original.key)))
        pool.nextkey += 1
        return pool
    end
end

"""
    add_trait!(pool::SpeciesPool, table; column = :species, authorship = nothing,
               merge = RefuseMerge(), replace = false, role = :species,
               provenance = ProvenanceRecord[], verbosity = :normal)
    add_trait!(f, pool::SpeciesPool, name::Symbol; replace = false)

Add values to the species of `pool`, and return the pool: the columns of a table keyed by species,
or a column computed from each species' row. Nothing is taken by position: a table's rows are
matched to species by name, through the pool's authority as the pool's own table was, and a
species with no row has `missing` values.

A table's names are reported as the pool's are: a name the authority does not recognise with a
warning, a recognised name of no species in the pool with a notice. Several rows reaching one
species are merged by `merge`. The table's names join the [`name_table`](@ref), and its record
the pool's inputs.

# Arguments

  - `pool`: the pool to add to.
  - `table`: a table of values, one row per species: a path to a CSV file, a
    [`ProvenanceRecord`](@ref) of one, or any Tables.jl table; every column but `column` and
    `authorship` is added.
  - `f`: a function of a species' row of the pool - a row of its table, identity and values -
    returning that species' value, called once per species; with a do-block,
    `add_trait!(pool, :range) do row ... end`.
  - `name`: the column `f`'s values go in.
  - `column`, `authorship`, `merge`, `provenance`, `verbosity`: as for [`SpeciesPool`](@ref).
  - `replace`: `true` lets a column replace one the pool already has; refused otherwise.
  - `role`: the part of the model the table's values are for, which its record is kept under:
    `:species`, or `:dispersal` or `:demography` for those traits.
"""
function add_trait!(pool::SpeciesPool, table; column::Symbol = :species,
                    authorship::Union{Nothing, Symbol} = nothing,
                    merge::AbstractMergeRule = RefuseMerge(),
                    replace::Bool = false, role::Symbol = :species,
                    provenance = ProvenanceRecord[],
                    verbosity::Symbol = :normal)
    _checkverbosity(verbosity)
    return _poolshared!(pool) do pool
        source = _readtable(table, role, column = column)
        input = source.table
        _checkcolumns(input, column, authorship)
        traitcols = _traitcolumns(input, column, authorship)
        _checknewcolumns(pool, traitcols, replace)
        given = _namecolumn(input, column)
        auths = _authorshipcolumn(input, authorship)
        index = _poolindex(pool)
        found = Dict{Tuple{String, Union{Nothing, String}}, Any}()
        groups = Dict{Int, Vector{Int}}()
        for i in eachindex(given)
            hit = get!(found, (given[i], auths[i])) do
                return _lookuprow(pool, index, given[i], auths[i])
            end
            isnothing(hit.row) || push!(get!(Vector{Int}, groups, hit.row), i)
        end
        problems = Problem[]
        _noteunmatched!(problems, found, verbosity)
        merged = _mergegroups!(problems, input, groups, traitcols, merge,
                               r -> pool.species.name[r], verbosity)
        # Everything is said, and a refused merge raised, before the pool changes.
        _reportproblems(problems, verbosity)
        for c in traitcols
            pool.species[!, c] = identity.(Any[haskey(merged, r) ?
                                               merged[r][c] :
                                               missing for r in 1:length(pool)])
        end
        append!(pool.names,
                _namerows(found, source.label, groups,
                          r -> (key = pool.species.key[r],
                                species = pool.species.name[r],
                                id = pool.species.id[r])))
        append!(pool.inputs,
                _recordvector(vcat(source.records,
                                   _recordvector(provenance))))
        return pool
    end
end

function add_trait!(f, pool::SpeciesPool, name::Symbol; replace::Bool = false)
    return _poolshared!(pool) do pool
        _checknewcolumns(pool, [name], replace)
        pool.species[!, name] = identity.(Any[f(row)
                                              for row in eachrow(_pooltable(pool))])
        return pool
    end
end

"""
    name_table(pool::SpeciesPool)

Return how each name of `pool`'s data was matched, as a table, one row per name of the pool's own
table and of each [`add_trait!`](@ref) table: the `name` as given, its `authorship` where one was
given, the `species` and `id` it reached (`missing` where it reached none), how it `match`ed - the
authority's answer (`:accepted`, `:synonym`, `:homonym`, `:unchecked`, `:not_found`,
`:unchecked_only`), `:exact` without an authority, or `:not_in_pool` for a recognised name of no
species of the pool - whether its row was `merged` with others reaching the same species, and the
`source` it came from. A removed species' names stay, since its data used them.

# Arguments

  - `pool`: the pool to report on.
"""
function name_table(pool::SpeciesPool)
    return DataFrames.select(pool.names, DataFrames.Not(:key))
end

"""
    provenance(pool::SpeciesPool)

Return the [`ProvenanceRecord`](@ref)s of what `pool` was read from: its table, each
[`add_trait!`](@ref) table, the name authority's release, and the records given as `provenance`.

# Arguments

  - `pool`: the pool to report on.
"""
provenance(pool::SpeciesPool) = copy(pool.inputs)

# Make a value once for every rank - a pool, a filter's answers for its names: under MPI rank 0 makes
# it, reading a name authority's release there alone, and every other rank takes the one it made.
# Every rank must make the call, in the same order, as layers are built.
function _madeonroot(build)
    _usempi() || return build()
    return _sharedwhole(build)
end

# Change a pool once for every rank, as `_madeonroot` makes one: the root changes its own and every
# other rank's is brought to the same state.
function _poolshared!(change!, pool::SpeciesPool)
    _usempi() || (change!(pool); return pool)
    shared = _sharedwhole() do
        change!(pool)
        return pool
    end
    shared === pool || _takestate!(pool, shared)
    return pool
end

# Bring `pool` to the state of `other`, field by field.
function _takestate!(pool::SpeciesPool, other::SpeciesPool)
    for f in fieldnames(SpeciesPool)
        setfield!(pool, f, getfield(other, f))
    end
    return pool
end

# A new pool: `table`'s names, standardised through `authority` where there is one and grouped into
# species - rows reaching one final name are one species, whose identifier joins every identifier
# they reached - merged by `merge` where a species has several rows, with the problems found
# reported at their levels. A species' key is its first row's position in `table`.
function _newpool(authority, table::DataFrame, column::Symbol, authorship,
                  merge::AbstractMergeRule, label::AbstractString,
                  records::Vector{ProvenanceRecord}, verbosity::Symbol;
                  notifymerge::Bool = true)
    _checkverbosity(verbosity)
    _checkcolumns(table, column, authorship)
    traitcols = _traitcolumns(table, column, authorship)
    clash = intersect(traitcols, _OWN_COLUMNS)
    isempty(clash) ||
        error("the table's column$(length(clash) == 1 ? "" : "s") " *
              "$(join(("`$c`" for c in clash), ", ")) would hide the pool's own; rename " *
              "$(length(clash) == 1 ? "it" : "them") first.")
    given = _namecolumn(table, column)
    auths = _authorshipcolumn(table, authorship)
    resolve = _givenresolver(authority)
    found = Dict{Tuple{String, Union{Nothing, String}}, Any}()
    groups = Dict{String, Vector{Int}}()
    order = String[]
    for i in eachindex(given)
        hit = get!(found, (given[i], auths[i])) do
            return resolve(given[i], auths[i])
        end
        isnothing(hit.name) && continue
        haskey(groups, hit.name) || push!(order, hit.name)
        push!(get!(Vector{Int}, groups, hit.name), i)
    end
    problems = Problem[]
    _notediscarded!(problems, found, authority, verbosity)
    _noteresolved!(problems, found, verbosity)
    # A name given twice in a list is one species, which needs no merge reported.
    merged = _mergegroups!(notifymerge ? problems : Problem[], table, groups,
                           traitcols, merge, identity, verbosity)
    _reportproblems(problems, verbosity)
    n = length(order)
    species = DataFrame(key = Int[first(groups[s]) for s in order],
                        name = copy(order),
                        id = Union{Missing, String}[_groupid(authority, found,
                                                             given, auths,
                                                             groups[s])
                                                    for s in order],
                        active = trues(n),
                        copyof = Vector{Union{Missing, String}}(missing, n),
                        copyofid = Vector{Union{Missing, String}}(missing, n),
                        copykey = Vector{Union{Missing, Int}}(missing, n))
    for c in traitcols
        species[!, c] = identity.(Any[merged[s][c] for s in order])
    end
    position = Dict(s => k for (k, s) in enumerate(order))
    names = _namerows(found, label, groups,
                      s -> (key = species.key[position[s]], species = s,
                            id = species.id[position[s]]))
    inputs = copy(records)
    isnothing(authority) || append!(inputs, _authorityrecords(authority))
    return SpeciesPool(authority, species, names, _recordvector(inputs),
                       nrow(table) + 1)
end

# The columns of a pool's table that are its own, shown, and those it keeps for itself; no table's
# column may take either.
const _IDENTITY_COLUMNS = [:name, :id, :active, :copyof, :copyofid]

const _OWN_COLUMNS = vcat(_IDENTITY_COLUMNS, [:key, :copykey])

# The columns of a pool that are values.
function _traitnames(pool::SpeciesPool)
    return setdiff(propertynames(pool.species), _OWN_COLUMNS)
end

# A pool as the table it is read as: its identity columns, then its values.
function _pooltable(pool::SpeciesPool)
    return DataFrames.select(pool.species, _IDENTITY_COLUMNS, _traitnames(pool))
end

# A table given to a pool or a placement, read into a `DataFrame` with its label for the name table
# and the records of where it came from, under `role`: a CSV file by its path, with its sidecar's
# record where it has one and its name otherwise; a record's file, with the record; any Tables.jl table
# as it is, with none. A file's `column` of species names is read as text, so a name made of digits
# keeps its own spelling (`007`) rather than becoming a number.
function _readtable(path::AbstractString, role::Symbol;
                    column::Symbol = :species)
    isfile(path) ||
        error("`$path` is not a file, so no table can be read from it.")
    record = something(Provenancer.recordat(path),
                       ProvenanceRecord(role = role, source = "file",
                                        path = basename(path)))
    return (table = _readcsv(path, column), label = basename(path),
            records = [record])
end

function _readtable(record::ProvenanceRecord, ::Symbol;
                    column::Symbol = :species)
    path = _specfile(SpeciesPool, record).location
    return (table = _readcsv(path, column), label = basename(path),
            records = [record])
end

# Deliberately untyped: anything else is a table or a mistake, told apart by `Tables.istable` so
# the mistake is named rather than raised as a `MethodError` deep in DataFrames.
function _readtable(table, ::Symbol; column::Symbol = :species)
    Tables.istable(table) ||
        error("a table is a path to a CSV file, a `ProvenanceRecord` of one, or a Tables.jl " *
              "table; got a $(typeof(table)).")
    return (table = DataFrame(table), label = "table",
            records = ProvenanceRecord[])
end

# A CSV file as a `DataFrame`, its `column` of names read as text. Not validated, so a file without
# that column is reported where its columns are checked rather than by the reader.
function _readcsv(path::AbstractString, column::Symbol)
    return DataFrame(CSV.File(path, types = Dict(column => String),
                              validate = false))
end

# Refuse a name or authorship column the table does not have, listing those it has.
function _checkcolumns(table::DataFrame, column::Symbol, authorship)
    cols = propertynames(table)
    for c in filter(!isnothing, [column, authorship])
        c in cols ||
            error("the table has no column `$c`; its columns are " *
                  "$(join(("`$x`" for x in cols), ", ")).")
    end
    return nothing
end

# The columns of a table that are values: every one but its names and their authors.
function _traitcolumns(table::DataFrame, column::Symbol, authorship)
    return [c for c in propertynames(table) if c != column && c != authorship]
end

# Refuse a column a pool already has, or one of its own, unless replacing is asked for.
function _checknewcolumns(pool::SpeciesPool, cols, replace::Bool)
    own = intersect(cols, _OWN_COLUMNS)
    isempty(own) ||
        error("$(join(("`$c`" for c in own), ", ")) name$(length(own) == 1 ? "s" : "") a " *
              "column of the pool's own, which no value can replace; rename it first.")
    replace && return nothing
    have = intersect(cols, _traitnames(pool))
    isempty(have) ||
        error("the pool already has $(join(("`$c`" for c in have), ", ")); pass " *
              "`replace = true` to replace $(length(have) == 1 ? "it" : "them").")
    return nothing
end

# A table's names, refusing a row with none.
function _namecolumn(table::DataFrame, column::Symbol)
    values = table[!, column]
    blank = findall(ismissing, values)
    isempty(blank) ||
        error("row$(length(blank) == 1 ? "" : "s") $(join(first(blank, 10), ", ")) of the " *
              "table name$(length(blank) == 1 ? "s" : "") no species in `$column`.")
    return String[string(v) for v in values]
end

# A table's authorships, `nothing` for a row with none, or for every row where no column is named.
function _authorshipcolumn(table::DataFrame, authorship)
    isnothing(authorship) && return fill(nothing, nrow(table))
    return Union{Nothing, String}[ismissing(v) ? nothing : string(v)
                                  for v in table[!, authorship]]
end

# How the names of a pool's own table are resolved: through the authority where there is one, each
# as itself otherwise.
function _givenresolver(::Nothing)
    return (name, _) -> (name = String(name), id = nothing, match = :exact)
end

_givenresolver(authority::AbstractNameAuthority) = _resolver(authority)

# A species' identifier: the one its names reached, or every identifier they reached joined by `+`
# where they reached several; `missing` without an authority.
_groupid(::Nothing, found, given, auths, rows) = missing

function _groupid(::AbstractNameAuthority, found, given, auths, rows)
    ids = unique(found[(given[r], auths[r])].id for r in rows)
    length(ids) == 1 && return only(ids)
    return join(sort!(unique(vcat(split.(ids, '+')...))), "+")
end

# The species each group of a table's rows gives its values, by the merge rule where a group has
# several rows: a single row's own values, or those `merge` chooses. `groups` maps a species to its
# rows, `label` names a species for the report.
function _mergegroups!(problems, table::DataFrame, groups, traitcols,
                       merge::AbstractMergeRule, label, verbosity)
    several = sort!([s for (s, rows) in groups if length(rows) > 1],
                    by = s -> first(groups[s]))
    isempty(several) || _merge!(problems, merge, table, groups, several, label,
            verbosity)
    values = Dict{keytype(groups), Dict{Symbol, Any}}()
    for (s, rows) in groups
        row = length(rows) == 1 ? table[only(rows), :] :
              _mergedrow(merge, DataFrames.view(table, rows, :))
        values[s] = Dict{Symbol, Any}(c => _rowvalue(row, c) for c in traitcols)
    end
    return values
end

# Report the species several rows reached: refused, listing them, or merged, noted.
function _merge!(problems, ::RefuseMerge, table, groups, several, label,
                 verbosity)
    return push!(problems,
                 Problem(ProblemError(), :merge_refused,
                         "$(length(several)) species $(length(several) == 1 ? "is" : "are") " *
                         "reached by several rows: " *
                         _namelist([_mergephrase(table, groups, s, label)
                                    for s in several], verbosity) *
                         ". Choose which values to keep with `merge = MergeWith(f)`, or remove " *
                         "the rows."))
end

function _merge!(problems, ::MergeWith, table, groups, several, label,
                 verbosity)
    return push!(problems,
                 Problem(ProblemNotice(), :merged,
                         "$(length(several)) species $(length(several) == 1 ? "was" : "were") " *
                         "reached by several rows, merged: " *
                         _namelist([_mergephrase(table, groups, s, label)
                                    for s in several], verbosity) * "."))
end

# One merged species as a report names it: the species, and the rows that reached it by their row
# numbers.
function _mergephrase(table, groups, s, label)
    rows = groups[s]
    return "`$(label(s))` (rows $(join(rows, ", ")))"
end

# The row of values a merge keeps.
_mergedrow(::RefuseMerge, rows) = first(rows)

function _mergedrow(merge::MergeWith, rows)
    row = merge.f(rows)
    isnothing(row) &&
        error("the function of `MergeWith` returned `nothing` for rows " *
              "$(join(DataFrames.parentindices(rows)[1], ", ")); it returns the row of values to keep.")
    return row
end

# A value of a row - a `DataFrameRow`, a named tuple or any Tables.jl row - by its column.
function _rowvalue(row, c::Symbol)
    c in Tables.columnnames(row) ||
        error("the row `MergeWith`'s function returned has no `$c`; it returns a row with every " *
              "value column of the table.")
    return Tables.getcolumn(row, c)
end

# What a pool looks names up by, built once for each operation: its species' rows by name and by
# identifier, the names of its data by name and authorship, and by name alone where a name reached
# one species whatever its authorship, and - where `resolve` asks for it - its authority's resolver,
# which reads the authority's names.
function _poolindex(pool::SpeciesPool; resolve::Bool = true)
    species = pool.species
    rowof(v) = Dict(x => r for (r, x) in enumerate(v) if !ismissing(x))
    bykey = rowof(species.key)
    given = Dict{Tuple{String, Union{Missing, String}}, Int}()
    keys = Dict{String, Set{Int}}()
    for r in eachrow(pool.names)
        (ismissing(r.key) || !haskey(bykey, r.key)) && continue
        given[(r.name, r.authorship)] = bykey[r.key]
        push!(get!(Set{Int}, keys, r.name), bykey[r.key])
    end
    byname = Dict(n => only(rows) for (n, rows) in keys if length(rows) == 1)
    return (name = rowof(species.name), id = rowof(species.id), given = given,
            givenname = byname,
            resolve = (resolve && !isnothing(pool.authority)) ?
                      _resolver(pool.authority) : nothing)
end

# Each species' row by name, for the species list built from a pool: its own name, and every name of
# the pool's data that reached it alone.
function _speciesnamemap(pool::SpeciesPool)
    index = _poolindex(pool, resolve = false)
    return Dict{String, Int}(merge(index.givenname, index.name))
end

# A name looked up in a pool, as `add_trait!` and the verbs do: its species' `row` in the pool, or
# `nothing`, with what the authority said of it - a recognised name reaching no species is
# `:not_in_pool`. Found by a species' own name (a copy's label among them), by a name of the pool's
# data, then through the authority by the identifier it reaches, else by its accepted name, since a
# name reached by another route may carry another identifier.
function _lookuprow(pool::SpeciesPool, index, name::AbstractString,
                    authorship = nothing)
    own = get(index.name, name, nothing)
    if isnothing(pool.authority)
        return (name = isnothing(own) ? nothing : String(name), id = nothing,
                match = isnothing(own) ? :not_found : :exact, row = own)
    end
    resolved = index.resolve(name, authorship)
    row = something(_directrow(index, name, authorship),
                    _resolvedrow(index, resolved), Some(nothing))
    # A species' own name the authority does not know is a copy's label, matched exactly.
    !isnothing(own) && isnothing(resolved.name) &&
        return (name = String(name), id = nothing, match = :exact, row = own)
    match = isnothing(row) && !isnothing(resolved.name) ? :not_in_pool :
            resolved.match
    return (name = resolved.name, id = resolved.id, match = match, row = row)
end

# A name's species row without asking the authority: by a species' own name (a copy's label among
# them), by a name of the pool's data with its authorship, or by that name alone where it reached one
# species; `nothing` where none.
function _directrow(index, name::AbstractString, authorship = nothing)
    return something(get(index.name, name, nothing),
                     get(index.given, (name, something(authorship, missing)),
                         nothing),
                     get(index.givenname, name, nothing), Some(nothing))
end

# The species row an authority's answer reaches: by the identifier it gives, else by its accepted
# name, since a name reached by another route may carry another identifier; `nothing` where none.
function _resolvedrow(index, resolved::NamedTuple)
    return something(isnothing(resolved.id) ? nothing :
                     get(index.id, resolved.id, nothing),
                     isnothing(resolved.name) ? nothing :
                     get(index.name, resolved.name, nothing), Some(nothing))
end

# What a filter or a placement matches a table's names against, for a pool: its species' rows by name
# and identifier as the verbs find them, whether each is active and how one is activated, and its
# authority's resolver, made on the first name not found otherwise.
function _matchindex(pool::SpeciesPool)
    return merge(_poolindex(pool, resolve = false),
                 (resolve = _lazyresolver(pool.authority),
                  active = Vector{Bool}(pool.species.active),
                  speciesname = Vector{String}(pool.species.name),
                  authority = pool.authority, holder = "the pool",
                  activate = "`activate_species!` activates them."))
end

# The rows of the species `names` name, for a verb, reporting the names that reach none: a name the
# authority does not recognise with a warning, a recognised name of no species in the pool with a
# notice.
function _poolrows(pool::SpeciesPool, names, verb::AbstractString,
                   verbosity::Symbol)
    wanted = String[n for x in names for n in (x isa AbstractString ? (x,) : x)]
    index = _poolindex(pool)
    rows = Int[]
    unrecognised = String[]
    absent = String[]
    for n in wanted
        hit = _lookuprow(pool, index, n)
        if !isnothing(hit.row)
            push!(rows, hit.row)
        elseif isnothing(hit.name)
            push!(unrecognised, n)
        else
            push!(absent, n)
        end
    end
    problems = Problem[]
    isempty(unrecognised) ||
        push!(problems,
              Problem(ProblemWarning(), :not_recognised,
                      "could not $verb $(_namelist(unrecognised, :verbose)): " *
                      (isnothing(pool.authority) ?
                       "no species of the pool has that name" :
                       "$(_authorityphrase(pool.authority)) does not recognise " *
                       (length(unrecognised) == 1 ? "it" : "them")) * "."))
    isempty(absent) ||
        push!(problems,
              Problem(ProblemNotice(), :not_in_pool,
                      "could not $verb $(_namelist(absent, :verbose)): recognised, but of no " *
                      "species in the pool."))
    _reportproblems(problems, verbosity)
    return rows
end

# The rows of the species the groups `names` name at `rank` hold, for a verb: each name resolved
# through the pool's authority, the parts of its identifier at that rank kept - so a name written
# alike at several ranks is read at the one asked for - and each species chosen where its
# classification holds one, a copy by its original's. A name not recognised, or found only at
# other ranks, is warned of. `"species"` reads names as the pool does.
function _selectrows(pool::SpeciesPool, names, rank::AbstractString,
                     verb::AbstractString, verbosity::Symbol)
    rank == "species" && return _poolrows(pool, names, verb, verbosity)
    authority = pool.authority
    isnothing(authority) &&
        error("`rank = \"$rank\"` needs the pool's authority to classify its species, and this " *
              "pool has none; name the species themselves.")
    wanted = String[n for x in names for n in (x isa AbstractString ? (x,) : x)]
    groups = Set{String}()
    unrecognised = String[]
    elsewhere = String[]
    for n in wanted
        resolved = resolve_name(authority, n)
        if isnothing(resolved.id)
            push!(unrecognised, n)
            continue
        end
        parts = split(resolved.id, '+')
        ranks = [first(lineage(authority, p)).rank for p in parts]
        atrank = parts[ranks .== rank]
        isempty(atrank) &&
            push!(elsewhere, "`$n` (a $(join(unique(ranks), " and a ")))")
        union!(groups, atrank)
    end
    problems = Problem[]
    isempty(unrecognised) ||
        push!(problems,
              Problem(ProblemWarning(), :not_recognised,
                      "could not $verb $(_namelist(unrecognised, :verbose)): " *
                      "$(_authorityphrase(authority)) does not recognise " *
                      (length(unrecognised) == 1 ? "it" : "them") * "."))
    isempty(elsewhere) ||
        push!(problems,
              Problem(ProblemWarning(), :other_rank,
                      "could not $verb $(join(elsewhere, ", ")): not a $rank in " *
                      "$(_authorityphrase(authority))."))
    _reportproblems(problems, verbosity)
    rows = Int[]
    for (r, s) in enumerate(eachrow(pool.species))
        id = ismissing(s.copyofid) ? s.id : s.copyofid
        ismissing(id) && continue
        any(l -> l.id in groups, lineage(authority, id)) && push!(rows, r)
    end
    return rows
end

# The rows of the species `f` chooses, given each species' row of the pool; anything but `true` or
# `false` is refused, naming the species, since a value that may be missing would choose silently.
function _chosenrows(f, pool::SpeciesPool)
    rows = Int[]
    for (r, row) in enumerate(eachrow(_pooltable(pool)))
        chosen = f(row)
        chosen isa Bool ||
            error("the function choosing species returned $(repr(chosen)) for `$(row.name)`; it " *
                  "returns `true` or `false`, so `coalesce` a value that may be missing.")
        chosen && push!(rows, r)
    end
    return rows
end

# Mark the species at `rows` active or not, noting how many changed.
function _setactive!(pool::SpeciesPool, rows, active::Bool, verbosity::Symbol)
    rows = unique(rows)
    already = count(r -> pool.species.active[r] == active, rows)
    pool.species.active[rows] .= active
    state = active ? "active" : "inactive"
    _reportproblems([Problem(ProblemNotice(),
                             active ? :activated : :deactivated,
                             "$(length(rows) - already) species $(active ? "activated" : "deactivated")" *
                             (already > 0 ? " ($already already $state)" : "") *
                             "; $(count(pool.species.active)) of $(length(pool)) active.")],
                    verbosity)
    return pool
end

# Remove the species at `rows` - keeping, and naming, the active ones where asked - refusing to leave
# the pool empty, and noting how many were removed.
function _remove!(pool::SpeciesPool, rows, keep_active::Bool, verbosity::Symbol)
    rows = sort!(unique(rows))
    problems = Problem[]
    if keep_active
        kept = filter(r -> pool.species.active[r], rows)
        filter!(r -> !pool.species.active[r], rows)
        isempty(kept) ||
            push!(problems,
                  Problem(ProblemNotice(), :kept_active,
                          "kept $(length(kept)) active species: " *
                          _namelist(pool.species.name[kept], verbosity) *
                          "; `keep_active = false` removes them."))
    end
    (!isempty(rows) && length(rows) == length(pool)) &&
        error("removing $(length(rows) == 1 ? "it" : "these $(length(rows)) species") would " *
              "leave the pool with no species.")
    deleteat!(pool.species, rows)
    push!(problems,
          Problem(ProblemNotice(), :removed,
                  "$(length(rows)) species removed; $(length(pool)) " *
                  (length(pool) == 1 ? "remains." : "remain.")))
    _reportproblems(problems, verbosity)
    return pool
end

# The name table's rows for the names one table used, in name order: each distinct name with its
# authorship, what it reached and how, whether its row was merged with others, and the table it
# came from. A name's hit is keyed into `groups` by its `row` where it has one, else by the species
# name it reached; `reached` gives that species' key, name and identifier.
function _namerows(found, source::AbstractString, groups, reached)
    entries = sort!(collect(found),
                    by = e -> (first(e)[1], something(first(e)[2], "")))
    target(hit) = haskey(hit, :row) ? hit.row : hit.name
    rows = [isnothing(target(h)) ? nothing : reached(target(h))
            for (_, h) in entries]
    field(r, f) = isnothing(r) ? missing : getfield(r, f)
    return DataFrame(name = String[first(k) for (k, _) in entries],
                     authorship = Union{Missing, String}[something(last(k),
                                                                   missing)
                                                         for (k, _) in entries],
                     species = Union{Missing, String}[field(r, :species)
                                                      for r in rows],
                     id = Union{Missing, String}[field(r, :id) for r in rows],
                     match = Symbol[h.match for (_, h) in entries],
                     merged = Bool[!isnothing(target(h)) &&
                                       length(get(groups, target(h), ())) > 1
                                   for (_, h) in entries],
                     source = fill(String(source), length(entries)),
                     key = Union{Missing, Int}[field(r, :key) for r in rows])
end

# The release records of a pool's authority, once its names have been read.
function _authorityrecords(authority::AbstractNameAuthority)
    return ProvenanceRecord[r for r in provenance(authority) if !isnothing(r)]
end

# Note the names an authority discarded, by why: not checked, or not in it at all.
_notediscarded!(problems, found, ::Nothing, verbosity) = problems

function _notediscarded!(problems, found, authority::AbstractNameAuthority,
                         verbosity)
    missed = [(n, h.match) for ((n, _), h) in found if isnothing(h.name)]
    isempty(missed) && return problems
    sort!(missed)
    unchecked = count(m -> m[2] === :unchecked_only, missed)
    return push!(problems,
                 Problem(ProblemWarning(), :names_discarded,
                         "$(length(missed)) name$(length(missed) == 1 ? "" : "s") discarded, " *
                         "not recognised by $(_authorityphrase(authority)) ($unchecked only " *
                         "unchecked, $(length(missed) - unchecked) not in it): " *
                         _namelist(first.(missed), verbosity) *
                         "; `name_table(pool)` lists every name and how it was matched."))
end

# Note the names an authority standardised: translated synonyms, homonyms kept under their joined
# identifiers, and unchecked names kept under their own.
function _noteresolved!(problems, found, verbosity)
    for (code, noun, one, several) in ((:synonym, "synonym", "translated to its accepted name",
        "translated to their accepted names"),
        (:homonym, "homonym",
        "kept under its identifiers joined by `+`",
        "kept under their identifiers joined by `+`"),
        (:unchecked, "unchecked name",
        "kept under its own identifier",
        "kept under their own identifiers"))
        hits = sort!([n for ((n, _), h) in found if h.match === code])
        isempty(hits) && continue
        n = length(hits)
        push!(problems,
              Problem(ProblemNotice(), code,
                      "$n $noun$(n == 1 ? "" : "s") $(n == 1 ? one : several): " *
                      _namelist(hits, verbosity) * "."))
    end
    return problems
end

# Note the names of an `add_trait!` table that reached no species of the pool: not recognised, with
# a warning, or recognised but of no species in it, with a notice.
function _noteunmatched!(problems, found, verbosity)
    unrecognised = sort!([n
                          for ((n, _), h) in found
                          if isnothing(h.row) && h.match !== :not_in_pool])
    absent = sort!([n for ((n, _), h) in found if h.match === :not_in_pool])
    isempty(unrecognised) ||
        push!(problems,
              Problem(ProblemWarning(), :not_recognised,
                      "$(length(unrecognised)) name$(length(unrecognised) == 1 ? "" : "s") " *
                      "of the table reach$(length(unrecognised) == 1 ? "es" : "") no species, " *
                      "not recognised: " *
                      _namelist(unrecognised, verbosity) * "."))
    isempty(absent) ||
        push!(problems,
              Problem(ProblemNotice(), :not_in_pool,
                      "$(length(absent)) name$(length(absent) == 1 ? "" : "s") of the table " *
                      "recognised, but of no species in the pool: " *
                      _namelist(absent, verbosity) * "."))
    return problems
end

# Names as a report lists them: every one at `:verbose`, otherwise the first ten and how many more.
function _namelist(names::AbstractVector, verbosity::Symbol)
    quoted = [startswith(n, '`') ? n : "`$n`" for n in names]
    (verbosity === :verbose || length(quoted) <= 10) &&
        return join(quoted, ", ")
    return join(first(quoted, 10), ", ") * " and $(length(quoted) - 10) more"
end

# Raise a pool's or a placement's problems at their levels, notices only where `verbosity` asks for
# them, warnings and errors always - an error last, once everything else is said.
function _reportproblems(problems, verbosity::Symbol)
    for p in sort(problems, by = p -> p.severity isa ProblemError)
        p.severity isa ProblemNotice && verbosity === :silent && continue
        _report(p.severity, p.message)
    end
    return nothing
end

# Refuse a `verbosity` a pool does not know.
function _checkverbosity(verbosity::Symbol)
    verbosity in (:silent, :normal, :verbose) ||
        error("`verbosity` must be `:silent`, `:normal` or `:verbose`; got `:$verbosity`.")
    return nothing
end

# The error a copy given no label raises.
function _requireas()
    return error("`copy_species!` needs the copy's label, `as = ...`: unique in the pool, and not a " *
                 "name the pool's authority recognises.")
end
