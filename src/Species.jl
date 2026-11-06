# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Everything the simulation knows about the species themselves: what they tolerate, what they
# demand, how they move, their demographic rates, and the phylogeny relating them.

using Diversity
using Missings
using Distributions
using DataFrames

"""
    SpeciesList{TL <: AbstractTolerance, DM <: AbstractDemand,
                MO <: AbstractVector{<:SpeciesDispersal}, T <: AbstractTypes,
                P <: AbstractParams} <: AbstractTypes

Everything the simulation knows about the species themselves - the species-side counterpart of a
[`GridHabitat`](@ref), and one of the two halves [`build_ecosystem`](@ref) checks against each other
name for name.

It is itself a `Diversity.AbstractTypes`, but it **holds** the real types in `types` rather than
being them, so every `AbstractTypes` question is forwarded there; see `src/DiversityInterface.jl`.

# Fields

  - `names`: one name per species.
  - `tolerance`: the niche of each species, an [`AbstractTolerance`](@ref) matching the habitat's
    regime axis for axis.
  - `demand`: what each species takes from the shared pool, an [`AbstractDemand`](@ref) matching the
    habitat's supply the same way.
  - `types`: the similarity structure - a `UniqueTypes` where species are maximally distinct, or a
    `PhyloBranches` derived from a phylogeny.
  - `movement`: one [`SpeciesDispersal`](@ref) per species, saying which individuals disperse and
    how far.
  - `params`: the demographic rates, an [`AbstractParams`](@ref).
  - `susceptible`: optional per-species disease susceptibility, `missing` where unset.
  - `inputs`: the [`ProvenanceRecord`](@ref)s of the published data the species were built from -
    occurrence records, trait tables, parameters from the literature - as given with `provenance`, and
    a [`SpeciesPool`](@ref)'s own.
  - `traits`: each species' identity and values, one row per species - its `key` in the pool it came
    from, `name`, `id`, whether it is `active`, what it is a copy of - and the values it was built
    from, filled ([`fill_traits`](@ref)); for species built otherwise, the identity alone.
  - `filled`: which of `traits`' values were filled or drawn rather than given, a column each.
  - `authority`: the name authority the species' names were standardised through, or `nothing`.
  - `namemap`: the species' rows by name - their own, and the names their data used.

# Type parameters

  - `TL`, `DM`, `MO`, `T`, `P`: the types of `tolerance`, `demand`, `movement`, `types` and `params`.
    `T` is what a method keys on to require a phylogeny, since `Diversity`'s own `Sim` slot cannot
    carry that constraint.
"""
mutable struct SpeciesList{TL <: AbstractTolerance,
                           DM <: AbstractDemand,
                           MO <: AbstractVector{<:SpeciesDispersal},
                           T <: AbstractTypes,
                           P <: AbstractParams} <: AbstractSpecies
    names::Vector{String}
    tolerance::TL
    demand::DM
    types::T
    movement::MO
    params::P
    susceptible::Vector{Union{Missing, Float64}}
    inputs::Vector{ProvenanceRecord}
    traits::DataFrame
    filled::DataFrame
    authority::Union{Nothing, AbstractNameAuthority}
    namemap::Dict{String, Int}

    function SpeciesList{TL, DM, MO, T, P}(names::Vector{String},
                                           tolerance::TL,
                                           demand::DM,
                                           types::T,
                                           movement::MO,
                                           params::P;
                                           provenance = ProvenanceRecord[],
                                           traits = _identityframe(names),
                                           filled = DataFrame(),
                                           authority = nothing,
                                           namemap = _namemap(names)) where {TL <:
                                                                             AbstractTolerance,
                                                                             DM <:
                                                                             AbstractDemand,
                                                                             MO <:
                                                                             AbstractVector{<:SpeciesDispersal},
                                                                             T <:
                                                                             AbstractTypes,
                                                                             P <:
                                                                             AbstractParams}
        # Check dimensions
        equal_param = equalpop(params, length(names))
        sus = Vector{Union{Missing, Float64}}(undef, length(names))
        return new{TL, DM, MO, T, typeof(equal_param)}(names,
                                                       tolerance,
                                                       demand,
                                                       types,
                                                       movement,
                                                       equal_param,
                                                       sus,
                                                       _recordvector(provenance),
                                                       traits, filled,
                                                       authority, namemap)
    end
end
# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# The species-side counterpart of `GridHabitat`'s display, and the asymmetry of having one without
# the other is the argument for it: an ecosystem is species in an environment, and both halves
# should answer "what have I got here?" the same way.
#
# The similarity structure is named because it is what decides whether a diversity measure means
# anything - `UniqueTypes` and a phylogeny give different answers to the same question.
function Base.show(io::IO, sl::SpeciesList)
    return print(io,
                 "SpeciesList($(length(sl.names)) species, ",
                 sprint(show, sl.tolerance), " | ", sprint(show, sl.demand),
                 ")")
end

function Base.show(io::IO, ::MIME"text/plain", sl::SpeciesList)
    println(io, "SpeciesList")
    println(io, "  species   ", length(sl.names))
    println(io, "  tolerance ", sprint(show, sl.tolerance))
    println(io, "  demand    ", sprint(show, sl.demand))
    println(io, "  movement  ", _movementphrase(sl.movement))
    println(io, "  params    ", nameof(typeof(sl.params)))
    print(io, "  types     ", nameof(typeof(sl.types)))
    return nothing
end

# The two tree-generating constructors live in `EcoSISTEMPhyloExt`. They do not merely accept a
# phylogeny - they build one and derive a `PhyloBranches` similarity from it - so they cannot work
# without `Phylo`. Their docstrings stay here, immediately above the stubs, because `api.md`'s
# `@autodocs` cannot see into an extension.
#
# Every other `SpeciesList` constructor is unaffected: they take the tolerances, and the similarity,
# explicitly, which is the path `build_species` uses.
"""
    SpeciesList(numspecies::Int64, numtraits::Int64, demand::DM, movement::MO, params::P;
      seed = nothing)
    SpeciesList(numspecies::Int64, numtraits::Int64, demand::DM, movement::MO, params::P,
      switch::Vector{Float64}; seed = nothing)
    SpeciesList(numspecies::Int64, numtraits::Int64, demand::DM, movement::MO, phy::T,
      params::P; seed = nothing)

Create a `SpeciesList` whose categorical niche tolerances are **evolved along a random ultrametric
phylogeny**, with a `PhyloBranches` similarity structure computed from that tree.

The third form takes a similarity structure explicitly instead of computing one from the tree; the
tolerances are still evolved along a random ultrametric tree.

**All three forms need `Phylo`** - each builds the tree, or forwards to one that does - so
`using Phylo` is required and without it the call is a `MethodError`. The implementation is in
`EcoSISTEMPhyloExt`.

# Arguments

  - `numspecies`: how many species.
  - `numtraits`: how many categorical niche classes the tolerances are drawn from.
  - `demand`: what each species takes from the shared pool.
  - `movement`: one [`SpeciesDispersal`](@ref) per species, saying which individuals disperse and
    how far.
  - `params`: the demographic rates.
  - `switch`: the mean time between changes of the niche class along a root-to-tip path, as a
    fraction of its length, so a larger value changes it less often. Defaults to `[0.5]`.
  - `phy`: a similarity structure to use in place of the one the tree would give.
  - `seed`: an integer fixing the tree, and so the tolerances evolved along it. They are drawn
    from a stream of their own, so the same seed gives the same phylogeny whatever else has been
    drawn; left unset, they come from the global generator.
"""
function SpeciesList(numspecies::Int64,
                     numtraits::Int64,
                     demand::DM,
                     movement::MO,
                     params::P;
                     seed = nothing) where {DM <: AbstractDemand,
                                            MO <:
                                            AbstractVector{<:SpeciesDispersal},
                                            P <: AbstractParams}
    return SpeciesList(numspecies, numtraits, demand, movement, params, [0.5],
                       seed = seed)
end
# ---------------------------------------------------------------------------
# The size-based (`SizeDemand`) constructor is parked - still commented out - in
# `ext/EcoSISTEMPhyloExt/`. It evolves body mass along a phylogeny, so were it ever uncommented it
# would need `Phylo` like the two constructors above, and keeping it beside them keeps the open
# decision about `SizeDemand` in one place rather than two. See `src/Demand.jl`.

"""
    SpeciesList(numspecies::Int64, tolerance::TL, demand::DM, movement::MO, params::P;
      names = nothing, provenance = ProvenanceRecord[], pool = nothing)

Create a `SpeciesList` from requirements and rates already built, rather than from a
[`SpeciesPool`](@ref) by [`build_species`](@ref). Species are treated as maximally distinct, with a
`UniqueTypes` similarity structure. Where the individuals start is not the species list's: it is
the placement [`build_ecosystem`](@ref) is given.

# Arguments

  - `numspecies`: how many species.
  - `tolerance`: their niches, as an [`AbstractTolerance`](@ref).
  - `demand`: what each species takes from the shared pool.
  - `movement`: one [`SpeciesDispersal`](@ref) per species, saying which individuals disperse and
    how far. At most four kinds, mode and kernel together, so the step loop can split them.
  - `params`: the demographic rates.
  - `names`: one name per species, unique, as a vector of strings. Left unset, the species are
    named `"1"` to `"numspecies"`. The names label the species' Diversity types as well, so a
    diversity measure reports them.
  - `provenance`: an [`ProvenanceRecord`](@ref), or a vector of them, for the published data the
    species were built from; [`provenance`](@ref) of an ecosystem built from them lists these.
  - `pool`: what [`build_species`](@ref) read from a [`SpeciesPool`](@ref) - its filled values, which
    were filled, its authority and its names - or `nothing`, the default, for species built otherwise.
"""
function SpeciesList(numspecies::Int64,
                     tolerance::TL,
                     demand::DM,
                     movement::MO,
                     params::P;
                     names = nothing,
                     provenance = ProvenanceRecord[],
                     pool = nothing) where {TL <:
                                            AbstractTolerance,
                                            DM <:
                                            AbstractDemand,
                                            MO <:
                                            AbstractVector{<:SpeciesDispersal},
                                            P <:
                                            AbstractParams}
    names = _speciesnames(names, numspecies)
    # Species are maximally distinct, and the types carry the names so that nothing downstream
    # can lose them.
    ty = _uniquetypes(names)
    counttypes(demand) == numspecies ||
        throw(DimensionMismatch("$(counttypes(demand)) demands for $numspecies species: give " *
                                "one demand per species."))
    length(movement) == numspecies ||
        throw(DimensionMismatch("$(length(movement)) movements for $numspecies species: " *
                                "give one `SpeciesDispersal` per species."))
    movement = _dispersalvector(movement)
    return SpeciesList{typeof(tolerance),
                       typeof(demand),
                       typeof(movement),
                       typeof(ty),
                       typeof(params)}(names,
                                       tolerance,
                                       demand,
                                       ty,
                                       movement,
                                       params;
                                       provenance = provenance,
                                       _poolfields(pool, names)...)
end

# == Functions ==================================================================================

# The species-side half of the four-way accessor family; the docstring covering both methods lives
# with the family in `Ecosystem.jl`. This returns the **demand** itself. The total resource usage is
# `_getdemand(abun, demand)` in `Demand.jl`, which is what the hot loop reads.
getdemand(sppl::SpeciesList) = sppl.demand

# The only route from this package to a `UniqueTypes`, so that a species' name cannot be lost on the
# way. Diversity offers `UniqueTypes(count)` as well as `UniqueTypes(names)`, and the count form
# generates `"1"`, `"2"`, ... - right for a genuinely anonymous community, wrong wherever names are
# held, as they always are here. Leaving nothing that can reach the lossy form is what keeps that from
# happening twice.
#
# The renaming is invisible on the synthetic route, where `build_species` already yields numeric
# names. It bites a phylogeny-built `SpeciesList`, whose types carry tip labels, so a test of this has
# to use real names to prove anything.
function _uniquetypes(names::AbstractVector{<:AbstractString})
    return UniqueTypes(collect(String, names))
end

# --- Per-species fields from what a caller wrote ------------------------------

# The species' names as a `Vector{String}`: `"1"` to `"n"` when none were given, else the caller's
# vector checked for length and for uniqueness, since a duplicate name would make every lookup by
# name (`AddAbundance`, `AddSpecies`, the storage recorders) answer for the first of them only.
_speciesnames(::Nothing, n::Integer) = map(string, 1:n)

function _speciesnames(names::AbstractVector{<:AbstractString}, n::Integer)
    length(names) == n ||
        throw(DimensionMismatch("`names` has $(length(names)) entries but there are $n species"))
    allunique(names) ||
        throw(ArgumentError("`names` must be unique; " *
                            "$(join(repr.(unique(filter(x -> count(==(x), names) > 1, names))), ", ")) " *
                            "appear more than once"))
    return collect(String, names)
end

# Turn a per-species argument into a length-`n` vector: a scalar is filled to
# every species, a vector is passed through after checking it has exactly `n`
# entries (erroring clearly, and naming the argument, if it does not).
function _tofield(x::AbstractVector, n::Integer, name::AbstractString)
    length(x) == n ||
        throw(DimensionMismatch("`$name` has $(length(x)) entries but there are $n species"))
    return collect(x)
end

_tofield(x, n::Integer, name::AbstractString) = fill(x, n)

# The word standing for the abundance split in a stream's seed - "abund" in ASCII - so that the
# split's stream can never coincide with a species' `(seed, j)`, a niche layout's or an
# intervention's. A literal rather than `hash(:abundance)`, because `Base.hash` changes between
# Julia versions and the split must not.
const _ABUNDANCE_STREAM = 0x6162756e64

# The word standing for placement in a stream's seed - "place" in ASCII - distinct from the other
# tags for the same reason.
const _PLACEMENT_STREAM = 0x706c616365

# One generator per species for placing its individuals, keyed on its key in the pool it came from (a
# copy on its own), so that a species' placement depends neither on its position nor on any other
# species, and draws nothing from the streams its dynamics use.
function _placementrngs(seed::Integer, species::SpeciesList)
    return [Random.Xoshiro(UInt64[seed % UInt64, _PLACEMENT_STREAM,
                                  key % UInt64])
            for key in species.traits.key]
end

# The generator a placement's community `total` is split across the species from, a stream of its
# own made from the run's seed.
function _abundancerng(seed::Integer)
    return Random.Xoshiro(UInt64[seed % UInt64, _ABUNDANCE_STREAM])
end

# The word standing for a phylogeny in a stream's seed - "phylo" in ASCII - distinct from the other
# tags for the same reason. The tree-growing `SpeciesList` constructors in `EcoSISTEMPhyloExt` draw
# the tree, and the traits evolved along it, from `_phylogenyrng(seed)`.
const _PHYLOGENY_STREAM = 0x7068796c6f

# The generator a tree-growing constructor draws from: a stream of its own with a seed, the global
# one without.
_phylogenyrng(::Nothing) = Random.default_rng()

function _phylogenyrng(seed::Integer)
    return Random.Xoshiro(UInt64[seed % UInt64, _PHYLOGENY_STREAM])
end

# The identity of species built from no pool: their names, each its own key and row, active, and a
# copy of nothing.
function _identityframe(names::AbstractVector{<:AbstractString})
    n = length(names)
    return DataFrame(key = collect(1:n), name = collect(String, names),
                     id = Vector{Union{Missing, String}}(missing, n),
                     active = trues(n),
                     copyof = Vector{Union{Missing, String}}(missing, n),
                     copyofid = Vector{Union{Missing, String}}(missing, n),
                     copykey = Vector{Union{Missing, Int}}(missing, n))
end

# Each species' row by its own name.
function _namemap(names::AbstractVector{<:AbstractString})
    return Dict(n => i
                for (i, n) in enumerate(names))
end

# What a filter or a placement matches a table's names against, for a species list: its species' rows
# by their own names, by every name of their pool's data that reached one of them, and by identifier,
# whether each is active, and its authority's resolver, made on the first name not found otherwise.
function _matchindex(species::SpeciesList)
    traits = species.traits
    rowof(v) = Dict(x => r for (r, x) in enumerate(v) if !ismissing(x))
    return (name = rowof(traits.name), id = rowof(traits.id),
            given = Dict{Tuple{String, Union{Missing, String}}, Int}(),
            givenname = species.namemap,
            resolve = _lazyresolver(species.authority),
            active = Vector{Bool}(traits.active),
            speciesname = Vector{String}(traits.name),
            authority = species.authority, holder = "the species list")
end

# The keywords a species list takes from what was read of a pool, none without one.
_poolfields(::Nothing, names) = NamedTuple()

function _poolfields(pool::NamedTuple, names)
    return (traits = pool.traits, filled = pool.filled,
            authority = pool.authority, namemap = pool.namemap)
end

# Give a species just appended to `sl` its identity: a key after every other, active, and no values
# of its own (`missing`, none filled), found by its name.
function _addidentity!(sl::SpeciesList, name::AbstractString)
    traits = sl.traits
    DataFrames.allowmissing!(traits)
    row = Dict{Symbol, Any}(c => missing for c in propertynames(traits))
    row[:key] = maximum(traits.key, init = 0) + 1
    row[:name] = String(name)
    row[:active] = true
    push!(traits, row)
    isempty(propertynames(sl.filled)) || push!(sl.filled, Dict(c => false
                          for c in propertynames(sl.filled)))
    sl.namemap[String(name)] = nrow(traits)
    return sl
end
