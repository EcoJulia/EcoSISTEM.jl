# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The simulation state itself - species, environment and abundances assembled - and the scratch
# space the hot loop writes into rather than allocating.

using Diversity

using Unitful

using Missings

using Random

using DimensionalData

using DimensionalData: Ti, At

using HCubature

using Distributions: ContinuousUnivariateDistribution, InverseGaussian, pdf

using SpecialFunctions: erfc

using EcoSISTEM.Units

using RecipesBase

using Dates: Dates

using Diversity.API: _calcabundance

"""
    Cache

Cache houses an integer array of moves made by all species in a timestep for the
[`update!`](@ref) function, `netmigration`.
"""
mutable struct Cache
    netmigration::Matrix{Int64}
    totaldemand::Matrix{Float64}
    # **Scratch for the normalised abundances behind the ordinariness, the weights, the
    # metacommunity abundances and the scale** (`src/DiversityInterface.jl`). Without it each of
    # those allocates a `Float64` matrix of the abundances' shape on every call - 7 GiB at 44,000
    # species on 21,000 cells - and every diversity measure asks for two of them.
    # **It is scratch, NOT a cache**: refilled on every call and never trusted, since nothing would
    # invalidate a cached value when the ecosystem changes. `_calcabundance`/`_calcordinariness`
    # **pass a `UniqueTypes` array straight through**, so the ordinariness handed back aliases this
    # buffer, and a caller wanting to keep it across calls must copy it.
    relativeabun::Matrix{Float64}
    valid::Bool

    # The scratch is an implementation detail of the cache, not something every construction site
    # should have to know about - so it is derived from `netmigration`'s shape here.
    function Cache(netmigration::Matrix{Int64}, totaldemand::Matrix{Float64},
                   valid::Bool)
        return new(netmigration, totaldemand, similar(netmigration, Float64),
                   valid)
    end
end

# A one-liner is the whole of what a cache needs: the default prints all four scratch arrays. Their
# contents are meaningless between timesteps - `relativeabun` is refilled on every call and never
# trusted - so its shape is the only fact worth showing.
function Base.show(io::IO, c::Cache)
    nsp, ncells = size(c.netmigration)
    return print(io, "Cache($(nsp) species × $(ncells) cells)")
end

# What a saved state recorded of each layer's change rule, by role and then layer name. A state
# holds a layer's values but not its rule, which is code. An ecosystem started from the state keeps
# this until its first step has installed the rules its interventions give, and then compares the
# two. `changing` names the layers the run was told would change rule. A name covers a regime and a
# supply layer of that name alike.
struct _SavedChanges
    path::String
    regime::Dict{String, @NamedTuple{change::String, phrase::String}}
    supply::Dict{String, @NamedTuple{change::String, phrase::String}}
    changing::Vector{Symbol}
end

# A layer whose change rule changed when a run continued a saved state, as the run's provenance
# records it: the layer's role and name, the rule it was saved with, the rule it changes by now, and
# whether `changing` named it.
const _RuleChange = @NamedTuple{role::String, layer::String, saved::String,
                                now::String, declared::Bool}

"""
    Ecosystem{Part <: AbstractHabitat} <:
       AbstractEcosystem{Part, SL, NF}

Ecosystem houses information on species and their interaction with their
environment. For species, it holds abundances and locations, as well as
properties such as trait information, `spplist`, and dispersal neighbourhoods, `lookup`.
For environments, it provides information on environmental conditions and
available resources,`habitat`. Finally, there is a slot for the nichefit
between the environment and the characteristics of the species, `nichefit`.
`elapsed` is the simulation clock (see [`simulationtime`](@ref)) and `seed` the
seed the per-species RNG streams were derived from. `epoch` is the real date elapsed time zero
corresponds to, or `nothing` for a run with no calendar - see [`simulationdate`](@ref) - and
`calendar` the [`AbstractRunCalendar`](@ref) its dates are placed by. `inputs`
holds the [`ProvenanceRecord`](@ref)s of published data belonging to the run as a whole - those given to
[`build_ecosystem`](@ref), those its placement read, and those of every intervention that has acted -
which [`provenance`](@ref) lists. `placement` is the placement its individuals were placed by, as it
was given (see [`AbstractPlacement`](@ref)), which [`repopulate!`](@ref) places by again.
`savedchanges` is what the saved state it was started from recorded of each layer's change rule, held
until its first step, and otherwise `nothing`. `rulechanges` lists each layer whose change rule
that first step found changed from the saved one, which [`provenance`](@ref) reports. It is empty
otherwise.
"""
mutable struct Ecosystem{Part <: AbstractHabitat,
                         SL <: SpeciesList,
                         NF <: AbstractNicheFit} <:
               AbstractEcosystem{Part, SL, NF}
    abundances::GridLandscape
    spplist::SL
    habitat::Part
    nichefit::NF
    lookup::Vector{Lookup}
    cache::Cache
    rngs::Vector{Random.Xoshiro}
    elapsed::typeof(1.0s)
    seed::UInt64
    # Typed on the *supertype*, never a concrete `DateTime`: that is what lets a 360-day or
    # noleap calendar be substituted later without changing a signature. Abstract, but never read in
    # the hot path - the epoch is resolved once at build and thereafter only reported.
    epoch::Union{Nothing, Dates.TimeType}
    # How dates map to elapsed time - set by `build_ecosystem`, and like the epoch never read in the
    # hot path.
    calendar::AbstractRunCalendar
    inputs::Vector{ProvenanceRecord}
    # Abstract on purpose: a placement may be any function, and nothing in the step reads it.
    placement::Any
    savedchanges::Union{Nothing, _SavedChanges}
    rulechanges::Vector{_RuleChange}

    function Ecosystem{Part, SL, NF}(abundances::GridLandscape,
                                     spplist::SL,
                                     habitat::Part,
                                     nichefit::NF,
                                     lookup::Vector{Lookup},
                                     cache::Cache,
                                     rngs::Vector{Random.Xoshiro},
                                     elapsed::typeof(1.0s),
                                     seed::UInt64,
                                     epoch::Union{Nothing, Dates.TimeType} = nothing) where {Part <:
                                                                                             AbstractHabitat,
                                                                                             SL <:
                                                                                             SpeciesList,
                                                                                             NF <:
                                                                                             AbstractNicheFit}
        _checksimulatable(habitat)
        # Tolerances, regimes and the nichefit between them must line up member for member - same
        # arity, same names in the same order, and each member measuring the same thing in the same
        # unit and of the same continuous/categorical kind.
        _checkaligned("species tolerances, environment regime and trait nichefit",
                      _toleranceside(spplist.tolerance),
                      _regimeside(habitat.regime), _nichefitside(nichefit))
        # The species' demands must likewise line up with the resources on offer.
        _checkaligned("species demand and environment supply",
                      _demandside(spplist.demand),
                      _supplyside(habitat.supply))
        _mcmatch(abundances.matrix, spplist, habitat) ||
            error("Dimension mismatch: the abundance matrix " *
                  "($(size(abundances.matrix, 1)) × $(size(abundances.matrix, 2))) " *
                  "does not match the $(_counttypes(spplist, true)) species and " *
                  "$(_countsubcommunities(habitat)) subcommunities of the species " *
                  "list and environment.")
        return new{Part, SL, NF}(abundances,
                                 spplist,
                                 habitat,
                                 nichefit,
                                 lookup,
                                 cache,
                                 rngs,
                                 elapsed,
                                 seed,
                                 epoch,
                                 ExactDates(),
                                 ProvenanceRecord[],
                                 nothing,
                                 nothing,
                                 _RuleChange[])
    end
end

"""
    CachedEcosystem{Part <: AbstractHabitat, SL <: SpeciesList,
        NF <: AbstractNicheFit} <: AbstractEcosystem{Part, SL, NF}

CachedEcosystem houses the same information as [`Ecosystem`](@ref) (see
?Ecosystem), but holds the time period abundances as a
[`CachedGridLandscape`](@ref), so that they may be present or missing.
"""
mutable struct CachedEcosystem{Part <: AbstractHabitat,
                               SL <: SpeciesList,
                               NF <: AbstractNicheFit} <:
               AbstractEcosystem{Part, SL, NF}
    abundances::CachedGridLandscape
    spplist::SL
    habitat::Part
    nichefit::NF
    lookup::Vector{Lookup}
    cache::Cache
    rngs::Vector{Random.Xoshiro}
    elapsed::typeof(1.0s)
    seed::UInt64
    epoch::Union{Nothing, Dates.TimeType}
    calendar::AbstractRunCalendar
    inputs::Vector{ProvenanceRecord}
    placement::Any
end

"""
    CachedEcosystem(eco::Ecosystem, outputfile::String, times::StepRangeLen;
                    saveinterval::Unitful.Time = step(times))

Create a CachedEcosystem given an existing [`Ecosystem`](@ref), `eco`, an output
folder to which the simulations are saved, `outputfile`, and a range of times
over which to simulate, `times`. The step of `times` is the simulation timestep;
`saveinterval` controls how often checkpoints are written to disk (a multiple of
the timestep, defaulting to every step). Because the simulation always advances
by the timestep, results are independent of `saveinterval`.
"""
function CachedEcosystem(eco::Ecosystem, outputfile::String,
                         times::StepRangeLen;
                         saveinterval::Unitful.Time = step(times))
    # Nothing here checks slice counts against the timestep. A layer holds one slice and its change
    # decides which, from elapsed time, so a series need not line up with the timestep at all;
    # running past the end of one is the series' own business, and it cycles, holds, or says so at
    # the step it happens - which is more precise than a length comparison made before the run.
    abundances = CachedGridLandscape(outputfile, times,
                                     saveinterval = saveinterval)
    abundances.matrix[1] = eco.abundances
    return CachedEcosystem{typeof(eco.habitat), typeof(eco.spplist),
                           typeof(eco.nichefit)}(abundances,
                                                 eco.spplist,
                                                 eco.habitat,
                                                 eco.nichefit,
                                                 eco.lookup,
                                                 eco.cache,
                                                 eco.rngs,
                                                 eco.elapsed,
                                                 eco.seed,
                                                 eco.epoch,
                                                 eco.calendar,
                                                 eco.inputs,
                                                 eco.placement)
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# Without these an ecosystem falls through to Diversity's `AbstractMetacommunity` generic, which
# prints its full type signature, too long to read.
#
# **Both sides of the model are named, not just the species count.** An ecosystem is species in an
# environment, so what a reader needs is what the environment offers and what the species are matched
# against it on: the regime axes it experiences and the supply axes it competes for. `_axisnames`
# reads them off the layer types, so nothing here touches a value.

function Base.show(io::IO, eco::Ecosystem)
    nsp, ny, nx = size(eco.abundances.grid)
    return print(io,
                 "Ecosystem($(nsp) species, $(ny) × $(nx), ",
                 _axisnames(getregime(eco.habitat)), " | ",
                 _axisnames(getsupply(eco.habitat)),
                 ", t = $(eco.elapsed))")
end

# The cached ecosystem mirrors the plain one, except that its abundances are a series rather than a
# single landscape -- so the grid shape comes from the habitat, which always has it, rather than
# from a landscape that may not currently be held.
function Base.show(io::IO, cache::CachedEcosystem)
    ny, nx = getgridshape(cache.habitat)
    return print(io,
                 "CachedEcosystem($(counttypes(cache.spplist, true)) species, $(ny) × $(nx), ",
                 _axisnames(getregime(cache.habitat)), " | ",
                 _axisnames(getsupply(cache.habitat)),
                 ", $(count(!ismissing, cache.abundances.matrix))/",
                 "$(length(cache.abundances.matrix)) held)")
end

function Base.show(io::IO, ::MIME"text/plain", cache::CachedEcosystem)
    ny, nx = getgridshape(cache.habitat)
    println(io, "CachedEcosystem")
    println(io, "  species   ", counttypes(cache.spplist, true))
    println(io, "  grid      $(ny) × $(nx) cells, ",
            count(cache.habitat.active), " active")
    println(io, "  regime    ", sprint(show, getregime(cache.habitat)))
    println(io, "  supply    ", sprint(show, getsupply(cache.habitat)))
    println(io, "  cached    ", sprint(show, cache.abundances))
    print(io, "  clock     t = $(cache.elapsed), seed = $(cache.seed)")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", eco::Ecosystem)
    nsp, ny, nx = size(eco.abundances.grid)
    println(io, "Ecosystem")
    println(io, "  species   $(nsp), $(sum(eco.abundances.matrix)) individuals")
    println(io, "  grid      $(ny) × $(nx) cells, ",
            count(eco.habitat.active), " active")
    println(io, "  regime    ", sprint(show, getregime(eco.habitat)))
    println(io, "  supply    ", sprint(show, getsupply(eco.habitat)))
    println(io, "  tolerance ", sprint(show, gettolerance(eco.spplist)))
    println(io, "  demand    ", sprint(show, getdemand(eco.spplist)))
    println(io, "  nichefit  ", nameof(typeof(eco.nichefit)))
    print(io, "  clock     t = $(eco.elapsed), seed = $(eco.seed)")
    return nothing
end

# == Functions ==================================================================================

"""
    abundances(cache::CachedEcosystem, tm::Unitful.Time)

Extract abundances for an ecosystem, `cache`, at a certain point in time, `tm`.
If the abundances for that time are missing from the ecosystem, then the
function checks on disk for the last saved version and simulates forward, saving as it goes - so a
cache built again over the same folder, by the same script run again, goes on from the last state
saved there.

# Arguments

  - `cache`: the cached ecosystem.
  - `tm`: the time to give its abundances at, one of the times it was built with.
"""
function abundances(cache::CachedEcosystem, tm::Unitful.Time)
    return _abundances(cache, tm)[2]
end

"""
    resettime!(eco::AbstractEcosystem)

Reset the simulation clock to zero, so that the next timestep is the first.

Layers that vary in time are indexed by elapsed time, so this is what puts a stored series back to
its first slice.
"""
function resettime!(eco::AbstractEcosystem)
    return eco.elapsed = 0.0s
end

"""
    getnichefit(eco::Ecosystem)

Extract niche fits.
"""
function getnichefit(eco::AbstractEcosystem)
    return eco.nichefit
end

"""
    clearcache!(cache::CachedEcosystem)

Delete all JLD2 cache files from the output folder of a
[`CachedEcosystem`](@ref). Returns a string reporting how many files were
removed.

The trailing `!` is not decoration: this **destroys** the recorded abundances on disk, and a
`CachedEcosystem` that has been cleared can no longer answer for the timepoints it held.

# Arguments

  - `cache`: the [`CachedEcosystem`](@ref) whose output folder to empty.
"""
function clearcache!(cache::CachedEcosystem)
    files = _searchdir(cache.abundances.outputfolder, ".jld2")
    rm.(joinpath.(cache.abundances.outputfolder, files))
    len = length(files)
    return "$len files cleared"
end

# **The four-way accessor family, and it is deliberately symmetric.** Each of `getregime`,
# `getsupply`, `gettolerance` and `getdemand` takes **one** argument and returns **the thing itself**,
# never a number derived from it - and each is defined twice: on an `AbstractEcosystem` (so it covers
# `Ecosystem` *and* `MPIEcosystem`) and on whichever half actually owns it, `GridHabitat` for the two
# environment sides and `SpeciesList` for the two species sides.
# **Totals are a separate family** - `totalsupply` - because a total is a different question from
# the object. Returning a total would make `getdemand` the odd
# one out; it now returns the demand, matching its three siblings.

"""
    getregime(eco::Ecosystem)
    getregime(habitat::GridHabitat)

Return the regime - what each cell *is like* - as the [`AbstractRegime`](@ref) layer (or collection
of them) the environment holds.

#### Arguments:

  - `eco`/`habitat`: the ecosystem to read the regime from, or the habitat directly.

#### Returns:

  - the regime layer or collection. Not its values: reach those through `.matrix`.
"""
function getregime(eco::AbstractEcosystem)
    return getregime(eco.habitat)
end

getregime(habitat::GridHabitat) = habitat.regime

"""
    getsupply(eco::Ecosystem)
    getsupply(habitat::GridHabitat)

Return the supply - what each cell *provides* - as the [`AbstractSupply`](@ref) layer (or collection
of them) the environment holds. The mirror of [`getregime`](@ref).

Returns the layer, not its values; a caller who wants the numbers writes
`getsupply(eco).matrix`. A layer's array type is an implementation detail, and the package does not
put one across its public boundary.

#### Arguments:

  - `eco`/`habitat`: the ecosystem to read the supply from, or the habitat directly.

#### Returns:

  - the supply layer or collection. For the *total* resource instead, see [`totalsupply`](@ref).
"""
function getsupply(eco::AbstractEcosystem)
    return getsupply(eco.habitat)
end

getsupply(habitat::GridHabitat) = habitat.supply

"""
    gettolerance(eco::Ecosystem)
    gettolerance(sppl::SpeciesList)

Return the species' tolerances - their [`Condition`](@ref)-role requirements, pairing member for
member with the regime. The species-side mirror of [`getregime`](@ref).

#### Arguments:

  - `eco`/`sppl`: the ecosystem to read the tolerances from, or the species list directly.

#### Returns:

  - the tolerance or collection of them.
"""
function gettolerance(eco::AbstractEcosystem)
    return gettolerance(eco.spplist)
end

gettolerance(sppl::SpeciesList) = sppl.tolerance

"""
    getdemand(eco::Ecosystem)
    getdemand(sppl::SpeciesList)

Return the species' demands - their [`Resource`](@ref)-role requirements, pairing member for member
with the supply. The species-side mirror of [`getsupply`](@ref).

Returns the demand itself rather than a total, which would make it the only member of this family
answering with a number rather than with the thing asked for.

#### Arguments:

  - `eco`/`sppl`: the ecosystem to read the demands from, or the species list directly.

#### Returns:

  - the demand or collection of them.
"""
function getdemand(eco::AbstractEcosystem)
    return getdemand(eco.spplist)
end

# --- The accessor 2x2: what a CELL offers against what a SPECIES brings ------
# Beside the four-way `get*` family above, which these complete: `getregime`/`getsupply` answer for
# the whole grid, `cellregime`/`cellsupply` for one cell, and the species pair mirrors both.

# ---------------------------------------------------------------------------
# The public four, and the 2×2 they complete
# ---------------------------------------------------------------------------
# Each returns a `NamedTuple` keyed by the corresponding `*_names` accessor, so a single layer
# gives a one-entry tuple and a collection one entry per member - the same normalisation
# `totalsupply` uses, and the reason one shape reads the same however many layers there are.
#
# Each takes the **assembled ecosystem** or the **object one of the three builders returns** - a
# `GridHabitat` from `GridHabitat`, a `SpeciesList` from `build_species`. Both exist before any
# ecosystem does, which is the workflow the builders encourage: build it, look at it, then assemble.

"""
    cellregime(eco::AbstractEcosystem, y::Int64, x::Int64)
    cellregime(habitat::AbstractHabitat, y::Int64, x::Int64)

Return what the environment **is like** at cell `(y, x)` - a `NamedTuple` keyed by
`keys`, with one entry per regime layer.

The **Condition** half of the cell side of the 2×2 this vocabulary is built on; see
[`speciestolerance`](@ref) for what a species brings to it.
"""
function cellregime(habitat::AbstractHabitat, y::Int64, x::Int64)
    return map(l -> _cellvalue(l, y, x), NamedTuple(habitat.regime))
end

function cellregime(eco::AbstractEcosystem, y::Int64, x::Int64)
    return cellregime(eco.habitat, y, x)
end

"""
    cellsupply(eco::AbstractEcosystem, y::Int64, x::Int64)
    cellsupply(habitat::AbstractHabitat, y::Int64, x::Int64)

Return what the environment **provides** at cell `(y, x)` - a `NamedTuple` keyed by
`keys`, with one entry per supply layer.

The **Resource** half of the cell side; see [`speciesdemand`](@ref) for what a species asks of it.
These are the supply's own values, not what is *available*: an inactive cell's resource is
included here, whereas [`totalsupply`](@ref) totals only the active ones.
"""
function cellsupply(habitat::AbstractHabitat, y::Int64, x::Int64)
    return map(l -> _cellvalue(l, y, x), NamedTuple(habitat.supply))
end

function cellsupply(eco::AbstractEcosystem, y::Int64, x::Int64)
    return cellsupply(eco.habitat, y, x)
end

"""
    speciestolerance(eco::AbstractEcosystem, sp::Int64)
    speciestolerance(spplist::SpeciesList, sp::Int64)

Return what species `sp` **can tolerate** - a `NamedTuple` keyed by `tolerance_names`, with one entry
per tolerance. Each value is whatever that tolerance kind holds: a response *distribution* for a
[`NicheTolerance`](@ref), a value for a categorical one, a vector of codes for land cover.

The **Condition** half of the species side; pairs with [`cellregime`](@ref).
"""
function speciestolerance(spplist::SpeciesList, sp::Int64)
    return map(t -> _speciestolerance(t, sp),
               NamedTuple(spplist.tolerance))
end

function speciestolerance(eco::AbstractEcosystem, sp::Int64)
    return speciestolerance(eco.spplist, sp)
end

"""
    speciesdemand(eco::AbstractEcosystem, sp::Int64)
    speciesdemand(spplist::SpeciesList, sp::Int64)

Return what species `sp` **requires** - a `NamedTuple` keyed by `demand_names`, with one entry per
resource: `(sunlight = 2.0kJ/day, water = 5.0L/day)`.

The **Resource** half of the species side; pairs with [`cellsupply`](@ref), and the values are
directly comparable with it.

These are the **raw unitful rates**, never `resource[sp] * exchange_rate`. The exchange rate is
arithmetic machinery - it defaults to `1/mean(resource)`, so the product is a dimensionless number on
a *community-relative* scale that shifts when a species is added or removed. `2.0 kJ/day` is what a
species actually needs; `0.5` is only "half the mean of whoever happens to be here".
"""
function speciesdemand(spplist::SpeciesList, sp::Int64)
    return map(d -> _speciesdemand(d, sp), NamedTuple(spplist.demand))
end

function speciesdemand(eco::AbstractEcosystem, sp::Int64)
    return speciesdemand(eco.spplist, sp)
end

"""
    speciesdispersal(eco::AbstractEcosystem, sp)
    speciesdispersal(spplist::SpeciesList, sp)

Return species `sp`'s dispersal **kernel**, an [`AbstractKernel`](@ref), or `nothing` for a species
that does not disperse. `sp` may be an
index or a name.

**The exception to the `NamedTuple` shape** of the other four: a species has one dispersal kernel,
not one per layer or per resource, so there is nothing to key a tuple by.

It returns the kernel itself rather than a parameterisation, which is both more honest - the
distance distribution is Rayleigh rather than Normal, and a `LongTailKernel` is a 2Dt with a `shape`
a Normal cannot express - and directly reusable: a dispersal mode called with it gives a species'
movement, so another species can be given this one's dispersal.

```julia
InsectDispersal()(speciesdispersal(eco, 3))
```
"""
function speciesdispersal(spplist::SpeciesList, sp::Int64)
    return spplist.movement[sp].kernel
end

function speciesdispersal(spplist::SpeciesList, sp::AbstractString)
    return speciesdispersal(spplist, _speciesindex(spplist, sp))
end

function speciesdispersal(eco::AbstractEcosystem, sp)
    return speciesdispersal(eco.spplist,
                            sp)
end

"""
    simulationtime(eco::AbstractEcosystem)

Return the simulation time elapsed since the ecosystem was built - zero at
construction, advanced by one `timestep` per [`update!`](@ref).
"""
function simulationtime(eco::AbstractEcosystem)
    return eco.elapsed
end

"""
    simulationdate(eco::AbstractEcosystem)

Return the real date the simulation has reached, or `nothing` for a run with no
epoch.

The epoch is the date [`simulationtime`](@ref) counts from - resolved at
[`build_ecosystem`](@ref) from the environment's dated series, or given there
explicitly. A run whose environment mentions no dates has no epoch, and so no
date: elapsed time is all there is to say about when it is. The date is reached
through the run's calendar, so under [`MeanMonths`](@ref) every month is one
`month_mean_duration` long.
"""
function simulationdate(eco::AbstractEcosystem)
    return _shiftdate(eco.epoch, simulationtime(eco), eco.calendar)
end

"""
    makeunique(eco::Ecosystem)
Convert type of similarity in [`SpeciesList`](@ref) to UniqueTypes, i.e. an
identity matrix.

The species keep their names, their identity and values from the pool they were built from, and
their name authority: only the similarity structure is dropped.
"""
function makeunique(eco::Ecosystem)
    sppl = eco.spplist
    EcoSISTEM.invalidatecaches!(eco)
    # Every field but `types` must reach the new list: a pool-built list's identity (`traits`,
    # `namemap`) and its `authority` cannot be rebuilt from the rest.
    newsppl = SpeciesList{typeof(sppl.tolerance),
                          typeof(sppl.demand),
                          typeof(sppl.movement),
                          UniqueTypes,
                          typeof(sppl.params)}(sppl.names,
                                               sppl.tolerance,
                                               sppl.demand,
                                               _uniquetypes(sppl.names),
                                               sppl.movement,
                                               sppl.params,
                                               provenance = sppl.inputs,
                                               traits = sppl.traits,
                                               filled = sppl.filled,
                                               authority = sppl.authority,
                                               namemap = sppl.namemap)
    newsppl.susceptible = sppl.susceptible
    stripped = Ecosystem{typeof(eco.habitat), typeof(newsppl),
                         typeof(eco.nichefit)}(eco.abundances,
                                               newsppl,
                                               eco.habitat,
                                               eco.nichefit,
                                               eco.lookup,
                                               eco.cache,
                                               eco.rngs,
                                               eco.elapsed,
                                               eco.seed,
                                               eco.epoch)
    stripped.calendar = eco.calendar
    stripped.inputs = eco.inputs
    stripped.placement = eco.placement
    stripped.savedchanges = eco.savedchanges
    stripped.rulechanges = eco.rulechanges
    return stripped
end

# The axes a layer or collection is on, as a `+`-joined string - `Temperature + Precipitation`.
function _axisnames(l::AbstractLayer)
    return l isa LayerCollection ?
           join(map(m -> nameof(axisof(m)), values(l)), " + ") :
           string(nameof(axisof(l)))
end

@recipe function f(::AbstractMovement, eco::AbstractEcosystem, sp::Int64)
    l = eco.lookup[sp]
    maxX = maximum(l.x)
    maxY = maximum(l.y)
    x, y = round(Int64, maxX / 2), round(Int64, maxY / 2)
    # Can't go over maximum dimension
    valid = findall((l.x .> -x) .& (l.y .> -y) .& (l.x .<= (maxX - x)) .&
                    (l.y .<= (maxY - y)))
    probs = l.p[valid]
    probs ./= sum(probs)
    xs = (l.x[valid] .+ x)
    ys = (l.y[valid] .+ y)
    A = zeros(maxX, maxY)
    # All three are sliced by the same `valid`, so they walk together and no index is needed.
    for (x, y, prob) in zip(xs, ys, probs)
        A[x, y] .= prob
    end
    seriestype := :heatmap
    grid --> false
    aspect_ratio --> 1
    title --> "Movement kernel (km)"
    ax = _plotaxes(getregime(eco))
    return ax.x, ax.y, A
end

"""
    makerngs(seed::Integer, n::Integer)

Build a vector of `n` independent, deterministically-seeded random number
generators, one per species. Species `j` is seeded from the pair `(seed, j)`
through the generator's own seeding, `Xoshiro(UInt64[seed, j])`, so its random
stream is a pure function of `(seed, j)` - independent of how species are
distributed across threads or MPI processes, and the same on every Julia version,
since `Base.hash`, whose values change between versions, is never involved. This
is what makes simulation results reproducible across different thread and process
counts (each species is always processed by exactly one task on one rank, drawing
in a fixed cell order) and across Julia releases. See [`getrng`](@ref).

Note: this per-species scheme is sufficient only because no single species' draws
are ever split across ranks/tasks. If a species' cells were ever partitioned
across ranks, a per-`(species, cell)` counter-based generator would be needed
instead.
"""
function makerngs(seed::Integer, n::Integer)
    return [Random.Xoshiro(UInt64[seed % UInt64, j]) for j in 1:n]
end

# An ecosystem of `spplist` in `habitat` with no individuals yet, its streams made from `seed`: what
# `build_ecosystem` places into.
function _emptyecosystem(spplist::SpeciesList{TL, DM}, habitat::GridHabitat,
                         nichefit::AbstractNicheFit;
                         seed::Integer = rand(UInt64)) where {TL, DM}
    # **First, and before `genlookups` below.** The inner constructor checks this too, but it runs
    # last, and `genlookups` divides the regime's cell size by a dispersal distance, so on a
    # geographic grid (`size` in `°`) it would raise `DimensionError: ° km^-1` before the tailored
    # message is reached.
    _checksimulatable(habitat)
    ml = empty_landscape(habitat, spplist)
    # One deterministically-seeded RNG per species, so births, deaths, dispersal and the placement
    # are reproducible across thread counts.
    rngs = makerngs(seed, size(ml.matrix, 1))
    lookup_tab = collect(map(k -> genlookups(habitat.regime, k),
                             getkernels(spplist.movement)))
    nm = zeros(Int64, size(ml.matrix))
    totaldemand = zeros(Float64, (size(ml.matrix, 2), numdemands(DM)))
    return Ecosystem{typeof(habitat), typeof(spplist), typeof(nichefit)}(ml,
                                                                         spplist,
                                                                         habitat,
                                                                         nichefit,
                                                                         lookup_tab,
                                                                         Cache(nm,
                                                                               totaldemand,
                                                                               false),
                                                                         rngs,
                                                                         0.0s,
                                                                         _storedseed(seed))
end

# Narrow a user-supplied seed of any integer type to the `UInt64` the ecosystem stores. `%` rather
# than `UInt64(...)`, so a negative seed reinterprets rather than throwing - `seed` only ever seeds
# a stream, it is never used as a magnitude. `makerngs` narrows the same way, so the stored seed
# and the one the species streams derive from are one value.
_storedseed(seed::Integer) = seed % UInt64

# Check the abundance matrix `m` (species × subcommunities) has one row per species in `sppl` and one
# column per subcommunity in `part`. A pure dimension check.
function _mcmatch(m::AbstractMatrix, sppl::SpeciesList, part::AbstractHabitat)
    return _counttypes(sppl, true) == size(m, 1) &&
           _countsubcommunities(part) == size(m, 2)
end

# The abundances at `tm`, simulating forward from the last checkpoint if that time was never saved.
# The gap exists because checkpoints land only on multiples of `saveinterval` while the run advances
# by `timestep`, so an arbitrary time usually falls between two of them.
function _abundances(cache::CachedEcosystem, tm::Unitful.Time)
    timestep = cache.abundances.timestep
    saveinterval = cache.abundances.saveinterval
    # Checkpoints are written to disk only at multiples of `saveinterval`,
    # indexed by how many save intervals have elapsed; the simulation itself
    # always advances by `timestep`, so results are independent of `saveinterval`.
    issave = iszero(mod(tm, saveinterval))
    idx = issave ? Int(round(uconvert(NoUnits, tm / saveinterval))) : missing
    if ismissing(cache.abundances.matrix[Ti(At(tm))])
        if checkfile(cache.abundances.outputfolder, idx)
            cache.abundances.matrix[Ti(At(tm))] = loadfile(cache,
                                                           cache.abundances.outputfolder,
                                                           idx,
                                                           cache.spplist.names,
                                                           getcoords(cache.habitat))
            return tm, cache.abundances.matrix[Ti(At(tm))]
        else
            newtm, abun = _abundances(cache, tm - timestep)
            if (newtm > 2 * timestep)
                cache.abundances.matrix[Ti(At(newtm -
                                              2 *
                                              timestep))] = missing
            end
        end
    else
        return tm, cache.abundances.matrix[Ti(At(tm))]
    end
    simulate!(cache, newtm, timestep)
    if issave
        _writestate(joinpath(cache.abundances.outputfolder,
                             string(idx, ".jld2")),
                    cache.abundances.matrix[Ti(At(tm))].matrix,
                    _streamwords(cache.rngs), cache.spplist, cache.habitat,
                    tm, cache.seed, cache.epoch, cache.calendar)
    end
    return _abundances(cache, newtm + timestep)
end

"""
    getlookup(eco::Ecosystem)

Extract movement lookup table of species from [`Ecosystem`](@ref) object.
"""
function getlookup(eco::AbstractEcosystem, sp::Int64)
    return eco.lookup[sp]
end

function getlookup(eco::AbstractEcosystem, sp::String)
    num = findall(eco.spplist.names .== sp)[1]
    return getlookup(eco, num)
end

@doc (@doc getlookup) getlookup(::AbstractEcosystem, ::String)

"""
    getrng(eco::AbstractEcosystem, sp::Int64)

Return the per-species random number generator for global species index `sp`.
Because each species has its own stream and is processed by exactly one task per
timestep, random draws are both thread-safe and reproducible independent of the
number of threads or MPI processes (see [`makerngs`](@ref)).
"""
function getrng(eco::AbstractEcosystem, sp::Int64)
    return eco.rngs[sp]
end

# The landing density per unit area of the long-tailed kernel, a 2-D Student-t, at a displacement
# `(u, v)` in units of its scale. Its tail falls off as a power of distance rather than exponentially,
# so a rare long jump stays possible. `b` is its shape.
function _2Dtdensity(u, v, b)
    return ((b - 1) / π) * (1 + u^2 + v^2)^-b
end

# The landing density per unit area, at a displacement `(u, v)`, of a kernel stated as the
# distribution of the distance travelled: the distance density spread round the circle at that
# distance, `pdf(dist, d) / (2pi d)`. The point `d = 0` has measure zero and is given zero.
function _distancedensity(u, v, dist::ContinuousUnivariateDistribution)
    d = sqrt(u^2 + v^2)
    return d > 0 ? pdf(dist, d) / (2π * d) : 0.0
end

# The landing density per unit area of the Gaussian kernel at a displacement `(u, v)` in units of its
# scale. Its lookup uses the closed form `_gaussianshare`; this is what that form integrates.
_gaussiandensity(u, v) = exp(-(u^2 + v^2)) / π

# The share of a Gaussian kernel's seeds that move `offset` cells along one axis from a cell `R`
# scales wide, the source position uniform in the cell. The kernel's two axes are independent, so a
# cell `(dx, dy)` away gets `_gaussianshare(dx, R) * _gaussianshare(dy, R)`. The exact form is a
# second difference of `|u| + _gaussiantail(u)`, the antiderivative of `erf`. The `|u|` terms cancel
# exactly, so they are left out, and nothing large cancels for a far cell.
function _gaussianshare(offset::Integer, R::Float64)
    offset == 0 && return 1 + (_gaussiantail(R) - _gaussiantail(0.0)) / R
    offset = abs(offset)
    return (_gaussiantail((offset + 1) * R) - 2 * _gaussiantail(offset * R) +
            _gaussiantail((offset - 1) * R)) / (2R)
end

# The part of `u erf(u) + exp(-u^2) / sqrt(pi)` that is not `|u|`, written with `erfc` so that it stays
# accurate where it is small.
_gaussiantail(u::Float64) = exp(-u^2) / sqrt(π) - abs(u) * erfc(abs(u))

# The share of a kernel's seeds that land in the cell `(dx, dy)` cells away, for a cell `R` scales
# wide and a landing density `density(u, v)` per unit area in the same units. The source position is
# uniform in its cell. A displacement `u` along one axis then lands `offset` cells over with the
# weight of a triangle, rising from zero at `(offset - 1)R` to one at `offset * R` and back to zero
# at `(offset + 1)R`. So the share is a 2-D integral of the density times one triangle per axis.
# Most of the density lies within a few scales of the origin, which for a wide cell is a thin strip
# of the triangles' range. So each axis is split where its triangle kinks and at powers of two
# either side of the origin, and each piece is integrated on its own.
function _cellshare(density, R::Float64, dx::Integer, dy::Integer)
    us, vs = _splitpoints(dx, R), _splitpoints(dy, R)
    weighted(r) = _triangle(r[1], dx, R) * _triangle(r[2], dy, R) *
                  density(r[1], r[2])
    total = 0.0
    for (u0, u1) in zip(us[1:(end - 1)], us[2:end]),
        (v0, v1) in zip(vs[1:(end - 1)], vs[2:end])
        total += hcubature(weighted, [u0, v0], [u1, v1], rtol = 1.0e-8,
                           atol = 1.0e-15)[1]
    end
    return total
end

# The weight with which a displacement `u` lands `offset` cells over, from a cell `R` wide.
_triangle(u, offset, R) = max(0.0, 1 - abs(u / R - offset))

# Where one axis's range is split for `_cellshare`: the ends and peak of its triangle, the origin, and
# the powers of two from one eighth of a scale outwards on each side of the origin, wherever they fall
# inside the range.
function _splitpoints(offset::Integer, R::Float64)
    lo, hi = (offset - 1) * R, (offset + 1) * R
    points = [lo, offset * R, hi]
    lo < 0.0 < hi && push!(points, 0.0)
    step = 0.125
    while step < max(abs(lo), abs(hi))
        lo < step < hi && push!(points, step)
        lo < -step < hi && push!(points, -step)
        step *= 2
    end
    return sort!(unique!(points))
end

"""
    genlookups(regime::AbstractRegime, kernel::AbstractKernel)

Generate lookup tables, which hold information on the probability of moving to
neighbouring squares.
"""
function genlookups(regime::AbstractRegime, kernel::GaussianKernel)
    R = uconvert(NoUnits, _uniformcellside(regime) / kernel.scale)
    return Lookup(maximum(getgridshape(regime)), kernel.thresh) do dx, dy
        return _gaussianshare(dx, R) * _gaussianshare(dy, R)
    end
end

function genlookups(regime::AbstractRegime, kernel::LongTailKernel)
    R = uconvert(NoUnits, _uniformcellside(regime) / kernel.scale)
    density(u, v) = _2Dtdensity(u, v, kernel.shape)
    return Lookup(maximum(getgridshape(regime)), kernel.thresh) do dx, dy
        return _cellshare(density, R, dx, dy)
    end
end

# Distances in units of the mean, so the inverse Gaussian is `InverseGaussian(1, scale / meandist)`.
function genlookups(regime::AbstractRegime, kernel::WaldKernel)
    R = uconvert(NoUnits, _uniformcellside(regime) / kernel.meandist)
    distance = InverseGaussian(1.0,
                               uconvert(NoUnits,
                                        kernel.scale / kernel.meandist))
    density(u, v) = _distancedensity(u, v, distance)
    return Lookup(maximum(getgridshape(regime)), kernel.thresh) do dx, dy
        return _cellshare(density, R, dx, dy)
    end
end

# A species that does not disperse has no neighbourhood, and its movement never reads one.
genlookups(::AbstractRegime, ::Nothing) = Lookup()

# Advance a date by an elapsed `Unitful.Time`. Milliseconds because that is the finest unit a
# `Dates.TimeType` generally shares, and rounded because elapsed time is a float while a date is not.
_shiftdate(::Nothing, _) = nothing

function _shiftdate(epoch::Dates.TimeType, elapsed::Unitful.Time)
    return epoch +
           Dates.Millisecond(round(Int, ustrip(uconvert(Unitful.ms, elapsed))))
end

# A `Date` has day resolution and rejects `+ Millisecond` outright, so a run given a `Date` epoch
# would throw here rather than at the line that set it. Promoted instead: a date plus an elapsed time
# is a *moment*, so the answer is a `DateTime` at midnight on that date plus however long has passed.
# Only `Date` is promoted - a calendar type from elsewhere (a 360-day or noleap one) must keep its
# own arithmetic, which is why the epoch is typed on the supertype in the first place.
function _shiftdate(epoch::Dates.Date, elapsed::Unitful.Time)
    return _shiftdate(Dates.DateTime(epoch), elapsed)
end

# The same under a run calendar: real elapsed time under `ExactDates`, and under `MeanMonths` the
# epoch's position in mean months moved on by elapsed time counted in them.
_shiftdate(epoch, elapsed, ::ExactDates) = _shiftdate(epoch, elapsed)

_shiftdate(::Nothing, _, ::MeanMonths) = nothing

function _shiftdate(epoch::Dates.TimeType, elapsed::Unitful.Time, ::MeanMonths)
    return _datefrommeanmonths(_meanmonthposition(epoch) +
                               ustrip(NoUnits, elapsed / month_mean_duration))
end

# The elapsed time a date names in a run from `epoch`, under the run calendar - the inverse of
# `_shiftdate`.
function _elapsedat(epoch::Dates.TimeType, date::Dates.TimeType, ::ExactDates)
    gap = Dates.DateTime(date) - Dates.DateTime(epoch)
    return uconvert(s, float(Dates.value(Dates.Millisecond(gap)) * Unitful.ms))
end

function _elapsedat(epoch::Dates.TimeType, date::Dates.TimeType, ::MeanMonths)
    months = _meanmonthposition(date) - _meanmonthposition(epoch)
    return uconvert(s, months * month_mean_duration)
end

# Advance the simulation clock by one timestep. Called at the end of `update!` (after the
# population dynamics, before the layer update), so a layer change sees the time it is changing
# *to*. Any `Unitful.Time` step is converted to the stored seconds.
function _advanceclock!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return eco.elapsed += uconvert(s, float(timestep))
end

"""
    checkfile(::String, ::Missing)

Check whether a cache file exists for a given timepoint. Always returns `false`
when the timepoint is `missing`.
"""
function checkfile(::String, ::Missing)
    return false
end

"""
    checkfile(file::String, idx::Int)

Check whether a JLD2 checkpoint file exists in the folder `file` for checkpoint
index `idx`. Returns `true` if `<idx>.jld2` is present.
"""
function checkfile(file::String, idx::Int)
    return isfile(joinpath(file, string(idx, ".jld2")))
end

"""
    loadfile(cache::CachedEcosystem, file::String, idx::Int, names::Vector{String},
             yx::Tuple{<:Y, <:X})

Load a cached [`GridLandscape`](@ref) from the folder `file` for checkpoint index
`idx`, a state in [`write_state`](@ref)'s layout, putting its abundances onto `names` and `grid` (the
ecosystem's own species names and grid). The checkpoint is checked against the cache's grid, species,
epoch and calendar, so a folder written by another configuration is refused rather than continued
from; and as a run started from a state does, the cache takes the saved species' activity, the
active cells, every layer's values and each species' random stream, so the resumed run goes on
exactly as the saved one would have.
"""
function loadfile(cache::CachedEcosystem, file::String, idx::Int,
                  names::Vector{String}, grid::StudyGrid)
    path = joinpath(file, string(idx, ".jld2"))
    state = _readstate(path)
    _checkstategrid(state, cache.habitat, path)
    columns = _statecolumns(state, cache.spplist, path)
    _stateepoch(state, cache.epoch, path)
    _statecalendar(state, cache.calendar, path)
    _checkstreamfields(state, path)
    _restoreground!(cache.spplist, cache.habitat, state, columns, path)
    _setstreams!(cache.rngs, state.streams, columns)
    saved = zeros(Int64, length(names), size(state.abundances, 2))
    for (row, column) in enumerate(columns)
        iszero(column) || (saved[row, :] .= view(state.abundances, column, :))
    end
    # The saved state the run resumes from is among what it was built from.
    record = ProvenanceRecord(role = :state, source = "saved state",
                              path = string(idx, ".jld2"))
    record in cache.inputs || push!(cache.inputs, record)
    return GridLandscape(saved, names, grid)
end

# **Ordinariness is not stored, so there is none to invalidate.** It is recomputed on demand into
# `cache.relativeabun`, so no stale value can survive a change to the ecosystem; a stored copy would
# have to be invalidated here on every change.
function invalidatecaches!(eco::AbstractEcosystem)
    eco.cache.netmigration .= 0
    return eco.cache.valid = false
end

# Write the current abundances into an abundance recording. Indexed by the species there are rather
# than by the whole first dimension, so storage of the caller's own with rows to spare can be given;
# rows beyond the species are left as they were.
function _record!(storage::AbstractArray, eco::AbstractEcosystem,
                  counting::Integer)
    nspp = size(eco.abundances.matrix, 1)
    nspp <= size(storage, 1) ||
        error("this run has $nspp species but the recording was sized for $(size(storage, 1)). " *
              "`generate_storage` sizes it from the ecosystem.")
    storage[1:nspp, :, counting] = eco.abundances.matrix
    return storage
end
