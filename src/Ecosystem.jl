# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The simulation state itself — species, environment and abundances assembled — and the scratch
# space the hot loop writes into rather than allocating.

using Diversity

using Unitful

using Missings

using Random

using JLD2

using DimensionalData

using DimensionalData: Ti, At

using HCubature

using DataFrames

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
    # **Scratch for `_getordinariness!`, and it is the difference between 32 bytes and ~1.9 GiB
    # a call.** The ordinariness chain allocates about one `Float64`
    # matrix** of the abundances' shape — on a UK 1 km grid with 1000 species that is **1.86 GiB
    # every time a diversity measure is asked for**.
    # **It is scratch, NOT a cache** — it is refilled on every call and never trusted. The cache
    # A cached one goes stale — nothing invalidates it when the ecosystem changes — and is safe only
    # by *position*. Filling in place gives the saving without the hazard.
    # `_calcabundance`/`_calcordinariness` **pass a `UniqueTypes` array straight through** (
    # measured: 32 B and 0 B, returning the same array), so the value handed back aliases this
    # buffer. A caller wanting to keep a result across calls must copy it.
    relativeabun::Matrix{Float64}
    valid::Bool

    # The scratch is an implementation detail of the cache, not something every construction site
    # should have to know about — so it is derived from `netmigration`'s shape here.
    function Cache(netmigration::Matrix{Int64}, totaldemand::Matrix{Float64},
                   valid::Bool)
        return new(netmigration, totaldemand, similar(netmigration, Float64),
                   valid)
    end
end

# A one-liner is the whole of what a cache needs: the default prints all four scratch arrays,
# measured at 462 054 characters. Its contents are meaningless between timesteps — `relativeabun` is
# refilled on every call and never trusted — so its shape is the only fact worth showing.
function Base.show(io::IO, c::Cache)
    nsp, ncells = size(c.netmigration)
    return print(io, "Cache($(nsp) species × $(ncells) cells)")
end

"""
    Ecosystem{Part <: AbstractHabitat} <:
       AbstractEcosystem{Part, SL, NF}

Ecosystem houses information on species and their interaction with their
environment. For species, it holds abundances and locations, as well as
properties such as trait information, `spplist`, and movement types, `lookup`.
For environments, it provides information on environmental conditions and
available resources,`habitat`. Finally, there is a slot for the nichefit
between the environment and the characteristics of the species, `nichefit`.
`elapsed` is the simulation clock (see [`simulationtime`](@ref)) and `seed` the
seed the per-species RNG streams were derived from. `epoch` is the real date elapsed time zero
corresponds to, or `nothing` for a run with no calendar — see [`simulationdate`](@ref).
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
    # the hot path — the epoch is resolved once at build and thereafter only reported.
    epoch::Union{Nothing, Dates.TimeType}

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
        # Tolerances, regimes and the nichefit between them must line up member for member — same
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
                                 epoch)
    end
end

"""
    Ecosystem(spplist::SpeciesList, habitat::GridHabitat,
        nichefit::AbstractNicheFit)

Create an `Ecosystem` given a species list, an abiotic environment and trait
nichefit. An optional population function can be added, `popfun`, which
defaults to generic random filling of the ecosystem. A `seed` may be supplied to
make the run reproducible: it deterministically seeds one random number
generator per species (see [`makerngs`](@ref)), so results are identical
regardless of the number of threads used. If no `seed` is given, one is drawn at
random.
"""
function Ecosystem(popfun::F,
                   spplist::SpeciesList{TL, DM},
                   habitat::GridHabitat,
                   nichefit::AbstractNicheFit;
                   seed::Integer = rand(UInt64)) where {F <: Function, TL, DM}
    # **First, and before `genlookups` below.** The inner constructor checks this too, but it runs
    # last — and `genlookups` divides the regime's cell size by a dispersal distance, so on a
    # geographic grid (`size` in `°`) it raises `DimensionError: ° km⁻¹` before the tailored message
    # is ever reached. That was invisible while the geographic cell size was silently converted to
    # kilometres; making it honest is what exposed the ordering.
    _checksimulatable(habitat)

    # Check there is enough resource to support number of individuals at set up
    # Commented out since before this rewrite, and left so deliberately — but note it would no
    # longer typecheck: `getdemand` returns the demand now, not a total. The total is `_getdemand`.
    #all(_getdemand(spplist.abun, spplist.demand) .<= totalsupply(habitat)) ||
    #error("Environment does not have enough resource to support species")
    # Create matrix landscape of zero abundances
    ml = emptygridlandscape(habitat, spplist)
    # One deterministically-seeded RNG per species, so births/deaths/dispersal
    # and the initial population draw are reproducible across thread counts
    rngs = makerngs(seed, size(ml.matrix, 1))
    # Populate this matrix with species abundances
    popfun(ml, spplist, habitat, nichefit, rngs)
    # Create lookup table of all moves and their probabilities
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

function Ecosystem(spplist::SpeciesList,
                   habitat::GridHabitat,
                   nichefit::AbstractNicheFit;
                   seed::Integer = rand(UInt64))
    return Ecosystem(populate!, spplist, habitat, nichefit, seed = seed)
end

@doc (@doc Ecosystem) Ecosystem(::SpeciesList,
                                ::GridHabitat,
                                ::AbstractNicheFit)

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
    # the step it happens — which is more precise than a length comparison made before the run.
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
                                                 eco.epoch)
end

# The only place `Lookup`'s fields are filled positionally, so the only place their `(y, x)` order
# has to be got right by hand.
function Lookup(df::DataFrame)
    return Lookup(df[!, :Y],
                  df[!, :X],
                  df[!, :Prob],
                  zeros(Float64, nrow(df)),
                  zeros(Int64, nrow(df)))
end

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------
# Without these an ecosystem falls through to Diversity's `AbstractMetacommunity` generic, which
# prints its full type signature — measured at 6 169 characters, constant but unreadable.
#
# **Both sides of the model are named, not just the species count.** An ecosystem is species in an
# environment, so what a reader needs is what the environment offers and what the species are matched
# against it on: the regime axes it experiences and the supply axes it competes for. `_axisnames`
# reads them off the layer types, so nothing here touches a value.

# The axes a layer or collection is on, as a `+`-joined string — `Temperature + Precipitation`.
function _axisnames(l::AbstractLayer)
    return l isa LayerCollection ?
           join(map(m -> nameof(axisof(m)), values(l)), " + ") :
           string(nameof(axisof(l)))
end

function Base.show(io::IO, eco::Ecosystem)
    nsp, ny, nx = size(eco.abundances.grid)
    return print(io,
                 "Ecosystem($(nsp) species, $(ny) × $(nx), ",
                 _axisnames(getregime(eco.habitat)), " | ",
                 _axisnames(getsupply(eco.habitat)),
                 ", t = $(eco.elapsed))")
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

"""
    abundances(cache::CachedEcosystem, tm::Unitful.Time)

Extract abundances for an ecosystem, `cache`, at a certain point in time, `tm`.
If the abundances for that time are missing from the ecosystem, then the
function checks on disk for the last saved version and simulates forward.
"""
function abundances(cache::CachedEcosystem, tm::Unitful.Time)
    return _abundances(cache, tm)[2]
end

"""
    getnichefit(eco::Ecosystem)

Extract niche fits.
"""
function getnichefit(eco::AbstractEcosystem)
    return eco.nichefit
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
# never a number derived from it — and each is defined twice: on an `AbstractEcosystem` (so it covers
# `Ecosystem` *and* `MPIEcosystem`) and on whichever half actually owns it, `GridHabitat` for the two
# environment sides and `SpeciesList` for the two species sides.
# **Totals are a separate family** — `totalsupply` — because a total is a different question from
# the object. Returning a total would make `getdemand` the odd
# one out; it now returns the demand, matching its three siblings.

"""
    getregime(eco::Ecosystem)
    getregime(habitat::GridHabitat)

Return the regime — what each cell *is like* — as the [`AbstractRegime`](@ref) layer (or collection
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

Return the supply — what each cell *provides* — as the [`AbstractSupply`](@ref) layer (or collection
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

Return the species' tolerances — their [`Condition`](@ref)-role requirements, pairing member for
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

Return the species' demands — their [`Resource`](@ref)-role requirements, pairing member for member
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
# gives a one-entry tuple and a collection one entry per member — the same normalisation
# `totalsupply` uses, and the reason one shape reads the same however many layers there are.
#
# Each takes the **assembled ecosystem** or the **object one of the three builders returns** — a
# `GridHabitat` from `GridHabitat`, a `SpeciesList` from `build_species`. Both exist before any
# ecosystem does, which is the workflow the builders encourage: build it, look at it, then assemble.

"""
    cellregime(eco::AbstractEcosystem, y::Int64, x::Int64)
    cellregime(habitat::AbstractHabitat, y::Int64, x::Int64)

Return what the environment **is like** at cell `(y, x)` — a `NamedTuple` keyed by
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

Return what the environment **provides** at cell `(y, x)` — a `NamedTuple` keyed by
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

Return what species `sp` **can tolerate** — a `NamedTuple` keyed by `tolerance_names`, with one entry
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

Return what species `sp` **requires** — a `NamedTuple` keyed by `demand_names`, with one entry per
resource: `(sunlight = 2.0kJ/day, water = 5.0L/day)`.

The **Resource** half of the species side; pairs with [`cellsupply`](@ref), and the values are
directly comparable with it.

These are the **raw unitful rates**, never `resource[sp] * exchange_rate`. The exchange rate is
arithmetic machinery — it defaults to `1/mean(resource)`, so the product is a dimensionless number on
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

Return species `sp`'s dispersal **kernel** — a `GaussianKernel` or `LongTailKernel`. `sp` may be an
index or a name.

**The exception to the `NamedTuple` shape** of the other four: a species has one dispersal kernel,
not one per layer or per resource, so there is nothing to key a tuple by.

It returns the kernel itself rather than a parameterisation, which is both more honest — the
distance distribution is Rayleigh rather than Normal, and a `LongTailKernel` is a 2Dt with a `shape`
a Normal cannot express — and directly reusable: it is exactly what [`AddSpecies`](@ref)'s
`dispersal` keyword takes, so an arrival can be given an existing species' dispersal.

```julia
AddSpecies(dispersal = speciesdispersal(eco, 3), abundance = 500)
```
"""
function speciesdispersal(spplist::SpeciesList, sp::Int64)
    return getkernels(spplist.movement)[sp]
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

Return the simulation time elapsed since the ecosystem was built — zero at
construction, advanced by one `timestep` per [`update!`](@ref).
"""
function simulationtime(eco::AbstractEcosystem)
    return eco.elapsed
end

"""
    simulationdate(eco::AbstractEcosystem)

Return the real date the simulation has reached, or `nothing` for a run with no
epoch.

The epoch is the date [`simulationtime`](@ref) counts from — resolved at
[`build_ecosystem`](@ref) from the environment's dated series, or given there
explicitly. A run whose environment mentions no dates has no epoch, and so no
date: elapsed time is all there is to say about when it is.
"""
function simulationdate(eco::AbstractEcosystem)
    return _shiftdate(eco.epoch, simulationtime(eco))
end

"""
    makeunique(eco::Ecosystem)
Convert type of similarity in [`SpeciesList`](@ref) to UniqueTypes, i.e. an
identity matrix.

The species keep their names: only the similarity structure is dropped.
"""
function makeunique(eco::Ecosystem)
    sppl = eco.spplist
    EcoSISTEM.invalidatecaches!(eco)
    newsppl = SpeciesList{typeof(sppl.tolerance),
                          typeof(sppl.demand),
                          typeof(sppl.movement),
                          UniqueTypes,
                          typeof(sppl.params)}(sppl.names,
                                               sppl.tolerance,
                                               sppl.abun,
                                               sppl.demand,
                                               _uniquetypes(sppl.names),
                                               sppl.movement,
                                               sppl.params,
                                               sppl.native)
    newsppl.susceptible = sppl.susceptible
    return Ecosystem{typeof(eco.habitat), typeof(newsppl),
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
generators, one per species. Species `j` is seeded as `Xoshiro(hash((seed, j)))`
so its random stream is a pure function of `(seed, j)` — independent of how
species are distributed across threads or MPI processes. This is what makes
simulation results reproducible across different thread and process counts (each
species is always processed by exactly one task on one rank, drawing in a fixed
cell order). See [`getrng`](@ref).

Note: this per-species scheme is sufficient only because no single species' draws
are ever split across ranks/tasks. If a species' cells were ever partitioned
across ranks, a per-`(species, cell)` counter-based generator would be needed
instead.
"""
function makerngs(seed::Integer, n::Integer)
    return [Random.Xoshiro(hash((seed, j))) for j in 1:n]
end

# Narrow a user-supplied seed of any integer type to the `UInt64` the ecosystem stores. `%` rather
# than `UInt64(...)`, so a negative seed reinterprets rather than throwing — `seed` is only ever
# hashed, never used as a magnitude. NB the *stored* value is not what `makerngs` hashes (that keeps
# the caller's own value and type), so narrowing here cannot perturb the species streams.
_storedseed(seed::Integer) = seed % UInt64

# Check the abundance matrix `m` (species × subcommunities) has one row per species in `sppl` and one
# column per subcommunity in `part`. A pure dimension check.
function _mcmatch(m::AbstractMatrix, sppl::SpeciesList, part::AbstractHabitat)
    return _counttypes(sppl, true) == size(m, 1) &&
           _countsubcommunities(part) == size(m, 2)
end

# The one implementation, reached only through the [`AddSpecies`](@ref) operation. Every trait
# defaults to `nothing`, which clones the last species — what a plain reintroduction wants.
function _addspecies!(eco::AbstractEcosystem, abun::Integer;
                      tolerance = nothing, demand = nothing,
                      dispersal = nothing, birth = nothing, death = nothing,
                      name = nothing)
    n = counttypes(eco.spplist, true)
    newnames = vcat(eco.spplist.names, isnothing(name) ? string(n + 1) : name)
    newmat = vcat(parent(eco.abundances.matrix),
                  zeros(Int64, 1, size(eco.abundances.matrix, 2)))
    yx = dims(eco.abundances.grid, (Y, X))
    # `GridLandscape` is immutable — the only way to change shape is to construct a whole new
    # one (via the constructor, which reshape-pairs `.matrix`/`.grid` correctly) and reassign the
    # `Ecosystem` field holding it.
    eco.abundances = GridLandscape(newmat, newnames, yx)
    # Give the new species its own RNG stream, derived from `hash((seed, j))` like every other
    # species. Drawing one from an existing species' stream instead would break the scheme twice
    # over: the new stream would depend on how many species had been added rather than on the seed,
    # and the draw would re-phase the demography of the species it came from.
    push!(eco.rngs, Random.Xoshiro(hash((eco.seed, n + 1))))
    repopulate!(eco, abun)
    push!(eco.spplist.names, isnothing(name) ? string(n + 1) : name)
    append!(eco.spplist.abun, abun)
    # An arrival carrying its own traits is by definition not native.
    append!(eco.spplist.native, all(isnothing, (tolerance, demand, dispersal)))
    push!(eco.spplist.susceptible, missing)
    _addtolerance!(eco.spplist.tolerance, tolerance)
    _addmovement!(eco.spplist.movement, dispersal)
    _addparams!(eco.spplist.params, birth, death)
    _adddemand!(eco.spplist.demand, demand)
    # Reassigned, not discarded: `_addtypes!` returns a new `UniqueTypes` because it is immutable.
    # Handed the name just pushed onto `names`, so the two name sources cannot drift apart.
    eco.spplist.types = _addtypes!(eco.spplist.types, eco.spplist.names[end])
    # A `Lookup` for the new species, or the first `move!` is a `BoundsError`. Built the same way
    # the constructor builds them, from the species' own movement kernel.
    push!(eco.lookup,
          genlookups(eco.habitat.regime, getkernels(eco.spplist.movement)[end]))
    # …and the migration cache has to match the new landscape, or applying net migration is a
    # `DimensionMismatch`.
    eco.cache.netmigration = zeros(Int64, size(eco.abundances.matrix))
    eco.cache.valid = false
    return eco
end

# **Each of these takes the new species' own value, and falls back to cloning the last species only
# when none is given** — which is what a plain reintroduction wants.
function _addtolerance!(tolerance::NicheTolerance, dist = nothing)
    return push!(tolerance.dists,
                 isnothing(dist) ? tolerance.dists[end] : dist)
end

# Takes the new species' **set** of tolerated classes, and clones an existing species' set where
# none is given. `push!` of a normalised set, never `append!` of a drawn element: a species'
# tolerance is one set per species, so appending into the set of the species being cloned from would
# alter that species as well as building the new one.
function _addtolerance!(tolerance::SimpleCategoricalTolerance, val = nothing)
    return push!(tolerance.vals,
                 isnothing(val) ? rand(tolerance.vals) : _categoryset(val))
end

# Takes a dispersal **kernel**, not an `AbstractMovement`: the movement *type* belongs to the whole
# assemblage, while the kernel and `disperse_safely` are per-species. That is why `AddSpecies`'
# keyword is `dispersal`, matching `build_species`, rather than `movement`.
#
# **`disperse_safely` must grow with `kernels`, or the first `move!` after `addspecies!` is a
# `BoundsError`** — the same trap the `Lookup` push above already documents, and it bit for the same
# reason: two per-species vectors that must stay the same length, only one of which was being
# extended. Caught by `extras_examples` (the invasion intervention adds a species mid-run), not by
# any unit test. The new species inherits the last one's setting, as its kernel and rates do.
function _addmovement!(mv::AbstractMovement, kernel = nothing)
    push!(mv.kernels, isnothing(kernel) ? mv.kernels[end] : kernel)
    hasproperty(mv, :disperse_safely) &&
        push!(mv.disperse_safely, mv.disperse_safely[end])
    return mv.kernels
end

# **Two vectors, and both must grow** — `birth` and `death` are separate per-species vectors on the
# same params object, so extending one and not the other leaves them different lengths and the
# demographic loop reads off the end. Cloning the last species' rates is the `nothing` default the
# family shares.
function _addparams!(params::AbstractParams, birth = nothing, death = nothing)
    append!(params.birth, isnothing(birth) ? params.birth[end] : birth)
    return append!(params.death, isnothing(death) ? params.death[end] : death)
end

# One entry per species — `counttypes(::Demand)` is `length(demand.resource)` — so this grows the
# vector by exactly one, cloning the last species' demand when none is given, as elsewhere in the
# family. A **multi-resource** demand is a `SpeciesRequirementCollection` of several `Demand`s, not
# a longer `resource`, so this method only ever sees one resource's worth.
function _adddemand!(d::AbstractDemand, resource = nothing)
    return append!(d.resource,
                   isnothing(resource) ? d.resource[end] : resource)
end

# `UniqueTypes` is immutable, so this **returns** a new one and the caller must reassign; assigning
# to the argument would throw the result away.
# Takes the new species' `name` and goes through `_uniquetypes`, so the existing species keep
# theirs — the count constructor used here before renamed **every** species to its index, silently
# and on every `addspecies!`.
function _addtypes!(ut::UniqueTypes, name::AbstractString)
    return _uniquetypes([gettypenames(ut, true); name])
end

# Adding a species to anything else is refused, with the reason rather than a `MethodError` off the
# end of dispatch. A phylogeny-backed `SpeciesList` is the real case: a new species has no position on
# the tree, and inventing one would change what every diversity measure computes.
function _addtypes!(types::AbstractTypes, ::AbstractString)
    return error("`addspecies!` cannot extend a $(typeof(types)): only `UniqueTypes` can gain a " *
                 "species without saying anything more. A phylogeny-backed species list needs the " *
                 "new species' position on the tree, which nothing here can supply — rebuild the " *
                 "`SpeciesList` from a tree that already contains it.")
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
                                                           dims(cache.habitat.active,
                                                                (Y, X)))
            return tm, cache.abundances.matrix[Ti(At(tm))]
        else
            newtm, abun = _abundances(cache, tm - timestep)
            if (newtm > 2 * timestep)
                cache.abundances.matrix[Ti(At(newtm - 2 * timestep))] = missing
            end
        end
    else
        return tm, cache.abundances.matrix[Ti(At(tm))]
    end
    simulate!(cache, newtm, timestep)
    if issave
        @save joinpath(cache.abundances.outputfolder, string(idx, ".jld2")) abuns=SavedLandscape(cache.abundances.matrix[Ti(At(tm))],
                                                                                                 cache.rngs)
    end
    return _abundances(cache, newtm + timestep)
end

# `getdispersaldist`/`getdispersalvar` moved to `deprecations.jl` — superseded by
# `speciesdispersal`, which returns the kernel itself.

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

# Complete a dispersal lookup so every offset has its mirror. The generator emits one of each pair,
# since the kernel is radially symmetric, and the hot loop needs both directions present.
function _symmetricgrid(grid::DataFrame)
    for x in 1:nrow(grid)
        if grid[x, 1] != grid[x, 2]
            push!(grid, hcat(grid[x, 2], grid[x, 1], grid[x, 3]))
        end
    end
    for x in 1:nrow(grid)
        if (grid[x, 1] > 0)
            push!(grid, hcat(-grid[x, 1], grid[x, 2], grid[x, 3]))
        end
        if (grid[x, 2] > 0)
            push!(grid, hcat(grid[x, 1], -grid[x, 2], grid[x, 3]))
        end
        if (grid[x, 1] > 0 && grid[x, 2] > 0)
            push!(grid, hcat(-grid[x, 1], -grid[x, 2], grid[x, 3]))
        end
    end
    return grid
end

# Define gaussian kernel function
function _gaussiandisperse(r)
    return exp(-((r[3] - r[1])^2 + (r[4] - r[2])^2)) / π
end

# The long-tailed dispersal kernel: a 2-D Student-t, whose tail falls off as a power of distance
# rather than exponentially, so a rare long jump stays possible. `b` is the shape parameter that
# `_gaussiandisperse` has no equivalent of, which is why the two take different argument lists and
# `_lookup` has a method for each.
function _2Dtdisperse(r, b)
    return ((b - 1) / (π)) * (1 + ((r[3] - r[1])^2 + (r[4] - r[2])^2))^-b
end

"""
    genlookups(regime::AbstractRegime, kernel::GaussianMovement)

Generate lookup tables, which hold information on the probability of moving to
neighbouring squares.
"""
function genlookups(regime::AbstractRegime, kernel::GaussianKernel)
    sd = (2 * kernel.dist) / sqrt(pi)
    relsize = uconvert(NoUnits, _uniformcellside(regime) / sd)
    m = maximum(getgridshape(regime))
    p = kernel.thresh
    return Lookup(_lookup(relsize, m, p, _gaussiandisperse))
end

function genlookups(regime::AbstractRegime, kernel::LongTailKernel)
    sd = (2 * kernel.dist) / sqrt(pi)
    relsize = uconvert(NoUnits, _uniformcellside(regime) / sd)
    m = maximum(getgridshape(regime))
    p = kernel.thresh
    b = kernel.shape
    return EcoSISTEM.Lookup(EcoSISTEM._lookup(relsize, m, p, b,
                                              EcoSISTEM._2Dtdisperse))
end

# Build one species' dispersal neighbourhood: every offset the kernel reaches with probability above
# `pThresh`, as a `(X, Y, Prob)` table normalised to sum to one. Bounded by the threshold rather than
# by a radius, which is what keeps a long-tailed kernel finite.
function _lookup(relSquareSize::Float64,
                 maxGridSize::Int64,
                 pThresh::Float64,
                 dispersalfn::F) where {F <: Function}
    # Create empty array
    lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

    # Loop through directions until probability is below threshold
    k = 0
    m = 0
    count = 0
    while (k <= maxGridSize && m <= maxGridSize)
        count = count + 1
        calc_prob = hcubature(r -> dispersalfn(r),
                              [0, 0, k * relSquareSize, m * relSquareSize],
                              [
                                  relSquareSize,
                                  relSquareSize,
                                  (k + 1) * relSquareSize,
                                  (m + 1) * relSquareSize
                              ],
                              maxevals = 10000)[1] / relSquareSize^2
        if m == 0 && calc_prob < pThresh
            break
        end
        if count == 1
            push!(lookup_tab, [k m calc_prob])
            k = k + 1
        elseif (calc_prob > pThresh && m <= k)
            push!(lookup_tab, [k m calc_prob])
            m = m + 1
        else
            m = 0
            k = k + 1
        end
    end
    # If no probabilities can be calculated, threshold is too high
    nrow(lookup_tab) != 0 || error("probability threshold too high")
    # Find all other directions
    lookup_tab = _symmetricgrid(lookup_tab)
    #info(sum(lookup_tab[:, 3]))
    # Normalise
    lookup_tab[!, :Prob] = lookup_tab[!, :Prob] / sum(lookup_tab[!, :Prob])
    return lookup_tab
end

# The same, for a kernel that takes a shape parameter `b` as well as a distance — the long-tailed
# family. Kept as a separate method rather than a default `b`, so a kernel that has no shape
# parameter cannot be handed one.
function _lookup(relSquareSize::Float64,
                 maxGridSize::Int64,
                 pThresh::Float64,
                 b::Float64,
                 dispersalfn::F) where {F <: Function}
    # Create empty array
    lookup_tab = DataFrame(X = Int64[], Y = Int64[], Prob = Float64[])

    # Loop through directions until probability is below threshold
    k = 0
    m = 0
    count = 0
    while (k <= maxGridSize && m <= maxGridSize)
        count = count + 1
        calc_prob = hcubature(r -> dispersalfn(r, b),
                              [0, 0, k * relSquareSize, m * relSquareSize],
                              [
                                  relSquareSize,
                                  relSquareSize,
                                  (k + 1) * relSquareSize,
                                  (m + 1) * relSquareSize
                              ],
                              maxevals = 10000)[1] / relSquareSize^2
        if m == 0 && calc_prob < pThresh
            break
        end
        if count == 1
            push!(lookup_tab, [k m calc_prob])
            k = k + 1
        elseif (calc_prob > pThresh && m <= k)
            push!(lookup_tab, [k m calc_prob])
            m = m + 1
        else
            m = 0
            k = k + 1
        end
    end
    # If no probabilities can be calculated, threshold is too high
    nrow(lookup_tab) != 0 || error("probability threshold too high")
    # Find all other directions
    lookup_tab = _symmetricgrid(lookup_tab)
    #info(sum(lookup_tab[:, 3]))
    # Normalise
    lookup_tab[!, :Prob] = lookup_tab[!, :Prob] / sum(lookup_tab[!, :Prob])
    return lookup_tab
end

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
# Only `Date` is promoted — a calendar type from elsewhere (a 360-day or noleap one) must keep its
# own arithmetic, which is why the epoch is typed on the supertype in the first place.
function _shiftdate(epoch::Dates.Date, elapsed::Unitful.Time)
    return _shiftdate(Dates.DateTime(epoch), elapsed)
end

# Advance the simulation clock by one timestep. Called at the end of `update!` (after the
# population dynamics, before the layer update), so a layer change sees the time it is changing
# *to*. Any `Unitful.Time` step is converted to the stored seconds.
function _advanceclock!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return eco.elapsed += uconvert(s, float(timestep))
end

# ══ Functions ══════════════════════════════════════════════════════════════════════════════════

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
`idx`, reshaping the abundance matrix onto `names`/`yx` (the ecosystem's own
species names and `(Y, X)` grid). The saved per-species RNG streams are
restored into `cache.rngs` so the resumed run continues a reproducible random
stream.
"""
function loadfile(cache::CachedEcosystem, file::String, idx::Int,
                  names::Vector{String}, yx::Tuple{<:Y, <:X})
    @load joinpath(file, string(idx, ".jld2")) abuns
    # Restore the per-species RNG streams for a reproducible resumed run
    cache.rngs .= copy.(abuns.rngs)
    return GridLandscape(abuns, names, yx)
end

# **There is no ordinariness to invalidate any more.** It is recomputed on demand into
# `cache.relativeabun`, so no stale value can survive a change to the ecosystem — which is
# exactly the hazard the old `eco.ordinariness` field carried.
function invalidatecaches!(eco::AbstractEcosystem)
    eco.cache.netmigration .= 0
    return eco.cache.valid = false
end

# Write the current abundances into an abundance recording. Indexed by *the species there are now*
# rather than by the whole first dimension, because `generate_storage`'s `maxspecies` may leave room
# for species that have not arrived yet — and `storage[:, :, i] = matrix` would be a
# `DimensionMismatch` the moment it does.
#
# Rows beyond the current count are left as they were. They are `undef` until a species occupies
# them, which is correct: "this species did not exist yet" is not the same as "it had zero
# abundance", and writing zeros would assert the latter.
function _record!(storage::AbstractArray, eco::AbstractEcosystem,
                  counting::Integer)
    nspp = size(eco.abundances.matrix, 1)
    nspp <= size(storage, 1) ||
        error("this run now has $nspp species but the recording was sized for $(size(storage, 1)). " *
              "Pass `maxspecies` to `generate_storage` to leave room for species that arrive " *
              "during the run.")
    storage[1:nspp, :, counting] = eco.abundances.matrix
    return storage
end
