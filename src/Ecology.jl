# SPDX-License-Identifier: LGPL-3.0-or-later
#
# THE MODEL, IN FOURTEEN DECLARATIONS. This file is the conceptual framework EcoSISTEM simulates,
# and nothing else: every type you can dispatch on without committing to a representation. It is the
# first file the package loads, and it should be readable on its own by someone who never opens
# another.
#
# The whole design turns on one distinction. A cell of the world is described in two independent
# ways: what it is LIKE (a `Condition` - temperature, rainfall, land cover) and what it PROVIDES
# (a `Resource` - the pool species compete for). Everything else is that distinction, seen from the
# environment's side or the species' side:
#
#     Role -+- Condition -- environment: a regime  -- species: a tolerance
#           +- Resource  -- environment: a supply  -- species: a demand
#
# The two sides mirror each other member for member. `AbstractLayer` is what the environment holds
# and `AbstractSpeciesRequirement` is what a species brings; each carries a role and a niche axis, so
# a regime pairs with the tolerance on the same axis and a supply with the demand on the same axis.
# The four aliases below name those four faces, and pinning the role in each is what makes a demand
# refuse to stand in for a tolerance.
#
# `AbstractNicheFit` closes the triangle: given a cell's layer and a species' requirement on the same
# axis, how well suited is that species to that cell?
#
# Then the three that assemble it: `AbstractSpecies` (what an ecosystem's species must be able to
# answer), `AbstractHabitat` (an environment over a set of cells) and `AbstractEcosystem` (species in
# an environment, with the fit between them).
#
# What is NOT here: anything that says *how*. `ContinuousLayer` commits to a numeric grid,
# `NicheTolerance` to a response distribution, `StudyGrid` to a regular lattice - all are choices of
# representation and live with their own machinery. The test is simple: if removing it would change
# what the model IS rather than how it is computed, it belongs here.

using Diversity
using Diversity.API
using EcoBase

"""
    Role

The two ways a cell of the world bears on the species living in it: as a
[`Condition`](@ref) it must tolerate, or as a [`Resource`](@ref) it competes for.

This distinction is the axis the whole model is built along, and it applies to **both sides at
once**. The environment states each role as a layer ([`AbstractLayer`](@ref)) - a regime and a
supply - and each species brings a matching requirement
([`AbstractSpeciesRequirement`](@ref)) - a tolerance and a demand. A role therefore names a
*relationship*, not a kind of data: the same temperature grid is a condition because of how species
meet it, not because of what it contains.

The two roles behave quite differently, which is why they are distinguished rather than merged: a
condition is experienced identically by everything in the cell, while a resource is finite and
shared, so what one species takes bears on the rest.
"""
abstract type Role end

"""
    Condition <: Role

What a cell **is like** - temperature, rainfall, land cover, elevation.

A condition is a state, not a stock. Every species in the cell experiences the same value; nothing
is divided between them and nothing is used up. What varies is how well each species copes, which is
what its [`AbstractTolerance`](@ref) says: a species with a narrow tolerance does well over a small
range of the condition and poorly outside it.

Conditions therefore set **where** a species can persist. They do not, by themselves, limit how many
individuals a cell can hold - that is what a [`Resource`](@ref) does.
"""
struct Condition <: Role end

"""
    Resource <: Role

What a cell **provides** - solar radiation, available water, land.

A resource is a shared pool that species compete for, and it is what regulates abundance. Each
species states how much it needs as an [`AbstractDemand`](@ref); the demands of everything present
are summed, and births fall and deaths rise as that total approaches what the cell supplies. Species
therefore interact only through this total. Carrying capacity is not a parameter anywhere in the
model - it emerges from supply divided by demand.
"""
struct Resource <: Role end

"""
    NicheAxis

The quantity being measured - [`Temperature`](@ref), [`Precipitation`](@ref),
[`SolarRadiation`](@ref) - independent of the unit it arrives in and orthogonal to [`Role`](@ref).

An axis is what makes a pairing meaningful: a species' tolerance for temperature is matched against
the *temperature* regime and nothing else. Both sides of the model carry one, so an
[`AbstractLayer`](@ref) and the [`AbstractSpeciesRequirement`](@ref) it is matched against name the
same axis, and a mismatch is a type error rather than a silent comparison of unrelated quantities.

Paired axes must be identical, not merely compatible: a `SoilTemperature` tolerance does not meet a
`Temperature` regime.

Extend it with [`@nicheaxis`](@ref). Related axes are grouped under an abstract intermediate named
`...Axis` ([`TemperatureAxis`](@ref), [`WaterAxis`](@ref)) that carries their shared interface.
"""
abstract type NicheAxis end

"""
    AbstractLayer{R <: Role, A}

What the environment holds: one measured quantity, stated for every cell, in one of the two
[`Role`](@ref)s - a **regime** if it is a [`Condition`](@ref), a **supply** if it is a
[`Resource`](@ref).

This is the environment's half of the model. Its counterpart is
[`AbstractSpeciesRequirement`](@ref), which is what a species brings to the same axis; the two are
matched role for role and axis for axis. Concrete kinds are [`ContinuousLayer`](@ref),
[`CategoricalLayer`](@ref) and the collection [`LayerCollection`](@ref), which holds several layers
over one set of cells.

The axis is carried in the type, not merely checked when an ecosystem is built, so that a signature
can require a regime and the tolerance matched to it to be on the same axis.

Note that `A` is the axis *structure* rather than always a single axis: a single layer carries its
own axis (`Temperature`), a collection a `Tuple` of its members' (`Tuple{Temperature,
Precipitation}`). That is why the parameter carries no `<: NicheAxis` bound, and it is what lets a
matched tolerance and regime be compared as one type whatever their arity. A collection's `A` is
always read off its members rather than chosen.
"""
abstract type AbstractLayer{R <: Role, A} end

"""
    AbstractRegime{A}

A **regime**: what the environment is like, on one axis, in every cell - the
[`Condition`](@ref) half of [`AbstractLayer`](@ref).

A species meets a regime with its [`AbstractTolerance`](@ref) on the same axis, and how well the two
agree is its suitability there. A regime is not divided between the species present: all of them
experience the same value.
"""
const AbstractRegime{A} = AbstractLayer{Condition, A}

"""
    AbstractSupply{A}

A **supply**: what the environment provides, on one axis, in every cell - the
[`Resource`](@ref) half of [`AbstractLayer`](@ref).

A species meets a supply with its [`AbstractDemand`](@ref) on the same axis. Unlike a regime, a
supply is shared: the demands of everything in the cell are summed against it, and that ratio is what
regulates births and deaths.
"""
const AbstractSupply{A} = AbstractLayer{Resource, A}

"""
    AbstractSpeciesRequirement{R <: Role, A, V}

What a species brings to its environment on one axis: a **tolerance** for a
[`Condition`](@ref), or a **demand** for a [`Resource`](@ref).

This is the species' half of the model, and the exact mirror of [`AbstractLayer`](@ref). A species
list holds the two requirements; a habitat holds the two layers; and they are matched role for role
and axis for axis:

| role | the environment holds | each species brings |
|---|---|---|
| [`Condition`](@ref) | a regime - [`AbstractRegime`](@ref) | a tolerance - [`AbstractTolerance`](@ref) |
| [`Resource`](@ref) | a supply - [`AbstractSupply`](@ref) | a demand - [`AbstractDemand`](@ref) |

`V` is the species' own response type: the type of the values the requirement is stated in, so that
it can be compared with the layer's. A temperature tolerance carries kelvin, a solar demand
kilojoules per day, a land-cover tolerance the class codes it accepts. For a collection of
requirements `V` is instead the members' `NamedTuple`, since there is no single one.

`A` follows [`AbstractLayer`](@ref) exactly, including that it is the axis *structure* rather than
always a single axis. That symmetry is the point: a tolerance and the regime it is matched against
carry the **same** `A`.

Write [`AbstractTolerance`](@ref) or [`AbstractDemand`](@ref) in preference to this. Each pins the
role, so a demand cannot be used where a tolerance is meant. This name is for the rarer case of
saying something true of both.
"""
abstract type AbstractSpeciesRequirement{R <: Role, A, V} end

"""
    AbstractTolerance{A, V}

A **tolerance**: how well a species copes with a [`Condition`](@ref), across the range
that condition can take - the species-side mirror of [`AbstractRegime`](@ref).

A tolerance is what decides where a species can persist. Matched against the regime on the same
axis, it yields that species' suitability in each cell, which raises its death rate where conditions
suit it poorly and lowers it where they suit it well.

This is a parametric alias for [`AbstractSpeciesRequirement`](@ref) with the role pinned to
[`Condition`](@ref), and pinning it is what keeps the constraint real: a demand is **not** an
`AbstractTolerance`, so every signature written against this name refuses one. `eltype` is inherited
from the shared supertype.
"""
const AbstractTolerance{A, V} = AbstractSpeciesRequirement{Condition, A, V}

"""
    AbstractDemand{A, V}

A **demand**: how much of a [`Resource`](@ref) one individual of a species needs - the
species-side mirror of [`AbstractSupply`](@ref).

Demands are what make species compete. Summed over everything present in a cell and set against what
that cell supplies, they regulate births and deaths, and so decide how many individuals the cell can
hold.

This is a parametric alias for [`AbstractSpeciesRequirement`](@ref) with the role pinned to
[`Resource`](@ref), and pinning it is what keeps the constraint real: a tolerance is **not** an
`AbstractDemand`, so every signature written against this name refuses one. `V` is the type the
amount is stated in, and `eltype` is inherited from the shared supertype.
"""
const AbstractDemand{A, V} = AbstractSpeciesRequirement{Resource, A, V}

"""
    AbstractNicheFit{A, V}

How well a species is suited to a cell: the rule that scores a species'
[`AbstractSpeciesRequirement`](@ref) against the [`AbstractLayer`](@ref) it is matched to.

This closes the triangle. The environment states a value, the species states what it needs, and the
fit says how well the two agree - a number that raises or lowers that species' birth and death rates
in that cell. Where several axes are in play, a [`CombiningFit`](@ref) reduces their scores to one.

`A` sits in the same slot here as in [`AbstractLayer`](@ref) and
[`AbstractSpeciesRequirement`](@ref), so a single signature can require a layer, the requirement
matched to it and the fit between them to be on the same axis. As there, it is the axis
*structure* - a single axis for one fit, a `Tuple` of them for a combining one.

`V` is the value type the two sides are compared in - the same `V` the requirement carries - so the
fit knows the frame it is working in and a layer's values can be brought into it before scoring.
"""
abstract type AbstractNicheFit{A, V} end

"""
    AbstractSpecies

The species living in an ecosystem, and everything the model needs to know about
them.

Species are what the simulation counts and moves. Each one brings a tolerance for every condition
and a demand for every resource, so that it can be matched against the environment axis by axis; it
disperses its offspring in its own way, and breeds and dies at its own rate. Those four things are
what make one species behave differently from another in the same place.

Species do not currently interact directly. Two of them affect one another only by drawing on the
same resource: their demands are summed against what a cell supplies, and everything present feels
the result equally. There are currently no pairwise terms - no predation, no interference - so an
assemblage is regulated entirely by what it collectively needs against what its surroundings provide.
"""
abstract type AbstractSpecies <: AbstractTypes end

"""
    AbstractHabitat{H <: AbstractRegime, B <: AbstractSupply,
                    L <: EcoBase.AbstractLocationData} <: AbstractPartition

The physical world a simulation runs in, with nothing living in it yet: what each
place is like (`H`, an [`AbstractRegime`](@ref)), what it provides (`B`, an
[`AbstractSupply`](@ref)), and where those places are (`L`).

A habitat is the setting rather than the story. It says that this place is this warm and this wet
and offers this much light, and that place is different - but it says nothing about who lives there
or how they fare. An [`AbstractEcosystem`](@ref) is what results from putting species into one.

Its conditions and supplies need not be fixed. A habitat can warm, dry out or lose land cover as the
simulation proceeds, which is how climate and land-use change enter the model; a species then
experiences a different environment from one year to the next without itself having changed.

A habitat also knows where its places are, so that neighbours can be identified for dispersal and
results can be put back on a map. Only that much is required: this package uses a regular grid of
cells, but irregular areas or scattered points would serve equally well.
"""
abstract type AbstractHabitat{H <: AbstractRegime, B <: AbstractSupply,
                              L <: EcoBase.AbstractLocationData} <:
              AbstractPartition{L} end

"""
    AbstractEcosystem{Part <: AbstractHabitat, SL <: AbstractSpecies,
        NF <: AbstractNicheFit} <: AbstractMetacommunity{Float64,
            Matrix{Int64}, Matrix{Float64}, SL, Part}

A community of species in a place, changing through time: an environment (`Part`),
the species living in it (`SL`), and the rule scoring how well each species suits each place (`NF`).

This is the whole model. What it holds is an abundance for every species in every place, and what it
does is let those abundances change: individuals are born, die and disperse to neighbouring places,
step by step. Conditions decide **where** a species can persist at all, by making it die faster
where its surroundings suit it poorly; resources decide **how many** can persist there, because the
demands of everything present are set against what the place supplies. Nothing sets a carrying
capacity - it is what those two pressures leave behind.

An ecosystem is a metacommunity: a set of local communities, one per place, drawn from a shared pool
of species. Diversity can therefore be asked of it at either level - of a single place, or of the
landscape as a whole, accounting for the fact that the same species appears in many of them.
"""
abstract type AbstractEcosystem{Part <: AbstractHabitat,
                                SL <: AbstractSpecies,
                                NF <: AbstractNicheFit} <:
              AbstractMetacommunity{Float64, Matrix{Int64}, Matrix{Float64}, SL,
                                    Part} end
