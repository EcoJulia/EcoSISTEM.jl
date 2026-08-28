# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The species side of the environment, mirroring `AbstractLayer` member for member: a `Condition`
# requirement is a tolerance, a `Resource` one a demand, and several of either are one collection.

"""
    SpeciesRequirementCollection{R <: Role, A, C <: NamedTuple}

Several species requirements of the same role `R` — tolerances over several regime layers, or
demands over several supplies — over the **axis structure** `A`, held by name in a `NamedTuple`,
`(Temperature = …, Precipitation = …)`. The species side's counterpart to
[`LayerCollection`](@ref): its members must line up, name for name and in the same order, with the
layer collection they are matched against.

**Members are named by their axis**, so a `Tuple` of them is accepted and named for you. Two members
on the *same* axis cannot be told apart that way, and are refused rather than silently numbered.

**The role is read off the members rather than chosen**, as [`LayerCollection`](@ref) does, so `R`
cannot be got wrong and a collection mixing a tolerance with a demand is refused instead of being
tagged with whichever role came first.

Members are reached through the standard container interface — `sr.Precipitation`, `sr[1]`,
`sr[:Precipitation]`, `keys`, `values`, `pairs`, `iterate`, `length`, `merge` and `NamedTuple(sr)`.
A **single requirement answers identically**, as a one-member container.

# Fields

  - `nt`: the members, by name. Reach them through the interface above rather than directly, so that
    no member name is reserved.
"""
struct SpeciesRequirementCollection{R <: Role, A, C <: NamedTuple} <:
       AbstractSpeciesRequirement{R, A, C}
    nt::C

    # The sole constructor, so there is exactly one spelling and `R` cannot be got wrong — see the
    # docstring, and `LayerCollection`, which this deliberately mirrors line for line.
    function SpeciesRequirementCollection(members::Union{Tuple, NamedTuple})
        nt = _asnamedtuple(members)
        return new{_sharedrole(AbstractSpeciesRequirement, values(nt),
                               "species requirement"),
                   _axisstructure(values(nt)), typeof(nt)}(nt)
    end
end

# --- Display ------------------------------------------------------------------
# The mirror of `LayerCollection`'s display in `Layer.jl`, and it shows the stored name only where
# it differs from the axis, for the same reason: a derived name merely repeats it.
function _requirementsummary(name::Symbol, r::AbstractSpeciesRequirement)
    axis = nameof(axisof(r))
    u = r isa SimpleCategoricalTolerance ? "" : " $(unit(eltype(r)))"
    return name === axis ? "$(axis)$(u)" : "$(name): $(axis)$(u)"
end

function Base.show(io::IO, c::SpeciesRequirementCollection{R}) where {R}
    members = join(map(_requirementsummary, keys(c), values(c)), " + ")
    return print(io, "SpeciesRequirementCollection{$(nameof(R))}($(members))")
end

function Base.show(io::IO, ::MIME"text/plain", c::SpeciesRequirementCollection)
    println(io, sprint(show, c))
    for (name, r) in pairs(c)
        println(io, "  ", rpad(name, 14), sprint(show, r))
    end
    return nothing
end

# --- Reading a requirement ----------------------------------------------------

# A requirement's axis structure, straight off the type parameter — the mirror of
# `axisof(::AbstractLayer{R, A})`.
axisof(::AbstractSpeciesRequirement{R, A}) where {R, A} = A

# The role a requirement is matched on, for `_sharedrole`. Mirrors `_roleof(::AbstractLayer{R})`.
_roleof(::AbstractSpeciesRequirement{R}) where {R} = R

# `eltype` is defined on the leaf types (`ContinuousTolerance`, `AbstractCategoricalTolerance`,
# `Demand`) and never on the abstract. A collection *is* an `AbstractTolerance`, so a method on the
# abstract would catch one and hand back its whole backing `NamedTuple` as an "element type" —
# the third parameter means the element type for a leaf and the members for a collection. A
# collection has no element type, and asking is a `MethodError` rather than a guess.
#
# The container interface itself is forwarded in `BaseInterface.jl`, with the rest of it.
