# SPDX-License-Identifier: LGPL-3.0-or-later

# ---------------------------------------------------------------------------
# Collections - the shared backing for every "several of these together" family
# ---------------------------------------------------------------------------
# `LayerCollection` (`Layer.jl`), `SpeciesRequirementCollection` (`SpeciesRequirement.jl`)
# and `CombiningFit` (`NicheFit.jl`) are three families of one shape: a
# single type per family holding its members in one field, `nt`, **always** a `NamedTuple` -
# carrying either the names the caller wrote or, for a `Tuple` of members, the names derived from
# each member's own niche axis (`_asnamedtuple`/`_derivenames` below). The arity, the names
# and every member's concrete type therefore all live in the backing's own type, which is what
# makes agreement between two collections a compile-time comparison rather than a runtime check.
#
# These families implement the **standard container interface** - `keys`, `values`, `iterate`,
# `getindex`, `pairs`, `merge` and the rest - by forwarding to that backing; the methods are in
# `BaseInterface.jl`, which also gives a **leaf** the same interface as a one-member container.
# `collect` and `eltype` are the two deliberate omissions, for the reasons given there.

# ---------------------------------------------------------------------------
# Deriving a positional collection's names
# ---------------------------------------------------------------------------
# One backing, always a `NamedTuple`: a `Tuple` is converted at construction and the type carries
# names either way, which makes **duplicate names unrepresentable for free**, since Julia refuses a
# repeated `NamedTuple` key with no code from us at all.
#
# Names come from the members' **axes** wherever those are distinguishable, so two structures that
# pair correctly - a tolerance against the regime it matches - derive the *same* names and check out
# against one another without the caller writing anything.

# The axis structure a collection carries in its type parameter: a `Tuple` of its members' axes. A
# `Tuple{...}` **type** rather than a tuple value, so that it can be a type parameter at all - and so
# that a tolerance collection and a regime collection over the same axes are the same type expression,
# which is what lets `===` compare them whatever their arity.
_axisstructure(members::Tuple) = Tuple{map(axisof, members)...}

# **Not the same question as `_axesagree`, and the difference is deliberate - do not unify them.**
# This asks whether two members *of one collection* can be told apart **by name**; `_axesagree` asks
# whether two *sides* describe the same quantity, where only identity will do. They sit close together
# and read alike, so a sweep that made either match the other would break the check it was tidying.
#
# A `nameof` test is the honest one here, because names **are** `nameof(axis)`: the only pair that
# cannot be told apart is one yielding the same symbol - the same axis twice, or
# `EcoSISTEM.Temperature` against another module's `Temperature`. A `typeintersect` test would be
# conservative rather than wrong, but it also calls `SoilTemperature` and `Temperature` a clash, when
# those name apart perfectly well.
_axesclash(a, b) = nameof(a) === nameof(b)

# Whether every member's axis is distinguishable from every other's.
function _distinctaxes(axes::Tuple)
    for i in eachindex(axes), j in eachindex(axes)
        i < j && _axesclash(axes[i], axes[j]) && return false
    end
    return true
end

# The names a positional collection's members get: their axes, always.
#
# **Two members on one axis is an error asking for names, never a refusal of the model.** Two
# temperature layers are legitimate - `CombiningFit`'s own docstring pairs a summer and a winter one,
# and two identical tolerances correctly give `suit^2` - so what cannot be done is *naming them apart*,
# and the caller is the only one who can say which is which. Numbering them silently would produce a
# collection whose members no other side could refer to.
function _derivenames(members::Tuple)
    axes = map(axisof, members)
    _distinctaxes(axes) && return map(nameof, axes)
    dup = _firstduplicate(axes)
    return error("two members of this collection are on the `$(nameof(dup))` niche axis, so they " *
                 "cannot be told apart by name. Name them yourself with a named tuple - for example " *
                 "`(summer = ..., winter = ...)` - which is also what lets the tolerance, regime and " *
                 "nichefit sides refer to the same layer.")
end

# The first axis that appears twice (by name), for `_derivenames`' error. Only ever called when
# `_distinctaxes` has already said there is one.
function _firstduplicate(axes::Tuple)
    for i in eachindex(axes), j in eachindex(axes)
        i < j && _axesclash(axes[i], axes[j]) && return axes[i]
    end
    return first(axes)
end

# The backing every collection stores: a `NamedTuple` either way. A caller's `NamedTuple` is kept
# exactly as written - naming is theirs to decide, and a derived name must never overwrite one.
_asnamedtuple(nt::NamedTuple) = nt
_asnamedtuple(t::Tuple) = NamedTuple{_derivenames(t)}(t)

# Rebuild `vals` under collection `c`'s own names. Used by a collection derived from another - a
# default nichefit built from a tolerance - so that the two carry identical names and check out
# against one another. Copying the names is what stops the derived collection re-deriving different
# ones from its own members.
_likebacking(c, vals::Tuple) = NamedTuple{keys(c)}(vals)

# The same idea one step earlier, and **not** the same function. `_likebacking` takes a built
# **collection**; this takes the caller's raw `tolerance = ...` or `demand = ...` spec, which is a
# plain `Tuple` or `NamedTuple` and not a collection at all. A caller who named their spec keeps those
# names; a caller who wrote a tuple gets a tuple back, and the collection's own constructor derives
# names from the built members' axes.
_asbacking(::Tuple, vals::Tuple) = vals
_asbacking(nt::NamedTuple, vals::Tuple) = NamedTuple{keys(nt)}(vals)

# Every collection's members must be of its own family's abstract type, and there must be at least
# one. Called from each family's *sole inner constructor*, so a collection of the wrong thing cannot
# be built at all, rather than failing later from somewhere that assumed otherwise.
function _checkbacking(::Type{T}, members::Tuple,
                       what::AbstractString) where {T}
    isempty(members) &&
        error("a $what collection must hold at least one $what.")
    all(m -> m isa T, members) ||
        error("a $what collection holds `$T`s; got $(map(typeof, members)).")
    return nothing
end

# The single role shared by every member of a collection, checked as the collection is built. A
# collection mixing condition and resource members has no meaning - and no `Role` to be tagged with -
# so a disagreement is an error rather than a widening to `Role`.
#
# Serves **both** role-parameterised families: `LayerCollection` and `SpeciesRequirementCollection`
# each pass their own abstract type and each supplies the one `_roleof` method it needs. That is what
# lets the role be *read off* the members rather than chosen, so neither collection needs a role
# parameter at its call site and neither can be tagged with a role its members contradict.
function _sharedrole(::Type{T}, members::Tuple,
                     what::AbstractString) where {T}
    _checkbacking(T, members, what)
    roles = map(_roleof, members)
    allequal(roles) ||
        error("a $what collection's members must all have the same role; got $roles.")
    return first(roles)
end

# ---------------------------------------------------------------------------
# Folds - how members are combined on any path that can be hot
# ---------------------------------------------------------------------------
# Recursive rather than `map` plus `reduce`: measured, `map` stops being unrolled at around seven
# members and starts allocating, while the recursion stays allocation-free and fully inferred at every
# arity. `map` is still fine for construction-time work - `eltype`, `iscontinuous`, coverage.

# Combine `f` applied to each member with the binary `op`, right to left. The two-argument form
# folds the members themselves, which is what most callers want - and it is what lets `op` be
# written as a do-block, since a do-block always binds the *first* argument.
_fold(op, f, t::Tuple{Any}) = f(t[1])
_fold(op, f, t::Tuple) = op(f(t[1]), _fold(op, f, Base.tail(t)))
_fold(op, t::Tuple) = _fold(op, identity, t)

# ---------------------------------------------------------------------------
# Pairing - what it means for two structures to line up member for member
# ---------------------------------------------------------------------------
# A species' tolerances, the environment's regimes and the nichefit between them are three parallel
# structures, as are a demand and a supply. Each pairing must agree in arity, in the names the caller
# wrote, and in what each member actually measures - and this is one check over all of that, which
# reports **which** member disagrees rather than only that something did.
#
# Each family builds its own side, beside its own accessors: a label for the error, the member names,
# each member's `eltype`, each member's niche axis, and each member's `iscontinuous` where that means
# anything - `nothing` for the resource side, which has no continuous or categorical distinction.

# One side of a pairing check, for any member family. One function rather than one per family,
# because the container interface gives every family the same `keys` and `values`, leaving only two
# things to vary: the error label, and whether the family distinguishes continuous from categorical.
#
# `kinds` is **passed** rather than derived from the `Role`, because `AbstractNicheFit` has no role at
# all - which is also why the three families cannot sit under one `{Role, Axis}` supertype.
function _side(x, label::AbstractString, kinds::Bool)
    members = values(x)
    return (label = label, names = keys(x),
            types = map(eltype, members), axes = map(axisof, members),
            kinds = kinds ? map(iscontinuous, members) : nothing)
end

# Check that structures which must line up member for member do. Nothing is ever silently
# reordered, since a different order is a mistake about which layer is which, so each error names
# both sides.
function _checkaligned(what::AbstractString, ref, others...)
    for other in others
        _checkarity(what, ref, other)
        _checknaming(what, ref, other)
        _checkmembers(what, ref, other)
    end
    return nothing
end

# Two paired structures must have the same number of members. Checked before names or axes, because
# a count mismatch makes every later comparison compare the wrong pairs and report a confusing fault
# about the first of them rather than the real one.
function _checkarity(what, ref::NamedTuple, other::NamedTuple)
    length(ref.names) == length(other.names) ||
        error("$what: the $(ref.label) has $(length(ref.names)) layer(s) but the " *
              "$(other.label) has $(length(other.names)) - pass one per environment layer.")
    return nothing
end

# "Named on one side, positional on the other" is not a state that can arise: every collection is
# `NamedTuple`-backed, a collection built from a tuple derives names from its members' axes, and two
# structures that pair correctly derive the **same** ones without the caller writing anything.
#
# So a mismatch here is a real one - different axes, or one side named by hand in a way the other's
# axes do not produce - and the message says where derived names come from, since a caller who wrote
# neither set needs to know what they are being compared against.
function _checknaming(what, ref::NamedTuple, other::NamedTuple)
    ref.names == other.names ||
        error("$what: the $(ref.label)'s layers are named $(ref.names) but the " *
              "$(other.label)'s are named $(other.names) - paired structures must carry the same " *
              "names in the same order; nothing is reordered for you. Names not written by hand " *
              "are each member's niche axis, so a mismatch here means the two sides are on " *
              "different axes or in a different order.")
    return nothing
end

# **Whether two niche axes agree, and it is identity.** Not `typeintersect`, and not `<:`: names are
# derived from axes, so two non-identical axes derive non-matching names and there would be no single
# name for the pair without picking a side or rewrapping a layer, and the design forbids both.
#
# The root is not special. `NicheAxis` is an ordinary axis meaning "I claim nothing", so it pairs with
# itself and with nothing else. That is what makes it a statement rather than a hole: a layer that
# declines to say what it measures cannot be silently read as saying whatever the species needs.
#
# Deliberately not `<:`, even though `SoilTemperature <: Temperature` reads like agreement. Accepting
# it would mean an analysis pairing two `Temperature` sides could be **invalidated** by someone
# declaring a subtype later.
_axesagree(a, b) = a === b

# What each member measures: its niche axis and its `eltype` always, and its continuous/categorical
# kind where both sides have one. Reported per member and by name, so a three-layer mismatch says
# which layer is wrong.
function _checkmembers(what, ref::NamedTuple, other::NamedTuple)
    # `i` also addresses `ref.types` and `other.types`, which are parallel to `names`.
    for (i, name) in enumerate(ref.names)
        # The axis comes **first**, because it is the question the unit cannot answer:
        # `Temperature` and `TemperatureRange` are both `K`, so an `eltype` comparison alone passes a
        # pair that measures two different things. `axes === nothing` on a side means it carries no
        # axis at all, not that it has nothing to say.
        (isnothing(ref.axes) || isnothing(other.axes) ||
         _axesagree(ref.axes[i], other.axes[i])) ||
            error("$what: layer `$name` is on the `$(nameof(ref.axes[i]))` niche axis in the " *
                  "$(ref.label) but `$(nameof(other.axes[i]))` in the $(other.label) - paired " *
                  "layers must be on the *same* axis, not merely a compatible unit. Declare " *
                  "the same axis on both sides (a grouping axis such as " *
                  "`EcoSISTEM.TemperatureAxis` is fine, so long as both say it), or check the " *
                  "order of the species niches against the environment layers.")
        ref.types[i] === other.types[i] ||
            error("$what: layer `$name` is $(ref.types[i]) in the $(ref.label) but " *
                  "$(other.types[i]) in the $(other.label) - the two must measure the same thing " *
                  "in the same unit.")
        (isnothing(ref.kinds) || isnothing(other.kinds) ||
         ref.kinds[i] == other.kinds[i]) ||
            error("$what: layer `$(ref.names[i])` is $(_kindlabel(ref.kinds[i])) in the " *
                  "$(ref.label) but $(_kindlabel(other.kinds[i])) in the $(other.label) - pair a " *
                  "continuous trait (e.g. NicheTolerance) with a continuous regime (UniformSpec / " *
                  "GradientSpec / PeakedSpec) and a continuous fit (NicheSuitability), or a " *
                  "categorical trait (SimpleCategoricalTolerance) with a categorical regime " *
                  "(NicheSpec) and a categorical fit (CategoricalSuitability).")
    end
    return nothing
end

# "continuous"/"categorical" for an `iscontinuous` answer.
_kindlabel(b::Bool) = b ? "continuous" : "categorical"
_kindlabel(bs::AbstractVector{Bool}) = "[" * join(_kindlabel.(bs), ", ") * "]"

# Reject a mismatch in arity between collections that must line up member for member, before it
# becomes a `BoundsError` from somewhere further in. Every length here is known from the backing's
# type, so the check folds away.
function _samearity(ts::Tuple...)
    allequal(map(length, ts)) ||
        error("collections that must line up member for member have differing arities: " *
              "$(map(length, ts)).")
    return nothing
end

# Apply `f` to the corresponding members of parallel tuples, giving a tuple of results. The shared
# arity is a *signature* constraint (`NTuple{N, Any}`, not `Tuple`), so a mismatch cannot reach the
# recursion and the less specific methods below only run to report it.
_zipmap(f, a::Tuple) = _zipped(f, a)
_zipmap(f, a::NTuple{N, Any}, b::NTuple{N, Any}) where {N} = _zipped(f, a, b)
function _zipmap(f, a::NTuple{N, Any}, b::NTuple{N, Any},
                 c::NTuple{N, Any}) where {N}
    return _zipped(f, a, b, c)
end
_zipmap(f, a::Tuple, b::Tuple) = _samearity(a, b)
_zipmap(f, a::Tuple, b::Tuple, c::Tuple) = _samearity(a, b, c)

# The recursion itself - arities already checked by `_zipmap`, so it can index freely.
#
# Recursive rather than `ntuple(..., Val(N))`, and the difference is not cosmetic: measured, the
# `ntuple` form runs into inference's recursion limiter on the path that matters most - `_suitability`
# over a collection assembling per-layer results - and widens the whole tuple to `Tuple{Any, Any}`, at
# 288 bytes a call. The recursion infers exactly, at zero.
_zipped(f, ::Tuple{}) = ()
_zipped(f, a::Tuple) = (f(a[1]), _zipped(f, Base.tail(a))...)
_zipped(f, ::Tuple{}, ::Tuple{}) = ()
function _zipped(f, a::Tuple, b::Tuple)
    return (f(a[1], b[1]), _zipped(f, Base.tail(a), Base.tail(b))...)
end
_zipped(f, ::Tuple{}, ::Tuple{}, ::Tuple{}) = ()
function _zipped(f, a::Tuple, b::Tuple, c::Tuple)
    return (f(a[1], b[1], c[1]),
            _zipped(f, Base.tail(a), Base.tail(b), Base.tail(c))...)
end
