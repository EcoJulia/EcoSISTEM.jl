# SPDX-License-Identifier: LGPL-3.0-or-later
#
# How offspring disperse: the dispersal kernels, the movement rules that draw from them, and the
# per-species lookup table of dispersal probabilities that the hot loop actually reads.

using Unitful

using SpecialFunctions: loggamma

"""
    AbstractKernel

The shape of the probability distribution over where an individual's offspring lands, relative to its
own cell. [`GaussianKernel`](@ref), [`LongTailKernel`](@ref) (2Dt) and [`WaldKernel`](@ref) (Inverse
Gaussian) are the concrete kernels.

A kernel says only *how far and how likely*. What happens when a draw crosses an edge is the grid's
[`AbstractTopology`](@ref), and what becomes of one aimed at a dead cell is `disperse_safely` on the
movement - three separate questions, answered by three separate things.
"""
abstract type AbstractKernel end

"""
    GaussianKernel(; meandist, sd, thresh = 1.0e-9)

A dispersal kernel whose displacement is Gaussian along each axis independently, so that the
distance travelled is Rayleigh distributed. The density of landing at distance `r` is
`exp(-(r / scale)^2) / (pi * scale^2)` per unit area.

Give exactly one of `meandist` and `sd`. They describe the same kernel: `meandist = sd *
sqrt(pi / 2)`.

# Keywords

  - `meandist`: the mean distance travelled.
  - `sd`: the standard deviation of the displacement along each axis.
  - `thresh`: the probability below which a destination is dropped, which is what bounds the lookup
    table to a finite neighbourhood.

# Fields

  - `scale`: the distance scale in the density above, in kilometres: `sqrt(2) * sd`, the
    root-mean-square distance travelled.
  - `thresh`: as above.
"""
struct GaussianKernel <: AbstractKernel
    scale::typeof(1.0km)
    thresh::Float64

    function GaussianKernel(;
                            meandist::Union{Nothing, Unitful.Length} = nothing,
                            sd::Union{Nothing, Unitful.Length} = nothing,
                            thresh::Real = 1.0e-9)
        _onedistance(GaussianKernel, (meandist = meandist, sd = sd))
        # Do not rearrange the mean's conversion: a kernel given its mean distance must build the
        # positional form's lookup to the bit, or the canonical references move.
        scale = isnothing(sd) ?
                (2 * _kilometres(meandist, :meandist)) / sqrt(pi) :
                sqrt(2) * _kilometres(sd, :sd)
        return new(scale, Float64(thresh))
    end
end

function Base.show(io::IO, k::GaussianKernel)
    return print(io, "GaussianKernel(meandist = ", _kernelphrase(_meandist(k)),
                 ", thresh = ", k.thresh, ")")
end

"""
    LongTailKernel(; shape, scale, meandist, thresh = 1.0e-9)

A dispersal kernel following the two-dimensional Student's t distribution, the "2Dt" of Clark et
al. (1999, Ecology 80: 1475-1494), whose tail falls off as a power of distance rather than
exponentially, so that rare long-distance dispersal is not effectively impossible. The density of
landing at distance `r` is `(shape - 1) / (pi * scale^2) * (1 + (r / scale)^2)^-shape` per unit
area.

Give `shape` and exactly one of `scale` and `meandist`, which are not equivalent: the mean distance
travelled is `scale * sqrt(pi) / 2 * gamma(shape - 3/2) / gamma(shape - 1)`, which is finite only for
`shape > 3/2`, and its variance is finite only for `shape > 2`. The same `scale` with a smaller
`shape` reaches further.

# Keywords

  - `shape`: how heavy the tail is, greater than 1, without which the density does not integrate.
  - `scale`: the distance scale in the density above.
  - `meandist`: the mean distance travelled, refused for `shape <= 3/2`, where it is infinite.
  - `thresh`: the probability below which a destination is dropped, which is what bounds the lookup
    table to a finite neighbourhood.

# Fields

  - `scale`: as above, in kilometres.
  - `shape`, `thresh`: as above.
"""
struct LongTailKernel <: AbstractKernel
    scale::typeof(1.0km)
    shape::Float64
    thresh::Float64

    function LongTailKernel(; shape::Real,
                            scale::Union{Nothing, Unitful.Length} = nothing,
                            meandist::Union{Nothing, Unitful.Length} = nothing,
                            thresh::Real = 1.0e-9)
        shape > 1 ||
            throw(ArgumentError("a `LongTailKernel`'s `shape` must exceed 1, or its density " *
                                "does not integrate; got $shape."))
        _onedistance(LongTailKernel, (scale = scale, meandist = meandist))
        if isnothing(meandist)
            a = _kilometres(scale, :scale)
        else
            shape > 3 / 2 ||
                throw(ArgumentError("a `LongTailKernel` of `shape` $shape has no finite mean " *
                                    "distance, which needs `shape` above 3/2; give its " *
                                    "`scale` instead."))
            a = 2 * _kilometres(meandist, :meandist) /
                (sqrt(pi) * _gammaratio(shape))
        end
        return new(a, Float64(shape), Float64(thresh))
    end
end

function Base.show(io::IO, k::LongTailKernel)
    return print(io, "LongTailKernel(shape = ", k.shape, ", scale = ",
                 _kernelphrase(k.scale), ", thresh = ", k.thresh, ")")
end

"""
    WaldKernel(; meandist, scale, thresh = 1.0e-9)

A dispersal kernel whose distance travelled follows the inverse Gaussian distribution, the "WALD"
kernel that Katul et al. (2005, American Naturalist 166: 368-381) derived for seeds dispersed by
wind. The density of the distance travelled is `sqrt(scale / (2 * pi * r^3)) * exp(-scale * (r -
meandist)^2 / (2 * meandist^2 * r))`, and its variance is `meandist^3 / scale`, so the larger
`scale` is beside `meandist`, the more tightly the distances gather round their mean.

# Keywords

  - `meandist`: the mean distance travelled.
  - `scale`: the distance scale in the density above, Katul et al.'s lambda.
  - `thresh`: the probability below which a destination is dropped, which is what bounds the lookup
    table to a finite neighbourhood.

# Fields

  - `meandist`, `scale`: as above, in kilometres.
  - `thresh`: as above.
"""
struct WaldKernel <: AbstractKernel
    meandist::typeof(1.0km)
    scale::typeof(1.0km)
    thresh::Float64

    function WaldKernel(; meandist::Unitful.Length, scale::Unitful.Length,
                        thresh::Real = 1.0e-9)
        return new(_kilometres(meandist, :meandist), _kilometres(scale, :scale),
                   Float64(thresh))
    end
end

function Base.show(io::IO, k::WaldKernel)
    return print(io, "WaldKernel(meandist = ", _kernelphrase(k.meandist),
                 ", scale = ", _kernelphrase(k.scale), ", thresh = ", k.thresh,
                 ")")
end

"""
    AbstractMovement

A dispersal mode: which individuals disperse at all, the family of kernel saying how far, and what
becomes of an individual dispersing into a dead cell. Calling a mode with its kernel's keywords gives
one species' movement, a [`SpeciesDispersal`](@ref), so species sharing a mode differ only in their
parameters; called with a kernel of its family, it gives that species' movement too:

```julia
insect = InsectDispersal()
insect(meandist = 2.0km)
insect(GaussianKernel(meandist = 2.0km))
```

[`SeedDispersal`](@ref), [`WindDispersal`](@ref), [`InsectDispersal`](@ref) and
[`VertebrateDispersal`](@ref) disperse only the newly born: the first states its family and
`disperse_safely`, and the other three default them for their agent. [`AlwaysDispersal`](@ref)
disperses any individual, and [`NoDispersal`](@ref) none. The grid's own edges are not part of a mode: they belong to the [`GridHabitat`](@ref), because
two species on one grid cannot be on different topologies.

`disperse_safely` says what becomes of an individual dispersing into a **dead cell** - one off the
grid, or an inactive one. `true` redistributes it among the reachable destinations as though it had
never aimed there; `false` loses it, so the share of the kernel pointing at dead cells becomes
mortality.

**It belongs to the mode because it is a fact about the disperser**, not about the map. A
wind-dispersed seed blown out to sea is gone; an animal-dispersed one is carried to somewhere the
animal can actually reach, so two species on one grid may legitimately differ. It cannot empty a cell
outright: staying put is one of the kernel's own destinations, and a cell only disperses while it is
active, so some share always survives.
"""
abstract type AbstractMovement end

"""
    SeedDispersal(family; disperse_safely)

A dispersal mode in which only the newly born disperse, as seeds do, with no particular agent named.
The generic mode, so both its family and its dead-cell fate are stated rather than defaulted.

# Arguments

  - `family`: the type of kernel each species' movement is built from, a concrete
    [`AbstractKernel`](@ref).
  - `disperse_safely`: whether an individual dispersing into a dead cell is redistributed (`true`)
    or lost, as described on [`AbstractMovement`](@ref).

# Fields

  - `disperse_safely`: as above.
"""
struct SeedDispersal{K <: AbstractKernel} <: AbstractMovement
    disperse_safely::Bool

    function SeedDispersal(family::Type{K};
                           disperse_safely::Bool) where {K <: AbstractKernel}
        return new{_concretefamily(family)}(disperse_safely)
    end
end

"""
    WindDispersal(family = WaldKernel; disperse_safely = false)

A dispersal mode in which only the newly born disperse, carried by the wind. A seed blown into a
dead cell is lost by default, as one blown out to sea is.

# Arguments

  - `family`: the type of kernel each species' movement is built from, a concrete
    [`AbstractKernel`](@ref); by default the [`WaldKernel`](@ref) derived for wind dispersal.
  - `disperse_safely`: whether an individual dispersing into a dead cell is redistributed (`true`)
    or lost, as described on [`AbstractMovement`](@ref).

# Fields

  - `disperse_safely`: as above.
"""
struct WindDispersal{K <: AbstractKernel} <: AbstractMovement
    disperse_safely::Bool

    function WindDispersal(family::Type{K} = WaldKernel;
                           disperse_safely::Bool = false) where {K <:
                                                                 AbstractKernel}
        return new{_concretefamily(family)}(disperse_safely)
    end
end

"""
    InsectDispersal(family = GaussianKernel; disperse_safely = true)

A dispersal mode in which only the newly born disperse, carried by insects.

# Arguments

  - `family`: the type of kernel each species' movement is built from, a concrete
    [`AbstractKernel`](@ref).
  - `disperse_safely`: whether an individual dispersing into a dead cell is redistributed (`true`)
    or lost, as described on [`AbstractMovement`](@ref).

# Fields

  - `disperse_safely`: as above.
"""
struct InsectDispersal{K <: AbstractKernel} <: AbstractMovement
    disperse_safely::Bool

    function InsectDispersal(family::Type{K} = GaussianKernel;
                             disperse_safely::Bool = true) where {K <:
                                                                  AbstractKernel}
        return new{_concretefamily(family)}(disperse_safely)
    end
end

"""
    VertebrateDispersal(family = LongTailKernel; disperse_safely = true)

A dispersal mode in which only the newly born disperse, carried by vertebrates, whose occasional long
journeys the long-tailed kernel's power-law tail allows for.

# Arguments

  - `family`: the type of kernel each species' movement is built from, a concrete
    [`AbstractKernel`](@ref).
  - `disperse_safely`: whether an individual dispersing into a dead cell is redistributed (`true`)
    or lost, as described on [`AbstractMovement`](@ref).

# Fields

  - `disperse_safely`: as above.
"""
struct VertebrateDispersal{K <: AbstractKernel} <: AbstractMovement
    disperse_safely::Bool

    function VertebrateDispersal(family::Type{K} = LongTailKernel;
                                 disperse_safely::Bool = true) where {K <:
                                                                      AbstractKernel}
        return new{_concretefamily(family)}(disperse_safely)
    end
end

"""
    AlwaysDispersal(family = GaussianKernel; disperse_safely = true)

A dispersal mode in which any individual may disperse, not only the newly born, as animals do.

# Arguments

  - `family`: the type of kernel each species' movement is built from, a concrete
    [`AbstractKernel`](@ref).
  - `disperse_safely`: whether an individual dispersing into a dead cell is redistributed (`true`)
    or lost, as described on [`AbstractMovement`](@ref).

# Fields

  - `disperse_safely`: as above.
"""
struct AlwaysDispersal{K <: AbstractKernel} <: AbstractMovement
    disperse_safely::Bool

    function AlwaysDispersal(family::Type{K} = GaussianKernel;
                             disperse_safely::Bool = true) where {K <:
                                                                  AbstractKernel}
        return new{_concretefamily(family)}(disperse_safely)
    end
end

"""
    NoDispersal()

A dispersal mode in which nothing disperses: every individual stays in the cell it was born in, so a
species' movement in this mode has no kernel.
"""
struct NoDispersal <: AbstractMovement end

# The modes whose species' movements carry a kernel of the mode's family `K`, and the ones among them
# that disperse only the newly born.
const _KernelDispersal{K} = Union{SeedDispersal{K}, WindDispersal{K},
                                  InsectDispersal{K}, VertebrateDispersal{K},
                                  AlwaysDispersal{K}}
const _BirthDispersal = Union{SeedDispersal, WindDispersal, InsectDispersal,
                              VertebrateDispersal}

# One species' movement: the mode, and a kernel of its family, given built or as its keywords.
function (mode::_KernelDispersal{K})(kernel::K) where {K}
    return SpeciesDispersal(mode, kernel)
end
function (mode::_KernelDispersal{K})(; kernel...) where {K}
    return SpeciesDispersal(mode, K(; kernel...))
end
(mode::NoDispersal)() = SpeciesDispersal(mode, nothing)

# A mode prints its name unqualified, as the modes with a family do.
Base.show(io::IO, ::NoDispersal) = print(io, nameof(NoDispersal), "()")

function Base.show(io::IO, m::_KernelDispersal{K}) where {K}
    return print(io, nameof(typeof(m)), "(", nameof(K), ", disperse_safely = ",
                 m.disperse_safely, ")")
end

"""
    SpeciesDispersal(mode::AbstractMovement, kernel)

One species' movement: its dispersal mode and its kernel, of the mode's family. Calling a mode with
the kernel's keywords builds one, as `InsectDispersal()(meandist = 2.0km)`, and a species list holds
one for each species.

# Arguments

  - `mode`: the [`AbstractMovement`](@ref), saying which individuals disperse and what becomes of one
    aimed at a dead cell.
  - `kernel`: the species' [`AbstractKernel`](@ref), of the mode's family, saying how far; `nothing`
    for [`NoDispersal`](@ref).

# Fields

  - `mode`, `kernel`: as above.
"""
struct SpeciesDispersal{M <: AbstractMovement,
                        K <: Union{Nothing, AbstractKernel}}
    mode::M
    kernel::K

    function SpeciesDispersal(mode::M,
                              kernel::K) where {M <: AbstractMovement,
                                                K <: Union{Nothing,
                                                      AbstractKernel}}
        _family(mode) === K ||
            throw(ArgumentError("a species dispersing by `$(sprint(show, mode))` takes " *
                                "a kernel of its family, `$(_family(mode))`; got a " *
                                "`$K`."))
        return new{M, K}(mode, kernel)
    end
end

function Base.show(io::IO, d::SpeciesDispersal)
    return print(io, "SpeciesDispersal(", sprint(show, d.mode), ", ",
                 sprint(show, d.kernel), ")")
end

"""
    Lookup

One species' dispersal neighbourhood, precomputed once and read by the hot loop: every destination
its kernel can reach, as an offset from the source cell, with the probability of landing there.

The last two fields are scratch, rewritten on every cell of every timestep. Holding them here rather
than allocating them per call is what keeps the movement step allocation-free.

# Fields

  - `y`, `x`: the destination offsets, relative to the cell dispersing. In that order, as everywhere
    else in the package.
  - `p`: the kernel's probability for each, fixed for the run.
  - `pnew`: those probabilities renormalised over the destinations actually available from this cell.
  - `moves`: how many individuals go to each destination this step.

`Lookup(maxsize, thresh) do dx, dy ... end` builds one from a function giving the probability of
landing `dx` cells along x and `dy` along y. The kernel is radially symmetric, so the function is
asked only for `0 <= dy <= dx`. `Lookup()` is the empty neighbourhood of a species that does not
disperse.
"""
struct Lookup
    y::Vector{Int64}
    x::Vector{Int64}
    p::Vector{Float64}
    pnew::Vector{Float64}
    moves::Vector{Int64}

    Lookup() = new(Int64[], Int64[], Float64[], Float64[], Int64[])

    # Every offset `share` reaches with a probability above `thresh`, normalised to sum to one. The
    # kernel is radially symmetric, so `share` is asked only for `0 <= dy <= dx`. The search along
    # `dx` stops at the first `(dx, 0)` below `thresh`, and along `dy` at the first `(dx, dy)` not
    # above it. Each offset found is then mirrored into the other seven directions. The order of the
    # offsets matters: the movement step draws over them in order, so another order would change
    # every run.
    function Lookup(share::Function, maxsize::Integer, thresh::Real)
        x, y, p = Int64[], Int64[], Float64[]
        function add!(dx, dy, q)
            push!(x, dx)
            push!(y, dy)
            return push!(p, q)
        end
        for dx in 0:maxsize
            onaxis = share(dx, 0)
            onaxis < thresh && break
            add!(dx, 0, onaxis)
            for dy in 1:dx
                offaxis = share(dx, dy)
                offaxis > thresh || break
                add!(dx, dy, offaxis)
            end
        end
        isempty(p) && error("probability threshold too high")
        for i in eachindex(p)
            x[i] == y[i] || add!(y[i], x[i], p[i])
        end
        for i in eachindex(p)
            x[i] > 0 && add!(-x[i], y[i], p[i])
            y[i] > 0 && add!(x[i], -y[i], p[i])
            x[i] > 0 && y[i] > 0 && add!(-x[i], -y[i], p[i])
        end
        normalised = p / sum(p)
        return new(y, x, normalised, zeros(Float64, length(normalised)),
                   zeros(Int64, length(normalised)))
    end
end

# One line, because the default prints all five vectors, far too much even for a single species'
# neighbourhood. What identifies a lookup is how far it reaches, not the numbers in it: the last
# three fields are scratch rewritten on every cell of every timestep.
function Base.show(io::IO, l::Lookup)
    return print(io, "Lookup($(length(l.y)) destinations)")
end

# == Functions ==================================================================================

"""
    getkernels(movement::AbstractVector{<:SpeciesDispersal})

Return the dispersal kernel of each species' movement, in species order: `nothing` for a species
that does not disperse.
"""
function getkernels(movement::AbstractVector{<:SpeciesDispersal})
    return [d.kernel for d in movement]
end

"""
    dispersesafely(movement::AbstractVector{<:SpeciesDispersal}, sp::Integer)

Return whether species `sp` is redistributed (`true`) or lost (`false`) when its dispersal is aimed
at a dead cell. A species that does not disperse answers `true`, since nothing can be lost.
"""
function dispersesafely(movement::AbstractVector{<:SpeciesDispersal},
                        sp::Integer)
    return _safely(movement[sp].mode)
end

# A vector's values as one phrase: the single value where all agree, else `lo to hi`. Quantities
# and floats are rounded to three digits, so a 2.4 km kernel reads as `2.4 km`.
function _rangephrase(v)
    isempty(v) && return "none"
    lo, hi = extrema(v)
    lo == hi && return _roundphrase(lo)
    return _roundphrase(lo) * " to " * _roundphrase(hi)
end
function _roundphrase(x::Unitful.Quantity)
    return string(round(typeof(1.0 * unit(x)), x, digits = 3))
end
_roundphrase(x::Real) = string(round(x, digits = 3))

# How far a movement's kernels reach, for its display: a Gaussian or WALD kernel's mean distance,
# and a long-tailed kernel's scale, since its mean may be infinite.
function _dispersalphrase(kernels::AbstractVector)
    return isempty(kernels) ? "none" : _reachphrase(kernels)
end
function _reachphrase(kernels::AbstractVector{<:Union{GaussianKernel,
                                                      WaldKernel}})
    return "mean " * _rangephrase(map(_meandist, kernels))
end
function _reachphrase(kernels::AbstractVector)
    return "scale " * _rangephrase(map(k -> k.scale, kernels))
end

# A kernel's distance in its own display, to six significant figures, since the mean is recovered
# from the stored scale and may differ from the distance given in its last bit.
function _kernelphrase(d::Unitful.Length)
    return string(round(typeof(1.0km), d, sigdigits = 6))
end

# A kernel's mean distance: a Gaussian kernel's from its stored scale, a WALD kernel's as stored.
_meandist(k::GaussianKernel) = k.scale * sqrt(pi) / 2
_meandist(k::WaldKernel) = k.meandist

# The ratio `gamma(shape - 3/2) / gamma(shape - 1)` in a long-tailed kernel's mean distance, through
# the log-gamma function so that a large shape cannot overflow.
_gammaratio(shape::Real) = exp(loggamma(shape - 3 / 2) - loggamma(shape - 1))

# Refuse a kernel given other than exactly one of its alternative distances, naming them.
function _onedistance(::Type{K}, given::NamedTuple) where {K <: AbstractKernel}
    count(!isnothing, values(given)) == 1 && return nothing
    names = join(("`$k`" for k in keys(given)), " and ")
    got = all(isnothing, values(given)) ? "neither" : "both"
    return throw(ArgumentError("a `$(nameof(K))` takes exactly one of $names; got $got."))
end

# The keywords a kernel family must be given to build one: every one of `all`, and exactly one of
# `oneof` where it names any. Kept in step with each constructor by `test_Movement.jl`.
_kernelrequires(::Type{GaussianKernel}) = (all = (), oneof = (:meandist, :sd))

function _kernelrequires(::Type{LongTailKernel})
    return (all = (:shape,),
            oneof = (:scale, :meandist))
end

_kernelrequires(::Type{WaldKernel}) = (all = (:meandist, :scale), oneof = ())

# A kernel's distance keyword in kilometres, refused unless positive.
function _kilometres(d::Unitful.Length, name::Symbol)
    d > zero(d) ||
        throw(ArgumentError("a kernel's `$name` must be positive; got $d."))
    return uconvert(km, d)
end

# A mode's dead-cell fate; one that disperses nothing loses nothing.
_safely(mode::AbstractMovement) = mode.disperse_safely
_safely(::NoDispersal) = true

# The kernel family a mode's species' movements carry: `Nothing` for one that disperses nothing.
_family(::_KernelDispersal{K}) where {K} = K
_family(::NoDispersal) = Nothing

# A mode's kernel family, refused unless concrete, since every species' kernel is built from it.
function _concretefamily(family::Type{K}) where {K <: AbstractKernel}
    isconcretetype(family) ||
        throw(ArgumentError("a dispersal mode's family must be a concrete kernel type, " *
                            "such as `GaussianKernel`; got `$family`."))
    return family
end

# A species list's movements as a vector the step loop reads without allocating: its element type
# the union of the movements' own types, of which there may be at most four, the most Julia splits
# at a call rather than dispatching dynamically.
function _dispersalvector(movement::AbstractVector{<:SpeciesDispersal})
    isempty(movement) && return movement
    types = unique(map(typeof, movement))
    length(types) <= 4 ||
        throw(ArgumentError("a species list's movements may be of at most four kinds, mode " *
                            "and kernel together, so the step loop can split them; got " *
                            "$(length(types)): " *
                            join(("`$(sprint(show, t))`" for t in types),
                                 ", ") * "."))
    return convert(Vector{Union{types...}}, movement)
end

# A species list's movement in one phrase, for its display: the mode and how far its species reach
# where they share one, else each mode and how many species take it.
function _movementphrase(movement::AbstractVector{<:SpeciesDispersal})
    isempty(movement) && return "none"
    modes = unique(d.mode for d in movement)
    if length(modes) == 1
        mode = only(modes)
        mode isa NoDispersal && return sprint(show, mode)
        kernels = [d.kernel for d in movement]
        return sprint(show, mode) * ", " * string(nameof(eltype(kernels))) *
               " " *
               _dispersalphrase(kernels)
    end
    return join(("$(sprint(show, m)) for $(count(d -> d.mode == m, movement))"
                 for m in modes), ", ")
end
