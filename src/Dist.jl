# SPDX-License-Identifier: LGPL-3.0-or-later

using Distributions
import Distributions: @check_args, ContinuousUnivariateDistribution, rand,
                      params, pdf

using Unitful
using Unitful: absoluteunit, NoUnits

using Random
import Random: rand, GLOBAL_RNG

"""
    Trapezoid{T<:Real} <: ContinuousUnivariateDistribution

Trapezoidal distribution as described at
https://en.wikipedia.org/wiki/Trapezoidal_distribution.
"""
struct Trapezoid{T <: Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
    c::T
    d::T

    function Trapezoid(a::T, b::T, c::T, d::T;
                       check_args::Bool = true) where {T <: Real}
        @check_args Trapezoid (a,
                               a ≤ b ≤ c ≤ d && a < d,
                               "Arguments must satisfy a ≤ b ≤ c ≤ d and a < d in Trapezoid distribution.")
        return new{T}(a, b, c, d)
    end
end

# Convenience constructors: `promote` mixed `Real` arguments to a common numeric type, and widen
# all-`Integer` arguments to `Float64`, before the checked inner constructor.
function Trapezoid(a::Real, b::Real, c::Real, d::Real; check_args::Bool = true)
    return Trapezoid(promote(a, b, c, d)...; check_args = check_args)
end

function Trapezoid(a::Integer, b::Integer, c::Integer, d::Integer;
                   check_args::Bool = true)
    return Trapezoid(Float64(a),
                     Float64(b),
                     Float64(c),
                     Float64(d);
                     check_args = check_args)
end

# Zero-argument default: `Trapezoid(0, 0, 1, 1)`. With `a == b` and `c == d` the rising and falling regions
# vanish, so it is the uniform distribution on the unit interval `[0, 1]`.
Trapezoid() = Trapezoid(0.0, 0.0, 1.0, 1.0)

# `Distributions` interface for `Trapezoid`: the `(a, b, c, d)` parameter tuple, and a single draw via the
# global RNG.
params(d::Trapezoid) = (d.a, d.b, d.c, d.d)
rand(d::Trapezoid) = rand(GLOBAL_RNG, d)

# Draw a sample by inverse-transform sampling: split `[0, 1]` into the three regions by probability mass —
# rising `[a, b]`, flat `[b, c]`, falling `[c, d]` — pick the region the uniform draw `u` falls in, then
# invert that region's CDF (√-shaped in the triangular tails, linear in the flat middle).
function rand(rng::AbstractRNG, T::Trapezoid)
    (a, b, c, d) = params(T)
    b_m_a = b - a
    c_m_b = c - b
    d_m_c = d - c
    Cϕ = 4 / ((b_m_a * 2) + (c_m_b * 4) + (d_m_c * 2))
    pi1 = Cϕ * b_m_a / 2
    pi2 = Cϕ * c_m_b
    pi3 = Cϕ * d_m_c / 2
    u = rand(rng)
    if (u >= 0 && u <= pi1)
        return a + (u / pi1)^(1 / 2) * b_m_a
    elseif (u > pi1 && u <= (1 - pi3))
        return b + ((u - pi1) / pi2) * c_m_b
    else
        return d - ((1 - u) / pi3)^(1 / 2) * d_m_c
    end
end

# Trapezoidal density: linear rise on `[a, b]`, the constant `2 / ((d + c) - (a + b))` on `[b, c]`, linear
# fall on `[c, d]`, and zero outside `[a, d]`.
function pdf(T::Trapezoid, x::Real)
    (a, b, c, d) = params(T)
    d_p_c = d + c
    a_p_b = a + b
    b_m_a = b - a
    d_m_c = d - c
    if (x >= a && x < b)
        return (2 / (d_p_c - a_p_b)) * ((x - a) / b_m_a)
    elseif (x >= b && x < c)
        return (2 / (d_p_c - a_p_b))
    elseif (x >= c && x <= d)
        return (2 / (d_p_c - a_p_b)) * ((d - x) / d_m_c)
    else
        return 0.0
    end
end

# Support bounds — a Trapezoid has support [a, d] (so `param_units` can recognise its parameters as
# positions on the support axis).
Base.minimum(d::Trapezoid) = d.a
Base.maximum(d::Trapezoid) = d.d

# The `Distributions` accessor function for each conversion role, used by `param_roles` to test which
# parameters a distribution family exposes as its location / scale / rate / shape.
const ROLE_ACCESSORS = (location = location, scale = scale,
                        rate = rate, shape = shape)

# Reconstruct a representative instance of distribution family `D`, purely for role/type introspection (not
# for modelling): try constructing from the leading `k` probe values for decreasing `k` until one succeeds.
# If float probes produced integer-typed fields, retry with integer probes (over every int/float combination
# of positions) to recover the clean parameter types, keeping the float instance if none construct.
# `probes = nothing` draws random, increasing probe values in `(0, 1)`.
function guess_params(D::Type{<:UnivariateDistribution}; probes = nothing)
    if isnothing(probes)
        len = length(fieldnames(D))
        probes = cumsum(rand(Uniform(0.0, 1.0), len))
        if len > 1
            probes ./= rand(Uniform(sum(probes), 2 * sum(probes)))
        end
    end

    d = nothing
    for k in length(probes):-1:0
        d = try
            D(probes[1:k]...)
        catch
            nothing
        end
        d === nothing || break
    end

    d === nothing &&
        error("guess_params: could not construct $D from the probe values $probes; " *
              "pass `probes` explicitly.")

    # If the float probes produced integer-typed fields, retry with integer probes to recover clean
    # parameters — trying every integer/float combination of the parameters (bitmask `j` over `k`
    # positions). Keep the good float `d` if no integer combination constructs.
    if any(typeof.(getfield.(Ref(d), fieldnames(D))) .<: Ref(Integer)) &&
       !any(typeof.(probes) .<: Ref(Integer))
        probeints = rand(1:10, length(fieldnames(D)))
        for k in length(probes):-1:0, j in 0:(2 ^ k - 1)
            newprobes = []
            for i in 1:k
                if (j >> (i - 1)) & 1 == 1
                    push!(newprobes, probeints[i])
                else
                    push!(newprobes, probes[i])
                end
            end
            dint = try
                D(newprobes...)
            catch
                nothing
            end
            if dint !== nothing
                d = dint
                break
            end
        end
    end

    return d
end

# Tag each parameter of `D` with the conversion role(s) it plays. Build a representative instance
# (`guess_params`), then for each role compare its `Distributions` accessor value against every parameter
# value: a match tags that parameter with the role; a role whose value matches no parameter is `derived`
# (carried by the distribution but not a stored parameter). Returns a named tuple
# `(dist, params, types, roles, derived)` consumed by `param_roles_resolved` / `all_positions`.
function param_roles(D::Type{<:UnivariateDistribution}; probes = nothing)
    d = guess_params(D, probes = probes)
    p = params(d)
    roles = [Symbol[] for _ in p]
    derived = Symbol[]
    for (role, acc) in pairs(ROLE_ACCESSORS)
        hasmethod(acc, Tuple{typeof(d)}) || continue
        v = try
            acc(d)
        catch
            continue
        end
        hit = false
        for (i, pv) in enumerate(p)
            pv == v && (push!(roles[i], role); hit = true)
        end
        hit || push!(derived, role)
    end
    return (dist = d, params = p, types = typeof.(p), roles = roles,
            derived = derived)
end

# Deterministic, well-separated probe values (distinct, increasing, in `(0,1)`) that construct valid
# instances of the target families — so role detection is reproducible by default. `probes = nothing`
# falls back to the random probes in `guess_params`.
function _default_probes(D::Type{<:UnivariateDistribution})
    return [k / (length(fieldnames(D)) + 1) for k in 1:length(fieldnames(D))]
end

# Resolve each parameter to a single *conversion role*. Accessor-matched roles (`param_roles`) plus
# `all_positions` (a support bound is a location). An unmatched, non-bound parameter is a dimensionless
# `:shape`; conflicting matched roles error. Transformed-scale families are special-cased below.
function param_roles_resolved(D::Type{<:UnivariateDistribution};
                              probes = _default_probes(D))
    info = param_roles(D; probes = probes)
    poss = all_positions(D, info)
    return map(enumerate(info.roles)) do (i, rs)
        poss && return :location
        :location in rs && return :location
        :scale in rs && return :scale
        :rate in rs && return :rate
        isempty(rs) && return :shape
        length(rs) == 1 && return only(rs)
        return error("parameter $i of $D matched conflicting roles $rs")
    end
end
# LogNormal/LogitNormal expose location/scale accessors, but their parameters live on the log/logit scale
# (dimensionless) — force `:shape` so they take the shape-only (offset/scale) placement path.
param_roles_resolved(::Type{<:LogNormal}; probes = nothing) = [:shape, :shape]
param_roles_resolved(::Type{<:LogitNormal}; probes = nothing) = [:shape, :shape]

# Convert one bare parameter value `x`, interpreted in support unit `us`, to a bare number in the target
# frame `uc`, per role: a location/bound is a *position* (proper affine conversion into `uc` — `0/°C→273.15`
# on a K frame; identity on a °C frame); a scale is an *interval* (a width, whose magnitude is frame-offset
# independent, so it is stripped in `absoluteunit(uc)` — `1/°C→1`, `1/°F→5/9` — never absolute-converting the
# width even when the frame `uc` is affine); a rate is `1/absolute(uc)`; a shape is dimensionless. For a
# non-affine `uc`, `absoluteunit(uc) == uc`; `uc` (not `absoluteunit(us)`) also handles `°F→K`/`cm→mm`.
_readparam(::Val{:location}, uc, us, x) = ustrip(uc, x * us)
_readparam(::Val{:scale}, uc, us, x) = ustrip(absoluteunit(uc), x * us - 0 * us)
_readparam(::Val{:rate}, uc, us, x) = ustrip(inv(absoluteunit(uc)), x * inv(us))
_readparam(::Val{:shape}, uc, us, x) = float(x)

# Reference position/width (bare, in the target frame `uc`) for a shape-only distribution's `offset`/
# `scale`, given either as a bare number (interpreted in the support unit `us`, exactly like a location/
# scale parameter) or as a `Quantity` (converted directly). `offset` is a *position* (proper affine into
# `uc`); `scale` is an *interval* (a width, stripped in `absoluteunit(uc)`). The `Quantity` methods must be
# listed since `Quantity <: Number`.
_refpos(uc, us, v::Unitful.AbstractQuantity) = ustrip(uc, v)
_refpos(uc, us, v::Number) = ustrip(uc, v * us)
function _refwidth(uc, us, v::Unitful.AbstractQuantity)
    return ustrip(absoluteunit(uc), v - 0 * unit(v))
end
_refwidth(uc, us, v::Number) = ustrip(absoluteunit(uc), v * us - 0 * us)

"""
    read_distribution(D, support_unit, bare; canonical = absoluteunit(support_unit),
                      offset = nothing, scale = nothing, probes = nothing)

Build one species' response distribution of type `D` in the `canonical` (absolute) frame from bare
parameter `values` interpreted in `support_unit`. Roles (`param_roles_resolved`) drive the per-parameter
conversion; a **shape-only** distribution (no location/scale/rate role — e.g. `Beta`, `LogNormal`) is
placed on the dimensioned axis via `LocationScale`, requiring `offset` (reference position) + `scale`
(reference width) to produce `dist * scale + offset`.
"""
function read_distribution(D::Type{<:UnivariateDistribution}, us, bare;
                           canonical = absoluteunit(us),
                           offset = nothing, scale = nothing,
                           probes = _default_probes(D),
                           roles = param_roles_resolved(D; probes = probes))
    us === NoUnits && return D(float.(bare)...)
    length(roles) == length(bare) ||
        error("$D expects $(length(roles)) parameters, got $(length(bare)).")
    if all(==(:shape), roles)
        (isnothing(offset) || isnothing(scale)) &&
            error("$D has no location/scale/rate parameter to carry the dimension of $us; pass " *
                  "`offset` (reference position) and `scale` (reference width), both in $us, to place it.")
        return _refpos(canonical, us, offset) +
               _refwidth(canonical, us, scale) * D(float.(bare)...)
    end
    return D((_readparam(Val(roles[i]), canonical, us, bare[i])
              for i in eachindex(bare))...)
end

# The ABSOLUTE unit each role's parameter carries (used for validation): location/scale → the absolute
# support unit (K for °C); rate → its inverse; shape → dimensionless. Everything is absolute — the affine
# handling lives in the value conversion (`_readparam`), not here.
function role_units(role::Symbol, u)
    return role === :shape ? NoUnits :
           role === :location ? absoluteunit(u) :
           role === :scale ? absoluteunit(u) :
           role === :rate ? inv(absoluteunit(u)) :
           error("unknown role :$role")
end

# Is `D` a "bounds" family — every parameter a *position* on the support axis (Uniform, Trapezoid,
# TriangularDist, Arcsine, …), so all parameters carry the support's own unit and affine support units are
# fine (positions have well-defined differences)? Detected by **perturbation**: a bounds family's finite
# support *moves* when a parameter changes, whereas a fixed-support shape family's (Beta, …) does not — so
# Beta's shape params (which merely fall inside `[0,1]`) are not mistaken for bounds.
function all_positions(D::Type{<:UnivariateDistribution}, info)
    lo, hi = try
        extrema(support(info.dist))
    catch
        return false
    end
    (isfinite(lo) && isfinite(hi)) || return false
    p = collect(float.(info.params))
    span = hi - lo
    for i in eachindex(p)
        for δ in (0.1span + 0.05, -(0.1span + 0.05), 0.37span + 0.11)
            q = copy(p)
            q[i] += δ
            d2 = try
                D(q...)
            catch
                nothing
            end
            d2 === nothing && continue
            supp2 = try
                extrema(support(d2))
            catch
                (lo, hi)
            end
            supp2 == (lo, hi) || return true   # a parameter moved the support ⇒ bounds family
        end
    end
    return false
end

"""
    param_units(D, support_units; probes = _default_probes(D))

The absolute unit of each parameter of distribution `D` given a support of `support_units`, inferred
by parameter role (see `param_roles_resolved`): a location or scale parameter carries the
**absolute** support unit (K for °C), a rate carries its inverse, and a shape parameter is
dimensionless — so a shape-only family (`Beta`, `LogNormal`, …) is all-`NoUnits`. This is the
introspection companion to [`read_distribution`](@ref), which performs the actual (affine-aware)
value conversion when a distribution is read from bare parameters.
"""
function param_units(D::Type{<:UnivariateDistribution},
                     support_units::Unitful.Units;
                     probes = _default_probes(D))
    support_units === NoUnits &&
        return fill(NoUnits, length(param_roles_resolved(D; probes = probes)))
    return [role_units(r, support_units)
            for r in param_roles_resolved(D; probes = probes)]
end
