# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoBase
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

# ---------------------------------------------------------------------------
# Layer dynamics
# ---------------------------------------------------------------------------
# (Kept named `HabitatUpdate` for now; renamed `LayerDynamics` when `update!` is unified.)

"""
    HabitatUpdate{F <: Function, DT}

Stores the update rule for a layer. `changefun` is applied to the layer matrix at each
timestep, and `rate` is the rate of change with appropriate units.
"""
struct HabitatUpdate{F <: Function, DT}
    changefun::F
    rate::DT
end

"""
    HabitatUpdate(changefun::F, rate::DT, ::Type{D})

Construct a `HabitatUpdate` with a type-checked rate. Errors if `rate * 1month` does not
have dimensions `D`.
"""
function HabitatUpdate(changefun::F, rate::DT,
                       ::Type{D}) where {F <: Function, DT,
                                         D <: Unitful.Dimensions}
    dimension(rate * 1month) isa D ||
        error("Failed to match types $(rate * 1month) vs $D")
    return HabitatUpdate(changefun, rate)
end

# ---------------------------------------------------------------------------
# AbstractLayer — the materialised, hot-loop grid-layer family
# ---------------------------------------------------------------------------
# A materialised layer is `matrix (+ time) + size + change`, tagged with a `Role`
# (`Habitat` condition vs `Budget` resource) and a `NicheAxis` (what it measures). It
# unifies today's `ContinuousHab`/`ContinuousTimeHab`/`DiscreteHab` (and, from sub-step 3,
# the budget structs). During the fold the old `*Hab` names below are **aliases** over these
# types, so existing methods/constructors keep working (defaulting to the `Unclassified`
# axis) until each path is threaded with its real axis.

abstract type AbstractLayer{R <: Role} <: EcoBase.AbstractGrid end

"""
    Unclassified <: NicheAxis

The default niche axis for a layer built without a declared one (an unclassified layer —
it still has data and a unit, but is not modelled on a named axis).
"""
struct Unclassified <: NicheAxis end

"""
    ContinuousLayer{R <: Role, A <: NicheAxis, V <: Number, Arr <: AbstractArray{V}}

A continuous (numeric) grid layer of role `R` on niche axis `A`, holding a value `matrix`
(eltype `V`). The array type `Arr` carries the **rank**: a `Matrix{V}` is a static layer
(`time` unused, `= 1`), an `AbstractArray{V, 3}` a monthly time series indexed by `time`.
"""
mutable struct ContinuousLayer{R <: Role, A <: NicheAxis, V <: Number,
                               Arr <: AbstractArray{V}} <: AbstractLayer{R}
    matrix::Arr
    time::Int64
    size::Unitful.Length
    change::HabitatUpdate
end

"""
    DiscreteLayer{A <: NicheAxis, V, Arr <: AbstractArray{V}}

A discrete (categorical) grid layer on niche axis `A` — always a `Habitat` (there is no
discrete budget).
"""
mutable struct DiscreteLayer{A <: NicheAxis, V, Arr <: AbstractArray{V}} <:
               AbstractLayer{Habitat}
    matrix::Arr
    size::Unitful.Length
    change::HabitatUpdate
end

"""
    LayerCollection2{R <: Role, L1, L2}
    LayerCollection3{R <: Role, L1, L2, L3}

A collection of 2 or 3 layers of the same role `R` over one grid (e.g. temperature +
rainfall). Sub-layer types are kept concrete in named fields for hot-loop type stability.
"""
struct LayerCollection2{R <: Role, L1 <: AbstractLayer{R},
                        L2 <: AbstractLayer{R}} <: AbstractLayer{R}
    one::L1
    two::L2
end
struct LayerCollection3{R <: Role, L1 <: AbstractLayer{R},
                        L2 <: AbstractLayer{R},
                        L3 <: AbstractLayer{R}} <: AbstractLayer{R}
    one::L1
    two::L2
    three::L3
end

# ---------------------------------------------------------------------------
# Back-compat aliases + constructors (habitat role)
# ---------------------------------------------------------------------------
# `AbstractBudget` is NOT aliased here — it stays Energy.jl's `AbstractBudget{Requirement}`
# until the budget fold (sub-step 3).

const AbstractHabitat = AbstractLayer{Habitat}

const ContinuousHab{C} = ContinuousLayer{Habitat, A, C, Matrix{C}} where {A}
const ContinuousTimeHab{C, M <: AbstractArray{C, 3}} = ContinuousLayer{Habitat,
                                                                       A, C,
                                                                       M} where {A}
const DiscreteHab{D} = DiscreteLayer{A, D, Matrix{D}} where {A}
const HabitatCollection2{H1, H2} = LayerCollection2{Habitat, H1, H2}
const HabitatCollection3{H1, H2, H3} = LayerCollection3{Habitat, H1, H2, H3}

# Old positional constructors → new layer types, defaulting to the `Unclassified` axis and
# time = 1 (the axis-aware `materialise` path constructs `ContinuousLayer{Habitat, A}` with
# the real axis directly).
function ContinuousHab(matrix::Matrix{C}, size::Unitful.Length,
                       change::HabitatUpdate) where {C}
    return ContinuousLayer{Habitat, Unclassified, C, Matrix{C}}(matrix, 1, size,
                                                                change)
end
function ContinuousTimeHab(matrix::AbstractArray{C, 3}, time::Integer,
                           size::Unitful.Length,
                           change::HabitatUpdate) where {C}
    return ContinuousLayer{Habitat, Unclassified, C, typeof(matrix)}(matrix,
                                                                     time, size,
                                                                     change)
end
function DiscreteHab(matrix::Matrix{D}, size::Unitful.Length,
                     change::HabitatUpdate) where {D}
    return DiscreteLayer{Unclassified, D, Matrix{D}}(matrix, size, change)
end

# The collection aliases fix `R = Habitat` but leave the sub-layer params free, so the
# auto-generated constructor doesn't cover `HabitatCollection2(l1, l2)`; forward to the bare
# `LayerCollection` constructor (which infers the role from the sub-layers).
function HabitatCollection2(l1::AbstractLayer{Habitat},
                            l2::AbstractLayer{Habitat})
    return LayerCollection2(l1, l2)
end
function HabitatCollection3(l1::AbstractLayer{Habitat},
                            l2::AbstractLayer{Habitat},
                            l3::AbstractLayer{Habitat})
    return LayerCollection3(l1, l2, l3)
end
