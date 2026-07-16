# SPDX-License-Identifier: LGPL-3.0-or-later

# ---------------------------------------------------------------------------
# AbstractLayer — the materialised, hot-loop grid-layer family
# ---------------------------------------------------------------------------
# A materialised layer is `matrix (+ time) + size + dynamics`, tagged with a `Role`
# (`Habitat` condition vs `Budget` resource) and a `NicheAxis` (what it measures). It
# unifies today's `ContinuousHab`/`ContinuousTimeHab`/`DiscreteHab` and the budget structs.
#
# ADDITIVE for now: these types are defined but nothing is rewired onto them yet (sub-step 1
# of the layer fold). The `AbstractHabitat`/`AbstractBudget` aliases and the retargeting of
# the sim core follow in later sub-steps. `dynamics::HabitatUpdate` reuses the existing
# dynamics struct (renamed to `LayerDynamics` when `update!` is unified).

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
    dynamics::HabitatUpdate
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
    dynamics::HabitatUpdate
end

"""
    LayerCollection2{R <: Role, L1, L2}
    LayerCollection3{R <: Role, L1, L2, L3}

A collection of 2 or 3 layers of the same role `R` over one grid (e.g. temperature +
rainfall). Sub-layer types are kept concrete in named fields for hot-loop type stability.
"""
struct LayerCollection2{R <: Role, L1 <: AbstractLayer{R},
                        L2 <: AbstractLayer{R}} <: AbstractLayer{R}
    l1::L1
    l2::L2
end
struct LayerCollection3{R <: Role, L1 <: AbstractLayer{R},
                        L2 <: AbstractLayer{R},
                        L3 <: AbstractLayer{R}} <: AbstractLayer{R}
    l1::L1
    l2::L2
    l3::L3
end
