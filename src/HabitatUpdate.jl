# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols

"""
    TempChange(eco::AbstractEcosystem, regime::ContinuousRegime, timestep::Unitful.Time)

Increase the temperature for one timestep of the ecosystem using
[`LayerUpdate`](@ref) information.
"""
function TempChange(eco::AbstractEcosystem, regime::ContinuousRegime,
                    timestep::Unitful.Time)
    val = regime.dynamics.rate
    v = uconvert(K / unit(timestep), val)
    return regime.matrix .+= v * timestep
end

"""
    RainfallChange(eco::AbstractEcosystem, regime::ContinuousRegime, timestep::Unitful.Time)

Change the rainfall for one timestep of the ecosystem using
[`LayerUpdate`](@ref) information.
"""
function RainfallChange(eco::AbstractEcosystem, regime::ContinuousRegime,
                        timestep::Unitful.Time)
    val = regime.dynamics.rate
    v = uconvert(mm / unit(timestep), val)
    return regime.matrix .+= v * timestep
end

"""
    TempFluct(eco::AbstractEcosystem, regime::ContinuousRegime, timestep::Unitful.Time)

Fluctuate the temperature for one timestep of the ecosystem using
[`LayerUpdate`](@ref) information.
"""
function TempFluct(eco::AbstractEcosystem, regime::ContinuousRegime,
                   timestep::Unitful.Time)
    val = regime.dynamics.rate
    v = uconvert(K / unit(timestep), val) * timestep
    offset = v / pi
    return regime.matrix .+= (sin.(regime.matrix ./ offset) .* v)
end

"""
    NoChange(eco::AbstractEcosystem, layer::AbstractLayer, timestep::Unitful.Time)

Keep a layer (regime or supply) the same for one timestep of the model.
"""
function NoChange(eco::AbstractEcosystem, layer::AbstractLayer,
                  timestep::Unitful.Time) end

# Axis-keyed layer dynamics: a regime's per-timestep change function is chosen from its niche axis via
# `dynamics(::NicheAxis)` (default `NoChange`, declared in `NicheInfo.jl`).
dynamics(::TemperatureAxis) = TempChange
# The dimensionless / degree-day temperature leaves are not live-changing temperatures, so they
# override the group default back to NoChange (matching their `canonicalunit` overrides).
dynamics(::CumulativeHeat) = NoChange
dynamics(::Isothermality) = NoChange
dynamics(::FrostChangeFrequency) = NoChange
# RainfallChange lives on the `Precipitation` leaf (not the topical `PrecipitationAxis`), so the
# dimensionless `PrecipitationSeasonality` inherits the `NoChange` default.
dynamics(::Precipitation) = RainfallChange

"""
    cyclicChange(eco::AbstractEcosystem, layer::ContinuousLayer, timestep::Unitful.Time)

Advance a time-varying (monthly) layer — a climate regime or a time supply — by one
timestep, wrapping back to the first month once the stored 3-D time series is exhausted
(with a warning). The per-source `eraChange`/`worldclimChange` names are aliases of this —
the time-step-and-wrap logic is identical regardless of source or role.
"""
function cyclicChange(eco::AbstractEcosystem,
                      layer::ContinuousLayer{R, A, V, Arr},
                      timestep::Unitful.Time) where {R, A, V,
                                                     Arr <: AbstractArray{V, 3}}
    monthstep = uconvert(month, timestep)
    layer.time += round(Int64, monthstep / month)
    if layer.time > size(layer.matrix, 3)
        layer.time = 1
        @warn "More timesteps than available, have repeated"
    end
end

const eraChange = cyclicChange
const worldclimChange = cyclicChange

"""
    HabitatLoss(eco::AbstractEcosystem, regime::AbstractRegime, timestep::Unitful.Time)

Destroy regime for one timestep of the ecosystem using [`LayerUpdate`](@ref)
information. The regime's `dynamics.rate` is a loss rate (per unit time); over
`timestep` it gives the per-cell probability that an active cell is lost. That many
active cells are drawn at random and have their supply and abundances zeroed.
"""
function HabitatLoss(eco::AbstractEcosystem, regime::AbstractRegime,
                     timestep::Unitful.Time)
    # Loss rate × timestep is the (dimensionless) probability of losing a cell.
    prob = uconvert(NoUnits, regime.dynamics.rate * timestep)
    pos = findall(vec(eco.habitat.active))
    smp = sample(pos, rand(Binomial(length(pos), prob)); replace = false)
    eco.habitat.supply.matrix[smp] .= zero(eltype(eco.habitat.supply.matrix))
    eco.abundances.matrix[:, smp] .= 0
    return eco
end

# One update rule for any layer regardless of role: apply its own `dynamics.changefun`
# (`NoChange` for a static layer, `TempChange`/`cyclicChange`/… for a dynamic one), and
# recurse into a collection. This unifies the old separate `regimeupdate!`/`supplyupdate!`
# (static supplies were no-ops; time supplies advanced identically to `cyclicChange`).
function _layerupdate!(eco::AbstractEcosystem, layer::AbstractLayer,
                       timestep::Unitful.Time)
    return layer.dynamics.changefun(eco, layer, timestep)
end
function _layerupdate!(eco::AbstractEcosystem, layer::LayerCollection2,
                       timestep::Unitful.Time)
    _layerupdate!(eco, layer.one, timestep)
    return _layerupdate!(eco, layer.two, timestep)
end
function _layerupdate!(eco::AbstractEcosystem, layer::LayerCollection3,
                       timestep::Unitful.Time)
    _layerupdate!(eco, layer.one, timestep)
    _layerupdate!(eco, layer.two, timestep)
    return _layerupdate!(eco, layer.three, timestep)
end

"""
    regimeupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the regime of an ecosystem for one timestep.
"""
function regimeupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _layerupdate!(eco, eco.habitat.regime, timestep)
end

"""
    supplyupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the supply of an ecosystem for one timestep.
"""
function supplyupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _layerupdate!(eco, eco.habitat.supply, timestep)
end
