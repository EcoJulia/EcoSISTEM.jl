# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols

"""
    TempChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Increase the temperature for one timestep of the ecosystem using
[`HabitatUpdate`](@ref) information.
"""
function TempChange(eco::AbstractEcosystem, hab::ContinuousHab,
                    timestep::Unitful.Time)
    val = hab.dynamics.rate
    v = uconvert(K / unit(timestep), val)
    return hab.matrix .+= v * timestep
end

"""
    RainfallChange(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Change the rainfall for one timestep of the ecosystem using
[`HabitatUpdate`](@ref) information.
"""
function RainfallChange(eco::AbstractEcosystem, hab::ContinuousHab,
                        timestep::Unitful.Time)
    val = hab.dynamics.rate
    v = uconvert(mm / unit(timestep), val)
    return hab.matrix .+= v * timestep
end

"""
    TempFluct(eco::AbstractEcosystem, hab::ContinuousHab, timestep::Unitful.Time)

Fluctuate the temperature for one timestep of the ecosystem using
[`HabitatUpdate`](@ref) information.
"""
function TempFluct(eco::AbstractEcosystem, hab::ContinuousHab,
                   timestep::Unitful.Time)
    val = hab.dynamics.rate
    v = uconvert(K / unit(timestep), val) * timestep
    offset = v / pi
    return hab.matrix .+= (sin.(hab.matrix ./ offset) .* v)
end

"""
    NoChange(eco::AbstractEcosystem, layer::AbstractLayer, timestep::Unitful.Time)

Keep a layer (habitat or budget) the same for one timestep of the model.
"""
function NoChange(eco::AbstractEcosystem, layer::AbstractLayer,
                  timestep::Unitful.Time) end

ChangeLookup = Dict(K => TempChange, NoUnits => NoChange)

"""
    cyclicChange(eco::AbstractEcosystem, layer::ContinuousLayer, timestep::Unitful.Time)

Advance a time-varying (monthly) layer — a climate habitat or a time budget — by one
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
    HabitatLoss(eco::AbstractEcosystem, hab::AbstractHabitat, timestep::Unitful.Time)

Destroy habitat for one timestep of the ecosystem using [`HabitatUpdate`](@ref)
information. The habitat's `dynamics.rate` is a loss rate (per unit time); over
`timestep` it gives the per-cell probability that an active cell is lost. That many
active cells are drawn at random and have their budget and abundances zeroed.
"""
function HabitatLoss(eco::AbstractEcosystem, hab::AbstractHabitat,
                     timestep::Unitful.Time)
    # Loss rate × timestep is the (dimensionless) probability of losing a cell.
    prob = uconvert(NoUnits, hab.dynamics.rate * timestep)
    pos = findall(vec(eco.abenv.active))
    smp = sample(pos, rand(Binomial(length(pos), prob)); replace = false)
    eco.abenv.budget.matrix[smp] .= zero(eltype(eco.abenv.budget.matrix))
    eco.abundances.matrix[:, smp] .= 0
    return eco
end

# One update rule for any layer regardless of role: apply its own `dynamics.changefun`
# (`NoChange` for a static layer, `TempChange`/`cyclicChange`/… for a dynamic one), and
# recurse into a collection. This unifies the old separate `habitatupdate!`/`budgetupdate!`
# (static budgets were no-ops; time budgets advanced identically to `cyclicChange`).
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
    habitatupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the habitat of an ecosystem for one timestep.
"""
function habitatupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _layerupdate!(eco, eco.abenv.habitat, timestep)
end

"""
    budgetupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)

Update the budget of an ecosystem for one timestep.
"""
function budgetupdate!(eco::AbstractEcosystem, timestep::Unitful.Time)
    return _layerupdate!(eco, eco.abenv.budget, timestep)
end
