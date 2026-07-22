# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols
using Diversity
using EcoSISTEM.Units
using Printf
using JLD2
import Diversity.Gamma

"""
    simulate!(eco::AbstractEcosystem, times::Unitful.Time, timestep::Unitful.Time)

Run an ecosystem, `eco` for a specified length of time, `times`, for a
particular timestep, `timestep`.
"""
function simulate!(eco::AbstractEcosystem, duration::Unitful.Time,
                   timestep::Unitful.Time)
    times = length((0s):timestep:duration)
    for i in 1:times
        update!(eco, timestep)
    end
end

"""
    simulate!(cache::CachedEcosystem, srt::Unitful.Time, timestep::Unitful.Time)

Run a cached ecosystem, `cache` at a specified timepoint, `srt`, for a
particular timestep, `timestep`.
"""
function simulate!(cache::CachedEcosystem, srt::Unitful.Time,
                   timestep::Unitful.Time)
    eco = Ecosystem{typeof(cache.habitat), typeof(cache.spplist),
                    typeof(cache.relationship)}(copy(cache.abundances.matrix[srt]),
                                                cache.spplist,
                                                cache.habitat,
                                                cache.ordinariness,
                                                cache.relationship,
                                                cache.lookup,
                                                cache.cache,
                                                cache.rngs)
    update!(eco, timestep)
    return cache.abundances.matrix[srt + timestep] = eco.abundances
end

"""
    simulate!(eco::Ecosystem, times::Unitful.Time, timestep::Unitful.Time,
              cacheInterval::Unitful.Time, cacheFolder::String,
              scenario_name::String)

Run an ecosystem, `eco` for specified length of times, `duration`, for a
particular timestep, 'timestep'. A cache interval and folder/file name are
specified for saving output.
"""
function simulate!(eco::Ecosystem,
                   times::Unitful.Time,
                   timestep::Unitful.Time,
                   cacheInterval::Unitful.Time,
                   cacheFolder::String,
                   scenario_name::String)
    time_seq = zero(times):timestep:times
    for i in eachindex(time_seq)
        update!(eco, timestep)
        # Save cache of abundances
        if mod(time_seq[i], cacheInterval) == zero(time_seq[i])
            @save joinpath(cacheFolder,
                           scenario_name *
                           (@sprintf "%02d.jld2" uconvert(NoUnits,
                                                          time_seq[i] /
                                                          cacheInterval))) abun=eco.abundances.matrix
        end
    end
end

"""
    generate_storage(eco::Ecosystem, times::Int64, reps::Int64)

Allocate an integer array of shape `(numSpecies, gridSize, times, reps)` for
recording species abundances across the ecosystem `eco` over multiple timesteps
and replicate runs.
"""
function generate_storage(eco::Ecosystem, times::Int64, reps::Int64)
    numSpecies = length(eco.spplist.abun)
    gridSize = _countsubcommunities(eco.habitat.regime)
    return abun = Array{Int64, 4}(undef, numSpecies, gridSize, times, reps)
end

"""
    generate_storage(eco::Ecosystem, qs::Int64, times::Int64, reps::Int64)

Allocate a float array of shape `(gridSize, qs, times, reps)` for recording
diversity values across the ecosystem `eco` for `qs` diversity orders over
multiple timesteps and replicate runs.
"""
function generate_storage(eco::Ecosystem, qs::Int64, times::Int64, reps::Int64)
    gridSize = _countsubcommunities(eco.habitat.regime)
    return abun = Array{Float64, 4}(undef, gridSize, qs, times, reps)
end

"""
    simulate_action!(action!::Function, eco::AbstractEcosystem, times::Unitful.Time,
                     interval::Unitful.Time, timestep::Unitful.Time;
                     scenario = nothing, offset = false)

Run an ecosystem `eco` up to time `times` in steps of `timestep`, calling the
user-supplied `action!` at regular intervals so that any periodic task can be
performed as the simulation proceeds — recording a quantity, logging progress,
applying a management intervention, checking a stopping condition, and so on.
This is the general engine behind the [`simulate_record!`](@ref) and
[`simulate_record_diversity!`](@ref) recorders; use it directly when you want to
do something they do not.

At each step the ecosystem is advanced with [`update!`](@ref) and, if a `scenario`
is given, modified with `runscenario!` (e.g. to remove regime or change climate).
Whenever the elapsed time falls on a multiple of `interval`, `action!(counting)`
is called, where `counting` is the 1-based index of that occurrence (handy as a
storage slot when the action is recording); `interval` must be a whole multiple of
`timestep`. Anything the action needs to read or update — the ecosystem, an output
array, an external counter — is captured by the closure, typically written as a
`do` block:

```julia
totals = zeros(Int, length((0s):interval:times))
simulate_action!(eco, times, interval, timestep) do counting
    totals[counting] = sum(eco.abundances.matrix)
end
```

`offset` shifts the action grid to start at `timestep` rather than `0`, which drops
the first occurrence (one fewer action in total). Use it to make the number of
actions match a pre-allocated array (such as one from [`generate_storage`](@ref));
the built-in diversity recorders pass `offset = iseven(size(storage, 3))`.

Returns the ecosystem `eco`, now advanced to `times`.
"""
function simulate_action!(action!::F,
                          eco::AbstractEcosystem,
                          times::Unitful.Time,
                          interval::Unitful.Time,
                          timestep::Unitful.Time;
                          scenario = nothing,
                          offset = false) where {F <: Function}
    iszero(mod(interval, timestep)) ||
        error("Interval must be a multiple of timestep")
    action_seq = offset ? (timestep:interval:times) : ((0s):interval:times)
    time_seq = offset ? (timestep:timestep:times) : ((0s):timestep:times)
    counting = 0
    for i in eachindex(time_seq)
        update!(eco, timestep)
        isnothing(scenario) ||
            runscenario!(eco, timestep, scenario, time_seq[i])
        if time_seq[i] in action_seq
            counting += 1
            action!(counting)
        end
    end
    return eco
end

"""
    simulate_record!(storage::AbstractArray, eco::Ecosystem, times::Unitful.Time,
         interval::Unitful.Time, timestep::Unitful.Time)
    simulate_record!(storage::AbstractArray, eco::Ecosystem, times::Unitful.Time,
         interval::Unitful.Time, timestep::Unitful.Time, scenario::AbstractScenario)

Run an ecosystem, `eco` for a specified length of time, `times`, for a
particular timestep, `timestep`, recording abundances into `storage` at each
time interval `interval`. If a `scenario` is given, it is also run at each
timestep to modify the ecosystem, allowing simulation of events such as regime
loss or climate change.

Pre-allocate `storage` with [`generate_storage`](@ref)`(eco, ntimes, reps)`,
where `ntimes = length((0s):interval:times)` is the number of recordings.

To record diversity rather than raw abundances, see
[`simulate_record_diversity!`](@ref); to perform an arbitrary action at regular
intervals via a callback, see [`simulate_action!`](@ref).
"""
function simulate_record!(storage::AbstractArray,
                          eco::Ecosystem,
                          times::Unitful.Time,
                          interval::Unitful.Time,
                          timestep::Unitful.Time,
                          scenario::Union{Nothing, AbstractScenario} = nothing)
    iszero(mod(interval, timestep)) ||
        error("Interval must be a multiple of timestep")
    record_seq = (0s):interval:times
    time_seq = (0s):timestep:times
    storage[:, :, 1] = eco.abundances.matrix
    counting = 1
    for i in 2:length(time_seq)
        update!(eco, timestep)
        isnothing(scenario) ||
            runscenario!(eco, timestep, scenario, time_seq[i])
        if time_seq[i] in record_seq
            counting = counting + 1
            storage[:, :, counting] = eco.abundances.matrix
        end
    end
    return storage
end

"""
    simulate_record_diversity!(storage, eco, times, interval, timestep,
                               divfun, qs::Vector{Float64})
    simulate_record_diversity!(storage, eco, times, interval, timestep,
                               scenario::SimpleScenario, divfun, qs::Vector{Float64})
    simulate_record_diversity!(storage, storage2, eco, times, interval, timestep,
                               qs::Vector{Float64})
    simulate_record_diversity!(storage, eco, times, interval, timestep,
                               divfuns::Array{Function}, q::Float64)
    simulate_record_diversity!(storage, eco, times, interval, timestep,
                               scenario::SimpleScenario, divfuns::Vector{Function},
                               q::Float64)

Run an ecosystem `eco` up to `times` in steps of `timestep`, recording diversity
into `storage` (and, for the alpha/beta/gamma form, `storage2`) every `interval`,
which must be a whole multiple of `timestep`. These are all thin wrappers over
[`simulate_action!`](@ref) — see it for the recording mechanics — and differ only
in what diversity they record:

  - `divfun, qs` — a single diversity function `divfun` (which returns a
    `DataFrame` with a `:diversity` column) evaluated over the diversity orders
    `qs`, reshaped into `storage`;
  - `storage, storage2, …, qs` — normalised alpha, normalised beta and gamma
    diversity over `qs`; subcommunity-level values are written to `storage`
    (gridSize × 3 × timepoints × qs) and metacommunity-level values to `storage2`
    (3 × timepoints × qs);
  - `divfuns, q` — several diversity functions at a single diversity order `q`,
    one per column of `storage`.

The `scenario::SimpleScenario` variants additionally apply `scenario` at each
timestep to modify the ecosystem (e.g. removal of regime patches).

For the `divfun`/`divfuns` forms, pre-allocate `storage` with
[`generate_storage`](@ref)`(eco, ncols, ntimes, reps)`, where `ncols` is
`length(qs)` (or `length(divfuns)`) and `ntimes = length((0s):interval:times)`.
"""
function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    scenario::SimpleScenario,
                                    divfun::F,
                                    qs::Vector{Float64}) where {F <: Function}
    simulate_action!(eco, times, interval, timestep;
                     scenario = scenario) do counting
        return storage[:, :, counting] = reshape(divfun(eco, qs)[!, :diversity],
                                                 countsubcommunities(eco),
                                                 length(qs))
    end
    return storage
end

function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    divfun::F,
                                    qs::Vector{Float64}) where {F <: Function}
    simulate_action!(eco, times, interval, timestep;
                     offset = iseven(size(storage, 3))) do counting
        diversity = divfun(eco, qs)[!, :diversity]
        return storage[:, :, counting] = reshape(diversity,
                                                 Int(length(diversity) /
                                                     length(qs)),
                                                 length(qs))
    end
    return storage
end

function simulate_record_diversity!(storage::AbstractArray,
                                    storage2::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    qs::Vector{Float64})
    simulate_action!(eco, times, interval, timestep;
                     offset = iseven(size(storage, 3))) do counting
        measures = [NormalisedAlpha, NormalisedBeta, Gamma]
        for (i, msr) in enumerate(measures)
            dm = msr(eco)
            diversity = subdiv(dm, qs)[!, :diversity]
            diversity2 = metadiv(dm, qs)[!, :diversity]
            storage[:, :, i, counting] = reshape(diversity,
                                                 Int(length(diversity) /
                                                     length(qs)),
                                                 length(qs))
            storage2[:, i, counting] = diversity2
        end
    end
    return storage, storage2
end

function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    divfuns::Array{Function},
                                    q::Float64)
    simulate_action!(eco, times, interval, timestep) do counting
        for j in eachindex(divfuns)
            storage[:, j, counting] .= divfuns[j](eco, q)[!, :diversity][1]
        end
    end
    return storage
end

function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    scenario::SimpleScenario,
                                    divfuns::Vector{Function},
                                    q::Float64)
    simulate_action!(eco, times, interval, timestep;
                     scenario = scenario) do counting
        for j in eachindex(divfuns)
            storage[:, j, counting] .= divfuns[j](eco, q)[!, :diversity][1]
        end
    end
    return storage
end
