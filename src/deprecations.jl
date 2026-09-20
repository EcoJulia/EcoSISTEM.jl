# SPDX-License-Identifier: LGPL-3.0-or-later

# ===========================================================================
# Deprecations - main `EcoSISTEM` module
#
# Every deprecated public API is collected here, in sections, and included late in
# `EcoSISTEM.jl`, after all the types it shims. Each shim warns - so downstream code gets a
# migration message rather than a silent `MethodError` - and forwards to the current API. Mirrored
# by `test/test_deprecations.jl`. The shims whose signatures name a `RasterDataSources` type are in
# `ext/EcoSISTEMRasterDataSourcesExt/deprecations.jl`, since they cannot be defined without it.
#
# ## Every section says which release deprecated it
#
# Each section header ends with a `Deprecated in vX.Y.Z` line saying which release the shim belongs
# to, and **the sections are ordered newest release first**, so the oldest deprecations sit at the
# foot of the file and dropping a release's worth is deleting from the bottom up to the first
# section carrying a newer label. `clean_Deprecations.jl` asserts every section has a label and
# that the labels never get newer going down the file.
#
# Two things a whole-section delete does not reach, and both fail quietly:
#
#   - an `export` or `public` declaration for a deprecated name that sits anywhere but beside its
#     shim. A stale `export` of a deleted binding is a precompile *warning*, not an error.
#   - `test/test_deprecations.jl`, which mirrors this file and must lose the same sections.
# ===========================================================================

# ---------------------------------------------------------------------------
# `simulate_action!`: a callback on `simulate!`
#
# `simulate!(f, eco, duration, timestep; every)` hands its callback the occurrence's count, elapsed
# time and date, and each occurrence observes the state at the time its schedule names.
# `simulate_action!` handed a bare count, one step after each multiple of the interval; the shim
# keeps that timing exactly.
#
# Deprecated in v0.8.0.
# ---------------------------------------------------------------------------
"""
    simulate_action!(action!::Function, eco::AbstractEcosystem, times::Unitful.Time,
                     interval::Unitful.Time, timestep::Unitful.Time;
                     intervention = nothing, offset = false)

Deprecated: give [`simulate!`](@ref) the callback instead, `simulate!(f, eco, duration, timestep;
every)`, whose callback is handed the occurrence's `count`, `elapsed` time and `date` and sees the
state at each time `every` names. This calls `action!(counting)` on the step after the clock stood
on a multiple of `interval`, which must be a whole multiple of `timestep`; `offset` starts that grid
at `timestep`, the run a step shorter. `intervention` is applied as for `simulate!`.
"""
function simulate_action!(action!::F, eco::AbstractEcosystem,
                          times::Unitful.Time, interval::Unitful.Time,
                          timestep::Unitful.Time; intervention = nothing,
                          offset = false) where {F <: Function}
    Base.depwarn("`simulate_action!` is deprecated: give `simulate!` the callback instead, as " *
                 "`simulate!(eco, duration, timestep, every = EveryInterval(interval)) do " *
                 "occurrence ... end`, whose callback is handed `(count, elapsed, date)`.",
                 :simulate_action!)
    return _simulateaction!(action!, eco, times, interval, timestep,
                            intervention = intervention, offset = offset)
end

export simulate_action!

# ---------------------------------------------------------------------------
# The recording functions: recorders passed to `simulate!`
#
# `RecordAbundance`, `RecordDiversity` and `SaveAbundance` are values handed to `simulate!` in place
# of a callback, each keeping the run's provenance. `simulate_record!` records into the same slots
# over the same steps through `RecordAbundance`; the diversity recorders and the caching `simulate!`
# keep their own loops and timing unchanged.
#
# Deprecated in v0.8.0.
# ---------------------------------------------------------------------------
"""
    simulate_record!(storage::AbstractArray, eco::Ecosystem, times::Unitful.Time,
                     interval::Unitful.Time, timestep::Unitful.Time; intervention = nothing)

Deprecated: use `simulate!(RecordAbundance(storage), eco, duration, timestep, every =
EveryInterval(interval))` with a [`RecordAbundance`](@ref). This records the starting state and each
multiple of `interval`, a whole multiple of `timestep`, into `storage[:, :, k]` over `times /
timestep` steps, applying `intervention` as `simulate!` does, and returns `storage`.
"""
function simulate_record!(storage::AbstractArray, eco::Ecosystem,
                          times::Unitful.Time, interval::Unitful.Time,
                          timestep::Unitful.Time; intervention = nothing)
    Base.depwarn("`simulate_record!` is deprecated: use `simulate!(RecordAbundance(storage), eco, " *
                 "duration, timestep, every = EveryInterval(interval))`, which keeps the run's " *
                 "provenance too.", :simulate_record!)
    iszero(mod(interval, timestep)) ||
        error("Interval must be a multiple of timestep")
    # The steps it always took, `times / timestep` rounded down, for any `times`.
    steps = length((0s):timestep:(times - timestep))
    simulate!(RecordAbundance(storage), eco, steps * timestep, timestep,
              every = EveryInterval(interval), intervention = intervention)
    return storage
end

"""
    simulate_record_diversity!(storage, eco, times, interval, timestep, divfun, qs::Vector{Float64})
    simulate_record_diversity!(substorage, metastorage, eco, times, interval, timestep,
                               qs::Vector{Float64})
    simulate_record_diversity!(storage, eco, times, interval, timestep,
                               divfuns::Array{Function}, q::Float64)

Deprecated: use `simulate!(RecordDiversity(storage, divfun, qs), eco, duration, timestep, every =
EveryInterval(interval))` with a [`RecordDiversity`](@ref) for the first form, and a callback on
[`simulate!`](@ref) computing the measures you want for the other two. These record on
[`simulate_action!`](@ref)'s timing, `interval` a whole multiple of `timestep`: `divfun` at the orders
`qs` into `storage`; normalised alpha, normalised beta and gamma at `qs` into `substorage` by cell and
`metastorage`, returned as `(subcommunity = substorage, metacommunity = metastorage)`; or each of
`divfuns` at the order `q` into a column of `storage`.
"""
function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    divfun::F,
                                    qs::Vector{Float64}) where {F <: Function}
    Base.depwarn("`simulate_record_diversity!` is deprecated: use " *
                 "`simulate!(RecordDiversity(storage, divfun, qs), eco, duration, timestep, " *
                 "every = EveryInterval(interval))`.",
                 :simulate_record_diversity!)
    _simulateaction!(eco, times, interval, timestep,
                     offset = iseven(size(storage, 3))) do counting
        diversity = divfun(eco, qs)[!, :diversity]
        return storage[:, :, counting] = reshape(diversity,
                                                 Int(length(diversity) /
                                                     length(qs)),
                                                 length(qs))
    end
    return storage
end

function simulate_record_diversity!(substorage::AbstractArray,
                                    metastorage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    qs::Vector{Float64})
    Base.depwarn("`simulate_record_diversity!` is deprecated: compute the measures you want in a " *
                 "callback on `simulate!`.", :simulate_record_diversity!)
    _simulateaction!(eco, times, interval, timestep,
                     offset = iseven(size(substorage, 3))) do counting
        measures = [NormalisedAlpha, NormalisedBeta, Gamma]
        for (i, msr) in enumerate(measures)
            dm = msr(eco)
            diversity = subdiv(dm, qs)[!, :diversity]
            diversity2 = metadiv(dm, qs)[!, :diversity]
            substorage[:, :, i, counting] = reshape(diversity,
                                                    Int(length(diversity) /
                                                        length(qs)),
                                                    length(qs))
            metastorage[:, i, counting] = diversity2
        end
    end
    return (subcommunity = substorage, metacommunity = metastorage)
end

function simulate_record_diversity!(storage::AbstractArray,
                                    eco::Ecosystem,
                                    times::Unitful.Time,
                                    interval::Unitful.Time,
                                    timestep::Unitful.Time,
                                    divfuns::Array{Function},
                                    q::Float64)
    Base.depwarn("`simulate_record_diversity!` is deprecated: compute the measures you want in a " *
                 "callback on `simulate!`.", :simulate_record_diversity!)
    _simulateaction!(eco, times, interval, timestep) do counting
        # `j` is a position: it addresses `storage`, allocated by `generate_storage`, as well as
        # picking the measure.
        for (j, divfun) in enumerate(divfuns)
            storage[:, j, counting] .= divfun(eco, q)[!, :diversity][1]
        end
    end
    return storage
end

"""
    simulate!(eco::Ecosystem, times::Unitful.Time, timestep::Unitful.Time,
              cacheInterval::Unitful.Time, cacheFolder::String, scenario_name::String)

Deprecated: use `simulate!(SaveAbundance(cacheFolder, scenario_name), eco, duration, timestep, every
= EveryInterval(cacheInterval))` with a [`SaveAbundance`](@ref), which writes the run's provenance
beside each file. This runs `eco` for `times` in steps of `timestep` and saves its abundances to
`<scenario_name>NN.jld2` in `cacheFolder` on the step after the clock stood on each multiple of
`cacheInterval`.
"""
function simulate!(eco::Ecosystem,
                   times::Unitful.Time,
                   timestep::Unitful.Time,
                   cacheInterval::Unitful.Time,
                   cacheFolder::String,
                   scenario_name::String)
    Base.depwarn("this six-argument `simulate!` is deprecated: use " *
                 "`simulate!(SaveAbundance(cacheFolder, scenario_name), eco, duration, timestep, " *
                 "every = EveryInterval(cacheInterval))`.", :simulate!)
    checkcoverage(eco, times, timestep)
    check_bounds(eco, times, timestep)
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

export simulate_record!, simulate_record_diversity!

# ---------------------------------------------------------------------------
# Demographic parameters: the `boost` field is gone
#
# The birth multiplier is `min(K/E, 1)`, as the model is written up: however plentiful the resource,
# a species reproduces no faster than its baseline rate. `EqualPop`, `PopGrowth` and `NoGrowth`
# lose their fifth field and `build_species` its keyword. The five-argument constructors and the
# keyword still work: each warns and discards the value. A value other than 1 changed results in the
# release that read it, and the warning says so.
#
# Deprecated in v0.8.0.
# ---------------------------------------------------------------------------
"""
    retrieve_era5(param, from_year, to_year, filename = "era5"; area = nothing)

Deprecated: fetch ERA5 monthly means through a [`EcoSISTEM.CDSRequest`](@ref) instead - one per
decade file, resolved by [`assetpath`](@ref) or on reading a `SourceSpec(ERA, code, files = ...)`
that names it. This forwards to exactly that, one request per decade of `from_year` to `to_year`,
each written to `<filename>_<decade>.nc`, and returns the paths. `param` is the variable's CDS name
(`"2m_temperature"`); the old MARS parameter codes no longer resolve.

# Arguments

  - `param`, `from_year`, `to_year`, `filename`: as above.
  - `area`: `[north, west, south, east]` in degrees, or `nothing` for the globe.
"""
function retrieve_era5(param::AbstractString, from_year::Integer,
                       to_year::Integer, filename::AbstractString = "era5";
                       area = nothing, kws...)
    Base.depwarn("`retrieve_era5` is deprecated; name the download as a `CDSRequest` entry in " *
                 "`SourceSpec(ERA, code, files = ...)`, or resolve one with `assetpath`.",
                 :retrieve_era5)
    isempty(kws) ||
        @warn "`retrieve_era5` ignores $(join(keys(kws), ", ")): the Climate Data Store's monthly " *
              "means take no such options."
    paths = String[]
    for decade in unique(fld.(from_year:to_year, 10) .* 10)
        years = filter(y -> fld(y, 10) * 10 == decade, from_year:to_year)
        request = CDSRequest(_ERA5_MONTHLY,
                             _era5request(param, years, area = area),
                             "$(filename)_$(decade).nc")
        push!(paths, assetpath(request))
    end
    return paths
end

public retrieve_era5

"""
    readfile(file::String; source = SyntheticData, unit = NoUnits, cut = nothing)

Deprecated: read a [`RasterFileSpec`](@ref) instead - `read(RasterFileSpec(file, axis =
NicheAxis, unit = unit, source = source, cut = cut))` - which is also what a layer built from the
file does. This forwards to it.

# Arguments

  - `file`: path to the raster, in any format GDAL reads.
  - `source`, `unit`, `cut`: as the spec's keywords of those names.
  - `xmin`, `xmax`, `ymin`, `ymax`: an older spelling of `cut`; all four together, or none.
"""
function readfile(file::String; source::Type = SyntheticData, unit = NoUnits,
                  cut = nothing, xmin = nothing, xmax = nothing, ymin = nothing,
                  ymax = nothing)
    Base.depwarn("`readfile(file; ...)` is deprecated; read " *
                 "`RasterFileSpec(file, axis = NicheAxis, unit = unit, source = source, " *
                 "cut = cut)` instead.", :readfile)
    n = count(!isnothing, (xmin, xmax, ymin, ymax))
    if n == 4
        isnothing(cut) ||
            error("`readfile`: pass either `cut` or the `xmin`/`xmax`/`ymin`/`ymax` extent, not both.")
        cut = Extents.Extent(Y = (ymin, ymax), X = (xmin, xmax))
    elseif n != 0
        error("`readfile` needs all four of `xmin`/`xmax`/`ymin`/`ymax` or none of them; got $n.")
    end
    return read(RasterFileSpec(file, axis = NicheAxis, unit = unit,
                               source = source, cut = cut))
end

# Warn that a `boost` was given and say what it would have done. `nothing` is the keyword's default
# and warns nothing.
_deprecatedboost(::Nothing, name::Symbol) = nothing
function _deprecatedboost(boost::Real, name::Symbol)
    tail = boost == 1 ? "" :
           " A release that read `boost` gave different results for `boost = $boost`, so " *
           "those runs will not reproduce."
    Base.depwarn("`boost` is deprecated and ignored: the birth multiplier is capped at 1, so a " *
                 "species reproduces no faster than its baseline rate however plentiful the " *
                 "resource. Drop the argument." * tail, name)
    return nothing
end

function EqualPop(birth, death, longevity, survival, boost::Real)
    _deprecatedboost(boost, :EqualPop)
    return EqualPop(birth, death, longevity, survival)
end

function PopGrowth{U}(birth::Vector{TimeUnitType{U}},
                      death::Vector{TimeUnitType{U}},
                      longevity::Float64,
                      survival::Float64,
                      boost::Real) where {U <: Unitful.Units}
    _deprecatedboost(boost, :PopGrowth)
    return PopGrowth{U}(birth, death, longevity, survival)
end

function NoGrowth{U}(birth::Vector{TimeUnitType{U}},
                     death::Vector{TimeUnitType{U}},
                     longevity::Float64,
                     survival::Float64,
                     boost::Real) where {U <: Unitful.Units}
    _deprecatedboost(boost, :NoGrowth)
    return NoGrowth{U}(birth, death, longevity, survival)
end

# ---------------------------------------------------------------------------
# `read` on the package's own source types -> `read(::RasterSpec)`
#
# `read(ERA, file, param)` and the rest named the files positionally; a `SourceSpec` names them
# with `file`, `files` or `directory` and reads the same way, so each is that spelling plus a
# warning. The CERA archive's decade time vectors, which the reader generated, are `_ceratimes`.
#
# Deprecated in v0.8.0.
# ---------------------------------------------------------------------------
function Base.read(::Type{CRUTS}, dir::AbstractString, var_name::AbstractString;
                   cut = nothing)
    Base.depwarn("`read(CRUTS, dir, var_name; cut)` is deprecated; read " *
                 "`SourceSpec(CRUTS, var_name, directory = dir)` instead.",
                 :read)
    return read(SourceSpec(CRUTS, var_name, directory = dir), cut = cut)
end

function Base.read(::Type{ERA}, file::AbstractString, param::AbstractString;
                   cut = nothing)
    Base.depwarn("`read(ERA, file, param; cut)` is deprecated; read " *
                 "`SourceSpec(ERA, param, file = file)` instead.", :read)
    return read(SourceSpec(ERA, param, file = file), cut = cut)
end

function Base.read(::Type{ERA}, file::AbstractString, param::AbstractString,
                   dim::Vector{<:Unitful.Time}; cut = nothing)
    Base.depwarn("`read(ERA, file, param, dim; cut)` is deprecated; read " *
                 "`SourceSpec(ERA, param, file = file, times = dim)` instead.",
                 :read)
    return read(SourceSpec(ERA, param, file = file, times = dim), cut = cut)
end

function Base.read(::Type{ERA}, dir::AbstractString, file::AbstractString,
                   param::AbstractString,
                   dim::Vector{<:AbstractVector{<:Unitful.Time}}; cut = nothing)
    Base.depwarn("`read(ERA, dir, file, param, dim; cut)` is deprecated; read " *
                 "`SourceSpec(ERA, param, files = matches, times = reduce(vcat, dim))` " *
                 "instead.", :read)
    files = joinpath.(dir, _searchdir(dir, file))
    return read(SourceSpec(ERA, param, files = files,
                           times = reduce(vcat, dim)),
                cut = cut)
end

function Base.read(::Type{CERA}, dir::AbstractString, file::AbstractString,
                   param::AbstractString; cut = nothing)
    Base.depwarn("`read(CERA, dir, file, param; cut)` is deprecated; read " *
                 "`SourceSpec(CERA, param, files = decades, times = ...)` instead.",
                 :read)
    files = joinpath.(dir, _searchdir(dir, file))
    return read(SourceSpec(CERA, param, files = files, times = _ceratimes()),
                cut = cut)
end

# The monthly elapsed-time coordinates the CERA reader labelled the archive with: 1901 to 2010,
# one decade of `month_mean_duration` steps per file.
function _ceratimes()
    times = collect((1901year + 1month_mean_duration):(1month_mean_duration):(1910year))
    for i in 2:12
        append!(times,
                1900year .+
                ifelse(i == 12,
                       collect(((i - 1) * 120month_mean_duration + 1month_mean_duration):(1month_mean_duration):((i - 1) * 120month_mean_duration + 1year)),
                       collect(((i - 1) * 120month_mean_duration + 1month_mean_duration):(1month_mean_duration):(i * 10year))))
    end
    return times
end

# ---------------------------------------------------------------------------
# Spec constructors: `ConstructedSpec` -> `ConstructedRasterSpec`
#
# The two compose the same way, one over rasters and one over geometry, and neither name said which
# it was, so the raster one was renamed when the vector mirror `ConstructedShapeSpec` was added.
#
# Deprecated in v0.7.0.
# ---------------------------------------------------------------------------
Base.@deprecate_binding ConstructedSpec ConstructedRasterSpec

# ---------------------------------------------------------------------------
# MPI landscape constructors: `emptyMPIgridlandscape` -> `empty_landscape`
#
# Deprecated in v0.6.0.
# ---------------------------------------------------------------------------
# `empty_mpi_gridlandscape` is replaced by `empty_landscape`, whose distributed method is chosen by
# the presence of the partition rather than by the name. This cannot be a redirecting `@deprecate`:
# the old signature took the partition alone, and the labelled views the landscape now carries need
# the habitat and species list, which it has no way to reach. So it errors, naming the replacement.
function emptyMPIgridlandscape(args...; kwargs...)
    return error("`emptyMPIgridlandscape` and `empty_mpi_gridlandscape` are replaced by " *
                 "`empty_landscape(habitat, spplist, sppcounts, sccounts)`, which takes the " *
                 "habitat and species list rather than the partition alone. The old form cannot " *
                 "be redirected: it has no way to reach the species names and grid coordinates " *
                 "the landscape now carries.")
end

function empty_mpi_gridlandscape(args...; kwargs...)
    return emptyMPIgridlandscape(args...; kwargs...)
end

export emptyMPIgridlandscape
