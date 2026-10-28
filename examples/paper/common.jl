# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Shared machinery for the paper's figures: the scale, the seed, the parameters, where the PDFs go,
# timing, Tables S1 and S2, and the provenance sidecar written beside them.
#
# Included by `examples/paper.jl`, never run directly - it reads `PAPER_ISSCRIPT` from there.

using EcoSISTEM
using EcoSISTEM.Units
using TOML: TOML
using Unitful
using Unitful.DefaultSymbols

# `small` under the test suite, which sets `ECOSISTEM_SCALE`; unset means the published
# configuration, exactly as `examples/interventions/configurations.jl` reads it.
paper_scale() = Symbol(get(ENV, "ECOSISTEM_SCALE", "large"))
paper_small() = paper_scale() === :small

# One seed for every figure; replicates use `PAPER_SEED + r`.
const PAPER_SEED = 1
const PAPER_TIMESTEP = 1month_mean_duration

# The parameters file: `parameters.toml` beside this one, or, in a data pipeline run, the copy
# `paper/pipeline.jl` read from the registry as `PAPER_PARAMETERS`.
function _parametersfile()
    return isdefined(Main, :PAPER_PARAMETERS) ? Main.PAPER_PARAMETERS :
           joinpath(@__DIR__, "parameters.toml")
end

# A parameter's value: a string is a quantity, parsed with the package's units as well as Unitful's;
# a table is a named tuple of its entries.
_parameter(x) = x
function _parameter(s::AbstractString)
    return uparse(s, unit_context = [Unitful, EcoSISTEM.Units])
end
_parameter(v::AbstractVector) = map(_parameter, v)
function _parameter(d::AbstractDict)
    return (; (Symbol(k) => _parameter(v) for (k, v) in d)...)
end

# The set of parameters for this scale, every figure's and both tables'.
function paper_parameters(path = _parametersfile())
    sets = TOML.parsefile(path)
    scale = string(paper_scale())
    haskey(sets, scale) ||
        error("`$path` has no `[$scale]` set of parameters for `ECOSISTEM_SCALE=$scale`")
    return _parameter(sets[scale])
end

const PAPER = paper_parameters()

# Launched as a script - `paper.jl`, or `paper/pipeline.jl` through the data pipeline - the first
# argument, else `outputs/` in this directory, which its `.gitignore` keeps out of the repository.
# Included by anything else, as the suite includes `paper.jl` with `ARGS` of its own, a temporary
# directory cleaned up at exit, so a test run never replaces a full run's figures.
function _outputdir(isscript::Bool)
    dir = !isscript ? mktempdir() :
          isempty(ARGS) ? joinpath(@__DIR__, "outputs") : abspath(first(ARGS))
    mkpath(dir)
    return dir
end

const PAPER_DIR = _outputdir(PAPER_ISSCRIPT)
println("Writing the paper's figures to ", PAPER_DIR, " at ", paper_scale(),
        " scale on ", Threads.nthreads(), " threads")

# Every output written, by its kind - `:figure`, `:table` or `:config` - and its name, with its
# path: what `paper/pipeline.jl` registers once everything is drawn.
const PAPER_OUTPUTS = @NamedTuple{kind::Symbol, name::String, path::String}[]

# Save a figure as a PDF under the name the paper's LaTeX includes, and say where.
function savefigure(plt, name)
    path = joinpath(PAPER_DIR, name * ".pdf")
    savefig(plt, path)
    push!(PAPER_OUTPUTS, (kind = :figure, name = name, path = path))
    println("    wrote ", path)
    return path
end

# Wall time per figure, so a full-scale run can be budgeted; written to the sidecar at the end.
const PAPER_TIMINGS = Pair{String, Float64}[]

# Which figures to draw: every one unless `ECOSISTEM_PAPER_FIGURES` names some, comma-separated by
# their file names, so a full-scale run can leave the continent for a machine that has the hours.
function paper_wanted(name)
    wanted = get(ENV, "ECOSISTEM_PAPER_FIGURES", "")
    return isempty(wanted) || name in strip.(split(wanted, ','))
end

# Run `f` as the work for figure `name`, reporting how long it took.
function timed(f, name)
    paper_wanted(name) || return nothing
    println("Figure ", name, " ...")
    elapsed = @elapsed f()
    push!(PAPER_TIMINGS, name => elapsed)
    println("    ", name, " took ", round(elapsed, digits = 1), " s")
    return nothing
end

# The commit the figures were made from, when the tree is a git checkout.
function _commit()
    dir = pkgdir(EcoSISTEM)
    try
        return readchomp(pipeline(`git -C $dir rev-parse HEAD`,
                                  stderr = devnull))
    catch
        return "unknown"
    end
end

# `figures.toml`: what produced the PDFs beside it.
function write_sidecar()
    path = joinpath(PAPER_DIR, "figures.toml")
    open(path, "w") do io
        println(io, "package = \"EcoSISTEM\"")
        println(io, "version = \"", pkgversion(EcoSISTEM), "\"")
        println(io, "commit = \"", _commit(), "\"")
        println(io, "julia = \"", VERSION, "\"")
        println(io, "scale = \"", paper_scale(), "\"")
        println(io, "seed = ", PAPER_SEED)
        println(io, "threads = ", Threads.nthreads())
        println(io)
        println(io, "[seconds]")
        for (name, elapsed) in PAPER_TIMINGS
            println(io, name, " = ", round(elapsed, digits = 1))
        end
    end
    push!(PAPER_OUTPUTS, (kind = :config, name = "Config", path = path))
    println("    wrote ", path)
    return path
end

# --- the island figures 5, 6 and 7 share -------------------------------------------------

# The island's temperature, in kelvin, since it is also the optimum of every species matched to it,
# and a tolerance's parameter vectors must share one unit with its widths.
island_temperature() = uconvert(K, PAPER.island.temperature)

# The island's environment, `cells` a side over `area`, uniform at `temperature` with the island's
# supplies.
function island_environment(cells = PAPER.island.cells,
                            area = PAPER.island.area;
                            temperature = island_temperature(),
                            topology = Torus())
    island = PAPER.island
    return uniform_environment(cells, area, temperature = temperature,
                               supply = (sunlight = island.sunlight_supply,
                                         water = island.water_supply),
                               topology = topology)
end

# A community on `env` with the island's demographic and demand values, each replaceable by a
# keyword of `community`'s.
function island_community(env, opts, widths; kw...)
    island = PAPER.island
    defaults = (dispersal = island.dispersal, individuals = island.individuals,
                birth = island.birth, death = island.death,
                sizes = island.size, longevity = island.longevity,
                survival = island.survival,
                demand = (sunlight = island.sunlight_demand,
                          water = island.water_demand),
                seed = PAPER_SEED)
    return community(env, opts, widths; merge(defaults, values(kw))...)
end

# A histogram of individuals over a trait: each species contributes its abundance at its trait
# value, binned by `edges` as `[edges[i], edges[i + 1])` with the last bin closed. A weighted count
# rather than one entry per individual, which at 100 million individuals would not fit.
function weighted_histogram(values, weights, edges)
    counts = zeros(Float64, length(edges) - 1)
    for (v, w) in zip(values, weights)
        i = v == edges[end] ? length(edges) - 1 : searchsortedlast(edges, v)
        1 <= i < length(edges) && (counts[i] += w)
    end
    return counts
end

# --- Tables S1 and S2, from the parameters the figures were built with ---------------------------

# A temperature's unit as the tables print it, for a temperature and for a difference of two (a
# niche width), which are the same size in degrees C as in kelvin.
const _DEGREES = raw"$^{\circ}$C"

# A number as the tables print it: a whole number below a billion in full with its thousands
# separated (`450,000`), anything else from 10^5 up as three significant figures times a power of
# ten in LaTeX (`5.03$\times 10^{7}$`), and the rest to six significant figures, trailing zeros
# dropped, so that a temperature converted to degrees C reads `25` rather than `24.999999999999972`.
function _latexnumber(x::Real)
    y = round(Float64(x), sigdigits = 6)
    if isinteger(x) && abs(x) < 1e9
        return replace(string(round(Int, x)), r"(?<=\d)(?=(\d{3})+$)" => ",")
    elseif abs(y) >= 1e5
        e = floor(Int, log10(abs(y)))
        return _trimmed(round(y / 10.0^e, sigdigits = 3)) * raw"$\times 10^{" *
               string(e) * raw"}$"
    end
    return _trimmed(y)
end

# A float without a trailing `.0`.
_trimmed(y) = isinteger(y) ? string(round(Int, y)) : string(y)

# `x` in `unit`, followed by `label`, its LaTeX.
_latex(x, unit, label) = _latexnumber(ustrip(unit, x)) * label

# The range of `values` in `unit` followed by `label`, or the one value where they all agree.
function _latexspan(values, unit, label)
    low, high = extrema(ustrip.(unit, values))
    span = low == high ? _latexnumber(low) :
           _latexnumber(low) * "-" * _latexnumber(high)
    return span * label
end

# The end of a table's row, ` \\`: in a raw string a backslash pair before the closing quote is one
# backslash, so it is written escaped.
const _ROWEND = " \\\\"

# A table in the form the paper prescribes for Tables S1 and S2: a complete `tabular`, its rows -
# pairs of a parameter's name and its value, both raw LaTeX - in a species group and an abiotic
# environment group.
function _tabular(species, environment)
    io = IOBuffer()
    println(io, raw"\begin{tabular}{p{20em} p{12em}}")
    println(io, raw"\toprule")
    println(io, raw"\textbf{Parameter name} & \textbf{Value}", _ROWEND)
    for (group, rows) in (("Species", species),
        ("Abiotic Environment", environment))
        println(io, raw"\midrule")
        println(io, raw"\multicolumn{2}{l}{\textbf{", group, "}}", _ROWEND)
        println(io, raw"\midrule")
        for (name, value) in rows
            println(io, name, " & ", value, _ROWEND)
        end
    end
    println(io, raw"\bottomrule")
    println(io, raw"\end{tabular}")
    return String(take!(io))
end

# Table S1, the island and patch dynamics: the island's species and environment, each range spanning
# every island experiment - figure 5's pools, figure 6's distances, figure 7's optima and widths.
function island_table(p = PAPER)
    island, pool, optima = p.island, p.testing, p.optima
    temperature = uconvert(K, island.temperature)
    furthest = pool.pool_spread * maximum(pool.pool_width)
    preferences = [temperature, temperature - optima.spread,
        temperature + optima.spread, temperature - furthest,
        temperature + furthest]
    widths = [island.width, optima.widths..., pool.pool_width...]
    dispersals = [island.dispersal, p.dispersal.distances...,
        pool.pool_dispersal...]
    month = month_mean_duration
    species = [
        "Number of species" => _latexnumber(island.species),
        "Initial number of individuals" => _latexnumber(island.individuals),
        raw"Sunlight demand per individual, $\epsilon_{1,s}$" =>
            _latex(island.sunlight_demand *
                   island.size,
                   kJ / month,
                   " kJ/month"),
        raw"Water demand per individual, $\epsilon_{2,s}$" =>
            _latex(island.water_demand *
                   island.size,
                   Unitful.L / month,
                   " L/month"),
        raw"Temperature preference, $\bar{T_s}$" =>
            _latexspan(preferences, °C,
                       _DEGREES),
        raw"Niche width, $\psi_s$" => _latexspan(widths, K, _DEGREES),
        raw"Average dispersal distance, $\mu_s$" => _latexspan(dispersals, km,
                   " km"),
        raw"Birth rate, $b_s$" => _latex(island.birth, year^-1, "/year"),
        raw"Death rate, $d_s$" => _latex(island.death, year^-1, "/year"),
        raw"Longevity, $\lambda$" => string(island.longevity),
        raw"Survival, $\tau$" => string(island.survival)]
    environment = [
        "Area" => _latex(island.area, km^2, raw" km$^2$"),
        raw"Grid cell width, $\eta$" =>
            _latex(sqrt(island.area) / island.cells,
                   km, " km"),
        raw"Total ecosystem sunlight supply, $K_1$" =>
            _latex(island.sunlight_supply,
                   kJ / (km^2 * month),
                   raw" kJ/km$^2$/month"),
        raw"Total ecosystem water supply, $K_2$" =>
            _latex(island.water_supply,
                   mm / month,
                   " mm/month"),
        raw"Total ecosystem temperature, $T$" => _latex(temperature, °C,
               _DEGREES)]
    return _tabular(species, environment)
end

# Table S2, the large scale dynamics: the continent's generalists and specialist and its
# environment, its area that of the circle the species live on.
function continent_table(p = PAPER)
    continent = p.continent
    temperature = uconvert(K, continent.temperature)
    widths = _latex(continent.generalist_width, K, _DEGREES) *
             " (generalists), " *
             _latex(continent.specialist_width, K, _DEGREES) * " (invasive; " *
             _latexspan(continent.widths, K, _DEGREES) *
             raw" in Figure \ref{fig:invasion})"
    species = [
        "Number of species" => _latexnumber(continent.species),
        "Initial number of individuals" => _latex(continent.supply /
               continent.demand, m^-2,
               raw" m$^{-2}$"),
        raw"Sunlight demand per individual, $\epsilon_{1,s}$" =>
            _latex(continent.demand,
                   MJ / day,
                   " MJ/day"),
        raw"Temperature preference, $\bar{T_s}$" => _latex(temperature, °C,
               _DEGREES),
        raw"Niche width, $\psi_s$" => widths,
        raw"Average dispersal distance, $\mu_s$" =>
            _latex(continent.dispersal,
                   km, " km"),
        raw"Birth rate, $b_s$" => _latex(continent.birth, year^-1, "/year"),
        raw"Death rate, $d_s$" => _latex(continent.death, year^-1, "/year"),
        raw"Longevity, $\lambda$" => string(continent.longevity),
        raw"Survival, $\tau$" => string(continent.survival)]
    environment = [
        "Area" => _latex(pi * continent.radius^2, km^2, raw" km$^2$"),
        raw"Grid cell width, $\eta$" => _latex(continent.cellsize, km, " km"),
        raw"Total ecosystem sunlight supply, $K_1$" =>
            _latex(continent.supply,
                   MJ / (m^2 * day),
                   raw" MJ/m$^2$/day"),
        raw"Total ecosystem temperature, $T$" => _latex(temperature, °C,
               _DEGREES)]
    return _tabular(species, environment)
end

# Write Tables S1 and S2 into the output directory, each a `tabular` the paper `\input`s, and
# record them.
function write_tables()
    for (name, table) in (("IslandParameters", island_table()),
                          ("ContinentParameters", continent_table()))
        path = joinpath(PAPER_DIR, name * ".tex")
        write(path, table)
        push!(PAPER_OUTPUTS, (kind = :table, name = name, path = path))
        println("    wrote ", path)
    end
    return nothing
end
