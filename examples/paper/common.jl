# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Shared machinery for the paper's figures: the scale, the seed, where the PDFs go, timing, and
# the provenance sidecar written beside them.
#
# Included by `examples/paper.jl`, never run directly - it reads `PAPER_ISSCRIPT` from there.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols

# `small` under the test suite, which sets `ECOSISTEM_SCALE`; unset means the published
# configuration, exactly as `examples/interventions/configurations.jl` reads it.
paper_scale() = Symbol(get(ENV, "ECOSISTEM_SCALE", "large"))
paper_small() = paper_scale() === :small

# One seed for every figure; replicates use `PAPER_SEED + r`.
const PAPER_SEED = 1
const PAPER_TIMESTEP = 1month_mean_duration

# The first argument when `paper.jl` is the program; under the suite `ARGS` holds the test
# arguments, so `ECOSISTEM_PAPER_DIR` or a temporary directory is used instead. A direct run's
# temporary directory survives the process, or it would delete its own figures; the suite's is
# cleaned up.
function _outputdir(isscript::Bool)
    dir = isscript && !isempty(ARGS) ? abspath(first(ARGS)) :
          get(ENV, "ECOSISTEM_PAPER_DIR", mktempdir(cleanup = !isscript))
    mkpath(dir)
    return dir
end

const PAPER_DIR = _outputdir(PAPER_ISSCRIPT)
println("Writing the paper's figures to ", PAPER_DIR, " at ", paper_scale(),
        " scale on ", Threads.nthreads(), " threads")

# Save a figure as a PDF under the name the paper's LaTeX includes, and say where.
function savefigure(plt, name)
    path = joinpath(PAPER_DIR, name * ".pdf")
    savefig(plt, path)
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
    println("    wrote ", path)
    return path
end

# --- the island configuration figures 5 and 7 share -------------------------------------

# The published island: 100 km^2 as a 10 x 10 grid, 100 species, 100,000,000 individuals, ten
# years at monthly steps. Small scale keeps the grid and shrinks the rest.
const PAPER_CELLS = 10
const PAPER_AREA = 100.0km^2
function island_configuration()
    return paper_small() ?
           (numspecies = 10, individuals = 10_000, years = 1year) :
           (numspecies = 100, individuals = 100_000_000, years = 10year)
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
