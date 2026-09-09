# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The paper's computational figures, regenerated from the package: figures 3, 5, 6, 7, 8 and 9 of
# the EcoSISTEM paper, written as PDFs under the filenames its LaTeX source already includes, so
# that what the paper shows is produced by code under test.
#
# **Two scales.** Run directly it reproduces the published experiments in full, which takes hours
# (figure 8 alone is 50,000 species on a 100 x 100 continent for 200 years; time one year first,
# below). Under the test suite, which sets `ECOSISTEM_SCALE=small`, every figure is drawn from a
# run that finishes in seconds - the same code, the same layout, a fraction of the species and the
# years - so the figures are exercised like every other top-level example. Unlike `models.jl`, this
# file never skips.
#
#     julia -t 12 --project=examples examples/paper.jl [output directory]
#
# The threads are the package's own: each simulation runs multithreaded, and the figures follow
# one another.
#
# The output directory defaults to `ECOSISTEM_PAPER_DIR`, else a temporary directory; it is printed
# at the start and the end. `ECOSISTEM_PAPER_FIGURES=NicheWidth2,Abundance` draws only the figures
# named, so the hours-long continent can be left for a machine that has them. Beside the PDFs it writes `figures.toml`, recording the package version,
# commit, scale, seed, thread count and the wall time of each figure, and `continent.jld2`, the
# per-cell maps figure 8 is drawn from, so that figure can be redrawn without rerunning.
#
# Everything it needs lives in `examples/paper/`:
#
#   - `paper/common.jl`     - scale, seed, output directory, saving, timing, the sidecar
#   - `paper/niche.jl`      - figure 3, the analytic niche and demand curves
#   - `paper/island.jl`     - figure 5, the island ecosystems
#   - `paper/dispersal.jl`  - figure 6, dispersal distance
#   - `paper/optima.jl`     - figure 7, temperature optima and niche widths
#   - `paper/continent.jl`  - figures 8 and 9, the circular continent and the invasions
#
# `models/ecosystems.jl` supplies the island builders (`uniform_environment`, `community`,
# `percell`, `perspecies`); the rest is here. Case study 2 (the GBIF and ECMWF Africa runs) is not
# part of this file.
#
# **A module, deliberately.** `test/extras_examples.jl` includes every top-level example into one
# module, and `examples/models.jl` defines the same builder names this file reuses. Without the
# wrapper the two would redefine each other whenever both ran.

module Paper

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Distributions
using Random
using Plots
using JLD2

gr()

# Whether this file is the program being run, as opposed to included by the suite - decided here,
# where `@__FILE__` is this file, and read by `common.jl` to know whether `ARGS` is ours.
const PAPER_ISSCRIPT = abspath(PROGRAM_FILE) == @__FILE__

include(joinpath(@__DIR__, "models", "ecosystems.jl"))
include(joinpath(@__DIR__, "paper", "common.jl"))
include(joinpath(@__DIR__, "paper", "niche.jl"))
include(joinpath(@__DIR__, "paper", "optima.jl"))
include(joinpath(@__DIR__, "paper", "island.jl"))
include(joinpath(@__DIR__, "paper", "dispersal.jl"))
include(joinpath(@__DIR__, "paper", "continent.jl"))

write_sidecar()
println("Done - figures are in ", PAPER_DIR)

end
