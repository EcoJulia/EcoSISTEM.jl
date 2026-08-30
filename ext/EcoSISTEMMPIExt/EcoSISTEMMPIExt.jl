# SPDX-License-Identifier: LGPL-3.0-or-later

# EcoSISTEMMPIExt
#
# The distributed `MPIEcosystem`, partitioned across ranks by species (rows) and grid cells
# (columns).
#
# **This module is a manifest** — the parts are ordinary files it includes, laid out the way
# `src/EcoSISTEM.jl` is.
#
# **The parts must live here rather than under `src/`.** `test/core_test.jl` walks `src/`
# and reports every file with no matching `test_*.jl`, and such a test cannot load a weak dependency
# — so anything MPI-specific placed under `src/` would sit permanently on its "potentially missing
# tests" list with no way off. Extension code lives with the extension;
# `test/ext_EcoSISTEMMPIExt.jl` is what covers these.
#
# **The filenames mirror `src/` deliberately.** `Landscape.jl`, `Ecosystem.jl`,
# `DiversityInterface.jl` and `dynamics.jl` are each named for the serial file they duplicate or
# conform alongside, so the two can be diffed by name. `dynamics.jl` was called `generate.jl` until
# 2026-08-30, after a `src/Generate.jl` that no longer exists -- which is the kind of drift that let
# a four-year-old divergence in the hot loop go unnoticed.

module EcoSISTEMMPIExt

include("Landscape.jl")

include("Ecosystem.jl")

# After `Ecosystem.jl`: these methods dispatch on the concrete `MPIEcosystem` declared there, and a
# `using` is module-scoped, so this file leans on the imports its neighbours have already run.
include("DiversityInterface.jl")

include("dynamics.jl")

end
