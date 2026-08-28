# SPDX-License-Identifier: LGPL-3.0-or-later

# EcoSISTEMMPIExt
#
# The distributed `MPIEcosystem`, partitioned across ranks by species (rows) and grid cells
# (columns).
#
# **This module is a manifest** — the parts are ordinary files it includes, laid out the way
# `src/EcoSISTEM.jl` is.
#
# **The three parts must live here rather than under `src/`.** `test/core_test.jl` walks `src/`
# and reports every file with no matching `test_*.jl`, and such a test cannot load a weak dependency
# — so anything MPI-specific placed under `src/` would sit permanently on its "potentially missing
# tests" list with no way off. Extension code lives with the extension;
# `test/ext_EcoSISTEMMPIExt.jl` is what covers these.

module EcoSISTEMMPIExt

include("landscape.jl")

include("ecosystem.jl")

include("generate.jl")

end
