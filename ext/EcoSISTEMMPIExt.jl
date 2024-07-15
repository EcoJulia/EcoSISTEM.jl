# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEMMPIExt

@info "Creating MPI interface for EcoSISTEM..."

include("../src/MPILandscape.jl")

include("../src/MPIEcosystem.jl")

include("../src/MPIGenerate.jl")

end
