# SPDX-License-Identifier: LGPL-3.0-or-later

using Pkg

# Update EcoSISTEM folder packages 
Pkg.activate(".")
Pkg.resolve()

# Update examples folder packages
if isdir("examples")
    if isfile("examples/Project.toml")
        Pkg.activate("examples")
        Pkg.resolve()
    end
end

# Update docs folder packages
Pkg.activate("docs")
Pkg.resolve()

# Reformat files in package
using JuliaFormatter
using EcoSISTEM
format(EcoSISTEM)

# Carry out crosswalk for metadata
using ResearchSoftwareMetadata
ResearchSoftwareMetadata.crosswalk()
