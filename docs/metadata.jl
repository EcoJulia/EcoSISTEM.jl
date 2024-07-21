# SPDX-License-Identifier: LGPL-3.0-or-later

using Pkg

# Update EcoSISTEM folder packages 
Pkg.activate(".")
Pkg.update()

# Update examples folder packages
if isdir("examples")
    if isfile("examples/Project.toml")
        Pkg.activate("examples")
        Pkg.update()
        "EcoSISTEM" ∈ [p.name for p in values(Pkg.dependencies())] &&
            Pkg.rm("EcoSISTEM")
        Pkg.develop("EcoSISTEM")
    end
end

# Update docs folder packages
Pkg.activate("docs")
Pkg.update()
"EcoSISTEM" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("EcoSISTEM")
Pkg.develop("EcoSISTEM")

# Reformat files in package
using JuliaFormatter
using EcoSISTEM
format(EcoSISTEM)

# Carry out crosswalk for metadata
using ResearchSoftwareMetadata
ResearchSoftwareMetadata.crosswalk()
