# SPDX-License-Identifier: LGPL-3.0-or-later

using Pkg
using EcoSISTEM

# 🔴 Every path here is resolved from the package root, NOT from the working directory. This script
# used to say `Pkg.activate(".")`, which was only ever correct when run from the root - and after it
# moved into `data/src`, running it from beside itself would have activated `data/src` and resolved
# the wrong environment entirely.
const ROOT = pkgdir(EcoSISTEM)

# Update EcoSISTEM folder packages
Pkg.activate(ROOT)
Pkg.resolve()

# Update examples folder packages
if isfile(joinpath(ROOT, "examples", "Project.toml"))
    Pkg.activate(joinpath(ROOT, "examples"))
    Pkg.resolve()
end

# Update docs folder packages
Pkg.activate(joinpath(ROOT, "docs"))
Pkg.resolve()

# Update the generator environment these scripts share
if isfile(joinpath(ROOT, "data", "src", "Project.toml"))
    Pkg.activate(joinpath(ROOT, "data", "src"))
    Pkg.resolve()
end

# Reformat files in package
using JuliaFormatter
format(EcoSISTEM)

# Carry out crosswalk for metadata
using ResearchSoftwareMetadata
ResearchSoftwareMetadata.crosswalk()
