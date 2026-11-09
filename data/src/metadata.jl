# SPDX-License-Identifier: LGPL-3.0-or-later

using Pkg
using EcoSISTEM

# 🔴 Every path here is resolved from the package root, NOT from the working directory:
# `Pkg.activate(".")` is right only when run from the root, and run from beside this script it
# activates `data/src` and resolves the wrong environment entirely.
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
