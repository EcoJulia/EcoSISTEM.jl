# SPDX-License-Identifier: LGPL-3.0-or-later

# EcoSISTEMPhyloExt
#
# The phylogenetics half of EcoSISTEM: everything that needs `Phylo`.
#
# **This module is a manifest** - the parts are ordinary files it includes, laid out the way
# `src/EcoSISTEM.jl` is, and in the shape `EcoSISTEMRasterDataSourcesExt` settled on first.
#
# **No docstring belongs in here.** `docs/src/api.md` is an `@autodocs` block over the package's
# own modules and cannot see inside an extension, so every public name whose implementation moved
# keeps its docstring on a method-less stub in the parent. `test/core_ext.jl` asserts it.

module EcoSISTEMPhyloExt

# **The phylogenetic half of EcoSISTEM, active only when `Phylo` is loaded.**
#
# **No docstrings live here, and that is not an oversight.** `docs/src/api.md` is one `@autodocs`
# block over `Modules = [EcoSISTEM, EcoSISTEM.ClimatePref, EcoSISTEM.Units]`; an extension module is
# not among them and cannot easily be, so a docstring written here would never reach the manual. Every
# function below is declared - and documented - method-less in the parent, and this file supplies the
# sole method. That is the same pattern `retrieve_era5` and `unziptemp` already use.
#
# An extension **adds** methods, never overwrites one: the parent declares `function foo end` and
# nothing more, so there is no method here to collide with and no precompilation warning.

using EcoSISTEM
# Named explicitly: every abstract type in this package is `public` rather than exported, so
# `using EcoSISTEM` does not bring them in - and the constructors moved here name several.
using EcoSISTEM: ClimatePref, AbstractDemand, AbstractMovement, AbstractParams,
                 SimpleCategoricalTolerance, SpeciesList
# The private helpers too: the moved bodies call each other by bare name, and a method defined as
# `EcoSISTEM.foo(...)` does not bind `foo` inside this module.
using EcoSISTEM: arenoderecordsempty, pair, root_to_tips, brownian_motion
# The phylogenetic trait set is `public` rather than exported (`[C7-VIS]`), so `using EcoSISTEM`
# above does not reach it - and several bodies here call these by bare name. The rest of the set
# (`resettraits!`, `reroot!`, `discrete_evolve`, `continuous_evolve`) is only ever *defined* here, with
# a qualified name, which binds nothing locally and so needs no import.
# `varcovar` is `public`, not exported, so `using EcoSISTEM` above does not reach it and the two
# fitters' unqualified calls to it would be an `UndefVarError`. Naming it here is the same fix the
# line above makes for `gettraits`/`assigntraits!`.
using EcoSISTEM: gettraits, assigntraits!, varcovar
using Phylo
using Phylo: NodeNameIterator
using Optim
using Calculus
using DataFrames
using Distributions
using LinearAlgebra
using Statistics
using Unitful
using Diversity
# `AbstractTypes` and `PhyloBranches` are `Diversity`'s and are not exported by it - the
# similarity structure a phylogenetic `SpeciesList` holds. `PhyloBranches`' concrete subtype lives
# in `DiversityPhyloExt`, which `using Phylo` above is what activates: the same abstract-parent
# pattern this extension is built on, one package down.
using Diversity: AbstractTypes, PhyloBranches

include("traits.jl")

include("species_list.jl")

include("models.jl")

include("diversity.jl")

end
