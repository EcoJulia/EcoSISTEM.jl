# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The distributed forms. Abstract here because their concrete subtypes live in `EcoSISTEMMPIExt`,
# and these stubs are where those subtypes are documented: a docstring written inside an extension
# does not reach the API reference, so the fields below describe what the extension actually builds.

"""
    MPIGridLandscape

An [`MPIEcosystem`](@ref)'s abundances, distributed across MPI ranks; `EcoSISTEMMPIExt` supplies the
concrete `MPIGridLandscape{RA, NT}`.

The same abundances are held **twice**, partitioned two different ways, because the simulation needs
a different split at different moments — demographics is per species, dispersal is per cell. The
synchronise steps copy between the two views.

# Fields

  - `rows_matrix`: this rank's species over *all* grid cells.
  - `cols_vector`: *all* species over this rank's grid cells, flattened.
  - `reshaped_cols`: `cols_vector` seen as one species-by-cell view per MPI block, which is what the
    collective operations address.
  - `rows_tuple`, `cols_tuple`: how each view is partitioned — `total`, `first`, `last`, and the
    per-rank `counts`.
"""
abstract type MPIGridLandscape end

"""
    MPIEcosystem{MPIGL, Part, SL, NF} <: AbstractEcosystem{Part, SL, NF}

An [`Ecosystem`](@ref) whose abundances are distributed across MPI ranks; `EcoSISTEMMPIExt` supplies
the concrete `MPIEcosystem`.

Only the abundances are split. The habitat, species list and nichefit are held whole on every rank,
and four extra fields record which slice of the work this rank owns.

# Fields

  - `abundances`: the distributed [`MPIGridLandscape`](@ref), in place of a `GridLandscape`.
  - `sppcounts`, `firstsp`: how many species each rank holds, and the first one this rank holds.
  - `sccounts`, `firstsc`: the same two, for grid cells.
  - `spplist`, `habitat`, `nichefit`, `lookup`, `cache`, `rngs`, `elapsed`, `seed`, `epoch`: as
    [`Ecosystem`](@ref), and identical on every rank. The `epoch` must be, or layers built
    redundantly per rank could phase their series differently.

# Type parameters

  - `MPIGL`: the [`MPIGridLandscape`](@ref) type holding the abundances.
  - `Part`, `SL`, `NF`: habitat, species list and nichefit, as [`AbstractEcosystem`](@ref).
"""
abstract type MPIEcosystem{MPIGL <: MPIGridLandscape,
                           Part <: AbstractHabitat,
                           SL <: SpeciesList,
                           NF <: AbstractNicheFit} <:
              AbstractEcosystem{Part, SL, NF} end
