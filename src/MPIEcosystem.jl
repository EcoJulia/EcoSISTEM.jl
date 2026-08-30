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
  - `dimrows`: `rows_matrix` labelled with the names of the species this rank owns, against every
    cell — the counterpart of a serial landscape's `dimmatrix`.
  - `dimgrid`: the same species against the habitat's real `Y` and `X`. `rows_matrix` covers every
    cell, so the map is complete and the coordinates mean what they say.
  - `dimcols`: every species, by name, against the global indices of the cells this rank owns.
    `cols_vector` is stored one block per contributing rank, because that is what the collective
    produces; this presents those blocks as a single ordinary matrix without copying. Those cells
    are a run of the flat ordering rather than a rectangle, so there is no `Y`/`X` view of them.

The labelled views matter more here than they do serially, where the matrix simply is the whole run:
a rank holds only a slice, and these say which one. They are views of the same memory, so a write
through either is seen by the raw fields and vice versa.

The raw fields are what the simulation uses. A scalar index into `dimcols` has to find its block
first, which is markedly slower than walking `reshaped_cols`, so it is for inspection rather than
for the loop.
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
