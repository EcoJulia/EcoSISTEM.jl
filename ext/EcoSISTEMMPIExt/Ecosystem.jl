# SPDX-License-Identifier: LGPL-3.0-or-later
#
# `MPIEcosystem` and the Diversity hooks it must supply.
#
# It mirrors `Ecosystem` field for field, but holds an `MPIGridLandscape` and the rank's own species
# and cell ranges. `gatherabundance`/`gatherdiversity` are what bring a distributed result back to
# one rank for recording.
#
# `_getmetaabundance`, `_getweight` and `_getordinariness!` are implementations of **Diversity's**
# generics, each imported immediately above its own definition. A call graph cannot see their callers,
# so they read as dead code and are not: deleting one silently breaks the interface for a
# distributed run.

import EcoSISTEM
using MPI
using Diversity
using HCubature
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using Missings
using Random
using Dates: Dates

"""
    MPIEcosystem{MPIGL <: MPIGridLandscape, Part <: AbstractHabitat,
                 SL <: SpeciesList, NF <: AbstractNicheFit} <: 
        AbstractEcosystem{Part, SL, NF}

MPIEcosystem houses information on species and their interaction with their
environment. It houses all information of a normal [`Ecosystem`](@ref) (see
documentation for more details), with additional fields to describe which
species are calculated on which machine. This includes: `sppcounts` - a vector
of number of species per node, `firstsp` - the identity of the first species
held by that particular node.
"""
mutable struct MPIEcosystem{MPIGL <: EcoSISTEM.MPIGridLandscape,
                            Part <: EcoSISTEM.AbstractHabitat,
                            SL <: EcoSISTEM.SpeciesList,
                            NF <: EcoSISTEM.AbstractNicheFit} <:
               EcoSISTEM.MPIEcosystem{MPIGL, Part, SL, NF}
    abundances::MPIGL
    spplist::SL
    habitat::Part
    nichefit::NF
    lookup::Vector{EcoSISTEM.Lookup}
    sppcounts::Vector{Int32}
    firstsp::Int64
    sccounts::Vector{Int32}
    firstsc::Int64
    cache::EcoSISTEM.Cache
    rngs::Vector{Random.Xoshiro}
    elapsed::typeof(1.0s)
    seed::UInt64
    # See the serial `Ecosystem` for why this is typed on the supertype. Under MPI the epoch must
    # be *identical on every rank*: layers are built redundantly per rank, so an epoch each rank
    # inferred for itself could phase their series differently. `build_ecosystem` resolves it
    # deterministically from the layers, which are the same everywhere, rather than from anything
    # rank-local.
    epoch::Union{Nothing, Dates.TimeType}

    function MPIEcosystem(abundances::MPIGL,
                          spplist::SL,
                          habitat::Part,
                          nichefit::NF,
                          lookup::Vector{EcoSISTEM.Lookup},
                          sppcounts::Vector,
                          firstsp::Int64,
                          sccounts::Vector,
                          firstsc::Int64,
                          cache::EcoSISTEM.Cache,
                          rngs::Vector{Random.Xoshiro},
                          elapsed::typeof(1.0s),
                          seed::UInt64,
                          epoch::Union{Nothing, Dates.TimeType} = nothing) where {MPIGL,
                                                                                  Part,
                                                                                  SL,
                                                                                  NF}
        EcoSISTEM._checkaligned("species tolerances, environment regime and trait nichefit",
                                EcoSISTEM._toleranceside(spplist.tolerance),
                                EcoSISTEM._regimeside(habitat.regime),
                                EcoSISTEM._nichefitside(nichefit))
        EcoSISTEM._checkaligned("species demand and environment supply",
                                EcoSISTEM._demandside(spplist.demand),
                                EcoSISTEM._supplyside(habitat.supply))
        return new{MPIGL, Part, SL, NF}(abundances,
                                        spplist,
                                        habitat,
                                        nichefit,
                                        lookup,
                                        sppcounts,
                                        firstsp,
                                        sccounts,
                                        firstsc,
                                        cache,
                                        rngs,
                                        elapsed,
                                        seed,
                                        epoch)
    end
end

EcoSISTEM.MPIEcosystem(args...; kwargs...) = MPIEcosystem(args...; kwargs...)

# With the MPI extension loaded, auto-selection (in `build_ecosystem`) treats the process as
# distributed only once MPI is initialised and there is more than one rank — so `mpirun -n 1` or a
# plain `using MPI` still builds a serial `Ecosystem`.
EcoSISTEM._should_mpi() = MPI.Initialized() && MPI.Comm_size(MPI.COMM_WORLD) > 1

using EcoSISTEM: getkernels, genlookups, numdemands, _checksimulatable
"""
    MPIEcosystem(spplist::SpeciesList, habitat::GridHabitat,
                 nichefit::AbstractNicheFit)

Create an `MPIEcosystem` given a species list, an abiotic environment and trait
nichefit.
"""
function MPIEcosystem(popfun::F,
                      spplist::EcoSISTEM.SpeciesList{T, DM},
                      habitat::EcoSISTEM.GridHabitat,
                      nichefit;
                      seed::Integer = rand(UInt64)) where {F <: Function, T,
                                                           DM}
    # **First, exactly as in the serial `Ecosystem` constructor and for the same reason.**
    # `genlookups` below divides the regime's cell size by a dispersal distance, so on a geographic
    # grid (`size` in `°`) it raises `DimensionError: ° km⁻¹` before any tailored message is reached.
    # **This guard was missing here.** It was added to the serial constructor on 2026-08-12, when
    # `[GEO-SIZE]` made a geographic cell size honest and exposed the ordering; the MPI constructor
    # has the identical fault and was overlooked. Nothing caught it because `test/SmallMPItest.jl`
    # builds a synthetic grid, so the geographic case is never exercised under MPI.
    _checksimulatable(habitat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    totalsize = MPI.Comm_size(comm)
    numspp = length(spplist.names)
    numsc = countsubcommunities(habitat.regime)

    # One deterministically-seeded RNG per global species, built identically on
    # every rank, so that species draws are reproducible regardless of how
    # species and cells are partitioned across processes and threads
    rngs = EcoSISTEM.makerngs(seed, numspp)

    count = div(numspp, totalsize)
    sppcounts = Int32.(fill(count, totalsize))
    sppcounts[1:(numspp - sum(sppcounts))] .+= 1
    sppindices = vcat([0], cumsum(sppcounts))
    firstsp = sppindices[rank + 1] + 1

    sccount = div(numsc, totalsize)
    sccounts = Int32.(fill(sccount, totalsize))
    sccounts[1:(numsc - sum(sccounts))] .+= 1
    scindices = vcat([0], cumsum(sccounts))
    firstsc = scindices[rank + 1] + 1

    # Create matrix landscape of zero abundances
    ml = EcoSISTEM.empty_mpi_gridlandscape(sppcounts, sccounts)

    # Populate this matrix with species abundances
    popfun(ml, spplist, habitat, nichefit, rngs)

    rankspp = firstsp:sppindices[rank + 2]
    lookup_tab = collect(map(k -> genlookups(habitat.regime, k),
                             @view getkernels(spplist.movement)[rankspp]))
    nm = zeros(Int64, (sppcounts[rank + 1], numsc))
    totaldemand = zeros(Float64, (numsc, numdemands(DM)))
    return MPIEcosystem(ml,
                        spplist,
                        habitat,
                        nichefit,
                        lookup_tab,
                        sppcounts,
                        firstsp,
                        sccounts,
                        firstsc,
                        EcoSISTEM.Cache(nm, totaldemand, false),
                        rngs,
                        0.0s,
                        EcoSISTEM._storedseed(seed))
end

function MPIEcosystem(spplist::EcoSISTEM.SpeciesList,
                      habitat::EcoSISTEM.GridHabitat, nichefit;
                      seed::Integer = rand(UInt64))
    return MPIEcosystem(EcoSISTEM.populate!, spplist, habitat, nichefit;
                        seed = seed)
end
@doc (@doc MPIEcosystem) MPIEcosystem(::EcoSISTEM.SpeciesList,
                                      ::EcoSISTEM.GridHabitat,
                                      ::Any)

function EcoSISTEM.gatherabundance(eco::MPIEcosystem)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    true_abuns = zeros(Int64, counttypes(eco), countsubcommunities(eco))
    if rank == 0
        # **Counted over the COLUMN partition, because that is what is being sent.**
        # `reshaped_cols` holds *every* species for *this rank's* cells, so a rank contributes
        # `counttypes(eco) * sccounts[rank]` values. The previous
        # `sppcounts .* sum(sccounts)` — this rank's *species* by *all* cells — is the row
        # partition, and the two agree only when species and cells both divide evenly across the
        # ranks. They always did in the old test fixture (8 species on a 4 × 4 grid), so the
        # mismatch never showed; with 7 species on 77 cells rank 0 sends 273 values while the old
        # expression asked for 308, and `MPI_Gatherv` fails outright.
        output_vbuf = VBuffer(true_abuns,
                              Int32.(counttypes(eco) .* eco.sccounts))
    else
        output_vbuf = VBuffer(nothing)
    end
    MPI.Gatherv!(vcat(eco.abundances.reshaped_cols...)[1:end], output_vbuf, 0,
                 comm)
    return true_abuns
end

function EcoSISTEM.gatherdiversity(eco::MPIEcosystem, divmeasure::F,
                                   q) where {F <: Function}
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    totalsize = MPI.Comm_size(comm)
    diversity = divmeasure(eco, q)
    totalabun = MPI.Gather(sum(eco.abundances.rows_matrix), 0, comm)
    mpidivs = MPI.Gather(diversity[!, :diversity], 0, comm)
    if rank == 0
        mpidivs = vcat(reshape(mpidivs, countsubcommunities(eco), totalsize),
                       diversity[!, :q])
        diversity[!, :diversity] .= mapslices(x -> Diversity.powermean(x[:,
                                                                         1:(end - 1)],
                                                                       1 .-
                                                                       x[:,
                                                                         end],
                                                                       totalabun .*
                                                                       1.0),
                                              mpidivs,
                                              dims = 2)[:,
                                                        1]
        return diversity
    else
        return diversity
    end
end
