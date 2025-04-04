# SPDX-License-Identifier: LGPL-3.0-or-later

using Diversity
using HCubature
using DataFrames
using Unitful
using EcoSISTEM.Units
using Missings
using RecipesBase

using Diversity.API: _calcabundance

using AMDGPU

"""
    Lookup

Lookup houses information on `x`, `y` grid locations and the probability of
occurrence at the location for the species in question `p`. `pnew` and `moves`
are initially empty storage and written over by the movement step in update!().
`pnew` is the recalculated probability based on which directions are available
and `moves` is the number of moves to that grid location in that step.
"""
mutable struct AMDLookup{IV, FV}
    x::IV
    y::IV
    p::FV
    pnew::FV
    moves::IV
end

"""
    AMDCache

AMDCache houses an integer array of moves made by all species in a timestep for the
update! function, `netmigration`.
"""
mutable struct AMDCache
    netmigration::ROCArray{Int64, 2, AMDGPU.Runtime.Mem.HIPBuffer}
    totalE::ROCArray{Float64, 2, AMDGPU.Runtime.Mem.HIPBuffer}
    exchanged_energy::ROCArray{Float64, 2, AMDGPU.Runtime.Mem.HIPBuffer}
    active_sc::ROCArray{Bool, 2, AMDGPU.Runtime.Mem.HIPBuffer}
    species_birth_rate::ROCArray{Float64, 1, AMDGPU.Runtime.Mem.HIPBuffer}
    species_death_rate::ROCArray{Float64, 1, AMDGPU.Runtime.Mem.HIPBuffer}
    width::Int64
    height::Int64
    valid::Bool
end

"""
    AMDEcosystem{Part <: AbstractAbiotic} <:
       AbstractEcosystem{Part, SL, TR}

AMDEcosystem houses information on species and their interaction with their
environment. For species, it holds abundances and locations, `abundances`,
as well as properties such as trait information, `spplist`, and movement types,
`lookup`. For environments, it provides information on environmental conditions
and available resources,`abenv`. Finally, there is a slot for the relationship
between the environment and the characteristics of the species, `relationship`.
"""
mutable struct AMDEcosystem{Part <: AbstractAbiotic, SL <: SpeciesList,
                         TR <: AbstractTraitRelationship} <:
               AbstractEcosystem{Part, SL, TR}
    abundances::AMDGridLandscape
    spplist::SL
    abenv::Part
    ordinariness::Union{Matrix{Float64}, Missing}
    relationship::TR
    lookup::Vector{Lookup}
    roclookup::Vector{AMDLookup{ROCVector{Int64}, ROCVector{Float64}}}
    cache::AMDCache

    function AMDEcosystem{Part, SL, TR}(abundances::AMDGridLandscape,
                                     spplist::SL, abenv::Part,
                                     ordinariness::Union{Matrix{Float64},
                                                         Missing},
                                     relationship::TR, lookup::Vector{Lookup},
                                     cache::AMDCache) where {
                                                          Part <:
                                                          AbstractAbiotic,
                                                          SL <: SpeciesList,
                                                          TR <:
                                                          AbstractTraitRelationship
                                                          }
        tematch(spplist, abenv) || error("Traits do not match habitats")
        trmatch(spplist, relationship) ||
            error("Traits do not match trait functions")
        #_mcmatch(abundances.matrix, spplist, abenv) ||
        #  error("Dimension mismatch")
        roclookup = Vector{AMDLookup{ROCVector{Int64}, ROCVector{Float64}}}(undef, length(lookup))
        # FIXME: write proper constructor
        for (i, lk) in enumerate(lookup)
            roclookup[i] = AMDLookup{ROCVector{Int64}, ROCVector{Float64}}(ROCVector(lk.x), ROCVector(lk.y), ROCVector(lk.p), ROCVector(lk.pnew), ROCVector(lk.moves))
        end
        return new{Part, SL, TR}(abundances, spplist, abenv, ordinariness,
                                 relationship, lookup, roclookup, cache)
    end
end

"""
    AMDEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
        rel::AbstractTraitRelationship)

Function to create an `AMDEcosystem` given a species list, an abiotic environment and trait relationship. An optional population function can be added, `popfun`, which defaults to generic random filling of the ecosystem.
"""
function AMDEcosystem(popfun::F, spplist::SpeciesList{T, Req},
                   abenv::GridAbioticEnv,
                   rel::AbstractTraitRelationship) where {F <: Function, T, Req}

    # Check there is enough energy to support number of individuals at set up
    #all(getenergyusage(spplist) .<= getavailableenergy(abenv)) ||
    #error("Environment does not have enough energy to support species")
    # Create matrix landscape of zero abundances
    ml = emptyAMDgridlandscape(abenv, spplist)
    # Populate this matrix with species abundances
    popfun(ml, spplist, abenv, rel)
    # Create lookup table of all moves and their probabilities
    lookup_tab = collect(map(k -> genlookups(abenv.habitat, k),
                             getkernels(spplist.movement)))
    nm = zeros(Int64, size(ml.matrix))
    totalE = ROCArray(zeros(Float64, (countsubcommunities(abenv), numrequirements(Req))))
    exchanged_energy = zeros(Float64, (counttypes(spplist), numrequirements(Req)))
    if numrequirements(Req) == 1
        exchanged_energy[:, 1] .= spplist.requirement.energy * spplist.requirement.exchange_rate .|> NoUnits
    else
        exchanged_energy[:, 1] .= spplist.requirement.r1.energy * spplist.requirement.r1.exchange_rate .|> NoUnits
        exchanged_energy[:, 2] .= spplist.requirement.r2.energy * spplist.requirement.r2.exchange_rate .|> NoUnits
    end
    sp_birth_rate = spplist.params.birth .* 1.0Unitful.d .|> NoUnits
    sp_death_rate = spplist.params.death .* 1.0Unitful.d .|> NoUnits
    width, height = _getdimension(abenv.habitat)
    cache = AMDCache(ROCArray(nm), totalE, ROCArray(exchanged_energy), ROCArray(abenv.active),
                     ROCArray(sp_birth_rate), ROCArray(sp_death_rate),
                     width, height, false)
    return AMDEcosystem{typeof(abenv), typeof(spplist), typeof(rel)}(ml, spplist,
                                                                  abenv,
                                                                  missing, rel,
                                                                  lookup_tab,
                                                                  cache)
end

function AMDEcosystem(spplist::SpeciesList, abenv::GridAbioticEnv,
                   rel::AbstractTraitRelationship)
    return AMDEcosystem(populate!, spplist, abenv, rel)
end

import Diversity.API: _getabundance
function _getabundance(eco::AMDEcosystem, raw::Bool)
    synchronise_from_gpu!(eco)
    if raw
        return eco.abundances.matrix
    else
        return _calcabundance(_gettypes(eco),
                              eco.abundances.matrix /
                              sum(eco.abundances.matrix))[1]
    end
end
