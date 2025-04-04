# SPDX-License-Identifier: LGPL-3.0-or-later

using StatsBase
using LinearAlgebra
using AMDGPU

function _inner_update!(i_sp, num_species, num_subcommunities, abundances, totalE, exchanged_energy, netmigration, active_sc, timestep, sp_birth_rate, sp_death_rate, width, height, bnd, lk_x, lk_y, lk_p, lk_pnew, lk_moves)
    # igpu is 0-indexed
    igpu = workitemIdx().x + (workgroupIdx().x - 1) * workgroupDim().x - 1
    if igpu >= num_subcommunities
        return
    end
    # i_sc, i_sp are 1-indexed
    i_sc = igpu + 1

    # Calculate how much birth and death should be adjusted
    #adjusted_birth, adjusted_death = energy_adjustment(width, totalE, energy, exchange_rate,
    #                                                   eco,
    #                                                   eco.abenv.budget,
    #                                                   i_sc, i_sp)

    adjusted_birth = 1.0
    adjusted_death = 1.0

    # Convert 1D dimension to 2D coordinates
    (x, y) = convert_coords(i_sc, width)
    # Check if grid cell currently active
    if active_sc[x, y] && (totalE[i_sc, 1] > 0)
        # Calculate effective rates
        birthrate = sp_birth_rate[i_sp] * timestep * adjusted_birth
        deathparam = sp_death_rate[i_sp] * timestep * adjusted_death

        # Turn deathparam into probability and cancel units of birthrate
        deathprob = 1.0 - exp(-deathparam)

        # if (birthrate < 0.0) || (deathprob < 0.0)
        #     @rocprintf "Birth: %f \n Death: %f \n \n i_sc: %d \n i_sp: %d" birthrate deathprob i_sc i_sp
        # end

        # Calculate how many births and deaths
        births = AMDGPU.rand_poisson(Cuint, 1, lambda=abundances[i_sp, i_sc] * birthrate)[1]
        # FIXME: this should follow a binomial distribution
        deaths = round(Int64, abundances[i_sp, i_sc] * deathprob)

        # Update population
        abundances[i_sp, i_sc] += (births - deaths)

        # Calculate moves and write to cache
        move!(width, height, lk_x, lk_y, lk_p, lk_pnew, lk_moves, bnd, i_sc, i_sp, netmigration,
              births, active_sc)
    end
    return
end

"""
    update!(eco::AMDEcosystem, time::Unitful.Time)

Function to update a ecosystem abundances and environment for one timestep.
"""
function update!(eco::AMDEcosystem, timestep::Unitful.Time)

    # Calculate dimenions of habitat and number of species
    n_sc = _countsubcommunities(eco.abenv.habitat)
    spp = size(eco.abundances.grid, 1)
    params = eco.spplist.params
    width = getdimension(eco)[1]

    timestep_days = uconvert(NoUnits, timestep / 1.0Unitful.d)

    # Set the overall energy budget of that square
    update_energy_usage!(eco)

    synchronise_to_gpu!(eco.abundances)

    group_size = 64
    grid_size = div(n_sc, group_size, RoundUp)

    AMDGPU.@sync begin
        for i_sp in 1:spp
            lk = eco.roclookup[i_sp]
            @roc groupsize=group_size gridsize=grid_size _inner_update!(
                i_sp,
                spp, n_sc, eco.abundances.rocmatrix, eco.cache.totalE,
                eco.cache.exchanged_energy,
                eco.cache.netmigration, eco.cache.active_sc,
                timestep_days, eco.cache.species_birth_rate,
                eco.cache.species_death_rate, eco.cache.width, eco.cache.height,
                getboundary(eco.spplist.movement),
                lk.x, lk.y, lk.p, lk.pnew, lk.moves)
        end
    end

    # Update abundances with all movements
    eco.abundances.rocmatrix .+= eco.cache.netmigration

    synchronise_from_gpu!(eco.abundances)

    # Invalidate all caches for next update
    invalidatecaches!(eco)

    # Update environment - habitat and energy budgets
    habitatupdate!(eco, timestep)
    return budgetupdate!(eco, timestep)
end

"""
    energy_adjustment(eco::Ecosystem, bud::AbstractBudget, i::Int64, sp::Int64)
Function to calculate how much birth and death rates should be adjusted by, according to how much energy is available, `bud`, in the grid square, `i`, and how much energy the species, `sp`, requires.
"""
function energy_adjustment(width, totalE, energy, exchange_rate, eco::AMDEcosystem, bud::AbstractBudget,
                           i::Int64, sp::Int64)
    if typeof(eco.spplist.params) <: NoGrowth
        return 0.0, 0.0
    else
        (x, y) = convert_coords(i, width)
        params = eco.spplist.params
        K = getbudget(eco)[x, y] * exchange_rate[1]
        # Get energy budgets of species in square
        ϵ̄ = energy[sp, 1] * exchange_rate[1]
        E = totalE[i, 1]
        # Traits
        ϵ̄real = 1 / traitfun(eco, i, sp)
        # Alter rates by energy available in current pop & own requirements
        birth_energy = ϵ̄^-params.longevity * ϵ̄real^-params.survival *
                       min(K / E, params.boost)
        death_energy = ϵ̄^-params.longevity * ϵ̄real^params.survival * (E / K)
    end
    return birth_energy, death_energy
end

# FIXME: adapt for AMD
function energy_adjustment(eco::AMDEcosystem, bud::BudgetCollection2,
                           i::Int64, sp::Int64)
    width = getdimension(eco)[1]
    (x, y) = convert_coords(eco, i, width)
    params = eco.spplist.params
    # FIXME: cache budget, check unit stuff on GPU
    # K1 = _getbudget(eco.abenv.budget.b1)[x, y] *
    #     eco.cache.exchange_rate[1]
    # K2 = _getbudget(eco.abenv.budget.b2)[x, y] *
    #     eco.cache.exchange_rate[2]
    K1 = 1.0
    K2 = 1.0
    # Get abundances of square we are interested in
    # Get energy budgets of species in square
    ϵ̄1 = eco.cache.exchanged_energy[sp, 1]
    ϵ̄2 = eco.cache.exchanged_energy[sp, 2]
    E1 = eco.cache.totalE[i, 1]
    E2 = eco.cache.totalE[i, 2]
    ϵ̄real1 = 1 / traitfun(eco, i, sp)
    ϵ̄real2 = 1 / traitfun(eco, i, sp)
    # Alter rates by energy available in current pop & own requirements
    birth_energy = (ϵ̄1 * ϵ̄2)^-params.longevity *
                   (ϵ̄real1 * ϵ̄real2)^-params.survival *
                   min(K1 / E1, K2 / E2, params.boost)
    death_energy = (ϵ̄1 * ϵ̄2)^-params.longevity *
                   (ϵ̄real1 * ϵ̄real2)^params.survival * max(E1 / K1, E2 / K2)
    return birth_energy, death_energy
end

"""
    calc_lookup_moves!(bound, x::Int64, y::Int64, sp::Int64, eco::Ecosystem, abun::Int64)

Function to calculate the number of moves taken by a species, `sp`, from a specific grid square location (`x`, `y`). There is a boundary condition, `bound`, which determines how the species can move across space (see AbstractBoundary). The total abundance of individuals is given in `abun`, which may be the number of births in the timestep, or total indiviuals.
"""
function calc_lookup_moves!(width, height, lk_x, lk_y, lk_p, lk_pnew, lk_moves, bound::NoBoundary, x::Int64, y::Int64, sp::Int64,
                            abun::Int64, active_sc)
    maxX = width - x
    maxY = height - y
    # Can't go over maximum dimension
    for i in eachindex(lk_x)
        valid = (-x < lk_x[i] <= maxX) && (-y < lk_y[i] <= maxY) &&
                (active_sc[lk_x[i] + x, lk_y[i] + y])

        lk_pnew[i] = valid ? lk_p[i] : 0.0
    end
    lk_pnew ./= sum(lk_pnew)
    lk_moves .= 1
    # dist = Multinomial(abun, lk_pnew)
    # return rand!(dist, lk_moves)
    return lk_moves
end

function calc_lookup_moves!(width, height, lk_x, lk_y, lk_p, lk_pnew, lk_moves, bound::Cylinder, x::Int64, y::Int64, sp::Int64,
                            abun::Int64, active_sc)
    maxX = width - x
    maxY = height - y
    # Can't go over maximum dimension
    for i in eachindex(lk_x)
        newx = -x < lk_x[i] <= maxX ? lk_x[i] + x :
               mod(lk_x[i] + x - 1, width) + 1

        valid = (-y < lk_y[i] <= maxY) &&
                (active_sc[newx, lk_y[i] + y])

        lk_pnew[i] = valid ? lk_p[i] : 0.0
    end
    lk_pnew ./= sum(lk_pnew)
    lk_moves .= 1
    # dist = Multinomial(abun, lk_pnew)
    # return rand!(dist, lk_moves)
    return lk_moves
end

function calc_lookup_moves!(width, height, lk_x, lk_y, lk_p, lk_pnew, lk_moves, bound::Torus, x::Int64, y::Int64, sp::Int64,
                            abun::Int64, active_sc)
    maxX = width - x
    maxY = height - y
    # Can't go over maximum dimension
    for i in eachindex(lk_x)
        newx = -x < lk_x[i] <= maxX ? lk_x[i] + x :
               mod(lk_x[i] + x - 1, width) + 1
        newy = -y < lk_y[i] <= maxY ? lk_y[i] + y :
               mod(lk_y[i] + y - 1, height) + 1
        valid = active_sc[newx, newy]

        lk_pnew[i] = valid ? lk_p[i] : 0.0
    end
    lk_pnew ./= sum(lk_pnew)
    lk_moves .= 1
    # dist = Multinomial(abun, lk_pnew)
    # return rand!(dist, lk_moves)
    return lk_moves
end

"""
    move!(width, height, lookup::Lookup, bnd::AbstractBoundary, i_sc::Int64, i_sp::Int64, grd::Matrix{Int64}, births::Int64, active_sc)

Function to calculate the movement of species `sp` from a given position in the
landscape `i`, using the lookup table found in the Ecosystem and updating the
movement patterns on a cached grid, `grd`. Optionally, a number of births can be
provided, so that movement only takes place as part of the birth process, instead
of the entire population
"""
function move!(width, height, lk_x, lk_y, lk_p, lk_pnew, lk_moves, bnd, i_sc::Int64, i_sp::Int64, grd::Matrix{Int64}, births::Int64, active_sc)
    (x, y) = convert_coords(i_sc, width)
    calc_lookup_moves!(width, height, lk_x, lk_y, lk_p, lk_pnew, lk_moves, bnd, x, y, i_sp, births, active_sc)
    # Lose moves from current grid square
    grd[i_sp, i_sc] -= births
    # Map moves to location in grid
    mov = lk_moves
    for i_sc in eachindex(lk_x)
        newx = mod(lk_x[i_sc] + x - 1, width) + 1
        newy = mod(lk_y[i_sc] + y - 1, height) + 1
        loc = convert_coords((newx, newy), width)
        grd[i_sp, loc] += mov[i_sc]
    end
end

"""
    populate!(ml::AMDGridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic, traits::Bool)

Function to populate a grid landscape given the abundances found in species list according to availability of resources.
"""

function populate!(ml::AMDGridLandscape, spplist::SpeciesList, abenv::AB,
                   rel::R) where {AB <: AbstractAbiotic,
                                  R <: AbstractTraitRelationship}
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b = reshape(ustrip.(_getbudget(abenv.budget)), size(grid))
    activity = reshape(copy(abenv.active), size(grid))
    units = unit(b[1])
    b[.!activity] .= 0.0 * units
    B = b ./ sum(b)
    # Loop through species
    for i in eachindex(spplist.abun)
        rand!(Multinomial(spplist.abun[i], B), (@view ml.matrix[i, :]))
    end
end

function populate!(ml::AMDGridLandscape, spplist::SpeciesList,
                   abenv::GridAbioticEnv{H, BudgetCollection2{B1, B2}},
                   rel::R) where {H <: AbstractHabitat, B1 <: AbstractBudget,
                                  B2 <: AbstractBudget,
                                  R <: AbstractTraitRelationship}
    # Calculate size of habitat
    dim = _getdimension(abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b1 = reshape(copy(_getbudget(abenv.budget, :b1)), size(grid))
    b2 = reshape(copy(_getbudget(abenv.budget, :b2)), size(grid))
    units1 = unit(b1[1])
    units2 = unit(b2[1])
    activity = reshape(copy(abenv.active), size(grid))
    b1[.!activity] .= 0.0 * units1
    b2[.!activity] .= 0.0 * units2
    B = (b1 ./ sum(b1)) .* (b2 ./ sum(b2))
    # Loop through species
    for i in eachindex(spplist.abun)
        rand!(Multinomial(spplist.abun[i], B ./ sum(B)),
              (@view ml.matrix[i, :]))
    end
end

"""
    repopulate!(eco::AMDEcosystem, abun::Int64)
Function to repopulate an ecosystem `eco`, with option for including trait preferences. An additional `abun` parameter can be included, in order to repopulate the ecosystem with a specified number of individuals.
"""
function repopulate!(eco::AMDEcosystem)
    eco.abundances = emptyAMDgridlandscape(eco.abenv, eco.spplist)
    eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun),
                                        length(eco.spplist.abun)))
    return populate!(eco.abundances, eco.spplist, eco.abenv, eco.relationship)
end
function repopulate!(eco::AMDEcosystem, abun::Int64)
    dim = _getdimension(eco.abenv.habitat)
    len = dim[1] * dim[2]
    grid = collect(1:len)
    # Set up copy of budget
    b = reshape(copy(_getbudget(eco.abenv.budget)), size(grid))
    units = unit(b[1])
    activity = reshape(copy(eco.abenv.active), size(grid))
    b[.!activity] .= 0.0 * units
    pos = sample(grid[b .> (0 * units)], abun)
    # Add individual to this location
    map(pos) do p
        return eco.abundances.matrix[end, p] += 1
    end
end

"""
    traitpopulate!(ml::AMDGridLandscape, spplist::SpeciesList,
                   abenv::AbstractAbiotic)

Function to populate a grid landscape given the abundances found in species list based upon how well the species traits match their environment.
"""
function traitpopulate!(ml::AMDGridLandscape, spplist::SpeciesList,
                        abenv::AB,
                        rel::R) where {AB <: AbstractAbiotic,
                                       R <: AbstractTraitRelationship}
    # Calculate size of habitat
    dim = _getdimension(abenv.habitat)
    numsquares = dim[1] * dim[2]
    numspp = length(spplist.names)
    maxrng = spplist.traits.mean .+ spplist.traits.var
    minrng = spplist.traits.mean .- spplist.traits.var
    hab = reshape(abenv.habitat.matrix, numsquares)
    probabilities = [_traitfun(abenv.habitat, spplist.traits, rel, i, sp)
                     for i in 1:numsquares, sp in 1:numspp]
    # Loop through species
    for i in eachindex(spplist.abun)
        if spplist.native[i]
            # Get abundance of species
            probs = probabilities[:, i] ./ sum(probabilities[:, i])
            probs[isnan.(probs)] .= 1 / numsquares
            abun = rand(Multinomial(spplist.abun[i], probs))
            # Add individual to this location
            ml.matrix[i, :] .+= abun
        end
    end
end

"""
    repopulate!(eco::Ecosystem, abun::Int64)
Function to repopulate an ecosystem `eco`, with option for including trait preferences. An additional `abun` parameter can be included, in order to repopulate the ecosystem with a specified number of individuals.
"""
function traitrepopulate!(eco::AMDEcosystem)
    eco.abundances = emptyAMDgridlandscape(eco.abenv, eco.spplist)
    eco.spplist.abun = rand(Multinomial(sum(eco.spplist.abun),
                                        length(eco.spplist.abun)))
    return traitpopulate!(eco.abundances, eco.spplist, eco.abenv,
                          eco.relationship)
end

"""
    emptypopulate!(ml::AMDGridLandscape, spplist::SpeciesList,
                   abenv::AB, rel::R) where {AB <: EcoSISTEM.AbstractAbiotic, R <: EcoSISTEM.AbstractTraitRelationship}
"""
function emptypopulate!(ml::AMDGridLandscape, spplist::SpeciesList,
                        abenv::AB,
                        rel::R) where {AB <: EcoSISTEM.AbstractAbiotic,
                                       R <: EcoSISTEM.AbstractTraitRelationship}
    @warn "Ecosystem not populated!"
end
"""
    reenergise!(eco::Ecosystem, budget::Union{Float64, Unitful.Quantity{Float64}}, grid::Tuple{Int64, Int64})
Function to refill an ecosystem `eco`, with energy from a budget value, `budget` and a grid size.
"""
function reenergise!(eco::AMDEcosystem,
                     budget::Union{Float64, Unitful.Quantity{Float64}},
                     grid::Tuple{Int64, Int64})
    return fill!(eco.abenv.budget.matrix, budget / (grid[1] * grid[2]))
end
function reenergise!(eco::AMDEcosystem,
                     budget::Tuple{Unitful.Quantity{Float64},
                                   Unitful.Quantity{Float64}},
                     grid::Tuple{Int64, Int64})
    fill!(eco.abenv.budget.b1.matrix, budget[1] / (grid[1] * grid[2]))
    return fill!(eco.abenv.budget.b2.matrix, budget[2] / (grid[1] * grid[2]))
end
