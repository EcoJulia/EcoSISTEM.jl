# SPDX-License-Identifier: LGPL-3.0-or-later

using EcoSISTEM
using Unitful
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Phylo
using DataFrames
using Diversity

function create_eco(paramDict::Dict, habitat::A; bound::B = Torus(),
                    size = (mean = 1.0m^2, std = 0.0001m^2)) where
    {A <: EcoSISTEM.AbstractHabitat,
     B <: EcoSISTEM.BoundaryCondition}
    # Set up initial parameters for ecosystem
    birth = haskey(paramDict, "birth") ? paramDict["birth"] : 0.6 / year
    death = haskey(paramDict, "death") ? paramDict["death"] : 0.6 / year
    l = haskey(paramDict, "l") ? paramDict["l"] : 1.0
    s = haskey(paramDict, "s") ? paramDict["s"] : 0.2
    boost = haskey(paramDict, "boost") ? paramDict["boost"] : 1.0

    numSpecies = paramDict["numSpecies"]
    numInvasive = paramDict["numInvasive"]
    individuals = paramDict["numIndiv"]
    demand = paramDict["demand"]
    opts = paramDict["opts"]
    vars = paramDict["vars"]
    kernel = paramDict["kernel"]

    # Set up how much resource each species consumes
    names = map(x -> "$x", 1:(numSpecies + numInvasive))
    tree = rand(Ultrametric{BinaryTree{DataFrame, DataFrame}}(names))
    units = unit(size.mean)
    tolerance = ContinuousEvolve(uconvert(NoUnits, size.mean / units),
                                 uconvert(NoUnits, size.std / units),
                                 tree)

    resource_vec1 = SolarDemand(abs.(tolerance.mean) .* (demand[1] * units))
    resource_vec2 = WaterDemand(abs.(tolerance.mean) .* (demand[2] * units))

    resource_vec = DemandCollection2(resource_vec1, resource_vec2)
    # Collect model parameters together (in this order!!)
    if length(birth) > 1
        param = PopGrowth{typeof(unit(birth[1]))}(birth, death, l, s, boost)
    else
        param = EqualPop(birth, death, l, s, boost)
    end

    # Create ecosystem

    movement = BirthOnlyMovement(kernel, bound)

    tolerance = GaussTrait(opts, vars)
    native = fill(true, numSpecies + numInvasive)
    native[(numSpecies + 1):end] .= false
    if length(individuals) > 1
        abun = [individuals; fill(0, numInvasive)]
    else
        abun = [rand(Multinomial(individuals, numSpecies))
                fill(0, numInvasive)]
    end
    sppl = SpeciesList(numSpecies + numInvasive,
                       tolerance,
                       abun,
                       resource_vec,
                       movement,
                       param,
                       native)
    nichefit = Gauss{typeof(first(paramDict["opts"]))}()
    return eco = Ecosystem(tolerancepopulate!, sppl, habitat, nichefit)
end

function recreate_eco!(eco::Ecosystem,
                       totalK::Tuple{Unitful.Quantity{Float64},
                                     Unitful.Quantity{Float64}},
                       individuals::Int64)
    numSpecies = sum(eco.spplist.native)
    numInvasive = sum(.!eco.spplist.native)
    eco.spplist.abun = [rand(Multinomial(individuals, numSpecies))
                        fill(0, numInvasive)]
    resupply!(eco, totalK, size(eco.habitat.regime.matrix))
    eco.habitat.active .= true
    return tolerancerepopulate!(eco)
end
