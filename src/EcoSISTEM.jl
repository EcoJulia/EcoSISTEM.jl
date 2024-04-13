module EcoSISTEM

# EcoSISTEM.Units sub-module
include("Units/Units.jl")

# EcoSISTEM.ClimatePref sub-module
include("ClimatePref/ClimatePref.jl")

# DataPipeline extension
"""
    EcoSISTEM.unziptemp(path::String)

Helper function for the FAIR Data Pipeline to unzip files that are stored as zips to a temporary folder.
"""
function unziptemp end

include("Dist.jl")
export Trapezoid

include("Phylo.jl")
export assign_traits!, get_traits, resettraits!, reroot!

include("TraitRelationship.jl")
export TraitRelationship, multiplicativeTR2, multiplicativeTR3, Gauss,
       Match, NoRelContinuous, NoRelDiscrete, Trapeze, Unif,
       additiveTR2, additiveTR3, LCmatch

include("Habitats.jl")
export ContinuousHab, ContinuousTimeHab, DiscreteHab, HabitatCollection2,
       HabitatCollection3, tempgrad, raingrad

include("Energy.jl")
export SimpleRequirement, SizeRequirement, SolarRequirement, WaterRequirement,
       VolWaterRequirement, SimpleBudget, SolarBudget, SolarTimeBudget,
       WaterBudget, VolWaterBudget, WaterTimeBudget, VolWaterTimeBudget,
       ReqCollection2, BudgetCollection2

include("AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, raingradAE, peakedgradAE,
       simplehabitatAE, degradedhabitatAE, eraAE, worldclimAE, bioclimAE, lcAE

include("Movement.jl")
export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement,
       NoMovement, getkernel, Torus, Cylinder, NoBoundary

include("Traits.jl")
export GaussTrait, DiscreteTrait, TempBin, RainBin,
       TraitCollection2, TraitCollection3, DiscreteEvolve,
       ContinuousEvolve, LCtrait

include("Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("SpeciesList.jl")
export SpeciesList

include("Landscape.jl")
export GridLandscape, CachedGridLandscape

include("Ecosystem.jl")
export Ecosystem, CachedEcosystem, getsize, gethabitat, gettraitrel,
       getgridsize,
       getdispersaldist, getdispersalvar, resetrate!, resettime!, getbudget,
       addspecies!

include("Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("HabitatUpdate.jl")
export getchangefun, TempChange, RainChange, TempFluct, eraChange,
       worldclimChange

include("Scenarios.jl")
export SimpleScenario, FluctScenario, MultiScenario

include("Generate.jl")
export populate!, repopulate!, traitpopulate!, traitrepopulate!, emptypopulate!,
       reenergise!, randomniches, update!, update_birth_move!, convert_coords,
       get_neighbours

include("Helper.jl")
export simulate!, simulate_record!, simulate_record_diversity!, expected_counts,
       generate_storage

include("Cache.jl")
export abundances, clearcache

include("DiversitySet.jl")
export DiversitySet, updatesimulation!, gettimes

include("AdditionalDiversity.jl")
export meta_simpson, meta_shannon, meta_speciesrichness, mean_abun,
       geom_mean_abun, sorenson, pd, makeunique

using Random

"""
    MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple}

MPIEcosystem abundances housed in the landscape, shared across multiple nodes.
"""
mutable struct MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple}
    rows_matrix::Matrix{Int64}
    cols_vector::Vector{Int64}
    reshaped_cols::Vector{RA}
    rows_tuple::NT
    cols_tuple::NT
    rngs::Vector{MersenneTwister}
end

"""
    MPIEcosystem{MPIGL <: MPIGridLandscape, Part <: AbstractAbiotic,
                 SL <: SpeciesList, TR <: AbstractTraitRelationship} <: AbstractEcosystem{Part, SL, TR}

MPIEcosystem houses information on species and their interaction with their environment. It houses all information of a normal `Ecosystem` (see documentation for more details), with additional fields to describe which species are calculated on which machine. This includes: `sppcounts` - a vector of number of species per node, `firstsp` - the identity of the first species held by that particular node.
"""
mutable struct MPIEcosystem{MPIGL <: MPIGridLandscape, Part <: AbstractAbiotic,
                            SL <: SpeciesList,
                            TR <: AbstractTraitRelationship} <:
               AbstractEcosystem{Part, SL, TR}
    abundances::MPIGL
    spplist::SL
    abenv::Part
    ordinariness::Union{Matrix{Float64}, Missing}
    relationship::TR
    lookup::Vector{Lookup}
    sppcounts::Vector{Int32}
    firstsp::Int64
    sccounts::Vector{Int32}
    firstsc::Int64
    cache::Cache
end

export MPIGridLandscape
export MPIEcosystem
function gather_abundance end
function gather_diversity end
export gather_abundance, gather_diversity
function emptyMPIgridlandscape end
function synchronise_from_rows! end
function synchronise_from_cols! end
export emptyMPIgridlandscape, synchronise_from_rows!, synchronise_from_cols!

end
