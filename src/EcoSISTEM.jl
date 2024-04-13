module EcoSISTEM

using Unitful
using Unitful.DefaultSymbols

include("units.jl")

module ClimatePref

"""
    env_bool(key)
Checks for an enviroment variable and fuzzy converts it to a bool
"""
env_bool(key, default=false) = haskey(ENV, key) ? lowercase(ENV[key]) ∉ ["0","","false", "no"] : default

include("ClimatePref/ClimatePref.jl")

function unzip end
export unzip

function processMet end
function writeMet end
function MetOfficeDownload end
function getMetparams end
function getMetdata end
export processMet, writeMet, MetOfficeDownload, getMetparams, getMetdata

end

"""
    enum: DiseaseState

    Disease state of a group, from: Susceptible Infectious Removed OtherDiseaseState
"""
@enum DiseaseState Susceptible Infectious Removed OtherDiseaseState
export Susceptible, Infectious, Removed, OtherDiseaseState

include("Biodiversity/Dist.jl")
export Trapezoid

include("Biodiversity/Phylo.jl")
export assign_traits!, get_traits, resettraits!, reroot!

include("Biodiversity/TraitRelationship.jl")
export TraitRelationship,multiplicativeTR2, multiplicativeTR3, Gauss,
 Match, NoRelContinuous, NoRelDiscrete, Trapeze, Unif,
 additiveTR2, additiveTR3, LCmatch

include("Biodiversity/Habitats.jl")
export ContinuousHab,ContinuousTimeHab, DiscreteHab, HabitatCollection2, HabitatCollection3, tempgrad, raingrad

include("Biodiversity/Energy.jl")
export SimpleRequirement, SizeRequirement, SolarRequirement, WaterRequirement, VolWaterRequirement, SimpleBudget, SolarBudget, SolarTimeBudget, WaterBudget, VolWaterBudget, WaterTimeBudget, VolWaterTimeBudget, ReqCollection2, ReqCollection3, BudgetCollection2, BudgetCollection3

include("Biodiversity/AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, raingradAE, peakedgradAE, simplehabitatAE, degradedhabitatAE, eraAE, worldclimAE, bioclimAE, lcAE

include("Biodiversity/Movement.jl")
export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernels, Torus, Cylinder, NoBoundary

include("Biodiversity/Traits.jl")
export GaussTrait, DiscreteTrait, TempBin,RainBin,
TraitCollection2, TraitCollection3,DiscreteEvolve,
ContinuousEvolve, LCtrait

include("Biodiversity/Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("Biodiversity/SpeciesList.jl")
export SpeciesList, SpeciesTypes, NoPathogen

include("Biodiversity/Landscape.jl")
export GridLandscape, CachedGridLandscape

include("Transitions/Transitions.jl")
export TransitionList, create_transition_list, addtransition!, specialise_transition_list

include("Epidemiology/MedianGenerator.jl")
export MedianGenerator

include("Epidemiology/EpiControl.jl")
export NoControl, Lockdown

include("Epidemiology/shrink.jl")
export shrink_to_active, convert_population

include("Epidemiology/EpiEnv.jl")
export GridEpiEnv, simplehabitatAE, ukclimateAE

include("Epidemiology/EpiParams.jl")
export transition

include("Epidemiology/EpiMove.jl")
export EpiMovement, LongDistance

include("Epidemiology/EpiList.jl")
export SpeciesList, SIS, SIR, SEIR, SEIRS, SEI2HRD

include("Epidemiology/EpiLandscape.jl")
export EpiLandscape, host, virus

include("Transitions/TransitionLookup.jl")

include("Transitions/TransitionSystem.jl")
export Ecosystem, CachedEcosystem

include("Transitions/TransitionSystemHelperFuns.jl")
export gettransitions, getsize, gethabitat, gettraitrel, getgridsize,
getdispersaldist, getdispersalvar, resetrate!, resettime!, getbudget, getlookup,
addspecies!

include("Biodiversity/Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("Biodiversity/HabitatUpdate.jl")
export getchangefun, TempChange, RainChange, TempFluct, eraChange, worldclimChange

include("Biodiversity/Scenarios.jl")
export SimpleScenario, FluctScenario, MultiScenario

include("Biodiversity/Generate.jl")
export populate!, repopulate!, traitpopulate!, traitrepopulate!, emptypopulate!,
reenergise!, randomniches, update!, update_birth_move!,
convert_coords, get_neighbours, update_energy_usage!, seedinfected!

include("Biodiversity/Cache.jl")
export abundances, clearcache

include("Biodiversity/DiversitySet.jl")
export DiversitySet, updatesimulation!, gettimes

include("Biodiversity/AdditionalDiversity.jl")
export meta_simpson, meta_shannon, meta_speciesrichness, mean_abun, geom_mean_abun, sorenson, pd, makeunique

include("Epidemiology/EpiGenerate.jl")
export populate!

include("Epidemiology/EpiPlots.jl")

include("Epidemiology/Inference.jl")
export SIR_wrapper, SIR_wrapper!, SEI3HRD_wrapper, SEI3HRD_wrapper!

include("Transitions/BiodiversityTransitions.jl")
export BirthProcess, DeathProcess, GenerateSeed, AllDisperse, SeedDisperse,
UpdateEnergy, UpdateEnvironment, update_environment!

include("Transitions/EpiTransitions.jl")
export ForceProduce, ForceDisperse, ViralLoad, Exposure, EnvExposure, Infection,
DevelopSymptoms, Hospitalise, DeathFromInfection,
Recovery, SeedInfection, UpdateEpiEnvironment, update_epi_environment!,
get_env

include("Transitions/BiodiversityRun.jl")
include("Transitions/EpiRun.jl")

include("Transitions/Run.jl")
export run_rule!, update!, simulate!, simulate_record!, generate_storage

# Path into package
path(paths...) = joinpath(@__DIR__, "..", paths...)

using Random

"""
    MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple}

MPIEcosystem abundances housed in the landscape, shared across multiple nodes.
"""
mutable struct MPIGridLandscape{RA <: Base.ReshapedArray, NT <: NamedTuple} <: AbstractLandscape
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
mutable struct MPIEcosystem{MPIGL <: MPIGridLandscape, Part <: AbstractAbiotic, SL <: SpeciesList, TR <: AbstractTraitRelationship} <: AbstractEcosystem{MPIGL, Part, SL, TR, SpeciesLookup, Cache}
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
