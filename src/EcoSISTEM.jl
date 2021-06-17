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
 additiveTR2, additiveTR3

include("Biodiversity/Habitats.jl")
export ContinuousHab,ContinuousTimeHab, DiscreteHab, HabitatCollection2, HabitatCollection3, tempgrad, raingrad

include("Biodiversity/Energy.jl")
export SimpleRequirement, SizeRequirement, SolarRequirement, WaterRequirement, VolWaterRequirement, SimpleBudget, SolarBudget, SolarTimeBudget, WaterBudget, VolWaterBudget, WaterTimeBudget, VolWaterTimeBudget, ReqCollection2, BudgetCollection2

include("Biodiversity/AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, raingradAE, peakedgradAE, simplehabitatAE, degradedhabitatAE, eraAE, worldclimAE

include("Biodiversity/Movement.jl")
export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernel, Torus, Cylinder, NoBoundary

include("Biodiversity/Traits.jl")
export GaussTrait, DiscreteTrait, TempBin, RainBin,
TraitCollection2, TraitCollection3, DiscreteEvolve,
ContinuousEvolve

include("Biodiversity/Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("Biodiversity/SpeciesList.jl")
export SpeciesList, SpeciesTypes, NoPathogen

include("Biodiversity/Landscape.jl")
export GridLandscape, CachedGridLandscape

include("Transitions/Transitions.jl")
export TransitionList, create_transition_list, addtransition!, getprob

include("Epidemiology/MedianGenerator.jl")
export MedianGenerator

include("Epidemiology/data_utils.jl")
export parse_hdf5, get_3d_km_grid_axis_array, get_bng, get_en

include("Epidemiology/EpiControl.jl")
export NoControl, Lockdown

include("Epidemiology/shrink.jl")
export shrink_to_active, convert_population

include("Epidemiology/EpiEnv.jl")
export GridEpiEnv, simplehabitatAE, ukclimateAE

include("Epidemiology/EpiParams.jl")
export transition

include("Epidemiology/EpiMove.jl")
export EpiMovement, Commuting

include("Epidemiology/EpiList.jl")
export SpeciesList, SIS, SIR, SEIR, SEIRS, SEI2HRD

include("Epidemiology/EpiLandscape.jl")
export EpiLandscape, human, virus

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

using Requires
function __init__()
    @require MPI="da04e1cc-30fd-572f-bb4f-1f8673147195" @eval begin
        include("Biodiversity/MPILandscape.jl")
        export MPIGridLandscape

        include("Biodiversity/MPIEcosystem.jl")
        export MPIEcosystem, gather_abundance!, gather_diversity

        include("Biodiversity/MPIGenerate.jl")
    end
end

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
export ForceProduce, ForceDisperse, ViralLoad, EnvViralLoad, Exposure, Infection,
DevelopSymptoms, Hospitalise, DeathFromInfection,
Recovery, SeedInfection, UpdateEpiEnvironment, update_epi_environment!,
get_env, deterministic_seed!

include("Transitions/BiodiversityRun.jl")
include("Transitions/EpiRun.jl")

include("Transitions/Run.jl")
export run_rule!, update!, simulate!, simulate_record!, generate_storage

# Path into package
path(paths...) = joinpath(@__DIR__, "..", paths...)

end
