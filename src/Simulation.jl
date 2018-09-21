module Simulation

include("Dist.jl")
export jnorm, jexp, jpois, jbinom, junif, jdir, jmulti, Trapezoid

include("Phylo.jl")
export jtree, jcoal, assign_traits!, get_traits

include("TraitRelationship.jl")
export TraitRelationship,multiplicativeTR2, multiplicativeTR3, Gauss,
 Match, NoRelContinuous, NoRelDiscrete, Trapeze, Unif,
 additiveTR2, additiveTR3

include("Habitats.jl")
export ContinuousHab,ContinuousTimeHab, DiscreteHab, HabitatCollection2, HabitatCollection3,
 tempgrad

include("Energy.jl")
export SimpleRequirement, SolarRequirement, WaterRequirement, SimpleBudget, SolarBudget,
    WaterBudget, ReqCollection2, BudgetCollection2

if VERSION >= v"0.7-"
    include("DataTypes.jl")
end

include("AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, simplehabitatAE, degradedhabitatAE,
    eraAE, worldclimAE

include("Movement.jl")
export GaussianKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernel

include("Traits.jl")
export GaussTrait, DiscreteTrait, TempBin,RainBin, TraitCollection2, TraitCollection3,
DiscreteEvolve, ContinuousEvolve

include("Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("SpeciesList.jl")
export SpeciesList

include("Landscape.jl")
export GridLandscape

include("Ecosystem.jl")
export Ecosystem, getsize, gethabitat, gettraitrel, getgridsize,
 getdispersaldist, getdispersalvar, resetrate!,resettime!, getbudget

include("Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("HabitatUpdate.jl")
export getchangefun, TempChange, eraChange, worldclimChange

include("Scenarios.jl")
export SimpleScenario, RandHabitatLoss!, ClustHabitatLoss!, DisturbanceScenario,
 HabitatDisturbance!, UniformDecline, ProportionalDecline, LargeDecline, RareDecline,
 CommonDecline, HabitatReplacement, Invasive

include("Generate.jl")
export populate!, repopulate!,reenergise!, randomniches, update!, update_birth_move!,
 convert_coords, get_neighbours

include("Helper.jl")
export simulate!, simulate_record!,simulate_record_diversity!, expected_counts, generate_storage

if VERSION < v"0.7-"
    include("ReadTOML.jl")
    export readTOML, runTOML, readoutput
end

include("plotting.jl")
export plot_move, plot_abun,plot_mean,plot_diversity, plot_divergence, freq_hist, plotdiv

include("AdditionalDiversity.jl")
export meta_simpson, meta_shannon, meta_speciesrichness

end
