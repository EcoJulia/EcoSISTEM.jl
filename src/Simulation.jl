module Simulation

struct BiDict{A, B}
    forward::Dict{A, B}
    backward::Dict{B, A}

    BiDict{A, B}() where {A, B} = new{A, B}(Dict{A, B}(), Dict{B, A}())
end
import Base: getindex, setindex!
getindex(dict::BiDict{A,B}, elt::A) where {A, B} = dict.forward[elt]
getindex(dict::BiDict{A,B}, elt::B) where {A, B} = dict.backward[elt]
function setindex!(dict::BiDict{A,B}, val::B, elt::A)  where {A, B}
    dict.forward[elt] = val
    dict.backward[val] = elt
end

GLOBAL_typedict = BiDict{String, Type}()
GLOBAL_funcdict = BiDict{String, Function}()

include("Dist.jl")
export jnorm, jexp, jpois, jbinom, junif, jdir, jmulti, Trapezoid

include("Phylo.jl")
export jtree, jcoal, assign_traits!, get_traits, resettraits!, reroot!

include("TraitRelationship.jl")
export TraitRelationship,multiplicativeTR2, multiplicativeTR3, Gauss,
 Match, NoRelContinuous, NoRelDiscrete, Trapeze, Unif,
 additiveTR2, additiveTR3

include("Habitats.jl")
export ContinuousHab,ContinuousTimeHab, DiscreteHab, HabitatCollection2, HabitatCollection3, tempgrad, raingrad

include("Energy.jl")
export SimpleRequirement, SizeRequirement, SolarRequirement, WaterRequirement, VolWaterRequirement, SimpleBudget, SolarBudget, SolarTimeBudget, WaterTimeBudget, VolWaterTimeBudget, ReqCollection2, BudgetCollection2

include("AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, raingradAE, peakedgradAE, simplehabitatAE, degradedhabitatAE, eraAE, worldclimAE

include("Movement.jl")
export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernel, Torus, Cylinder, NoBoundary

include("Traits.jl")
export GaussTrait, DiscreteTrait, TempBin,RainBin, TraitCollection2, TraitCollection3,
DiscreteEvolve, ContinuousEvolve

include("Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("SpeciesList.jl")
export SpeciesList

include("Landscape.jl")
export GridLandscape, CachedGridLandscape

include("MPILandscape.jl")
export MPIGridLandscape

include("Ecosystem.jl")
export Ecosystem, CachedEcosystem, getsize, gethabitat, gettraitrel, getgridsize,
 getdispersaldist, getdispersalvar, resetrate!,resettime!, getbudget, addspecies!

include("MPIEcosystem.jl")
export MPIEcosystem

include("Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("HabitatUpdate.jl")
export getchangefun, TempChange, RainChange, TempFluct, eraChange, worldclimChange

include("Scenarios.jl")
export SimpleScenario, TempIncrease, RandHabitatLoss!, ClustHabitatLoss!, DisturbanceScenario,
 HabitatDisturbance!, UniformDecline, ProportionalDecline, LargeDecline, RareDecline,
 CommonDecline, HabitatReplacement, Invasive, SusceptibleDecline

include("Generate.jl")
export populate!, repopulate!, traitpopulate!, traitrepopulate!, emptypopulate!, reenergise!, randomniches, update!, update_birth_move!,
 convert_coords, get_neighbours

include("MPIGenerate.jl")

include("SantiniScenarios.jl")
export trait_populate!, trait_repopulate!

include("Helper.jl")
export simulate!, simulate_record!,simulate_record_diversity!, expected_counts, generate_storage

if (VERSION <= v"0.6") || (VERSION >= v"1.0-")
    include("ReadTOML.jl")
    export readTOML, runTOML, readoutput
end

include("Cache.jl")
export abundances, clearcache

include("DiversitySet.jl")
export DiversitySet, updatesimulation!, gettimes

#include("plotting.jl")
#export plot_move, plot_abun,plot_mean,plot_diversity, plot_divergence, freq_hist, plotdiv

include("AdditionalDiversity.jl")
export meta_simpson, meta_shannon, meta_speciesrichness, mean_abun, geom_mean_abun,
    sorenson, pd, makeunique

end
