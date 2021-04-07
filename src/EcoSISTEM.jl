module EcoSISTEM

include("units.jl")

module ClimatePref

include("ClimatePref/ClimatePref.jl")

end

include("Dist.jl")
export Trapezoid

include("Phylo.jl")
export assign_traits!, get_traits, resettraits!, reroot!

include("TraitRelationship.jl")
export TraitRelationship,multiplicativeTR2, multiplicativeTR3, Gauss,
 Match, NoRelContinuous, NoRelDiscrete, Trapeze, Unif,
 additiveTR2, additiveTR3

include("Habitats.jl")
export ContinuousHab,ContinuousTimeHab, DiscreteHab, HabitatCollection2, HabitatCollection3, tempgrad, raingrad

include("Energy.jl")
export SimpleRequirement, SizeRequirement, SolarRequirement, WaterRequirement, VolWaterRequirement, SimpleBudget, SolarBudget, SolarTimeBudget, WaterBudget, VolWaterBudget, WaterTimeBudget, VolWaterTimeBudget, ReqCollection2, BudgetCollection2

include("AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, raingradAE, peakedgradAE, simplehabitatAE, degradedhabitatAE, eraAE, worldclimAE

include("Movement.jl")
export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernel, Torus, Cylinder, NoBoundary

include("Traits.jl")
export GaussTrait, DiscreteTrait, TempBin,RainBin,
TraitCollection2, TraitCollection3,DiscreteEvolve,
ContinuousEvolve

include("Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("SpeciesList.jl")
export SpeciesList

include("Landscape.jl")
export GridLandscape, CachedGridLandscape

include("Ecosystem.jl")
export Ecosystem, CachedEcosystem, getsize, gethabitat, gettraitrel, getgridsize,
 getdispersaldist, getdispersalvar, resetrate!,resettime!, getbudget, addspecies!

include("Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("HabitatUpdate.jl")
export getchangefun, TempChange, RainChange, TempFluct, eraChange, worldclimChange

include("Scenarios.jl")
export SimpleScenario, FluctScenario, MultiScenario

include("Generate.jl")
export populate!, repopulate!, traitpopulate!, traitrepopulate!, emptypopulate!, reenergise!, randomniches, update!, update_birth_move!, convert_coords, get_neighbours

using Requires
function __init__()
    @require MPI="da04e1cc-30fd-572f-bb4f-1f8673147195" @eval begin
        include("MPILandscape.jl")
        export MPIGridLandscape

        include("MPIEcosystem.jl")
        export MPIEcosystem, gather_abundance!, gather_diversity

        include("MPIGenerate.jl")
    end
end

include("Helper.jl")
export simulate!, simulate_record!,simulate_record_diversity!, expected_counts, generate_storage

include("Cache.jl")
export abundances, clearcache

include("DiversitySet.jl")
export DiversitySet, updatesimulation!, gettimes

include("AdditionalDiversity.jl")
export meta_simpson, meta_shannon, meta_speciesrichness, mean_abun, geom_mean_abun, sorenson, pd, makeunique

end
