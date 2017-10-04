module Simulation

include("Units.jl")
export day, week, month, year, Rates

include("Dist.jl")
export jnorm, jexp, jpois, jbinom, junif, jdir, jmulti

include("Phylo.jl")
export jtree, jcoal, assign_traits!, get_traits

include("Envtypes.jl")
export Temp

include("TraitRelationship.jl")
export TraitRelationship,multiplicativeTR2, multiplicativeTR3, GaussTemp, SimpleNiche, NoRel

include("Habitats.jl")
export ContinuousHab, DiscreteHab, HabitatCollection2, HabitatCollection3,
 tempgrad

include("Energy.jl")
export SimpleRequirement

include("AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, simplehabitatAE, degradedhabitatAE

include("Movement.jl")
export GaussianKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernel

include("Traits.jl")
export ContinuousTrait, DiscreteTrait, TraitCollection2, TraitCollection3,
DiscreteEvolve, ContinuousEvolve

include("Demographics.jl")
export PopGrowth, EqualPop

include("SpeciesList.jl")
export SpeciesList

include("Landscape.jl")
export GridLandscape

include("Ecosystem.jl")
export Ecosystem, getsize, getenvtype, gethabitat, gettraitrel, getgridsize,
 getdispersaldist, getdispersalvar, resetrate!

include("Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("HabitatUpdate.jl")
export getchangefun, TempChange

include("Scenarios.jl")
export SimpleScenario, RandHabitatLoss!, ClustHabitatLoss!

include("Generate.jl")
export populate!, repopulate!, randomniches, update!, update_birth_move!,
 convert_coords, get_neighbours

include("Helper.jl")
export simulate!, simulate_record!, expected_counts, generate_storage

include("plotting.jl")
export plot_move, plot_abun,plot_mean, plot_divergence, freq_hist, plotdiv


end
