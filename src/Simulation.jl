module Simulation

include("Dist.jl")
export jnorm, jexp, jpois, jbinom, junif, jdir, jmulti

include("Phylo.jl")
export jtree, jcoal, assign_traits!, get_traits

include("TraitRelationship.jl")

include("Habitats.jl")
export Habitats, Niches

include("Energy.jl")
export SimpleRequirement

include("AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE

include("Movement.jl")
export GaussianKernel, BirthOnlyMovement, AlwaysMovement, getkernel

include("Traits.jl")
export BasicTrait

include("Demographics.jl")
export PopGrowth, EqualPop

include("SpeciesList.jl")
export SpeciesList

include("Landscape.jl")
export GridLandscape

include("Ecosystem.jl")
export Ecosystem

include("Generate.jl")
export populate!, repopulate!, randomniches, update!, update_birth_move!

include("Helper.jl")
export simulate!, simulate_record!, expected_counts, generate_storage

include("plotting.jl")
export plot_move, plot_abun, plot_divergence, freq_hist, plotdiv


end
