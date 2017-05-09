module Simulation

include("Dist.jl")
export jnorm, jexp, jpois, jbinom, junif, jdir, jmulti

include("Phylo.jl")
export jtree, jcoal, assign_traits!, get_traits

include("TraitRelationship.jl")

include("Habitat.jl")
export Habitats, Niches

include("Energy.jl")
export SimpleRequirement

include("AbioticEnv.jl")
export GridAbioticEnv

include("Movement.jl")
export GaussianMovement

include("Landscape.jl")
export GridLandscape

include("Ecosystem.jl")
export Ecosystem, Habitats, Niches, StringTraits, RealTraits,
RealEnergy, GaussianMovement, SpeciesList, MatrixAbioticEnv,  MatrixLandscape,
Ecosystem, copy_eco, GridLandscape

include("Generate.jl")
export populate!, repopulate!, randomniches, update!, update_birth_move!

include("Helper.jl")
export run_sim, run_sim_spatial, expected_counts

include("plotting.jl")
export plot_move, plot_abun, plot_divergence, freq_hist

end
