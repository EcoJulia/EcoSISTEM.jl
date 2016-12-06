module Simulation

include("Ecosystem.jl")
export Ecosystem, Habitats, StringTraits, MatrixLandscape

include("Diversity_funcs.jl")

include("Generate.jl")
export populate, SR

include("Dist.jl")
export jnorm, jexp, jpois, jbinom, junif, jdir

include("Phylo.jl")
export jtree, jcoal, assign_trait, get_traits

end
