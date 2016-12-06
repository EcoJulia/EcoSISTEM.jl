module Simulation

include("Ecosystem.jl")
export Ecosystem, Habitats, StringTraits, MatrixLandscape

include("Diversity_funcs.jl")

include("Generate.jl")
export populate, SR

include("Phylo.jl")
export jtree, jcoal, assign_trait

include("Dist.jl")
export jexp, jpois, jbinom, junif, jdir, jnorm

end # module
