# Integration with Diversity.jl

Simulation is integrated with the [Diversity](https://github.com/richardreeve/Diversity.jl) package, so that diversity measures can be calculated directly on ecosystems.

See [Basics](https://boydorr.github.io/Simulation.jl/dev/basics/) for more information on setting up an Ecosystem.

```julia
using Diversity
# Subcommunity measures
norm_sub_alpha(eco, 1.0)
# Or metacommunity measures
norm_meta_alpha(eco, 1.0)
# Or multiple values of q
norm_sub_beta(eco, 0.0:3.0)
```
