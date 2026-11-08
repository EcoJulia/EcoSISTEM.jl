# SPDX-License-Identifier: LGPL-3.0-or-later
#
# `GeneralistInvasive` and `SpecialistInvasive` from `examples/models/scenarios.jl`, recreated.
#
# **An invasive species *arrives*** - it is not present from the start with zero abundance, and it
# does not acquire its niche on arrival. So it is declared in the pool with its own niche, left out
# of the run, and brought in by `AddSpecies`.
#
# **The originals mutated an existing species' traits mid-run**
# (`eco.spplist.tolerance.var[invasive] .= maximum(...)`), the only way they had to give an arrival
# a niche of its own; and they wrote fields `NicheTolerance` does not have, so neither runs.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
# Only `DataFrame` is wanted, and it is named rather than swept in. A blanket `using` of a large
# package lands every one of its exports in whichever module the example is included into, where
# any of them may collide with a name already defined there -- and the example harness checks stderr
# is empty, so such a collision fails the suite. 1.11 warns about such a clash and 1.12 does not, so
# one can pass locally and fail on an older runner.
using DataFrames: DataFrame

# The native assemblage, spread evenly along the temperature gradient, and the invader `kind`
# declared beside it but left out of the run: a **generalist** with the widest tolerance in the
# assemblage, or a **specialist** with the narrowest, centred on the hot end of the gradient.
# `numnative` is one short of the configuration's species count, so the pool holds `numspecies`
# however the scale is set.
function invasion_ecosystem(kind::Symbol = :generalist;
                            numnative = configuration().numspecies - 1,
                            individuals = configuration().individuals, seed = 1)
    env = GridHabitat(regime = GradientSpec(10.0°C, 20.0°C,
                                            axis = Temperature),
                      supply = UniformSpec(supply_density(individuals),
                                           axis = SolarRadiation),
                      area = study_area())
    optimum, width = kind === :generalist ? (15.0°C, 8.0K) : (19.5°C, 0.5K)
    pool = SpeciesPool(DataFrame(species = vcat(string.(1:numnative),
                                                string(kind)),
                                 optimum = vcat([t * °C
                                                 for t in range(11.0, 19.0,
                                                                length = numnative)],
                                                optimum),
                                 width = vcat(fill(3.0K, numnative), width)))
    deactivate_species!(pool, string(kind), verbosity = :silent)
    species = build_species(pool,
                            tolerance = ToleranceSpec(Temperature,
                                                      μ = TraitColumn(:optimum),
                                                      σ = TraitColumn(:width)),
                            demand = DemandSpec(SolarRadiation,
                                                intervention_demand()),
                            seed = seed, verbosity = :silent)
    return build_ecosystem(BySupply(total = individuals), species, env,
                           seed = seed)
end

# ---------------------------------------------------------------------------
# The arrival
# ---------------------------------------------------------------------------
# `abundance` individuals of the invader, scattered by resource availability as every species is
# placed, over the whole grid: `AllCells()`, so the region does not narrow where it lands.
function invasion(kind::Symbol; at = 2.0year, abundance = 1000)
    return Intervention(AtTime(at), AllCells(),
                        AddSpecies(string(kind),
                                   BySupply(total = abundance)))
end

generalist_invasion(; kw...) = invasion(:generalist; kw...)
specialist_invasion(; kw...) = invasion(:specialist; kw...)
