# SPDX-License-Identifier: LGPL-3.0-or-later
#
# `GeneralistInvasive` and `SpecialistInvasive` from `examples/models/scenarios.jl`, recreated.
#
# **An invasive species *arrives*** - it is not present from the start with zero abundance, and it
# does not acquire its niche on arrival. `AddSpecies` says exactly that, carrying the invader's own
# tolerance, so the recreation is a translation rather than a workaround.
#
# **The originals could not do this, and the reason is worth knowing.** Each *mutated an existing
# species' traits* mid-run - `eco.spplist.tolerance.var[invasive] .= maximum(...)` - because the only
# way to add a species cloned the last one. Those lines have also been
# dead since `NicheTolerance` replaced `GaussTrait`: there is no `.var` or `.mean` field to write, so
# both scenarios throw `FieldError` on their first line and have done for some time.
#
# A run that gains a species needs room to record it - see `generate_storage`'s `maxspecies`.

using EcoSISTEM
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
# Only `Normal` is wanted, and it is named rather than swept in. A blanket `using` of a large package
# lands every one of its exports in whichever module the example is included into, where any of them
# may collide with a name already defined there -- and the example harness checks stderr is empty, so
# such a collision fails the suite. Naming the import keeps the next export from finding the next
# name.
#
# Worth knowing that a clash here is invisible on a current Julia: 1.11 warns about it and 1.12 does
# not warn at all, so a collision can pass locally and fail on an older runner.
using Distributions: Normal

# The native assemblage the invader arrives into - spread evenly along the temperature gradient.
# `numnative` is one short of the configuration's species count: the invader arriving mid-run
# makes it up, so the assemblage ends at `numspecies` however the scale is set.
function invasion_ecosystem(; numnative = configuration().numspecies - 1,
                            individuals = configuration().individuals, seed = 1)
    env = GridHabitat(regime = GradientSpec(283.0K, 293.0K,
                                            axis = Temperature),
                      supply = UniformSpec(supply_density(individuals),
                                           axis = SolarRadiation),
                      area = study_area())
    means = collect(range(284.0K, 292.0K, length = numnative))
    widths = fill(3.0K, numnative)
    species = build_species(numnative, tolerance = (means, widths),
                            toleranceaxis = Temperature,
                            demand = intervention_demand(),
                            demandaxis = SolarRadiation,
                            abundance = individuals,
                            seed = seed)
    return build_ecosystem(species, env, seed = seed)
end

# ---------------------------------------------------------------------------
# The arrival
# ---------------------------------------------------------------------------
# A **generalist** invader carries the widest tolerance in the assemblage; a **specialist** the
# narrowest, centred on the hot end of the gradient. Both are properties of the arriving species, so
# both are stated here rather than imposed on the ecosystem afterwards.
function invasive_species(kind::Symbol; abundance = 1000)
    niche = kind === :generalist ? Normal(288.0, 8.0) : Normal(292.5, 0.5)
    return AddSpecies(tolerance = niche, abundance = abundance,
                      name = string(kind))
end

# `AllCells()`: the arrival is scattered by resource availability, the same way every species is
# placed, so the region is not the intervention's business. Saying `AllCells()` says so.
function invasion(kind::Symbol; at = 2.0year, abundance = 1000)
    return Intervention(AtTime(at), AllCells(),
                        invasive_species(kind, abundance = abundance))
end

generalist_invasion(; kw...) = invasion(:generalist; kw...)
specialist_invasion(; kw...) = invasion(:specialist; kw...)
