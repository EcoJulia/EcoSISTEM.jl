# SPDX-License-Identifier: LGPL-3.0-or-later

# **Why these methods exist at all, and why they cannot simply be inherited.**
# Diversity's own `faith_pd`/`generalisedfaith_pd` are keyed on the metacommunity's **`Sim` type
# parameter** being a `PhyloBranches`. An `AbstractEcosystem` puts its **`SpeciesList`** in that
# slot - a `SpeciesList` *is* an `AbstractTypes` (so the declaration is legal) but is never a
# `PhyloBranches`; the phylogeny sits one level further in, as `SpeciesList`'s fourth type
# parameter. So the constraint is expressible, just not where Diversity looks for it.
# This is the same *holds-versus-is* distinction that produced seven silent Diversity defects
# (`[TF-FORWARD]`): a `SpeciesList` **holds** an `AbstractTypes` and must forward everything about
# it. Faith's PD is the one hook that cannot be forwarded, because the requirement is on the
# concrete type rather than on a value the wrapper could pass along.

# **From `Diversity.Ecology`, not `Diversity`.** `Diversity` re-exports both names but does not
# own the bindings, and importing from a re-exporter cannot extend them - Julia refuses with
# *"exported function Diversity.generalisedfaith_pd does not exist"* at precompile. Diversity's own
# `DiversityPhyloExt` imports from the owning submodule for the same reason.
import Diversity.Ecology: generalisedfaith_pd, faith_pd
using Diversity: DiversityLevel, subcommunityDiversity, metacommunityDiversity
using Diversity: norm_sub_alpha, meta_gamma
using Diversity.API: _getscale

# An ecosystem whose species list holds a `PhyloBranches` - the only kind for which Faith's PD is
# defined. `SpeciesList`'s `T` parameter makes this a **static** constraint, so this dispatches
# rather than testing a value at run time.
const _PhyloEcosystem = EcoSISTEM.AbstractEcosystem{<:Any,
                                                    <:EcoSISTEM.SpeciesList{<:Any,
                                                                            <:Any,
                                                                            <:Any,
                                                                            <:PhyloBranches,
                                                                            <:Any}}

# **Deliberately a transcription of `DiversityPhyloExt`'s body, not a new formula.** The measure
# is Diversity's and must not drift from it; the multiplication by `_getscale` is what makes this
# Faith's PD rather than Chao, Chiu and Jost's `^0D̄(T)` (PD *per unit* branch length), which is a
# different measure and is not abundance-independent as PD must be.
# **Why transcribe rather than wrap.** The obvious alternative is to build a
# `Metacommunity(getabundance(eco, true), sppl.types, eco.habitat)` and let Diversity's own method
# run, keeping the body in one place. Measured, and it is far too expensive: that constructor
# allocates **~16× the abundance matrix** (0.29 MB -> 13.85 MB, 18.75 MB -> 296.83 MB, the ratio
# converging as fixed overhead washes out), because it builds the derived phylogenetic state rather
# than merely copying. `getabundance(eco, true)` itself allocates **nothing** - the cost is all in
# the constructor. At this package's target scale - a UK 1 km grid with ~1400 species, a **2 GB+**
# abundance matrix - that is **~32 GB+** for one call.
# If Diversity ever gains a hook for *"my types are phylogenetic"* - the natural forwarding
# question, and the one this file works around - these two methods should be deleted in favour of
# it rather than kept in step by hand. That, not wrapping, is the way out of the duplication.
function generalisedfaith_pd(level::DiversityLevel, eco::_PhyloEcosystem)
    if level == subcommunityDiversity
        gs = norm_sub_alpha(eco, 0)
    elseif level == metacommunityDiversity
        gs = meta_gamma(eco, 0)
    else
        error("Can't calculate Faith's PD for $level")
    end
    gs[!, :diversity] .*= _getscale(eco)
    gs[!, :measure] .= "Faith's PD"
    select!(gs, Not(:q))
    return gs
end

function faith_pd(eco::_PhyloEcosystem)
    return generalisedfaith_pd(subcommunityDiversity, eco)
end
