# SPDX-License-Identifier: LGPL-3.0-or-later

module TestSpecies

using EcoSISTEM
using Test
using Unitful
using Unitful.DefaultSymbols
using Distributions
using EcoSISTEM.Units
using Diversity
using EcoBase
# **Load-bearing, and it took a parallel runner to notice.** The `SpeciesList` forms this file
# builds are completed by `EcoSISTEMPhyloExt` — `src/Species.jl`'s 7-argument constructor forwards
# to an 8-argument one (`switch::Vector{Float64}`) that only exists once `Phylo` is loaded. This file
# never said so, and passed anyway: run serially, some *earlier* test file had already loaded `Phylo`
# into the shared process. Each parallel worker is a fresh process, so the order-dependency became
# a `MethodError`. The rule it broke is `CLAUDE.md`'s: a file needing a weak dependency must name it.
using Phylo

@testset "SpeciesList" begin
    ## Run simulation over a grid and plot
    numSpecies = 4
    numTraits = 2

    # Set up how much resource each species consumes
    resource_vec = Demand{SolarRadiation}(fill(2.0kJ / day, numSpecies))

    # Set probabilities
    birth = 6.0 / year
    death = 6.0 / year
    long = 1.0
    surv = 0.0
    boost = 1000.0
    timestep = 1.0month_mean_duration

    # Collect model parameters together (in this order!!)
    param = EqualPop(birth, death, long, surv, boost)

    individuals = 100

    # Create ecosystem
    kernel = GaussianKernel.(fill(2.0km, numSpecies), 10e-4)
    movement = AlwaysMovement(kernel)

    opts = fill(5.0°C, numSpecies)
    vars = rand(Uniform(0, 25 / 9), numSpecies) * °C
    tolerance = NicheTolerance(Temperature, Normal, opts, vars)
    abun = rand(Multinomial(individuals, numSpecies))
    native = fill(true, numSpecies)
    @test_nowarn sppl = SpeciesList(numSpecies, tolerance, abun, resource_vec,
                                    movement, param, native)
    @test_nowarn sppl = SpeciesList(numSpecies, numTraits, abun, resource_vec,
                                    movement, param, native)

    sppl = SpeciesList(numSpecies, tolerance, abun, resource_vec, movement,
                       param,
                       native)
    newsppl = SpeciesList{typeof(sppl.tolerance),
                          typeof(sppl.demand),
                          typeof(sppl.movement),
                          typeof(sppl.types),
                          typeof(sppl.params)}(sppl.tolerance, sppl.abun,
                                               sppl.demand, sppl.types,
                                               sppl.movement, sppl.params,
                                               sppl.native)
    @test newsppl.names == sppl.names

    # The mass-based `SpeciesList` constructor is **commented out** in `src/Species.jl`,
    # together with the `SizeDemand` it was the only caller of: with the free/dimensionless supply
    # gone there was nothing for that demand to pair with, so a species list built this way could
    # be constructed but never simulated. Both are parked against a metabolic reading of size that
    # nobody has implemented — see the note above `SizeDemand` in `src/Layer.jl`.
    @test_throws MethodError SpeciesList(numSpecies, numTraits, 10.0, 10.0, 0.5,
                                         100.0km^2, movement, param, native,
                                         [0.5, 0.5])

    # Test
    @test_nowarn sppl = SpeciesList(numSpecies, numTraits, abun, resource_vec,
                                    movement, UniqueTypes(numSpecies), param,
                                    native)
end

# **Diversity's API hooks drifted, and nothing noticed for a release.** `_calcsimilarity`'s second
# argument became a **scale** (`Real`) where ours still took an `AbstractArray`, so our method matched
# nothing and calls fell through to Diversity's *"not implemented"* fallback. Invisible in the
# suite, because the only thing that needs it is a species-restricted `view` — which nothing here
# exercised. These assert the delegation by the route that actually broke.
@testset "Diversity API hooks delegate to the wrapped types" begin
    species = build_species(EcoSISTEM.DefaultEcosystem(), verbosity = :silent)
    habitat = build_habitat(EcoSISTEM.DefaultEcosystem(), verbosity = :silent)
    eco = build_ecosystem(species, habitat, seed = 1)

    # A whole-assemblage view keeps the `SpeciesList` — nothing is restricted, so nothing is subset.
    whole = EcoBase.view(eco)
    @test whole isa EcoBase.AbstractAssemblage
    @test gettypes(whole) === eco.spplist

    # **Restricting SITES must not touch the types** — Diversity is explicit that subsetting them
    # is what would cost a phylogeny its tree, so a site-only view leaves them alone.
    bysite = EcoBase.view(eco, sites = 1:10)
    @test gettypes(bysite) === eco.spplist
    @test countsubcommunities(bysite) == 10
    @test counttypes(bysite) == counttypes(eco)

    # Restricting SPECIES is the route that failed: it needs `_subsettypes`, which needs
    # `_calcsimilarity` with a **scale**.
    bysp = EcoBase.view(eco, species = 1:3, sites = 1:10)
    @test counttypes(bysp) == 3
    @test countsubcommunities(bysp) == 10
    # Delegation preserves the *kind* Diversity can preserve: a `UniqueTypes` subsets to a
    # `UniqueTypes` rather than collapsing to a `GeneralTypes`.
    @test gettypes(bysp) isa UniqueTypes

    # And the hook itself, directly — a scale, not an array.
    @test EcoSISTEM._calcsimilarity(eco.spplist, 1.0) ==
          EcoSISTEM._calcsimilarity(eco.spplist.types, 1.0)

    # **The ordinariness path must not hardcode a scale of 1** — it is the phylogenetic case that
    # breaks, and the scale must come from the same `_calcabundance` call as the abundances.
    @test getordinariness!(eco) isa AbstractMatrix
    @test Diversity.API._getscale(eco) == 1          # `UniqueTypes`: no similarity, so unit scale

    # **Never cached — it must track the ecosystem.** A cache on `eco.ordinariness` is one nothing
    # invalidates when the ecosystem changes: after `_addspecies!` it returns a `(10, 100)` matrix
    # for an 11-species ecosystem. Such a cache is safe only by *position*, every mutator happening
    # to run after `update!`'s single `invalidatecaches!`.
    before = size(getordinariness!(eco), 1)
    EcoSISTEM._addspecies!(eco, 10)
    @test counttypes(eco) == before + 1
    @test size(getordinariness!(eco), 1) == before + 1
end

# **The forwarding rule, asserted by ENUMERATION rather than by a list.** A `SpeciesList` holds
# an `AbstractTypes`; every Diversity hook about the types must therefore give the same answer for
# the list as for what it wraps. Every defect found on 2026-08-18 was a breach of exactly this —
# hardcoded (`floattypes`, `_hassimilarity`), drifted (`_calcsimilarity` gained a **scale**), or
# absent (`_subsettypes`, `_calcordinariness`, `_addedoutputcols`, `_getaddedoutput`).
# Written as a sweep over `Diversity.API` so that a hook added upstream is caught here rather
# than by whatever breaks first — the failure mode was always silence, never an error.
@testset "SpeciesList forwards every AbstractTypes hook" begin
    species = build_species(EcoSISTEM.DefaultEcosystem(), verbosity = :silent)
    habitat = build_habitat(EcoSISTEM.DefaultEcosystem(), verbosity = :silent)
    sl = build_ecosystem(species, habitat, seed = 1).spplist
    wrapped = sl.types

    # The one-argument hooks: the list must answer exactly as the types it wraps.
    for f in (Diversity.API._hassimilarity, Diversity.API.floattypes,
        Diversity.API._getdiversityname, Diversity.API._getaddedoutput,
        Diversity.API._addedoutputcols)
        @test f(sl) == f(wrapped)
    end
    # …and those taking the raw/processed flag, or a scale.
    for raw in (true, false)
        @test Diversity.API._counttypes(sl, raw) ==
              Diversity.API._counttypes(wrapped, raw)
        @test Diversity.API._gettypenames(sl, raw) ==
              Diversity.API._gettypenames(wrapped, raw)
    end
    @test Diversity.API._calcsimilarity(sl, 1.0) ==
          Diversity.API._calcsimilarity(wrapped, 1.0)

    # **Not a list of names — a sweep.** Any `AbstractTypes` hook we do not forward is a silent
    # claim about a similarity structure the list does not own.
    AT = Diversity.API.AbstractTypes
    missing_hooks = Symbol[]
    for n in names(Diversity.API; all = true)
        isdefined(Diversity.API, n) || continue
        (n === :eval || n === :include) && continue
        f = getfield(Diversity.API, n)
        f isa Function || continue
        takes_types = any(methods(f)) do m
            ps = Base.unwrap_unionall(m.sig).parameters
            length(ps) >= 2 || return false
            t = Base.unwrap_unionall(ps[2])
            t isa Type && (t === AT || AT <: t)
        end
        takes_types || continue
        ours = any(m -> String(nameof(m.module)) == "EcoSISTEM" &&
                        occursin("SpeciesList", string(m.sig)), methods(f))
        ours || push!(missing_hooks, n)
    end
    @test missing_hooks == Symbol[]
end

end
