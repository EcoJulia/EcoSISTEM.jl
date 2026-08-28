```@meta
CurrentModule = EcoSISTEM
```

# How the model works

[The basics](@ref "The basics of EcoSISTEM.jl") shows how to assemble and run a simulation. This
page says what that simulation *is* — the population model underneath, and what each parameter
does to it. It is worth reading before interpreting any output, because two of the parameters
change quite different things and their names do not say which.

## A spatial metacommunity

EcoSISTEM simulates **abundances**, not individuals: a count per species per grid cell, changing
each timestep as individuals are born, die and disperse to neighbouring cells. There is one
community per cell, drawn from a shared pool of species, which is what makes an
[`Ecosystem`](@ref) a metacommunity in the sense
[Diversity.jl](@ref "Integration with Diversity.jl") means.

Every demographic parameter is a **rate**, so the answer does not depend on how finely you step
the clock — see [Time in EcoSISTEM](@ref).

## Two ways a cell matters

A cell bears on the species living in it in two independent ways, and the whole model is built
along that distinction:

| | the environment holds | each species brings | what it decides |
| --- | --- | --- | --- |
| [`Condition`](@ref) | a **regime** — temperature, rainfall, land cover | a **tolerance** | **where** a species can persist |
| [`Resource`](@ref) | a **supply** — light, water, space | a **demand** | **how many** can persist there |

A condition is a state: every species in the cell experiences the same value, and nothing is
divided between them. What differs is how well each one copes, which is what its tolerance says.

A resource is a shared pool that species compete for. Each states how much it needs, the demands
of everything present are summed, and that total is set against what the cell supplies.

## Regulation, and where carrying capacity comes from

In each cell the model holds two numbers: `K`, what the cell supplies, and `E`, the **total**
demand summed over every species present. Births and deaths are then scaled by their ratio:

  - births rise when resource is plentiful, by `min(K / E, boost)` — the `boost` ceiling stops an
    empty cell producing an unbounded birth rate;
  - deaths rise as demand approaches supply, by `E / K`.

**There is no carrying-capacity parameter anywhere in the package.** A cell fills until births and
deaths balance, and where that lands is a consequence of supply, of who is present, and of how well
each of them is suited to the cell. Adding a species to a cell raises `E` for everyone in it, which
is the only route by which species affect one another: there are currently no pairwise terms, no
predation and no interference, so an assemblage is regulated entirely by what it collectively needs
against what its surroundings provide.

Where a species draws on several resources at once, the scarcest binds — births take the `min` of
the availability ratios and deaths the `max` of the demand ratios, so every demand must be met for
a population to grow. That is Liebig's law of the minimum, and it means an abundant resource
correctly makes no difference.

Births are then drawn from a Poisson distribution and deaths from a Binomial, per species per cell,
from that species' own random stream.

## The two per-species exponents, and why they are not interchangeable

Two parameters weight those rates, and this is the part most worth knowing:

```julia
birth ∝ demand^-longevity  *  suitability^-survival  *  min(K/E, boost)
death ∝ demand^-longevity  *  suitability^+survival  *  (E/K)
```

Look at the signs.

**`longevity` carries the same exponent on both rates**, so it cancels from the birth/death ratio
entirely. It makes a species slow-and-long-lived or fast-and-short-lived — a body-size proxy — and
it does **not** change where that species can persist or what abundance it reaches. It sets the
*tempo* of turnover.

**`survival` carries opposite exponents**, lowering deaths and raising births where a species is
well suited and doing the reverse where it is not. That is what makes a tolerance mean anything, so
`survival` sets the *niche*. At `survival = 0` a species' tolerances stop affecting its demography
at all.

If a change moves the wrong thing in your results, that table is where to look first.

## Suitability

A species' suitability in a cell comes from matching its tolerance against the regime on the same
axis. For a continuous tolerance it is the density of the species' response distribution evaluated
at the cell's value, so a **narrow** niche peaks higher at its optimum than a broad one and falls
away faster. That specialist advantage is deliberate: peak-normalising would leave a specialist
strictly worse off than a generalist everywhere, which is not what a niche means.

Suitability is a relative weight rather than a probability, so it may exceed 1, and nothing clamps
it. Where several regime layers are in play their suitabilities are combined before anything else
sees them.

## What the model does not do

Stated so you can tell whether it suits your question:

  - **individuals are not tracked** — there is no age, size or genotype, only counts;
  - **species do not interact in pairs** — no predation, competition coefficients or interference,
    only shared demand on a common pool;
  - **dispersal is by kernel** — offspring are placed by distance, with no directed movement or
    habitat selection beyond what suitability already implies.

## Reproducibility

Each species has its own deterministic random stream, seeded from the run's `seed`. A result
therefore does not depend on how the work was divided: the same seed gives the same answer on one
thread or many, and serially or across MPI ranks. See
[Running at scale](@ref "Running a simulation at scale") for what that buys.

## Where to go next

  - [Layers, conditions and resources](@ref) — how a layer carries its meaning, and where climate
    data comes from.
  - [Axes, units and roles: how a layer is classified](@ref) — why a layer's axis, not its unit,
    says what it is.
  - [Time in EcoSISTEM](@ref) — the clock, and environments that change as a run proceeds.
