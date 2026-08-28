# Canonical tests

A canonical test records what the model **actually produced** and checks that it keeps producing it.
It answers a different question from the rest of the suite: a unit test asks *is this right?*, a
canonical test asks *has this changed?* Both matter, and neither substitutes for the other.

Run them:

```julia
julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
```

They are deliberately **not** part of the unit suite — they run whole simulations and read real
rasters — so they are a separate gate, alongside the examples, notebooks and hygiene checks.

## Re-blessing

When you change the model and the change is *intended*, the blessed numbers must be updated:

```julia
ECOSISTEM_BLESS=true julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_canonical.jl"])'
```

This rewrites `reference.toml`. **The diff to that file is the deliverable** — it is the machine-checked
statement of what your change did to the model's output, and it belongs in the pull request alongside
the code. Review it before committing:

- Did exactly the results you expected to move, move?
- Did anything move that you did not expect? That is the finding, not the noise.
- Are the new numbers still physically sensible?

⚠️ **Re-blessing is not a way to make a failing test pass.** A canonical failure means the output
changed; your job is to explain *why* before recording the new value. If you cannot say why, do not
bless it.

### ⚠️ A dependency can move the numbers without the model changing

**2026-08-15 — `simulated/total_abundance` and `simulated/abundance_by_species` were re-blessed for a
reason that was not a model change**, and the episode is recorded here because the next one will look
identical from the outside.

`Distributions` went **0.25.127 → 0.25.130**, which changed [which Poisson sampler is used
when](https://github.com/JuliaStats/Distributions.jl/commit/c6c9d824a6ebba82b47598cba76dcc5bf7fa0a6b).
Births are drawn once per species per cell per timestep by `rand(getrng(eco, sp), Poisson(...))`
(`src/dynamics.jl`), so a different sampler consumes the stream differently and every subsequent draw
shifts. The result was ~0.03 %: one individual per species after two simulated years.

⭐ **What made it diagnosable, and what to repeat:**

- ✅ **Construction was bit-identical** — habitat values, cell size, active mask, and an initial
  abundance of exactly 1 000 000. So the divergence was in the *draws*, not in what was built.
- ✅ **Only 2 of 58 blessed values moved.** `occupied_cells` and `species_surviving` were unchanged,
  as were all of `varying/*`, `intervention/*` and `realdata/*`. A model change would not be so
  selective; a change to the random stream is exactly that selective.
- ✅ **The decisive test was a pin.** Re-running the fixture in a throwaway environment with
  `Distributions` held at 0.25.127 reproduced the blessed numbers **exactly**, element for element.
  That is what turns "probably a dependency" into a finding.

🔴 **The general lesson: a blessed absolute count is only reproducible against a fixed sampler.** These
numbers are still the right thing to check — they caught this, which is the point — but a move in them
means "the output changed", not "the model changed", and the two have to be told apart before
blessing. ⭐ Derived, less brittle quantities (`occupied_cells`, `species_surviving`) rode the bump
untouched and are worth preferring where a test can be written either way.

## The five kinds, and why there are five

| file | input | what only it can catch |
|---|---|---|
| `test_simulated.jl` | synthetic, **uniform, square, static** | a change in the **model**: birth/death arithmetic, dispersal, the resource ratio |
| `test_varying.jl` | synthetic, **varying, non-square, warming** | anything **spatial** or **temporal**: grid orientation, gradient fields, the layer-change subsystem |
| `test_intervention.jl` | synthetic, varying, plus interventions | the **intervention** mechanism: the cell-selection stream, and the clock → intervention → layer ordering |
| `test_categorical.jl` | synthetic, **categorical**, non-square, two penalties | a change in the **categorical branch**: the tolerance's `penalty`, how it reaches suitability, and the soft-vs-hard distinction |
| `test_realdata.jl` | real WorldClim rasters | a change in the **data path**: readers, units, accumulation periods, resampling |

Keeping them apart is the point. When a number moves, which file moved tells you immediately whether
to look at the simulation, the grid, the intervention machinery, the categorical branch or the
reader — and a simulated-only suite cannot see a data-path regression at all.

🔴 **`test_categorical.jl` covers a branch that had NO canonical coverage at all**, and it became
worth having the moment the soft/hard exclusion distinction stopped being a *type* (`Match` vs
`LCmatch`, scoring 0.5 and 0) and became a `penalty` **number** on the tolerance. A number is far more
expressive and correspondingly easier to change by accident — a default, a comparison, the exponent it
enters the demographics through — and every one of those changes would previously have left
`reference.toml` untouched. ⭐ Two runs identical but for the penalty are what makes it detectable: a
change that collapses the distinction moves them *together*, which the file's `hard ≠ soft` assertions
refuse whatever the blessed numbers say.

### ⚠️ Two traps that file hit, and any new fixture can

Both produced numbers that were stable when run one way and different when run another — which is the
worst kind of canonical failure, because re-blessing makes it go away until next time.

1. 🔴 **Do not build two ecosystems on one habitat and simulate both.** An `Ecosystem` *shares* its
   habitat with the caller rather than copying it, so the second run inherits state the first left.
   ✅ Measured: it made the blessed totals depend on the **thread count** (agreeing at 2, 4 and 8 and
   differing at 1). ⭐ Build, simulate and finish one run before starting the next — which is what
   `test_simulated.jl`'s own reproducibility check already does. Passing a *built layer* to a fresh
   `GridHabitat` is fine, because it is copied in.
2. 🔴 **Bless the fixture's INPUT, not only its outputs**, whenever the input is generated rather than
   written down. `NicheSpec` draws its layout from the global stream, and a canonical file is
   `include`d into a shared process when re-blessing but runs in its own when checking — so a
   self-consistency check ("each run got the same map") cannot see the map itself moving between the
   two. Blessing `categorical/regime_map` is what makes every number conditional on a *recorded*
   input. ✅ The gate for this is cheap and worth running on any new file: three consecutive checking
   runs, plus a re-bless that leaves `reference.toml` byte-identical.

🔴 **`test_varying.jl` and `test_intervention.jl` exist because of the same class of hole
`test_realdata.jl` did.** Until they were written, *every* blessed number came from an environment
where **every cell was identical and nothing changed over time**, on a **square** grid. That is
structurally incapable of detecting a transposed index (see CLAUDE.md — a past bug had exactly that
shape), any part of the layer-change subsystem, or which cells an intervention selected. Four steps
of the axis/units work passed with `reference.toml` unmoved, which was *correct* for those steps —
but the experiment could not have told the difference.

⭐ They are **additive**: no key in them replaces one in `test_simulated.jl`, so the baseline those
original numbers provide is preserved. The shared fixture is
[`test/varyingcase.jl`](../varyingcase.jl), which the MPI cross-rank test now decomposes too — that
comparison ran on the same uniform square grid, where every decomposition looks alike.

## Heavy downloads are local-only

⚠️ **The GitHub runners already download close to their time budget**, so the canonical reads that
need CHELSA `BioClimPlus` layers — `pet_penman_mean`, `cmi_mean`, `gdd0`, `swb`, which together pin
the "when is an accumulated layer divided?" table from `docs/src/units.md` — are guarded by
`heavydata()`. It is `true` locally and `false` on a runner; `ECOSISTEM_HEAVY_DATA=true` forces them
on, `=false` suppresses them locally.

⭐ **Skipping them cannot lose their blessed values.** `writereference` *merges*, so a re-blessing run
that skipped them leaves their keys exactly as they were blessed locally. ✅ Verified by putting a
sentinel key in `reference.toml` and re-blessing with `RUNNER_OS` set: the key survived.

⭐ **And skipping does not leave the underlying fact unchecked.** Which of those layers divides is a
*catalogue* property, needing no download at all — it is asserted for all five in
`test_LayerCatalogue.jl`, and so runs on every CI build. Only the **read values** are gated. That split
matters: it is the catalogue half that a change to an axis's canonical unit moves.

⚠️ An unblessed key reports as `@test_broken`, not a failure, so the first local run after adding one
tells you to bless rather than failing the suite.

⚠️ **Both files must stay fast enough to run on demand.** A canonical test earns its keep by being
*run*, and one that takes an hour will not be.

The ecological property tests that used to sit in this directory — "commonest at the temperature
optimum", "abundance scales with area", and so on — now live in
[`examples/biodiversity.jl`](../../examples/biodiversity.jl). They are not canonical tests: they
assert *properties* the model must never stop satisfying, which no re-blessing can silence, rather
than recording what it happened to produce. They are also the slow ones (~20 whole simulations, about
3 minutes). Nothing runs that file automatically; run it directly.

⭐ **`test_realdata.jl` exists because of a real hole.** Until it was written, *nothing in the entire
suite read a monthly climate layer*. Monthly totals were being divided by a fixed 30.4375-day month —
wrong for every real month, and 7.7% wrong for February — and that error could have been introduced,
or corrected, without a single test moving. A "number-moving" change that leaves every gate green
means **missing coverage**, not safety.

## Writing one

```julia
include("canonical.jl")
using .Canonical

canonical("simulated/total_abundance", sum(abundances))
```

- **Strip units at the call site**: `ustrip(u"mm/d", x)`. A blessed `7.5398` means nothing on its own,
  and stripping explicitly is what pins which unit it is in. The helper refuses a `Quantity` and says
  so.
- **Name as `file/thing`**, so `reference.toml` sorts by file and a diff stays readable.
- **Prefer several specific numbers to one summary.** Per-species totals catch a change that
  redistributes abundance while preserving the grand total; a single sum does not. The monthly rate
  vector in `test_realdata.jl` is there for exactly this reason — an annual total would have hidden
  the February bug entirely, since dividing every month by the same figure conserves the year.
- **Keep ordinary assertions alongside the blessed ones.** A blessed number tells you *something
  changed*; a property tells you *the model is still right*. Re-blessing silences the first and must
  never be able to silence the second. Every file here asserts positivity, plausible magnitudes and
  the qualitative ecology as well as its blessed figures.
- **Seed anything stochastic.** `Ecosystem(..., seed = …)` gives one deterministic RNG stream per
  species, so results do not depend on thread or MPI-rank count. Draw initial abundances from a seeded
  generator too, or the starting state differs between runs and nothing downstream is blessable.
- **Choose data locations by coordinate, not array index.** An index is meaningless if the source
  changes resolution — and index `(400, 300)` turned out to be ocean, where every value is `NaN` and
  every assertion vacuously true.

## Notes

- A value with no blessed counterpart reports as **Broken**, not failed: adding a canonical test
  should not break the build before you have blessed it.
- `reference.toml` is **merged**, not replaced, on blessing. A partial run therefore cannot silently
  delete the blessed values of files it did not execute.
- The default tolerance is tight (`rtol = 1e-8`). Widen it only where a result genuinely is not
  reproducible to more digits, and say why at the call site.
