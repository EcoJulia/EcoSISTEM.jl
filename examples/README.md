# Examples

**Every `.jl` file in the top level of this directory is run as part of the test suite.** They are
tests as much as demonstrations: an example that is never executed stops working without anyone
noticing, and these are what show the package's interfaces actually do what the documentation says.

Run them on their own:

```julia
julia --project -e 'using Pkg; Pkg.test(test_args = ["extras_examples.jl"])'
```

⚠️ **Do not add an arbitrary example here.** Anything at this level must run, in reasonable time,
cleanly enough to pass `@test_nowarn` — which fails on *any* output to stderr, not just warnings —
and without needing data or hardware a test runner will not have. If your example does not meet all
of that, put it in a subdirectory instead.

## What is where

| | |
|---|---|
| `*.jl` (top level) | **run by the test suite on every full run** |
| `interventions/` | the scenario library the top-level `interventions.jl` exercises |
| `landscapes/` | the five landscape configurations the top-level `landscapes.jl` exercises |
| `models/` | the published experiments and diversity measures, included by the top-level `models.jl` |
| `other/` | examples that are *not* run — long runs, ones needing data or a cluster, work in progress |
| `HPC/`, `pipeline/`, `scripts/`, `benchmarks/` | tooling and one-off scripts, not examples |

⭐ A subdirectory is included by a top-level file when it should be tested — `interventions.jl`
includes each file in `interventions/`, so that library *is* covered while staying organised. That
is the pattern to copy: put the substance in a subdirectory and add one top-level file that runs it.

## Scale

⭐ **`ECOSISTEM_SCALE` decides how much work these do.** `test/runtests.jl` sets it to `small`, so
anything reached through `Pkg.test` runs small; a direct
`julia --project=examples examples/landscapes.jl` leaves it unset, which reads as full scale. An
explicit `ECOSISTEM_SCALE=large` wins even under the suite.

| | under the suite | run directly |
|---|---|---|
| `landscapes.jl` | small grids, tens of species | the originals' sizes, up to 50 × 50 × 10,000 species |
| `interventions.jl` | 10 species, 10,000 individuals | the published 100 species, 100,000,000 individuals |
| `biodiversity.jl` | **skipped entirely** | ~20 whole simulations, ~3 min |
| `models.jl` | **skipped entirely** | 10–50-year runs that draw figures |
| everything else | runs as written | runs as written |

⚠️ `biodiversity.jl` and `models.jl` skip rather than shrink because there is no smaller version of
them worth running — they are the most expensive things in the repository. They are still top-level
files, so they *are* found by the test suite; they simply decline to run.

⭐ **`biodiversity.jl` and `models.jl` investigate the same ecology to different ends**: the first
asserts the end state as property tests, the second follows the same experiments *through time* with
a full suite of Diversity.jl measures and draws the figures. `models.jl` also reuses
`interventions/` verbatim for its climate, habitat-loss and invasion scenarios.

⭐ **Supply is derived from the population**, in both `landscapes/` and `interventions/`, so every
configuration starts at its carrying capacity. 🔴 Scaling the population alone does not work: it
leaves capacity where it was and the run collapses straight back to it — measured, an island seeded
with 100,000,000 individuals against the small supply crashed to **124**. At full scale the derived
figure comes out at the published `4.5e11 kJ/km²/day`.

⚠️ **Every top-level example is a `module`.** `test/extras_examples.jl` includes them all into one
module, several define `runscale()`/`CONFIGURATIONS`, and Julia 1.12 lets a later `const` silently
overwrite an earlier one — so without the wrappers the examples would quietly reconfigure each other
depending on the order they ran in.

⚠️ `other/` is a holding area, not an archive. Its contents are not compiled, not run, and will rot;
if something there is worth keeping, give it a top-level runner.
