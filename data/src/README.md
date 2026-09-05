# `data/src` - the code that makes and checks what the repository keeps

Everything here is **run deliberately, by hand**, to produce or check an artefact, a measurement or
an invariant that the repository keeps. That is the whole membership rule. These scripts are not
tests, not examples, and **nothing runs them automatically** - no CI job, no `Pkg.test` set, no
`runtests.jl` discovery.

Outputs live in `data/`, so a generator and what it generates are one directory apart rather than
scattered across `docs`, `examples` and `test`. An audit that writes no file still belongs here: what
it produces is an answer, and it is run the same way for the same reason.

Four of these are wired to a `test/clean_*.jl` gate, which includes the script by explicit path and
asserts the audit is clean. That is the only thing here that runs automatically, and it runs the
*check*, never the regeneration.

## Do not run these as a batch

They sit together because they are the same *kind* of thing, not because they are a suite. Several
rewrite a tracked file, so running the folder would overwrite committed artefacts wholesale, and the
damage would look like ordinary uncommitted work.

A script here must not overwrite a tracked file unless it was asked to. `architecture.jl` is the
pattern to copy: it audits by default and writes only under `--fix`.

`metadata.jl` is the one to be most careful with. It reformats the entire source tree and rewrites
`codemeta.json` in a single unprompted pass, so it changes more tracked files than everything else
here combined.

## What is here

| script | makes | notes |
|---|---|---|
| `architecture.jl` | `data/architecture.md` | audits by default; `--fix` repairs inheritance edges. Also run by `test/clean_Architecture.jl` |
| `naturalearth_regions.jl` | `data/NaturalEarth/regions.csv` | 2444 rows from the Natural Earth 1:10m shapefiles; needs a download |
| `rate_table.jl` | `data/rate_table.csv` | birth and death rates against the demographic parameters |
| `bioclimplus_units.jl` | a measurement, not a file | the evidence establishing the units of two CHELSA-BIOCLIM+ layers |
| `benchmarks/` | timings, CSVs and plots | `run_benchmarks.jl` drives worker processes, each launched with this directory as its project |
| `typeorder.jl` | an answer | is every type declared after what its declaration needs? Gated by `test/clean_TypeOrder.jl` |
| `constructororder.jl` | an answer | is every constructor written after its type's own file? A constructor above it is **silently discarded**. Gated by `test/clean_ConstructorOrder.jl` |
| `api_surface.jl` | `data/api_surface.md` | every name the package defines, by visibility and kind. Regenerable, so not committed |
| `overloads.jl` | `data/overloads.md` | every method we define on a generic we do not own, by reflection. Regenerable, so not committed |
| `metadata.jl` | `codemeta.json`, and a formatted tree | resolves all four environments, formats the package, runs the metadata crosswalk |

## Running them

Every script here runs in this directory's own environment:

```
julia --project=data/src data/src/rate_table.jl
```

`data/src/Project.toml` carries what these scripts import, and also the three weak dependencies - MPI,
Phylo and RasterDataSources - so every EcoSISTEM extension loads and the audits can see the types
that live in them.

That matters beyond convenience. `architecture.jl` decides whether a name is ours or an upstream
package's by asking `Base.identify_package`, which sees only the **direct** dependencies of the
active project. Run under an environment where Diversity and EcoBase are merely transitive, it
reported seven differences that do not exist - and `--fix` would have written them into the
document, deleting real inheritance edges. So add a dependency to this file rather than reaching for
another project, and be suspicious of an audit that suddenly has a lot to say.

Because nothing runs these scripts, a path one computes is checked by nothing until a person runs
it. Two are relative to this directory - `rate_table.jl` writes to its parent, and
`benchmarks/run_benchmarks.jl` launches its workers with this directory as their project - so moving
a script means re-checking what each resolves to.
