# The paper's figures

`examples/paper.jl` regenerates the computational figures of the EcoSISTEM paper from the package
itself: figures 3, 5, 6 and 7 (the analytic niche curves, the island experiments, dispersal
distance, and temperature optima and niche widths), and figures 8 and 9 (the circular continent
and the invasions on it). Each is written as a PDF under the filename the paper's LaTeX source
includes, so the paper needs no change to pick them up. Every value the figures are built from is
in `parameters.toml` here, a set for each scale (`large`, the published experiments, and `small`),
and the paper's Tables S1 and S2 are written from the set that ran, as `IslandParameters.tex` and
`ContinentParameters.tex`, complete `tabular`s the paper `\input`s. Beside the PDFs it writes
`figures.toml`, recording the package version, commit, Julia version, scale, seed, thread count and
the wall time of every figure, and `continent.jld2`, the per-cell diversity maps figure 8 is drawn
from, so that figure can be redrawn without rerunning.

The pieces are in this directory: `parameters.toml`, `common.jl` (scale, seed, parameters, output
directory, timing, the tables, the sidecar), `niche.jl` (figure 3), `island.jl` (figure 5),
`dispersal.jl` (figure 6), `optima.jl` (figure 7) and `continent.jl` (figures 8 and 9), with
`pipeline.jl` and `paper.yaml` for a run through the data pipeline (below). The island builders come
from `examples/models/ecosystems.jl`.

## Running it in full

Full scale is the published configuration. Two of the figures dominate the time: figure 5's ten
replicates of a 5,000-species pool, about half an hour on twelve threads, and figure 8's 50,000
species on a 100 by 100 continent for 125 years, a matter of hours. On one machine:

```sh
julia -t 12 --project=examples examples/paper.jl /path/to/output
```

The first argument is the output directory; without it, `examples/paper/outputs/`, which this
directory's `.gitignore` keeps out of the repository. It is printed at the start and the end.

`ECOSISTEM_PAPER_FIGURES=NicheWidth,Abundance` draws only the figures named, by their file
names, so the continent can be left for a machine that has the hours. The continent run prints a
line every fifth simulated year with the projected finish and how many individuals the generalists
and the specialist hold; redirected to a file, those lines are flushed as they are written.

The output directory's `figures.toml` is overwritten by each run and lists only the figures that
run drew, so a run of some figures into a directory holding others leaves the sidecar describing
the latest run alone.

## Through the data pipeline

`paper.yaml` draws the same figures and tables as a FAIR data pipeline run, in namespace
`EcoSISTEM`. `fair pull` registers `parameters.toml`, from the repository's own file, as
`Examples/Parameters`, at a new version whenever it changes; `fair run` then reads it from the
registry and registers every output as derived from it - each figure as `Examples/Figures/<name>`,
each table as `Examples/Tables/<name>`, and `figures.toml` as `Examples/Config`:

```sh
julia --project=examples/paper -e 'using Pkg; Pkg.instantiate()'   # once
fair pull --local examples/paper/paper.yaml
fair run --local examples/paper/paper.yaml
```

It runs `pipeline.jl` here, in this directory's own environment (`Project.toml`: EcoSISTEM from this
checkout, DataPipeline, and what the figures load), and `pipeline.jl` includes `examples/paper.jl`,
so the two draw exactly the same; a plain run never loads DataPipeline.

The figures are drawn into `outputs/` here, as for a plain run with no argument, and each run
replaces the last's files there: the registered copies are in the data store, and `continent.jld2`,
which is not registered, is the latest run's. `ECOSISTEM_SCALE` and `ECOSISTEM_PAPER_FIGURES` work
as for a plain run; set them before `fair run`, which passes them on.

## Under the test suite

`test/extras_examples.jl` runs `examples/paper.jl` like every other top-level example, with
`ECOSISTEM_SCALE=small`: a few species and years, the same code and layout, every figure drawn, in
about ten seconds, into a temporary directory that is cleaned up. It never skips. To run just that:

```sh
ECOSISTEM_SCALE=small julia --project=examples examples/paper.jl
```
