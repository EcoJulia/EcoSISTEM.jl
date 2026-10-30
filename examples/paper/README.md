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
`pipeline.jl` and `paper.yaml` for a run through the data pipeline, and `hpc/` for that run as a
job on a cluster (both below). The island builders come from `examples/models/ecosystems.jl`.

## Running it in full

Full scale is the published configuration. Two of the figures dominate the time: figure 5's ten
replicates of a 5,000-species pool, about half an hour on twelve threads, and figure 8's 50,000
species on a 100 by 100 continent for 125 years, a matter of hours. On one machine:

```sh
julia -t 12 --project=examples examples/paper.jl /path/to/output
```

The first argument is the output directory; without it, `examples/paper/outputs/`, which this
directory's `.gitignore` keeps out of the repository. It is printed at the start and the end.

`ECOSISTEM_PAPER_FIGURES=NicheWidth,IslandTesting` draws only the figures named, by their file
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

## On an HPC cluster

`hpc/` runs that `fair run` as a Slurm job on one node, set up for the `smp` nodes of the University
of Glasgow's MARS cluster: 128 cores and 1 TB each. Its `#SBATCH` lines name MARS's partition and
project `project0000`, and `hpc/env.sh` MARS's juliaup module and project storage, so change those
for another machine. A full run spends nearly all its time on figure 8 (eighteen hours on twelve
threads of a laptop); the job asks for three days, and the continent's progress lines give its
projected finish early.

Everything large goes under `ECOSISTEM_HPC_ROOT`, by default `/mnt/data/project0000/$USER`: the
Julia depot, and with it EcoSISTEM's asset cache and DataPipeline's copy of `fair`, then the local
registry and its data store. `env.sh` also picks the juliaup channel (`JULIAUP_CHANNEL`, `1.13.3`)
and the registry's port (`ECOSISTEM_HPC_PORT`, 8030; a node is shared, so not 8000). Set any of the
three in `~/.bashrc` rather than editing `env.sh`: `fair run` refuses a clone with changes.

Once, on the login node - the clone, at the release to run, and who the registry records the runs
as:

```sh
module load apps/juliaup/1.18.9
juliaup add 1.13.3
export ECOSISTEM_HPC_ROOT=/mnt/data/project0000/$USER
git clone https://github.com/EcoJulia/EcoSISTEM.jl.git $ECOSISTEM_HPC_ROOT/EcoSISTEM.jl
cd $ECOSISTEM_HPC_ROOT/EcoSISTEM.jl
git checkout v0.9.1
mkdir -p $ECOSISTEM_HPC_ROOT/fair
cp examples/paper/hpc/fair-user.yaml $ECOSISTEM_HPC_ROOT/fair/user.yaml
```

Fill in `user.yaml` - your name, email, and your ORCID or GitHub page - then, from the root of the
clone, set the rest up as a job. It instantiates and precompiles this directory's environment on an
`smp` node, installs DataPipeline's copy of `fair` (`DataPipeline.install_cli()`) and the registry,
writes a CLI configuration from `user.yaml` and the paths above, and runs `fair init --using` with
it:

```sh
sbatch examples/paper/hpc/setup.bash
```

If it fails part way, with a download cut short say, submit it again: it keeps what it finished.

Then the small-scale run, which tests every piece in minutes, and the full one:

```sh
ECOSISTEM_SCALE=small sbatch --time=01:00:00 examples/paper/hpc/run.bash
sbatch examples/paper/hpc/run.bash
```

Each starts the registry, runs `fair pull` and `fair run` as above, and stops it. Both runs stay in
the registry under the same names, each `Examples/Config` recording its scale. `run` is given no
`--local`, so it refuses a commit the clone's remote does not have. Each job writes
`<name>-<job ID>.out` and `.err` in the clone, with the `Manifest.toml` it ran with beside them,
since the manifest is not part of the repository.

A run job that is killed, by its time limit or for running out of memory, cannot stop the registry,
which stays recorded as running on that job's node. The next run job, if it lands on another node,
refuses with "The local registry installed in '...' is recorded as running on the host '...'":
the killed job has gone, so delete the file the message names, then submit again.

## Under the test suite

`test/extras_examples.jl` runs `examples/paper.jl` like every other top-level example, with
`ECOSISTEM_SCALE=small`: a few species and years, the same code and layout, every figure drawn, in
about ten seconds, into a temporary directory that is cleaned up. It never skips. To run just that:

```sh
ECOSISTEM_SCALE=small julia --project=examples examples/paper.jl
```
