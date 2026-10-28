# Data pipeline

The [FAIR Data Pipeline](https://www.fairdatapipeline.org/) records where every input of a run came
from and what the run made of it, in a registry that can be pushed to the online
[data registry](https://data.fairdatapipeline.org/). EcoSISTEM takes part through
[DataPipeline.jl](https://github.com/FAIRDataPipeline/DataPipeline.jl), a weak dependency: with it
loaded, a spec given a [`EcoSISTEM.Track`](@ref) as its `track` keyword has its files go through a
pipeline run.

## Two runs, the same specs

A **fetch run** fetches the files its specs read and registers each one with its source, its
authors and its licence, from the catalogue; a **model run** then builds the same specs, and every
file is read from the pipeline's store rather than fetched. Which a run is depends on the function
called, never on what is there: [`fetchfiles`](@ref EcoSISTEM.fetchfiles) on a tracked spec is the
fetch run's, and building from one only reads.

```julia
using EcoSISTEM, DataPipeline
using EcoSISTEM: Track
temperature(handle) = SourceSpec(WorldClim{BioClim}, :bio1,
                                 track = Track(handle, "worldclim"))
```

A track names the run's `handle`, from `DataPipeline.initialise()`, and a root: each file is the data
product `<root>/data/<stem>`, the stem its name without its extension, and its provenance record
`<root>/provenance/<stem>`. A handle of `nothing` tracks nothing, so one script runs with the
pipeline or without it:

```julia
launched = haskey(ENV, "FDP_CONFIG_DIR")      # started by `fair run`
if launched
    nothing_to_fetch = DataPipeline.rehearse() do handle
        return isempty(EcoSISTEM.missingfiles(inputs(handle)))
    end
    nothing_to_fetch && exit()
end
handle = launched ? DataPipeline.initialise() : nothing
EcoSISTEM.fetchfiles(inputs(handle))
isnothing(handle) || DataPipeline.finalise(handle)
```

[`EcoSISTEM.missingfiles`](@ref) says what is still to do, and why, without fetching: a fetch run
with nothing to do stops in the rehearsal, before a code run exists.

## What a file is registered as

The registry derives every output of a run from every input of it, so a fetch run registers one
source's files. Its one input is the source's deposit, registered with no file and named by the
catalogue's identifier in the source's namespace (`datasets.csv`'s `Namespace`), so every user
registers one source's files alike. Its outputs, each derived from that deposit alone, are:

  - every file, whether downloaded as its provider publishes it, made by the Climate Data Store to a
    request (as an ERA5 extract is), or given to a spec by you (`file`, `files` or `directory`). A
    download is the cited item's own file (`primary`) where the identifier names the data - a CHELSA
    deposit, Natural Earth's web page - and not where it names a paper describing the data, as
    WorldClim's does; an extract or a file of yours never is. A netCDF extract's `history`
    attribute, which records when and where the CDS converted it, is replaced by a note of the same
    length, so fetches of the same data are identical. A file of yours is of the catalogue's release
    unless the track says otherwise:
    `Track(handle, "cru", release_version = "4.8", release_date = "2024-05")`;
  - every file's provenance record, `<file>.provenance.toml`.

A fetch run asked to register a second source's files refuses. A provider can change a file without
a new release: a file whose bytes differ from those already registered under its release is
registered under that release with build metadata, `2.1.0+20230104`, the day its provider changed
it, so both stay traceable.

## The example

`examples/pipeline/` is a model run over Scotland and the fetch runs it reads from, with an
environment of its own (`examples/pipeline/Project.toml`), since DataPipeline is a weak dependency:
`inputs.jl` names the specs, `fetch.jl` fetches one source's files, with a config for each source
(`fetch-worldclim.yaml`, `fetch-naturalearth.yaml`, `fetch-era5.yaml`), and `model.jl` and
`model.yaml` are the model run, one species built from WorldClim's annual temperature and rainfall
within Natural Earth's outline of Scotland. Its environment is instantiated once, with
`julia --project=examples/pipeline -e 'using Pkg; Pkg.instantiate()'`; then, in Julia started with
`julia --project=examples/pipeline` from the root of a clone of this repository:

```julia
using DataPipeline
DataPipeline.install_cli!()      # where there is no `fair` installed
fair"registry install"           # the first time: the local registry
fair"init"                       # the first time in this clone: it asks who you are, and more
fair"run examples/pipeline/fetch-worldclim.yaml"      # fetches and registers WorldClim's files
fair"run examples/pipeline/fetch-naturalearth.yaml"
fair"run examples/pipeline/fetch-worldclim.yaml"      # stops at once: everything is registered
fair"pull --local examples/pipeline/model.yaml"
fair"run examples/pipeline/model.yaml"                # reads them from the store, registers the figure
```

`fair"..."` runs the pipeline's command-line tool with its arguments as on the command line, so
from a shell, with `fair-cli` installed from PyPI, the same commands are `fair registry install`,
`fair init` and so on. `pull` takes `--local` because it would otherwise look for what the
config's patterns match in a remote registry, where none of this is; a run matches them in the
local one. A run records the commit it ran from, and refuses unless the clone's git remote has
that commit and nothing is uncommitted; `--dirty` lets it run regardless, recorded as made from
changed code.

`fetch-era5.yaml` fetches and registers a month of ERA5 over Scotland the same way, an extract the
Climate Data Store makes to a request; fetching it needs a key in `~/.cdsapirc`, and the model does
not read it.

The rest of this page runs those steps when this documentation is built with the environment
variable `ECOSISTEM_DOCS_RUN_PIPELINE` set to `true`, as the package's own documentation build does,
against a registry that starts empty; built without it, the blocks only say so.

### Fetching each source

A fetch config takes its source's namespace for the run's inputs and outputs, registers it, and
writes the files and their records:

```@example pipeline
using EcoSISTEM
example = joinpath(pkgdir(EcoSISTEM), "examples", "pipeline")
print(read(joinpath(example, "fetch-worldclim.yaml"), String))
```

Its run fetches WorldClim's two layers where they are not held, and registers each file and its
record as derived from WorldClim's deposit. The first fetch run in a registry also registers the
catalogue's authors, so that each citation's authors are found by their identifiers; the rehearsal
the run starts with says which it would.

```@example pipeline
using DataPipeline
run_example = get(ENV, "ECOSISTEM_DOCS_RUN_PIPELINE", "") == "true"
run_example || println("Not run: ECOSISTEM_DOCS_RUN_PIPELINE is not set.")
# `fair` runs from the root of the repository, which the configs' paths start from
run_example && cd(pkgdir(EcoSISTEM)) do
    fair"run examples/pipeline/fetch-worldclim.yaml"
end
nothing # hide
```

Natural Earth's config is the same in its own namespace, and its run registers the zip of outlines
that Scotland's is read from:

```@example pipeline
run_example && cd(pkgdir(EcoSISTEM)) do
    fair"run examples/pipeline/fetch-naturalearth.yaml"
end
nothing # hide
```

Run again, WorldClim's fetch run finds every file registered as it would register it now, and
stops in its rehearsal, before a code run exists, registering nothing:

```@example pipeline
run_example && cd(pkgdir(EcoSISTEM)) do
    fair"run examples/pipeline/fetch-worldclim.yaml"
end
nothing # hide
```

### The model run

The model config reads each source's files and their records in its namespace, and writes the
figure in EcoSISTEM's:

```@example pipeline
print(read(joinpath(example, "model.yaml"), String))
```

`fair pull` finds what the config reads in the registry, and `fair run` runs `model.jl`, which
builds the habitat from the files in the store, simulates five years and registers the figure:

```@example pipeline
run_example && cd(pkgdir(EcoSISTEM)) do
    fair"pull --local examples/pipeline/model.yaml"
    fair"run examples/pipeline/model.yaml"
end
nothing # hide
```

### Where the figure came from

DataPipeline's `provenance` traces a data product back through the runs that made it to what each
of them read, and each one's source:

```@example pipeline
trace = run_example ? DataPipeline.provenance("EcoSISTEM:Scotland-plot") : nothing
```

The figure, read back from the store where the registry says it is: the abundance at the start and
after five years, beside the annual temperature and rainfall the habitat was built from.

```@example pipeline
run_example && include(joinpath(pkgdir(EcoSISTEM), "docs", "pipeline.jl"))
run_example ? registeredimage(trace.products[only(trace.roots)]) : nothing
```

The same provenance as a graph, drawn by `docs/pipeline.jl` from the registry's own record two
code runs back, which reaches the providers: grey for each item a provider's data is cited by,
with its authors; yellow for each data product; blue for each code run, with the code it ran and
that code's authors. A dashed arrow joins a cited item to a file that is the item itself, as
Natural Earth's zip is; WorldClim's files are cited by the paper describing their data.

```@example pipeline
run_example ? provenancegraph(trace, depth = 2) : nothing
```
