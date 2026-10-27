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

`examples/pipeline/` has the runs over Scotland, with an environment of their own
(`examples/pipeline/Project.toml`), since DataPipeline is a weak dependency: `inputs.jl` names the
specs, `fetch.jl` fetches one source's files, with a config for each source (`fetch-worldclim.yaml`,
`fetch-naturalearth.yaml`, `fetch-era5.yaml`), and `model.jl` and `model.yaml` are the model run. A
fetch config takes its source's namespace for the run's inputs and outputs, registers it, and
writes the files and their records:

```yaml
run_metadata:
  default_input_namespace: ECMWF
  default_output_namespace: ECMWF
  script: |
    julia --project=examples/pipeline examples/pipeline/fetch.jl era5
register:
- namespace: ECMWF
  full_name: European Centre for Medium-Range Weather Forecasts
  website: https://www.ecmwf.int/
write:
- data_product: era5/data/*
  description: An ERA5 extract from the Climate Data Store
  file_type: nc
- data_product: era5/provenance/*
  description: EcoSISTEM's provenance record of a file under era5
  file_type: toml
```

The model config reads each source's files and their records in its namespace; a `*` matches one
segment of a name, so `era5/*/*` takes both and leaves the deposit, which has no file:

```yaml
read:
- data_product: era5/*/*
  use:
    namespace: ECMWF
```

With `fair-cli` installed (see DataPipeline.jl's documentation), and the environment instantiated
once with `julia --project=examples/pipeline -e 'using Pkg; Pkg.instantiate()'`:

```sh
fair init
fair run examples/pipeline/fetch-worldclim.yaml      # fetches and registers WorldClim's files
fair run examples/pipeline/fetch-naturalearth.yaml
fair run examples/pipeline/fetch-era5.yaml
fair run examples/pipeline/fetch-worldclim.yaml      # stops at once: everything is registered
fair pull examples/pipeline/model.yaml
fair run examples/pipeline/model.yaml                # reads them from the store, registers the figure
```

Fetching the ERA5 month needs a Climate Data Store key in `~/.cdsapirc`; a file already in the asset
cache needs none.
