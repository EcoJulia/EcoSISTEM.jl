# Using the data pipeline

This page summarises the basic usage of the data pipeline.
For more information see the [`data_pipeline_api`](https://github.com/ScottishCovidResponse/data_pipeline_api)
repository.

## Wrapper structure

This package uses `SimulationData.jl`, which provides a Julia wrapper around the Python package
`data_pipeline_api`. The functions in `SimulationData.jl` are mostly designed to have the same signatures
as those in `data_pipeline_api`, except one calls `f(api, ...)` rather than `api.f(...)`.

## Basic flow

The basic flow is as follows:

1. Write a `config.yaml` file, which specifies the data products which you need.
2. Download the specified data products from the registry.
3. From within your model script, use the `StandardAPI` from [`SimulationData.jl`](https://github.com/ScottishCovidResponse/SimulationData.jl) to read the data.
4. Use the same API to write model outputs.

Use a single instantiation of the API for the entire model run, for reading and writing.

These steps are outlined in more detail below.

### Config file

See `examples/Epidemiology/data_config.yaml` for an example.

### Downloading data

Use the `download_data_registry` function from `SimulationData.jl`. For example:
```julia
julia> using SimulationData
[ Info: Precompiling SimulationData [3d44aec0-db1b-416b-9784-2428c815ea7f]
[ Info: Found data-pipeline-api v0.7.0

julia> download_data_registry("examples/Epidemiology/data_config.yaml")
┌ Info:     Downloading from registry: https://data.scrc.uk/api/
│     Using config file: examples/Epidemiology/data_config.yaml
└     ...
[ Info: Registry download done
```

The above will download the required data to the path specified in the config file.
Note: this path is relative to the config file location, not to where you ran the download.

### Reading data

An example of reading a parameter:

```julia
api = StandardAPI("path/to/config.yaml", "test_uri", "test_git_sha")
read_estimate(
    api,
    "human/infection/SARS-CoV-2/asymptomatic-period",
    "asymptomatic-period"
)
```

Note: you must call `close(api)` at the end of the model run, in order to write out the access file.
Alternatively, use the following idiom which automatically closes out the API:

```julia
StandardAPI("path/to/config.yaml", "test_uri", "test_git_sha") do api
    # Do everything in here with `api`
end
```

### Writing model outputs

Example of writing an array:

```julia
write_array(api, "simulation-outputs", "final-abundances", abuns)
```

This will write outputs to the data folder specified in the API config file.
