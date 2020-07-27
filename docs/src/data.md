# Using the data pipeline

This page summarises the basic usage of the data pipeline.
For more information see the [`data_pipeline_api`](https://github.com/ScottishCovidResponse/data_pipeline_api)
repository.

## Basic flow

The basic flow is as follows:

1. Write a `config.yaml` file, which specifies the data products which you need.
2. Download the specified data products from the registry.
3. From within your model script, use the `StandardAPI` from [`SimulationData.jl`](https://github.com/ScottishCovidResponse/SimulationData.jl) to read the data.
4. Use the same API to write model outputs.

These steps are outlines in more detail below.

### Config file

See `examples/Epidemiology/data_config.yaml` for an example.

### Downloading data

There is an [open issue](https://github.com/ScottishCovidResponse/SimulationData.jl/issues/13) in
`SimulationData.jl` to wrap the download script. Until then, data can be downloaded manually by the
following procedure:

```bash
# Clone the data_pipeline_api repo (at a specified commit where this has been tested)
git clone -b d2efc99 https://github.com/ScottishCovidResponse/data_pipeline_api.git
cd data_pipeline_api

# Install all the requirements, here by creating a virtualenv and installing the package
# through pip.  Or you could use conda instead.
virtualenv -p python3 venv
source venv/bin/activate
pip install data_pipeline_api

# Run the download script on your config.yaml file.
# Note: this runs the one in the working directory, not the one installed through pip.
python -m data_pipeline_api.registry.download --config <path_to_config_file>
```

The above will download the required data to the path specified in the config file.
Note: this path is relative to the config file location, not to where you ran the download script.

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

Use a single instantiation of `StandardAPI` for the entire model run, for reading and writing.

### Writing model outputs

Example of writing an array:

```julia
write_array(api, "simulation-outputs", "final-abundances", abuns)
```

## Wrapper structure

This package uses `SimulationData.jl`, which provides a Julia wrapper around the Python package
`data_pipeline_api`.
