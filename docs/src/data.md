# Using the data pipeline

This page summarises the basic usage of the data pipeline.
For more information see the [`DataRegistryUtils.jl`](https://github.com/ScottishCovidResponse/DataRegistryUtils.jl)
repository and [documentation](https://scottishcovidresponse.github.io/DataRegistryUtils.jl/stable/).

## Overview

This package uses `DataRegistryUtils.jl`, which provides a Julia interaction with the [SCRC data registry](https://data.scrc.uk). An overview of the pipeline itself can be found [here](https://scottishcovidresponse.github.io).

## Basic flow

The basic flow is as follows:

1. Write a `config.yaml` file, which specifies the data products which you need.
2. Download the specified data products from the registry.
3. From within your model script, use the read functions provided by `DataRegistryUtils` or access directly through the downloaded database. An optional SQL file can be provided to describe how a data product should be displayed.
4. Write model outputs back to pipeline (WIP).

### Config file

See `examples/Epidemiology/data_config.yaml` for an example.

### Downloading data

Use the `download_data_registry` function from `SimulationData.jl`. For example:
```julia
julia> using DataRegistryUtils
[ Info: Precompiling DataRegistryUtils [1a5903ec-9658-4297-bb6d-314c615f2e02]

julia> db = fetch_data_per_yaml("examples/Epidemiology/data_config.yaml")
processing config file: examples/Epidemiology/data_config.yaml
 - hint: use the 'verbose' option to see more stuff
 - downloading human/demographics/population/scotland/1.0.1.h5, please wait...
WARNING - HASH DISCREPANCY DETECTED:
 server file := human/demographics/population/scotland/1.0.1.h5
 hash: 835fa268cb115510d3195b957fe8dd61665e5f6b
 downloaded: ./out/human/demographics/population/scotland/1.0.1.h5
 hash: b238974dd07fbaf4495ba51282f6ad44db486e9a
 - files refreshed, but issues were detected.
Dict{Any,Any} with 7 entries:
  "fixed-parameters/T_hos"                              => Dict{String,Any}("T_hos"=>Dict{String,Any}("value"=>5,"type"=>"point-estim…
  "geography/lookup_table/gridcell_admin_area/scotland" => Dict{Any,Any}("/conversiontable/scotland"=>NamedTuple{(:grid1km_id, :grid1…
  "human/demographics/population/scotland"              => Dict{Any,Any}("/travel_to_work_area/age/genders"=>[87.0 93.0 … 22.0 37.0; …
  "fixed-parameters/T_rec"                              => Dict{String,Any}("T_rec"=>Dict{String,Any}("value"=>11,"type"=>"point-esti…
  "human/infection/SARS-CoV-2/*"                        => Dict{String,Any}("symptom-probability"=>Dict{String,Any}("value"=>0.692,"t…
  "records/pollution"                                   => Dict{Any,Any}("/array"=>[NaN NaN; NaN NaN; … ; NaN NaN; NaN NaN])
  "prob_hosp_and_cfr/data_for_scotland"                 => Dict{Any,Any}("/cfr_byage"=>NamedTuple{(:p_h, :cfr, :p_d),Tuple{Float64,Fl…

```

The above will download the required data to a specified path or an automatic folder if none defined.

### Reading data

An example of reading a parameter:

```julia
db["human/infection/SARS-CoV-2/*"]["symptom-probability"]
```

More examples can be found on DataRegistryUtils [`README`](https://github.com/ScottishCovidResponse/DataRegistryUtils.jl#readme) or [`examples folder`](https://github.com/ScottishCovidResponse/DataRegistryUtils.jl/tree/main/examples).
