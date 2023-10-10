# Data pipeline

**DataPipeline.jl** is a [Julia](https://github.com/FAIRDataPipeline/DataPipeline.jl) package that provides functionality for the [FAIR Data Pipeline](https://www.fairdatapipeline.org/). The pipeline is intended to enable tracking of provenance of FAIR (findable, accessible, interoperable and reusable) data. We use it in examples under the `pipeline` folder for a run of a single species in Africa. 

See [here](https://www.fairdatapipeline.org/docs/data_registry/installation/) for details on what to install to set up the Data Pipeline. Once installed, the workflow is to initialise the pipeline in the repository, pull in the external data needed for the simulation (as described in the `AfricaRun.yaml`) and run the simulation. The output is also described in the `AfricaRun.yaml`, which is produced from the corresponding run file `AfricaRun.jl`. The output and provenance can then be pushed back to the online [data registry](https://data.fairdatapipeline.org/) to be inspected further.

```
## Initialise the data pipeline in the git repository
fair init
# Pull in any external data described in the yaml
fair pull .\examples\pipeline\AfricaRun.yaml
# Run the simulation described in the yaml
fair run .\examples\pipeline\AfricaRun.yaml
# Stage the code run using the unique identifier
fair add <code-run>
# Push the run and corresponding metadata back to the online registry
fair push
```
