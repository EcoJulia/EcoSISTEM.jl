# Model development

## Top priorities
* **Virus seasonality.** There will be at least some seasonality with climate and other environmental information could also be considered, such as pollution levels. Work underway to incorporate these types of data into the model.
* **Clearer designation of the different compartment types.** E.g. Disease state, age, gender etc. Work has begun on this, creating a higher dimensional storage array with labelled axes (through AxisArrays).
* **MPI.** In order to make best use of HPC resources, we should re-implement MPI in this section of the code, based upon the original in the biodiversity code.

## Further developments
* **Add in household structure.** This will include a more "individual based" model in which different ages and genders are tracked within a cell. However, timing of events is unlikely to be tracked as this is difficult to parallelise. This work was underway, but has been shelved in favour of other priorities.
* **Add in different location types.** E.g. schools, work, care homes, general social situations (proxy for pubs, restaurants etc).
* **Regionalisation.** Working towards an 'irregular' grid structure. However, in the first instance, we could track information like health board or local authorities for each grid cell, and explore lockdown strategies among these different areas on a regular grid.
