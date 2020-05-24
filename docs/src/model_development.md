# Model development

## Top priorities
* **Incorporate higher dimensionality in disease classes.** I.e. Single susceptible category needs to be easily broken up into higher dimensions for infection category, age, gender, job type etc. This may involve a mapping between different array dimensions.
* **Two types of infection process from direct and indirect contact.** Need to have an instantaneous infection pressure (from direct contact) that simulates aerosol transmission, which would be eliminated at the end of each step. Also have an environmental transmission which decays over time and has a lower transmission rate.  

## Further developments
* **Add in household structure.** This will include a more "individual based" model in which different ages and genders are tracked within a cell. However, timing of events is unlikely to be tracked as this is difficult to parallelise.
* **Add in different location types.** E.g. schools, work, carehomes, general social situations (proxy for pubs, restaurants etc).
* **Regionalisation.** Working towards an 'irregular' grid structure. However, in the first instance, we could track information like health board or local authorities for each grid cell, and explore lockdown strategies among these different areas on a regular grid.
* **Virus seasonality.** There will be at least some seasonality with climate and other environmental information could also be considered, such as pollution levels.
