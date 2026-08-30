# NEWS

- v0.5.1
  - Fixed
    - A distributed (MPI) run gave different results from a serial one on the same seed, breaking the
      reproducibility the design guarantees. This was a bugfix applied to the serial code several years ago forgotten on MPI.
    - `update!` for an abstract ecosystem was faulty but masked by equivalent functions for the
      concrete types. Now unified to a single function.
    - Fixing v0.5.0 hot loop allocation bug for serial code - revert to pre-v0.5.0 raw arrays for hot
      loop access to species counts, but keep DimArrays referencing the same memory to keep records
      of species and locations
  - Internal
    - The distributed code is renamed and rearranged to mirror the serial code file for file and name
      for name, so that the two can be read side by side.
    - More distributed and serial code redundancy elimination - update_resource_usage!() and move!().
    - Serial and distributed runs are now pinned to the same blessed results, so a divergence in
      the distributed code alone is caught at every rank count. The previous check compared
      distributed runs against each other, and all of them shared the duplicated code.
- v0.5.0
  - Removed and deprecated
    - The v0.4.0 vocabulary is renamed onto the ecological distinction between conditions (what a cell
      is like) and resources (what it provides): the `*Hab` habitat types become layers, the `*Budget`
      types supplies, the `*Requirement` types demands, and `GridAbioticEnv` becomes `GridHabitat`.
      Every name v0.4.0 exported survives as a deprecated shim.
    - The scenario callback mechanism is removed, along with `addspecies!`, `resupply!` and
      `reenergise!`; declarative interventions replace it.
    - Also removed, with no replacement: `Reference` and `create_reference`, six diversity wrappers
      that duplicated Diversity.jl, the fixed-length `month` and `quarter` units, and Pagel's-lambda
      fitting, which had been broken since Julia 0.7. Faith's PD is now Diversity.jl's `faith_pd`, and
      is the subcommunity measure rather than a mean-height-scaled richness.
    - Breaking: `ERA`, `CERA` and `CRUTS` are data sources now, rather than container types wrapping
      an array, and a reader returns a `ClimateRaster`. `ERA(array)` still works, deprecated.
    - `Phylo` and `RasterDataSources` become weak dependencies, and `AxisArrays` and `IndexedTables`
      are dropped, so installing EcoSISTEM no longer brings in a phylogenetics stack or the raster
      download machinery. Every abstract type is now `public` rather than exported.
  - A new interface
    - The grid is decided first: a `StudyArea` fixes extent, cell size and CRS before anything is
      built on it, and `investigate_study_area` reports what a run would cost before you commit to it.
    - An environment is assembled from spec recipes — `UniformSpec`, `GradientSpec`, `PeakedSpec`,
      `NicheSpec`, `SourceSpec`, `ConstructedSpec` — by `GridHabitat`, or by `build_habitat`, which
      supplies what you do not name and reports what it chose. Species come from `build_species`, and
      `build_ecosystem` pairs the two sides and checks them against each other.
    - One layer family covers both halves of the environment, parameterised by role (`Condition` or
      `Resource`) and by niche axis, and tolerances and demands mirror it exactly on the species side.
      Everything is named by its niche axis, and a layer pairs only with a requirement on the
      identical axis. Collections of any of them implement the standard container interface — `keys`,
      `values`, `pairs`, `iterate`, `getindex` — in place of the bespoke accessors.
    - Environmental change is declared rather than programmed: any layer may carry a layer change,
      which is a pure function of elapsed time. Change to the ecosystem itself — deactivating cells,
      adding or removing abundance, introducing a species — is a separate mechanism, `Intervention`,
      scheduled by time and aimed at a region.
    - A run may carry an epoch, the real date its elapsed time zero stands for, and calendar month
      durations replace the old fixed 30.4375-day month.
  - Improvements
    - Climate reading moved onto Rasters.jl: coordinates come from the file's own metadata, a read is
      windowed to the study area rather than pulled in whole and cropped, downloads are cached under
      `EcoSISTEM.assetdir()`, and a simulation grid can be genuinely projected. A read that asks for a
      `scale` without a window is the exception: it coarsens the whole file, and memoises the result
      so that only the first one is slow.
    - Cell size and area can be asked of anything that knows the grid, and a geographic grid answers
      with an angle and a true solid angle rather than a fabricated length — so a supply on a
      latitude/longitude grid is scaled by its own cell's area instead of the whole grid's.
    - Several fixes changed what the model computes: the suitability term was applied once per resource
      in a multi-resource environment; an axis's canonical unit silently moved every equilibrium; a
      generated niche read uninitialised memory, so it was not reproducible from its seed; the
      ordinariness cache went stale; a `SpeciesList` did not forward Diversity's type interface;
      `gatherabundance` was wrong whenever species or cells did not divide evenly across MPI ranks; and
      a monthly climatology was addressed from zero, so asking for March returned February.
    - Grid dimensions are (Y, X) everywhere, now enforced by the type, and sixteen types gained a
      `show`, so a value that holds a grid no longer prints the grid.
    - New documentation pages on layers, time, units, interventions, how the model works and running at
      scale. The code in the documentation and every file in `examples/` is run as a test, and the
      suite is now eight separately nameable sets.
    - The shipped layer catalogues record each layer's actual raster unit, the period it accumulated
      over and the sources that supply it, and several unit and scale errors in the shipped data are
      corrected on read.
- v0.4.0
  - Speed up the multithreaded update loop with cache-line-sized species blocks and greedy scheduling
  - Require Julia v1.11 for greedy scheduling
  - Simplify interface and clean up code
  - readfile calls now uses keyword bounds (xmin/xmax/ymin/ymax)
  - Modernise in-repo climate wrapping from the deprecated Worldclim_bioclim(...) to ClimateRaster(WorldClim{BioClim}, ...)
- v0.3.0
  - Tidying up and adding missing documentation
  - Fix race condition in multithreaded code post Julia 1.9
  - Update compats
  - Add tests for CachedEcosystem reproducibility
  - Add tests for MPIEcosystem reproducibility
  - Make MPIEcosystem reproducible
  - Add tests that CachedEcosystem and MPIEcosystem match Ecosystem
  - Refactor to condense some repeated code
  - Fix bug in NoGrowth energy use for multiple energy budgets
  - Fix bug in MPI use of multiple energy budgets
- v0.2.6
  - Compat fixes and resolve Pluto notebook error
- v0.2.5
  - Use ResearchSoftwareMetadata
  - Add in metadata and code hygene testing
- v0.2.4
  - Add metadata and crosswalk
- v0.2.3
  - Security fix on unzipping: #140
  - Minor bugfix on up- and down-scaling images: #139
- v0.2.2
  - Remove all manifests and associated code
  - Fix MPI example on HPC
- v0.2.1
  - Move MPI structs into extension
- v0.2.0
  - Require Julia v1.9 for extensions
  - Create package extensions
  - Update testing
  - Fix plotting and unzipping code
- v0.1.4
  - Add in Pluto example
  - Some Windows fixes
- v0.1.3
  - Code restructure
  - Move to EcoJulia
  - Add in MPI testing
- v0.1.2
  - Fix phylogenetic diversity management
- v0.1.1
  - Remove unnecessary Compat dependency
  - Fix incompatibility with latest Diversity release
  - Fix license recognition by GitHub
  - Fix MPI process allocation for high number of processes
  - CompatHelper updates
- v0.1.0
  - Initial release to Julia registry
- v0.0.1
  - First tagged release stored on Zenodo
