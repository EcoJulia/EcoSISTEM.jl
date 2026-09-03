# NEWS

- v0.7.0
  - Added
    - `AllTerritories` and `LargestLandmass`, which say how much of a named region to take. A name
      almost never denotes one connected piece of ground - "France" includes Guadeloupe, "Norway"
      includes Bouvet Island in the South Atlantic - so a selection either takes everything the name
      covers, or the largest connected pieces of ground it covers.
    - `NaturalEarthSpec`, which names a study area's active cells as a **region** - a country, a
      continent, an island - rather than as a file. `StudyArea(within = NaturalEarthSpec("Scotland"))`
      downloads and cuts the outline when it is built; `outline = false` takes the region's bounding
      box instead. The name is checked against the shipped table when the spec is written, by the
      same rule `boundingbox` uses, so the box that function reports is the box the shape has.
    - `ConstructedShapeSpec`, the vector mirror of `ConstructedRasterSpec`: it composes **geometry**
      where that one composes rasters, so the result is exact and carries no resolution of its own.
      Members are any shape specs - a `ShapeSpec` of your own study area, a `NaturalEarthSpec` named
      by country, or another `ConstructedShapeSpec`. The union of the United Kingdom, Ireland and the
      Isle of Man reaches Shetland at 60.85, which Natural Earth's own "BRITISH ISLES" polygon cuts
      off at 59.80.

      Operations are `ShapeUnion`, `ShapeIntersection` and `ShapeDifference` for combining, and
      `ShapeBuffer`, `ShapeSimplify` and `ShapeConvexHull` for transforming one shape -
      `ShapeBuffer(50km)` being how "within 50 km of this coastline" is said. As on the raster side
      an **arbitrary function** is accepted too, so anything ArchGDAL offers is reachable without a
      new operation type.
    - `LandmassesAbove`, a third coverage keeping every component at least a given area - the United
      Kingdom without Rockall, which is 0.031 km2 against a next-smallest of 2.536. `boundingbox`
      refuses it and says why: the shipped table records the sizes of only the largest few
      components, so only a built shape can answer it.
  - Changed
    - **`boundingbox` is breaking.** Its `islands::Bool` keyword becomes `coverage`, taking
      `LargestLandmass()` (the default, as `islands = false` was) or `AllTerritories()`; it gains a
      `level` keyword saying what kind of region a name means; and its values move, because they now
      come from Natural Earth rather than from a hand-made file. Mainland extents are within about
      0.02 degrees of the old ones, but island-inclusive extents can move much further - Scotland's
      western edge goes from -8.65 to -13.69, which is Rockall.

      Five names are gone, having had no Natural Earth equivalent: `UK` is `"United Kingdom"`,
      `NI` is `"Northern Ireland"`, `SouthAmerica` is `"South America"`, `GB` is
      `boundingbox("United Kingdom", coverage = LargestLandmass())` or the physical island
      `boundingbox("GREAT BRITAIN", level = "Physical Island")`, and `BritishIsles` has to be
      assembled, Natural Earth's own polygon of that name cutting off Shetland.

      **The default coverage is `AllTerritories()`**, which is what Natural Earth means by a name:
      its "France" is the one that includes Guadeloupe. Taking only the principal landmass is a real
      choice about which ground is wanted, so it is now written out rather than assumed - the old
      `islands = false` behaviour is `coverage = LargestLandmass()`. One consequence worth knowing:
      54 of the 2 444 rows cross the antimeridian, so names like "Russia" and "North America" now
      have no bounding box under the default and say so, naming `LargestLandmass()` as the remedy.

      A name meaning genuinely different ground at different levels is now refused rather than
      guessed at, and the error tabulates what each level would give so that the choice can be made
      from the message: as a continent "Africa" is 55 countries and as a UN region 62, whose full
      extents differ by 54 degrees of longitude. The comparison is against the coverage asked for,
      not the whole selection, so a name is only refused when the answer you actually requested is
      ambiguous. A name whose levels agree, such as "Scotland", still needs no level. A region
      crossing the antimeridian is refused too, since an `Extents.Extent` holds an interval and
      cannot express one.
  - Fixed
    - Hot loop allocation fix. `GridHabitat` now carries the topology's type as a parameter.
    - Speed-up for `ShapeSpec` mask building.
    - Remote reads through GDAL - a URL given to `ShapeSpec`, or any `/vsicurl/` path - work again
      on macOS under Julia 1.12, where they had failed with an error naming a null pointer rather
      than a certificate. That one release builds libcurl against OpenSSL with no certificate roots,
      where 1.11 uses the system keychain and 1.13 ships roots of its own, so the fix is confined to
      it and defers to `CURL_CA_BUNDLE`, `SSL_CERT_FILE` and `SSL_CERT_DIR` if you set them.
  - Internal
    - The hot loop's inference check is now a sweep over every field of the types it reaches, and
      the distributed loop gains the same allocation and inference checks, which it had never had.
- v0.6.1
  - Added
    - `AlwaysMovement` now works in a distributed (MPI) run, so all three movement types do. Each
      rank owns a block of species across the whole grid while the dynamics run, so dispersing an
      established individual stays rank-local exactly as dispersing a newborn does.
  - Changed
    - `getabundance` on a distributed ecosystem refuses rather than returning one rank's block, and
      names `gatherabundance` instead. Diversity's consumers each reduce that matrix over a different
      axis, so a block silently answered at least one of them wrongly. The measures built on it
      refuse too for now, rather than returning numbers that depended on the rank count.
  - Fixed
    - Another inconsistency between serial and distributed code - changing abundances through an
      intervention diverged from the serial run on the same seed because of synchronisation
      order, which is now switched to give a consistent result.
    - The diversity measures now run distributed, each rank computing the cells it owns against the
      full similarity matrix, and `gatherdiversity` assembles them into the serial answer on every
      rank. `_getordinariness!`, `_getmetaabundance`, `_getweight` and `_getscale` now work too.
- v0.6.0
  - Breaking
    - A landscape's `matrix` and `grid` are plain arrays again, as they were before v0.5.0. The
      labelled views are now `dimmatrix` and `dimgrid`, which share the same memory, so nothing is
      copied and both are always available.
    - A `GridLandscape` is built from the `StudyGrid` it sits on rather than from a `(Y, X)` tuple,
      which is what lets each cell name its own extent. The form taking only a size is gone, since
      it had no grid to take positions from.
    - `emptygridlandscape` and `empty_mpi_gridlandscape` are replaced by a single `empty_landscape`.
      Given a habitat and species list it builds a serial landscape; given the partition as well it
      builds a distributed one, so the signature says which rather than the name. The old MPI name
      errors, naming the replacement: it took the partition alone and cannot reach the species names
      or grid the labelled views need.
    - `copy` of a landscape is removed - not needed and hard to reimplement for new landscape fields.
  - Fixed
    - A distributed (MPI) run gave different results from a serial one on the same seed, breaking the
      reproducibility the design guarantees. This was a bugfix applied to the serial code but
      missed on MPI.
    - `update!` for an abstract ecosystem was faulty but masked by equivalent functions for the
      concrete types. Now unified to a single function.
    - Fixing v0.5.0 hot loop allocation bug for serial code - revert to pre-v0.5.0 raw arrays for hot
      loop access to species counts, but keep DimArrays referencing the same memory to keep records
      of species and locations.
  - Internal
    - The distributed code is renamed and rearranged to mirror the serial code file for file and name
      for name, so that the two can be read side by side.
    - More distributed and serial code redundancy elimination - update_resource_usage!() and move!().
    - Serial and distributed runs are now pinned to the same blessed results, so a divergence in
      the distributed code alone is caught at every rank count. The previous check compared
      distributed runs against each other, and all of them shared the duplicated code.
    - A landscape's flat `location` dimension now names each cell by its extent, computed on demand
      from the grid. The distributed landscape gains the same kind of labelled views: `dimgrid` for
      its own species against the real `Y` and `X`, and `dimcols` for every species against its own
      cells, presented using BlockArrays as a single ordinary matrix without copying.
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
    - An environment is assembled from spec recipes - `UniformSpec`, `GradientSpec`, `PeakedSpec`,
      `NicheSpec`, `SourceSpec`, `ConstructedSpec` - by `GridHabitat`, or by `build_habitat`, which
      supplies what you do not name and reports what it chose. Species come from `build_species`, and
      `build_ecosystem` pairs the two sides and checks them against each other.
    - One layer family covers both halves of the environment, parameterised by role (`Condition` or
      `Resource`) and by niche axis, and tolerances and demands mirror it exactly on the species side.
      Everything is named by its niche axis, and a layer pairs only with a requirement on the
      identical axis. Collections of any of them implement the standard container interface - `keys`,
      `values`, `pairs`, `iterate`, `getindex` - in place of the bespoke accessors.
    - Environmental change is declared rather than programmed: any layer may carry a layer change,
      which is a pure function of elapsed time. Change to the ecosystem itself - deactivating cells,
      adding or removing abundance, introducing a species - is a separate mechanism, `Intervention`,
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
      with an angle and a true solid angle rather than a fabricated length - so a supply on a
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
