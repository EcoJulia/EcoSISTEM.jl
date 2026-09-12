# NEWS

- v0.8.0
  - Added
    - `RasterSpec`, one lazy spec for raster data, written `SourceSpec(source, code)` for a
      catalogued layer and `RasterFileSpec(path; axis)` for a file that belongs to no dataset. The
      read options `cut`, `scale` and `fn` are its fields for either spelling; the read is windowed
      to the study area, cached, and coarsened on read by `scale`. `readfile` gains `unit`.
    - `show` methods for some over-long types.
    - `build_species` and the direct `SpeciesList` constructor take `names`, so species can carry
      real names; the names label their Diversity types too.
    - `class_fractions` and `dominant_class`, which nested as two `ConstructedRasterSpec`s regrid
      a layer of class codes by the plurality of the covering cells. `compress_landcover` is
      `dominant_class` on EarthEnv, and a cell with no data is now absent rather than class 1.
    - `examples/paper.jl`, which writes the paper's computational figures (3, 5, 6, 7, 8 and 9)
      as PDFs from the package's own examples, at the published scale when run directly and from
      a small run under the test suite.
  - Changed
    - Random streams are seeded through the generator's own seeding rather than `Base.hash`, whose
      values change between Julia versions, so a run reproduces from its seed on every Julia. A
      seed's results differ from earlier releases, once; the canonical references are re-blessed.
    - The aggregate cache is keyed by a SHA-256 digest, for the same reason, and each entry is
      written whole before it is visible. Every existing entry is re-primed on first use.
    - Tested on Julia 1.13; the continuous integration matrix runs 1.11, 1.12 and the latest release,
      and the type-order audit reads 1.13's parser as well as 1.12's.
    - Every layer reaches the study grid by aggregation of the source cells covering each grid
      cell, and nothing is interpolated: exact block aggregation where the grid is an aligned whole
      multiple of the layer's cells, which the report has always claimed, and nearest-neighbour
      sampling onto a finer lattice then aggregation otherwise, reprojection included. The
      read-time `scale` is the same computation, and a layer far finer than the grid is read
      pre-aggregated by it. Values on any grid that is not a layer's own change.
    - A coarsening reduces over the cells that carry data with a reducer chosen from the layer's
      axis - the mean, or the most frequent class for class codes. A grid cell is covered by a layer
      when the layer has data at its centre.
    - `readfile` returns a `ClimateRaster`, with `source` and `unit` keywords.
    - A categorical tolerance and its regime no longer need the same numeric type for their
      codes: an integer class list pairs with a layer of float codes.
    - EarthEnv land cover is no longer coarsened 10× by default. A read is at the file's own
      resolution unless it asks for a `scale`, and a study area chooses one from its cell size.
    - Building a layer onto a study area reads only the grid's own window, and reuses the read
      the area made when it was decided; it no longer reads the whole file a second time.
    - The `boost` parameter is gone and the birth multiplier is capped at 1, as the model is written
      up. The old constructors and keyword warn and discard it; results with `boost = 1` are
      unchanged, and runs with any other value will not reproduce.
    - `SpatialEcology`, `Proj`, `OnlineStats` and seven standard libraries are no longer
      dependencies; nothing in the package loaded them.
    - EcoBase 0.2, which answers the whole gridded interface for anything holding a grid, so
      `xmin`, `xrange`, `xedges`, `indices` and `cellanchor` now work on a `GridHabitat` and on an
      `Ecosystem` as well as on a `StudyGrid`. `StudyGrid` is an `EcoBase.AbstractRegularGrid`,
      EcoBase's new name for what it used to call `AbstractGrid`. It declares `YThenX`, so
      `indices` and `coordinates` now report `(y, x)` columns, the package's own order; ask
      `XThenY()` for the other.
  - Fixed
    - A `StudyGrid` declares that it labels cells by their lower corner. EcoBase assumed a centre,
      so every edge it derived - and every heatmap drawn from those edges - sat half a cell low.
    - An angular `cellsize` such as `30arcminute` is accepted on a geographic grid; a length there,
      and an angle on a projected grid, are refused.
    - `ShapeSpec` documents that a URL must name a self-contained file.
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
    - `LandmassesAbove`, a third coverage keeping every component that clears a threshold - the
      United Kingdom without Rockall, which is 0.031 km2 against a next-smallest of 2.536. The
      threshold is either an area (`1km^2`) or a share of the region's own total (`5percent`), the
      latter travelling better between regions of different sizes. A bare `0.05` is refused: it is
      the same quantity as `5percent` but does not say so when read beside `1km^2`. `boundingbox`
      refuses the whole coverage and says why: the shipped table records the sizes of only the
      largest few components, so only a built shape can answer it.
    - Each region now reports its `share` - what fraction of its area its largest component holds -
      which is what says whether `LargestLandmass()` suits it. New Zealand's is 56%, so asking for
      its principal landmass returns South Island alone. Derived from the shipped columns, so the
      table did not grow.
    - `investigate_regions`, which asks which named regions relate to something you have - a
      coordinate, a raster, a layer, a study area - as `investigate_study_area` reports on a grid
      before one is built. `Encloses(x)` is the default; `Overlaps(x)` gives regions your data
      reaches into, ordered by how much ground they share, and `Within(x)` those it covers entirely.
      A row of the report converts straight into a `NaturalEarthSpec`.

      By default regions are compared by **bounding box**, which costs no download and is loose:
      Norway's box encloses Edinburgh, running west to Jan Mayen and north to Svalbard. `exact = true`
      checks the survivors against the real outlines instead, which removes those false positives and
      reaches the 54 regions that cross the antimeridian and so have no comparable box. Refinement is
      lazy and stops as soon as the answer cannot change, so a query matching 361 regions by box
      typically fetches a few dozen.
    - `EcoSISTEM.naturalearth_levels()` and `naturalearth_regions(level)`, which list what exists -
      the latter with each region's bounding box, area and component count, so a level can be
      browsed rather than just enumerated.
    - A "Named regions" documentation page, whose recipes are executed by the test suite.
    - `examples/NamedRegions.jl`, which finds a study area by name, builds the British Isles from
      three countries, and simulates on the result.
    - `data/NaturalEarth/regions.csv`, a table of named regions generated from Natural Earth's
      1:10m polygons, and the levels that index it (`EcoSISTEM.NATURALEARTH_LEVELS`). `boundingbox`
      gains about 1 660 names, still answers offline and instantly, and now agrees with the shape the
      same name gives, both being derived from the same geometry.
  - Changed
    - **`ConstructedSpec` is now `ConstructedRasterSpec`**, paired with the new
      `ConstructedShapeSpec`. Neither old name said which medium it composed, and they compose the
      same way. The old name still resolves, as a deprecated binding.
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
      `NicheSpec`, `SourceSpec`, `ConstructedRasterSpec` - by `GridHabitat`, or by `build_habitat`, which
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
