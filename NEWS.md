# NEWS

- v0.5.0 (unreleased)
  - Added / changed
    - New `build_environment`/`build_species`/`build_ecosystem` convenience API for assembling ecosystems
    - New `(Role, NicheAxis)` layer abstraction (`AbstractLayer{R, A}`, `ContinuousLayer`/`DiscreteLayer`);
      a layer's meaning now comes from its niche axis, not its unit
    - **Core vocabulary renamed onto the ecological *resources vs conditions* distinction**
      (Begon–Townsend–Harper). The environment (now unambiguously `Habitat`) and each species meet on two
      sides. On the **condition** side an environmental `Regime` (was the `Habitat`-role layer) is matched to
      a species `Tolerance` (`NicheTolerance`, was `Bin`/`*Trait`), their fit a `Suitability`
      (`NicheSuitability`, was `DistRel`), under the umbrella `AbstractNicheFit` (was
      `AbstractTraitRelationship`). On the **resource** side an environmental `Supply` (was `Budget`) is
      consumed against a species `Demand` (was `Requirement`). The layer `Role` markers are now
      `{Condition, Resource}` (were `{Habitat, Budget}`), and "energy" quantities are now "resources" (they
      were never only energy — solar is `kJ/day`, water is `L/day`). Nearly every old public name keeps
      working via a warning shim — see "Deprecated" for the full map
    - **`Resource`-role Supply/Demand values are now genuinely rate-dimensioned**: `SolarSupply`/`SolarDemand`
      store `kJ/day` (were bare `kJ`) and `WaterSupply`/`WaterDemand` store `L/day`, a genuine volume, not a
      depth (were bare `mm`) — a resource is consumed and replenished *per timestep*, so its unit should say
      so. This also fixes a pre-existing, previously-untested dimensional bug where `WaterSupply` never
      actually multiplied a depth by cell area, so two grids of different sizes with the "same" rainfall
      silently held the same total water. `cancel()` (`src/AbioticEnv.jl`), the shared "arealRate × area →
      absolute rate" primitive, now dispatches purely on dimension (any native time unit — year/month/day —
      converts transparently) rather than on a fixed unit, and covers the Solar/Water/Simple cases uniformly.
      Per the water-density identity (1 kg water ≈ 1 L ≈ 1 mm over 1 m²), an areal water rate only needs
      `/day` (e.g. `mm/day`); an areal solar rate needs an explicit `/area` too (`kJ/m²/day`), since energy has
      no equivalent implicit density
    - `canonicalunit` gained a `(::Type{<:Role}, ::NicheAxis)` overload: a `Condition`-role reading of an axis
      (a niche tolerance / descriptive climatological normal) and a `Resource`-role reading of the same axis
      (a literal consumption rate) are legitimately different physical questions — `Precipitation` is `mm` as
      a `Condition` but `L/day` as a `Resource`, `SolarRadiation` is `kJ/day` as a `Resource` (no `Condition`
      meaning is defined). The existing 1-arg form is unchanged and remains the default for any role/axis
      without a specific override. Scope: `Precipitation` and `SolarRadiation` only for now — other `:rate`
      axes (`Evapotranspiration`, `ClimateMoisture`, `CarbonFlux`, etc.) don't yet have a dedicated `Resource`
      type and are deferred
    - The old Dict-based `SUPPLYDICT`/`checksupply` dispatch is replaced by genuine multiple dispatch:
      `supplytype`/`demandtype` now have direct `Unitful.Power`/`VolumeFlow`/`Real` quantity-type overloads (a
      new `VolumeFlow` `@derived_dimension`, `𝐋³·𝐓⁻¹`, since Unitful has no built-in name for it), with a
      generic fallback method raising a clear error for anything else — no boolean check-and-error chain. A
      new guard test walks every concrete `Abstract1Demand` subtype (via `subtypes`) and asserts its stored
      element dimension is genuinely substance × 𝐓⁻¹ (`_basedimension`)
    - New distribution-based continuous tolerance `NicheTolerance` (any `Distributions` response) and nichefit
      `NicheSuitability`, replacing the temperature/rainfall-specific `TempBin`/`RainBin` split (kept as the
      `TempTolerance`/`RainTolerance` aliases)
    - New symmetric environment accessors `condition(::AbstractRegime)` and `resource(::AbstractSupply)`
    - Environment constructors renamed to drop the legacy `AE` ("abiotic env") suffix and end in `habitat`
      (e.g. `simplehabitatAE`→`simplehabitat`, `worldclimAE`→`worldclimhabitat`, `lcAE`→`landcoverhabitat`);
      the four data constructors (`erahabitat`/`worldclimhabitat`/`bioclimhabitat`/`landcoverhabitat`) now take
      an `axis` keyword to tag the built regime with its `NicheAxis` (default `Unclassified`, i.e. unchanged)
    - The `lc` land-cover abbreviation is expanded throughout — `LCtolerance`→`LandCoverTolerance`,
      `LCsuitability`→`LandCoverSuitability`, `compressLC`→`compressLandCover` (CamelCase forms use `LandCover`,
      matching RasterDataSources' `EarthEnv{LandCover}`)
    - Climate readers migrated onto Rasters.jl (real coordinates from file metadata; faster, disk-cached
      aggregation)
    - The `NicheAxis` catalogue was reorganised into `XxxAxis` grouping supertypes (`TemperatureAxis`,
      `WaterAxis`→`PrecipitationAxis`/`HumidityAxis`/…, `SolarRadiationAxis`, `WindSpeedAxis`, `CloudCoverAxis`,
      `DayAxis`, `CarbonAxis`, `TypologyAxis`). Following the max/min/mean-collapse-with-separate-`…Range`
      policy, the many BioClim temperature/precipitation axes collapse into `Temperature`/`TemperatureRange`
      and `Precipitation`, and every fetchable RasterDataSources layer is now classified onto an axis
    - `layerinfo`/`LayerRecord` now report each layer's full list of supporting sources and its temporal
      resolution / number of layers (e.g. "Monthly (12 layers) — select with `month=`")
    - The shipped layer tables (`data/RasterDataSources/*.csv`) now record each layer's **actual** raster unit
      rather than the canonical unit — notably absolute temperatures are `°C` (as WorldClim/CHELSA store them),
      converted to `K` only when a layer is built. The readers no longer apply a hidden `°C→K` shift, so
      `read(WorldClim{Climate}, :tmin)` and `readCRUTS`/`readCHELSA_monthly` now return `°C`, and the spurious
      `ncdf` "layer" was dropped
  - Breaking changes
    - `SolarSupply`/`SolarDemand`/`SolarTimeSupply` now require `kJ/day` (were bare `kJ`) and
      `WaterSupply`/`WaterDemand`/`WaterTimeSupply` now require `L/day` (were bare `mm`); `SimpleSupply`/
      `SimpleDemand`/`SizeDemand` are now stored as `typeof(1.0/day)` (were bare `Float64`). This is a
      dimensional type change, so there is no deprecation shim — every construction site (`build_species`'s
      `resource =`/`build_environment`'s `supply =` defaults included) must be updated to the new units
    - Removed the volumetric-water resource types (`VolumetricWater` axis, `VolWaterSupply`/`VolWaterTimeSupply`/
      `VolWaterDemand`, and their v0.4.0 `VolWaterBudget`/`VolWaterTimeBudget`/`VolWaterRequirement` names),
      added in v0.4.0. There is no deprecation shim; an `m³` supply/resource now errors
    - `EcoSISTEM.AbstractHabitat` (unexported) **changed meaning**: it is now the environment container (was
      `AbstractAbiotic`), because the old condition-layer `AbstractHabitat` was renamed `AbstractRegime`.
      Code using the qualified name for the old condition layer must switch to `AbstractRegime`
    - `TempTolerance`/`RainTolerance` (and their `TempBin`/`RainBin` deprecated aliases) matrix constructors
      now require valid distribution parameters (4 columns for the `Trapezoid` response, 2 for `Uniform`);
      previously-invalid matrices that "worked" under lazy construction now error at build time
    - A regime's per-timestep change function now comes from its niche `axis` (`dynamics(::NicheAxis)`), not
      from the value's unit. A `simplehabitat`/`simpleregime` temperature regime that should change under a
      scenario must now pass `axis = Temperature` (previously any `K`-unit habitat was given `TempChange`
      automatically); a static regime is unaffected
  - Deprecated (still work, emit a warning) — the vocabulary rename, by side
    - Condition side, environment layer: `AbstractHabitat`(-role)→`AbstractRegime`, `ContinuousHab`→
      `ContinuousRegime`, `ContinuousTimeHab`→`ContinuousTimeRegime`, `DiscreteHab`→`DiscreteRegime`,
      `HabitatCollection2/3`→`RegimeCollection2/3`, `gethabitat`→`getregime`
    - Condition side, species response: `AbstractTraits`→`AbstractTolerance`, `DiscreteTrait`→
      `DiscreteTolerance`, `LCtrait`→`LandCoverTolerance`, `TraitCollection2/3`→`ToleranceCollection2/3`,
      `TempBin`/`RainBin`→`TempTolerance`/`RainTolerance`, `traitpopulate!`/`traitrepopulate!`→
      `tolerancepopulate!`/`tolerancerepopulate!`
    - Condition side, the matcher: `AbstractTraitRelationship`→`AbstractNicheFit`, `Match`→`MatchSuitability`,
      `LCmatch`→`LandCoverSuitability`, `NoRelContinuous`/`NoRelDiscrete`→`NoFitContinuous`/`NoFitDiscrete`,
      `multiplicativeTR2/3`→`multiplicativeFit2/3`, `additiveTR2/3`→`additiveFit2/3`, `gettraitrel`→
      `getnichefit`
    - Resource side, environment: `AbstractBudget`→`AbstractSupply`, the `*Budget` supply types →`*Supply`
      (`SimpleBudget`, `SolarBudget`/`SolarTimeBudget`, `WaterBudget`/`WaterTimeBudget`),
      `BudgetCollection2`→`SupplyCollection2`, `getbudget`→`getsupply`, `reenergise!`→`resupply!`
    - Resource side, species: `AbstractRequirement`→`AbstractDemand`, the `*Requirement` demand types →
      `*Demand` (`SimpleRequirement`, `SizeRequirement`, `SolarRequirement`, `WaterRequirement`),
      `ReqCollection2`→`DemandCollection2`
    - Environment container: `AbstractAbiotic`→`AbstractHabitat`, `GridAbioticEnv`→`GridHabitat`
    - Layer dynamics: `HabitatUpdate`→`LayerUpdate`
  - Deprecated (still work, emit a warning) — other
    - `GaussTrait(A, mean, sd)` → `NicheTolerance(A, Normal, mean, sd)`; a *unitful* axis-less
      `GaussTrait(mean, sd)` still works but now warns that inferring the niche axis from the unit is
      deprecated — name the axis
    - `Gauss`/`Trapeze`/`Unif` relationships → `NicheSuitability` (the `Gauss` 3-argument `(current, opt, sd)`
      call is retained)
    - `Worldclim_monthly`, `CHELSA_monthly` → `ClimateRaster(WorldClim{Climate}/CHELSA{Climate}, …)` (new
      this release); `Worldclim_bioclim`, `CHELSA_bioclim`, `Landcover` remain deprecated
    - `readworldclim` → `read(WorldClim{Climate}, layers; …)`; `readfile(file, xmin, xmax, ymin, ymax)` →
      `readfile(file; cut = …)`
    - Environment constructors (drop the `AE` suffix): `simplehabitatAE`→`simplehabitat`,
      `tempgradAE`→`tempgradhabitat`, `raingradAE`→`raingradhabitat`, `peakedgradAE`→`peakedgradhabitat`,
      `simplenicheAE`→`simplenichehabitat`, `eraAE`→`erahabitat`, `worldclimAE`→`worldclimhabitat`,
      `bioclimAE`→`bioclimhabitat`, `lcAE`→`landcoverhabitat`
    - `compressLC` → `compressLandCover`
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
