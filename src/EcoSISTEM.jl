# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEM

using Scratch: get_scratch!
using Downloads: Downloads

# **The model, before anything else.** Fourteen declarations saying what EcoSISTEM simulates —
# conditions and resources, the environment's layers and the species' requirements, the fit between
# them, and the habitat and ecosystem that assemble them. It depends on nothing of ours, so it comes
# first; everything below is a choice about *how*.
include("Ecology.jl")

public NicheAxis, AbstractLayer, Role, Condition, Resource, AbstractRegime,
       AbstractSupply,
       AbstractSpecies,
       AbstractSpeciesRequirement, AbstractTolerance, AbstractDemand,
       AbstractNicheFit,
       AbstractEcosystem, AbstractHabitat

# The asset cache's descriptor type, before the two functions below name it.
include("Asset.jl")

public CachedAsset

public assetdir

public assetpath

# The units submodule: the arcminute/arcsecond subdivisions of a degree, and the calendar-month
# durations. First, because other files depend on it.
include("Units/Units.jl")

# The coordinate vocabulary — the two-dimensional position/size family and the geographic point type
# — used across the `ClimatePref` submodule and the main module, and by `CircleMaskSpec`
# (`Spec.jl`), so it is defined here, before all of them.
include("Coordinates.jl")

export LatLong

# Members of the `get*` accessor family, which is flatcase by design and `public` throughout. Not
# `lat`/`long`: those are short enough to collide with a local variable and with the English word
# "long", both of which occur in this repo's own examples.
public getlat, getlong

# `public`, not exported. [`LatLong`](@ref) is the type a user names; these are the unit-general
# forms behind it, reached when writing code that works on projected grids as well as geographic.
public SpatialLocation, SpatialSize

# `public`, not exported: a named-region lookup is reached for deliberately, and `boundingbox` is a
# common enough English word that it should not land in every `using EcoSISTEM` namespace.
public boundingbox

public SpatialKind, Spatial2D

# The axis machinery, then the axes declared with it. `@nicheaxis` emits methods whose signatures
# name `Resource`, so `Ecology.jl` — which declares the roles — must already have been included.
include("nicheaxis_macro.jl")

include("NicheAxis.jl")

# The axis grouping supertypes, declared by `@nicheaxis` above. `public`, not exported: an axis
# group is what one dispatches on or subtypes; the leaves below are what a user names.
public TemperatureAxis, WaterAxis, PrecipitationAxis, HumidityAxis,
       VapourPressureDeficitAxis, RelativeHumidityAxis, EvapotranspirationAxis,
       ClimateMoistureAxis, SolarRadiationAxis, WindSpeedAxis, CloudCoverAxis,
       DayAxis,
       CarbonAxis, TypologyAxis, SpaceAxis

# The declaration macro is the supported way to add an axis, and the only one — the
# `canonicalunit`/`supplytype`/`demandtype`/`bounds` methods behind it stay internal. Exported rather
# than `public` because it is written, not merely dispatched on.
export @nicheaxis

# The response-distribution parameter roles, and `Trapezoid`.
include("Dist.jl")

export Trapezoid

public read_distribution, paramunits

public ParamRole

# Where a layer's data comes from, and which grid a combine runs on.
include("DataSource.jl")

# `public`, not exported: the abstract types (the house rule), the code-spelling union, and
# `in_memory_raster` — a supported pathway for data you already hold, but not one to encourage.
public EcoSISTEMSource, CODE_TYPE
include("CombineStage.jl")

export CombineOnTargetGrid, CombineOnSourceGrid

public AbstractCombineStage

# The layer recipes.
include("Spec.jl")

export UniformSpec, GradientSpec, PeakedSpec, NicheSpec, CircleMaskSpec,
       SurfaceSpec

public AbstractSpec, AbstractSyntheticSpec, AbstractSyntheticLayerSpec,
       AbstractSyntheticMaskSpec

# A raster of climate data, then the specs that read one.
# The shipped layer catalogue — dataset-agnostic parsing that `SourceSpec` consults for a layer's
# unit and axis, so it must precede it.
include("AccumulationPeriod.jl")

# The catalogue's own vocabulary, moved here with it from `ClimatePref`. `public`, not exported:
# a record is read off the shipped table, never written by a caller.
public AbstractAccumulationPeriod, ConstantAccumulationPeriod,
       PerSliceAccumulationPeriod,
       PerCellAccumulationPeriod

include("LayerCatalogue.jl")

public LayerRecord, AxisNode

public layerinfo, layersbyaxis, layerrate

# Exported rather than `public`, unlike their three siblings above -- that split is inherited from
# `ClimatePref`, which exported exactly these two. Kept as it was rather than changed in passing;
# whether the catalogue family should be uniform is a separate question.
export layeraxis, layerunit

# A raster of climate data, and the specs that lazily read one — including the three unions naming
# everything a caller may hand to `regime`/`supply`/`within`.
include("Climate.jl")

# The concrete data-source types, and the readers and sampler that go with them. All were exported by
# `EcoSISTEM.ClimatePref` before it was dissolved.
export ERA, CERA, CRUTS

public AbstractClimate

public in_memory_raster

export ClimateRaster

include("LazySpec.jl")

# **What a layer is, and how one is declared** — the raster type, what may name its source, and the
# specs. **Before `ClimatePref`, deliberately**: the readers in that submodule construct
# `ClimateRaster`s, so these must already exist. See the file's own header for the rule that decides
# what lives here against what stays with the climate data.
# Exported here, and **re-exported by `ClimatePref`**, so that `using EcoSISTEM.ClimatePref` reaches
# them too.
export SourceSpec, ConstructedSpec, ShapeSpec

public AbstractLazySpec

# How a layer changes in time: what a series does at its end, how a change value is interpreted,
# and the recipes a caller writes. `SeriesPolicy` first — `SeriesChange` holds an `AbstractSeriesEnd`.
include("SeriesPolicy.jl")

# what a series does once elapsed time runs past its last stored slice (the `atend` keyword)
export ErrorAtEnd, HoldAtEnd, RepeatAtEnd, RevertToLayer

# what a series' time coordinates mean, and so what a run's epoch can do with them (the `calendar`
# keyword); see `build_ecosystem`'s `epoch`
export DatedSeries, MonthOfYearSeries, UndatedSeries

public AbstractSeriesEnd, AbstractSeriesCalendar

include("ChangeMode.jl")

# the per-layer change rules a layer can hold (defined in Layer.jl; applied in LayerChange.jl)
# The change *modes* are derived from a spec, never named — `_changemode(::ReplaceWith)`
# gives `AbsoluteChange()`, and no public signature takes a mode. They are what a spec becomes.
public NoChange, AbsoluteChange, RelativeChange, RateChange

public AbstractChangeMode

include("ChangeSpec.jl")

# The layer-change unit contract and the apply/drive methods. Included late — it dispatches on
# both `AbstractLayer` (Layer.jl) and `AbstractEcosystem` (Ecosystem.jl), and calls `axisof`
# (Simplify.jl) to name the offending axis when a change's unit does not match its layer.
export ReplaceWith, OffsetBy, IncrementBy, PatternedChange, SeriesChange

# Adding changes together is spelled `A + B`, so the type that spelling builds is supported but
# rarely written by hand — `public`, not exported, to keep the name out of the way.
public CombinedChange

public AbstractChangeSpec

include("collections.jl")

# The materialised changes, and `Varying` — whose constructor accepts either a spec or one of these.
include("LayerChange.jl")

# `NoLayerChange` stays exported because **no spec produces it** — a static layer is
# written `NoLayerChange()` directly, alongside `IncrementBy(rate)`. The other four are what
# `_attachchange` materialises a spec into, so the spec is the name a user writes.
export NoLayerChange

public SteadyLayerChange, PatternedLayerChange, SeriesLayerChange,
       SumOfLayerChanges

# The default `shape` of a `PatternedChange`. `public` rather than exported: a caller writes a shape
# far more often than they name the default, and the name is common enough to be worth keeping out
# of the way of one — but it must be nameable, or the signature advertises something unreachable.
public sinusoidal

# Declaring a change alongside a layer spec, at the `GridHabitat` boundary
export Varying

# Installing a change on a layer directly is supported but not exported — the user-facing route is
# to declare the change alongside the layer when the environment is built.
public changeunit

public AbstractLayerChange

include("Layer.jl")

export ContinuousLayer, CategoricalLayer, LayerCollection

# Regime and supply aliases over the layer types, defined in `Layer.jl`. There are no `*Time*`
# counterparts, because a layer that varies in time is the **same** type as one that does not,
# differing only in the change it carries.
export ContinuousRegime, CategoricalRegime

# **One supply name, parameterised by axis** — `Supply{SolarRadiation}`, never a per-resource alias.
# The axis is the single declaration of what a supply measures, and its value type follows from
# `canonicalunit(Resource, axis)`, so an alias per resource would restate in the type what the axis
# already says.
export Supply

public setchange!

# Checking a run against the series driving it is supported API — `simulate!` calls it before the
# first step, but it is worth being able to ask before committing to a run.
public checkcoverage, check_bounds

# **The species side mirrors the layer side member for member**, so it follows straight on: a
# `Condition` requirement is a tolerance, a `Resource` one a demand, and the fits below say how a
# tolerance is scored against the regime it is paired with.
include("SpeciesRequirement.jl")

export SpeciesRequirementCollection

include("Tolerance.jl")

# Reaching inside a tolerance is how you *write* one, not how you run a model.
public getpref, getdist

export SimpleCategoricalTolerance, NicheTolerance, TempTolerance, RainTolerance

public ContinuousTolerance, AbstractCategoricalTolerance
include("Demand.jl")

export Demand

include("NicheFit.jl")

# Supported but unexported: `nichefitcombine` is how a `CombiningFit` folds its members, which callers
# need when writing their own, but it is not part of the everyday vocabulary.
public nichefitcombine

export suitability

export CombiningFit, MultiplicativeFit, CategoricalSuitability, NoFitContinuous,
       NoFitCategorical, NicheSuitability, AdditiveFit

# The species' own behaviour and rates, and the grid's edge rule — none of which need the
# environment, so they sit here rather than beside the habitat that carries the topology.
include("Movement.jl")

export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement,
       NoMovement

# `getkernels` reaches inside a movement type, which is what writing a new one needs.
public AbstractMovement, AbstractKernel, getkernels
include("Demographics.jl")

export PopGrowth, EqualPop, NoGrowth

public AbstractParams

include("Topology.jl")

# The three aliases are the API. They cover three of the four combinations —
# `EdgeTopology{Periodic, Bounded}` has no name — and the maintainer's rule is that a fourth would get
# an **exported alias of its own** rather than callers reaching for the underlying type. So the alias
# layer stays the surface, and what it is built from does not.
export Torus, Cylinder, Island

public EdgeTopology, Periodic, Bounded

public AbstractTopology, AbstractBoundaryCondition

# **The study area, and the report that says how its grid was decided.** Both must precede
# `GridHabitat`, which holds an `area::StudyArea{L}` — a struct field type is resolved when the
# struct is defined, so the area cannot be declared after the habitat that carries it.
include("DecisionSource.jl")

public AbstractDecisionSource, GivenByUser, AdoptedFromLayers,
       NoRealWorldPosition,
       TakenFromAlignedLayer, AgreedByAllLayers, MeasuredAcrossProjection,
       RoundedFromMeasurement

include("ReportTerms.jl")

# What a report says about each decision and each layer. `public`, not exported, to match
# `Problem`/`LayerPlan`/`StudyAreaReport` themselves — these are read off a report, not written.
public AbstractProblemSeverity, ProblemNotice, ProblemWarning

public AbstractLayerFate, LayerKeptExactly, LayerAggregated, LayerResampled

# `public`, not exported, for the same reason as the neighbours above: a stage is **read off** a
# report, never written by a caller. This diverges from the usual "abstract type `public`, concrete
# leaves exported" rule, and deliberately — the rule is about types you *construct*, and nothing
# outside the package should be constructing a claim that a report is as-built.
public AbstractReportStage, AsInvestigated, AsBuilt

include("StudyAreaReport.jl")

# The report and its parts are meant to be read programmatically (`report.problems`,
# `report.layers`), not only printed, so they are supported names rather than internals.
public StudyAreaReport, LayerPlan, Problem, LayerCache
include("StudyArea.jl")

# **Ask anything that knows the grid how big it is** — a `StudyArea`, a raster, a layer, a habitat
# or an ecosystem. One dispatch hub (`_gridyx`) with a method per kind, so there is one answer and
# one place it comes from. Here because `GridHabitat` is the last type that hub names, and because
# `Layer.jl`'s `_uniformcellside` calls `getcellsizes` — a forward reference until now.
# `public`, not exported: supported API reached deliberately, not names wanted in every namespace.
public getspeciesstorage

# **Same reason again, one level up**: `GridHabitat` holds an `area::StudyArea{L}`, so both the
# area and the `StudyGrid` its parameter names must exist before `GridHabitat.jl`. `StudyArea.jl`
# itself cannot move up here — it names `GridHabitat` — so the two type definitions live in this
# file and every function that decides an area stays below.
# Neither is ever constructed by a user: a `StudyGrid` is read off a habitat that already exists,
# and `CellNames` is what its `placenames` hands back. Supported and documented, but not in the
# `using EcoSISTEM` namespace.
public StudyGrid, CellNames

# The study area — the simulation grid, decided from the layers before anything is built on it.
export StudyArea

# **The two halves of a simulation, and then the simulation itself.** A `GridHabitat` holds the
# study area above; an `Ecosystem` holds a habitat, a species list and a landscape of abundances.
# **Putting layers on a decided grid** — the read cache, `materialise` (inspection) and
# `_buildonarea` (building), kept in one file because the last two implement the same thing and
# have drifted apart three times. Placed here, before `GridHabitat.jl` calls `_buildonarea`.
include("materialise.jl")

public materialise

include("GridHabitat.jl")

# Supported, but not exported: a short name on a package's own concept, of the kind a user is likely
# to have taken for something else in their own script.
public totalsupply

export GridHabitat

include("cellgeometry.jl")

public getcellareas, getcellsizes, getcellarea, getcellsize, getcellat,
       getgridarea,
       getgridshape, getcellcount
export SpeciesList

include("Landscape.jl")
export GridLandscape, CachedGridLandscape

include("Ecosystem.jl")
export Ecosystem,
       CachedEcosystem,
       getsize,
       getregime,
       getnichefit,
       getgridsize,
       getdispersaldist,
       getdispersalvar,
       resetrate!,
       resettime!,
       getsupply,
       addspecies!

include("Simplify.jl")
export build_environment
export datamask, landmask, circlemask, shapemask

include("Traitfuns.jl")
export suitability, getpref, getdist, getnichefit, getregime

# Deprecated public API (trait line): `GaussTrait` → `NicheTolerance`, `Gauss`/`Trapeze`/`Unif` → `NicheSuitability`. Included
# late, after every type it shims; the shim names stay exported (above). See also
# `src/ClimatePref/deprecations.jl` for the ClimatePref submodule's deprecations.
include("deprecations.jl")

include("HabitatUpdate.jl")
export TempChange, RainfallChange, TempFluct, eraChange,
       worldclimChange

include("Scenarios.jl")
export SimpleScenario, FluctScenario, MultiScenario

include("Generate.jl")
export populate!,
       repopulate!,
       tolerancepopulate!,
       tolerancerepopulate!,
       emptypopulate!,
       resupply!,
       randomniches,
       update!

include("Helper.jl")
export simulate!,
       simulate_action!, simulate_record!, simulate_record_diversity!,
       generate_storage

include("Cache.jl")
export abundances, clearcache

include("DiversitySet.jl")
export DiversitySet, updatesimulation!, gettimes

include("AdditionalDiversity.jl")
export meta_simpson,
       meta_shannon, meta_speciesrichness, mean_abun, geom_mean_abun, sorenson,
       pd, makeunique

using Random
using Hwloc

# Number of contiguous species processed as an inner block in `update!`, sized so
# a block spans one CPU cache line. The abundance matrix is (species, cells) and
# column-major, so blocking species this way makes each per-cell access use a
# whole cache line instead of one element of it. Detected at startup in
# `__init__`; falls back to a value covering 128-byte lines if detection fails.
const _SPECIES_BLOCK = Ref(16)

"""
    species_blocksize()

Number of species iterated together as a contiguous inner block in `update!`,
chosen so one block spans a CPU cache line (`cachelinesize ÷ sizeof(Int)`).
"""
species_blocksize() = _SPECIES_BLOCK[]

function __init__()
    _SPECIES_BLOCK[] = try
        max(1, Hwloc.cachelinesize() ÷ sizeof(Int64))
    catch
        16
    end
    # Point RasterDataSources at its own subdirectory of our scratch space (unless the user has
    # already set RASTERDATASOURCES_PATH), keeping downloads under EcoSISTEM's scratch lifecycle.
    get!(ENV, "RASTERDATASOURCES_PATH") do
        return assetdir(RasterDataSources)
    end
    return nothing
end

abstract type MPIGridLandscape end
export MPIGridLandscape

abstract type MPIEcosystem{MPIGL <: MPIGridLandscape,
                           Part <: AbstractHabitat,
                           SL <: SpeciesList,
                           NF <: AbstractNicheFit} <:
              AbstractEcosystem{Part, SL, NF} end
export MPIEcosystem

"""
    gather_abundance(eco::MPIEcosystem)

Gather the full abundances matrix onto the root node.
"""
function gather_abundance end

"""
    gather_diversity(eco::MPIEcosystem, divmeasure::F, q) where F <: Function

Gather diversity calculated by `divmeasure` at value `q` from all MPI nodes onto
the root node (rank 0), combining subcommunity diversity values using a power
mean weighted by total abundances across nodes.
"""
function gather_diversity end
export gather_abundance, gather_diversity

function emptyMPIgridlandscape end
function synchronise_from_rows! end
function synchronise_from_cols! end
export emptyMPIgridlandscape, synchronise_from_rows!, synchronise_from_cols!

function _use_mpi()
    return !isnothing(Base.get_extension(@__MODULE__, :EcoSISTEMMPIExt)) &&
           _should_mpi()
end
# **THE WORKFLOW, LAST — the counterpart to `Ecology.jl` at the front.** `Ecology.jl` declares
# what the model is; this declares what you do with it, in the order you do it: investigate a grid,
# build the pieces, assemble them, run it. Nothing else lives there — the machinery behind each verb
# stays with its own concept, and the file says so.
include("actions.jl")

# Whether this process should build a distributed `MPIEcosystem`: the MPI extension overrides this
# to `MPI.Initialized() && MPI.Comm_size(MPI.COMM_WORLD) > 1`. The default (extension not loaded) is
# always false, so the base package builds a serial `Ecosystem` and never references MPI symbols.
function _should_mpi end

end
