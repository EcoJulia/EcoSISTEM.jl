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

include("Layer.jl")
export condition, resource
export AbstractLayer, ContinuousLayer, DiscreteLayer,
       LayerCollection2, LayerCollection3, Unclassified
# back-compat regime + supply aliases over the layer types (defined in Layer.jl)
export ContinuousRegime, ContinuousTimeRegime, DiscreteRegime,
       RegimeCollection2, RegimeCollection3
export SimpleSupply,
       SolarSupply,
       SolarTimeSupply,
       WaterSupply,
       WaterTimeSupply,
       SupplyCollection2

include("Habitats.jl")
export tempgrad, raingrad

include("Energy.jl")
export SimpleDemand,
       SizeDemand,
       SolarDemand,
       WaterDemand,
       DemandCollection2

include("AbioticEnv.jl")
export GridHabitat,
       simplenichehabitat,
       tempgradhabitat,
       raingradhabitat,
       peakedgradhabitat,
       simplehabitat,
       erahabitat,
       worldclimhabitat,
       bioclimhabitat,
       landcoverhabitat

include("Movement.jl")
export GaussianKernel,
       LongTailKernel,
       BirthOnlyMovement,
       AlwaysMovement,
       NoMovement,
       getkernels,
       Torus,
       Cylinder,
       NoBoundary

include("Traits.jl")
export GaussTrait,
       DiscreteTolerance,
       NicheTolerance,
       TempTolerance,
       RainTolerance,
       ToleranceCollection2,
       ToleranceCollection3,
       DiscreteEvolve,
       ContinuousEvolve,
       LandCoverTolerance

include("Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("SpeciesList.jl")
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
