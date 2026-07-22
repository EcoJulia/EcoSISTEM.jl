# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEM

using Scratch: get_scratch!
import RasterDataSources

"""
    EcoSISTEM.assetdir(mod::Module = EcoSISTEM)

Path to the `mod` subdirectory of EcoSISTEM's single Scratch.jl space, for storing `mod`-related
files (e.g. downloaded data) outside the repository. Using one EcoSISTEM-owned space with
per-package subdirectories (rather than a separate space per package) keeps the whole cache under
EcoSISTEM's lifecycle — created on first use, reclaimed by `Pkg.gc()` when EcoSISTEM is removed.
`mod` defaults to EcoSISTEM's own subdirectory.

On load, EcoSISTEM sets `RASTERDATASOURCES_PATH` to `assetdir(RasterDataSources)` (in `__init__`).
"""
function assetdir(mod::Module = EcoSISTEM)
    return mkpath(joinpath(get_scratch!(EcoSISTEM, "assets"),
                           string(nameof(mod))))
end
public assetdir

# EcoSISTEM.Units sub-module
include("Units/Units.jl")

# Layer roles + specs (recipes) — defined early so later modules can build on them.
include("NicheInfo.jl")
export NicheAxis, AbstractTemperature, AbstractPrecipitation,
       MeanTemperature, MinTemperature, MaxTemperature, DiurnalTemperatureRange,
       TemperatureSeasonality, TemperatureAnnualRange,
       WettestQuarterTemperature,
       DriestQuarterTemperature, WarmestQuarterTemperature,
       ColdestQuarterTemperature,
       Precipitation, WettestMonthPrecipitation, DriestMonthPrecipitation,
       WettestQuarterPrecipitation, DriestQuarterPrecipitation,
       WarmestQuarterPrecipitation, ColdestQuarterPrecipitation,
       SolarRadiation, WindSpeed, VaporPressure, Isothermality,
       PrecipitationSeasonality, Heterogeneity, LandType, Altitude,
       VolumetricWater

# Geographic coordinate type — used across the ClimatePref sub-module and the main module, so it is
# defined here (before ClimatePref) rather than inside it.
include("Coordinates.jl")
export LatLong

# EcoSISTEM.ClimatePref sub-module
include("ClimatePref/ClimatePref.jl")

# DataPipeline extension
"""
    EcoSISTEM.unziptemp(path::String)

Helper function for the FAIR Data Pipeline to unzip files that are stored as
zips to a temporary folder.
"""
function unziptemp end

include("Dist.jl")
export Trapezoid
public read_distribution, param_units

include("Phylo.jl")
export assign_traits!, get_traits, resettraits!, reroot!

include("TraitRelationship.jl")
export multiplicativeFit2,
       multiplicativeFit3,
       Gauss,
       MatchSuitability,
       NoFitContinuous,
       NoFitDiscrete,
       NicheSuitability,
       Trapeze,
       Unif,
       additiveFit2,
       additiveFit3,
       LCsuitability

# Materialised layer family (Role × NicheAxis): the AbstractLayer types + HabitatUpdate,
# with the old *Hab names kept as aliases. Included BEFORE Habitats.jl, whose methods
# dispatch on those aliases.
include("Layer.jl")
export AbstractLayer, ContinuousLayer, DiscreteLayer,
       LayerCollection2, LayerCollection3, Unclassified
# back-compat regime + supply aliases over the layer types (defined in Layer.jl)
export ContinuousRegime, ContinuousTimeRegime, DiscreteRegime,
       RegimeCollection2, RegimeCollection3
export SimpleSupply,
       SolarSupply,
       SolarTimeSupply,
       WaterSupply,
       VolWaterSupply,
       WaterTimeSupply,
       VolWaterTimeSupply,
       SupplyCollection2

include("Habitats.jl")
export tempgrad, raingrad

include("Energy.jl")
export SimpleDemand,
       SizeDemand,
       SolarDemand,
       WaterDemand,
       VolWaterDemand,
       DemandCollection2

include("AbioticEnv.jl")
export GridHabitat,
       simplenicheAE,
       tempgradAE,
       raingradAE,
       peakedgradAE,
       simplehabitatAE,
       eraAE,
       worldclimAE,
       bioclimAE,
       lcAE

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
       LCtolerance

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
                           TR <: AbstractNicheFit} <:
              AbstractEcosystem{Part, SL, TR} end
export MPIEcosystem

function gather_abundance end
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

# Whether this process should build a distributed `MPIEcosystem`: the MPI extension overrides this
# to `MPI.Initialized() && MPI.Comm_size(MPI.COMM_WORLD) > 1`. The default (extension not loaded) is
# always false, so the base package builds a serial `Ecosystem` and never references MPI symbols.
function _should_mpi end

end
