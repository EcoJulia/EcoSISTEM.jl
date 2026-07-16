# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEM

# Path into package
path(path...; dir::String = "test") = joinpath(@__DIR__, "..", dir, path...)

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

include("Phylo.jl")
export assign_traits!, get_traits, resettraits!, reroot!

include("TraitRelationship.jl")
# `TraitRelationship` — dead export (likely meant `AbstractTraitRelationship`); revisit.
export multiplicativeTR2,
       multiplicativeTR3,
       Gauss,
       Match,
       NoRelContinuous,
       NoRelDiscrete,
       Trapeze,
       Unif,
       additiveTR2,
       additiveTR3,
       LCmatch

# Materialised layer family (Role × NicheAxis): the AbstractLayer types + HabitatUpdate,
# with the old *Hab names kept as aliases. Included BEFORE Habitats.jl, whose methods
# dispatch on those aliases.
include("Layer.jl")
export AbstractLayer, ContinuousLayer, DiscreteLayer,
       LayerCollection2, LayerCollection3, Unclassified
# back-compat habitat aliases over the layer types (defined in Layer.jl)
export ContinuousHab, ContinuousTimeHab, DiscreteHab,
       HabitatCollection2, HabitatCollection3

include("Habitats.jl")
export tempgrad, raingrad

include("Energy.jl")
export SimpleRequirement,
       SizeRequirement,
       SolarRequirement,
       WaterRequirement,
       VolWaterRequirement,
       SimpleBudget,
       SolarBudget,
       SolarTimeBudget,
       WaterBudget,
       VolWaterBudget,
       WaterTimeBudget,
       VolWaterTimeBudget,
       ReqCollection2,
       BudgetCollection2

include("AbioticEnv.jl")
export GridAbioticEnv,
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
       DiscreteTrait,
       TempBin,
       RainBin,
       TraitCollection2,
       TraitCollection3,
       DiscreteEvolve,
       ContinuousEvolve,
       LCtrait

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
       gethabitat,
       gettraitrel,
       getgridsize,
       getdispersaldist,
       getdispersalvar,
       resetrate!,
       resettime!,
       getbudget,
       addspecies!

include("Traitfuns.jl")
export traitfun, getpref, gettraitrel, gethabitat

include("HabitatUpdate.jl")
export TempChange, RainfallChange, TempFluct, eraChange,
       worldclimChange

include("Scenarios.jl")
export SimpleScenario, FluctScenario, MultiScenario

include("Generate.jl")
export populate!,
       repopulate!,
       traitpopulate!,
       traitrepopulate!,
       emptypopulate!,
       reenergise!,
       randomniches,
       update!,
# update_birth_move!,  # dead export (maybe meant `move!`?); revisit
       convert_coords,
       get_neighbours

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
    return _SPECIES_BLOCK[] = try
        max(1, Hwloc.cachelinesize() ÷ sizeof(Int64))
    catch
        16
    end
end

abstract type MPIGridLandscape end
export MPIGridLandscape

abstract type MPIEcosystem{MPIGL <: MPIGridLandscape,
                           Part <: AbstractAbiotic,
                           SL <: SpeciesList,
                           TR <: AbstractTraitRelationship} <:
              AbstractEcosystem{Part, SL, TR} end
export MPIEcosystem

function gather_abundance end
function gather_diversity end
export gather_abundance, gather_diversity

function emptyMPIgridlandscape end
function synchronise_from_rows! end
function synchronise_from_cols! end
export emptyMPIgridlandscape, synchronise_from_rows!, synchronise_from_cols!

end
