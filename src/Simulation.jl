module Simulation

module Units

import Unitful
using Unitful: @unit
day = Unitful.d
week = Unitful.wk
@unit month "month" Month 2.628e6 * Unitful.s false
@unit year "year" Year 31536000 * Unitful.s false

const days = day
const weeks = week
const months = month
const years = year
const January = 0month
const February = 1month
const March = 2months
const April = 3months
const May = 4months
const June = 5months
const July = 6months
const August = 7months
const September = 8months
const October = 9months
const November = 10months
const December = 11months

const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(Units)
end

export day, days, week, weeks, month, months, year, years, Rates,
January, February, March, April, May, June, July, August,
September, October, November, December

end

module ClimatePref

include("ClimatePref/ClimatePref.jl")

end

include("Biodiversity/Dist.jl")
export Trapezoid

include("Biodiversity/Phylo.jl")
export assign_traits!, get_traits, resettraits!, reroot!

include("Biodiversity/TraitRelationship.jl")
export TraitRelationship,multiplicativeTR2, multiplicativeTR3, Gauss,
 Match, NoRelContinuous, NoRelDiscrete, Trapeze, Unif,
 additiveTR2, additiveTR3

include("Biodiversity/Habitats.jl")
export ContinuousHab,ContinuousTimeHab, DiscreteHab, HabitatCollection2, HabitatCollection3, tempgrad, raingrad

include("Biodiversity/Energy.jl")
export SimpleRequirement, SizeRequirement, SolarRequirement, WaterRequirement, VolWaterRequirement, SimpleBudget, SolarBudget, SolarTimeBudget, WaterTimeBudget, VolWaterTimeBudget, ReqCollection2, BudgetCollection2

include("Biodiversity/AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, raingradAE, peakedgradAE, simplehabitatAE, degradedhabitatAE, eraAE, worldclimAE

include("Biodiversity/Movement.jl")
export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernel, Torus, Cylinder, NoBoundary

include("Biodiversity/Traits.jl")
export GaussTrait, DiscreteTrait, TempBin,RainBin,
TraitCollection2, TraitCollection3,DiscreteEvolve,
ContinuousEvolve

include("Biodiversity/Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("Biodiversity/SpeciesList.jl")
export SpeciesList

include("Biodiversity/Landscape.jl")
export GridLandscape, CachedGridLandscape

include("Biodiversity/MPILandscape.jl")
export MPIGridLandscape

include("Biodiversity/Ecosystem.jl")
export Ecosystem, CachedEcosystem, getsize, gethabitat, gettraitrel, getgridsize,
 getdispersaldist, getdispersalvar, resetrate!,resettime!, getbudget, addspecies!

include("Biodiversity/MPIEcosystem.jl")
export MPIEcosystem, gather_abundance!, gather_diversity

include("Biodiversity/Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("Biodiversity/HabitatUpdate.jl")
export getchangefun, TempChange, RainChange, TempFluct, eraChange, worldclimChange

include("Biodiversity/Scenarios.jl")
export SimpleScenario, FluctScenario, MultiScenario

include("Biodiversity/Generate.jl")
export populate!, repopulate!, traitpopulate!, traitrepopulate!, emptypopulate!, reenergise!, randomniches, update!, update_birth_move!, convert_coords, get_neighbours

include("Biodiversity/MPIGenerate.jl")

include("Biodiversity/Helper.jl")
export simulate!, simulate_record!,simulate_record_diversity!, expected_counts, generate_storage

include("Biodiversity/Cache.jl")
export abundances, clearcache

include("Biodiversity/DiversitySet.jl")
export DiversitySet, updatesimulation!, gettimes

include("Biodiversity/AdditionalDiversity.jl")
export meta_simpson, meta_shannon, meta_speciesrichness, mean_abun, geom_mean_abun, sorenson, pd, makeunique

include("Epidemiology/EpiControl.jl")
export NoControl

include("Epidemiology/EpiEnv.jl")
export GridEpiEnv, simplehabitatAE

include("Epidemiology/EpiParams.jl")
export SIRGrowth, SEI2HRDGrowth, transition

include("Epidemiology/EpiList.jl")
export EpiList, SIR, SEI2HRD

include("Epidemiology/EpiLandscape.jl")
export EpiLandscape

include("Epidemiology/EpiSystem.jl")
export EpiSystem

include("Epidemiology/EpiTraits.jl")

include("Epidemiology/EpiGenerate.jl")
export populate!

include("Epidemiology/EpiHelper.jl")
export simulate!, simulate_record!

end
