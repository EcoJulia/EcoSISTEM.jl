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

export day, days, week, weeks, month, months, year, years, Rates,
January, February, March, April, May, June, July, August,
September, October, November, December

end

module ClimatePref

include("ClimatePref/ClimatePref.jl")

end

struct BiDict{A, B}
    forward::Dict{A, B}
    backward::Dict{B, A}

    BiDict{A, B}() where {A, B} = new{A, B}(Dict{A, B}(), Dict{B, A}())
end
import Base: getindex, setindex!
getindex(dict::BiDict{A,B}, elt::A) where {A, B} = dict.forward[elt]
getindex(dict::BiDict{A,B}, elt::B) where {A, B} = dict.backward[elt]
function setindex!(dict::BiDict{A,B}, val::B, elt::A)  where {A, B}
    dict.forward[elt] = val
    dict.backward[val] = elt
end

GLOBAL_typedict = BiDict{String, Type}()
GLOBAL_funcdict = BiDict{String, Function}()

include("Dist.jl")
export Trapezoid

include("Phylo.jl")
export assign_traits!, get_traits, resettraits!, reroot!

include("TraitRelationship.jl")
export TraitRelationship,multiplicativeTR2, multiplicativeTR3, Gauss,
 Match, NoRelContinuous, NoRelDiscrete, Trapeze, Unif,
 additiveTR2, additiveTR3

include("Habitats.jl")
export ContinuousHab,ContinuousTimeHab, DiscreteHab, HabitatCollection2, HabitatCollection3, tempgrad, raingrad

include("Energy.jl")
export SimpleRequirement, SizeRequirement, SolarRequirement, WaterRequirement, VolWaterRequirement, SimpleBudget, SolarBudget, SolarTimeBudget, WaterTimeBudget, VolWaterTimeBudget, ReqCollection2, BudgetCollection2

include("AbioticEnv.jl")
export GridAbioticEnv, simplenicheAE, tempgradAE, raingradAE, peakedgradAE, simplehabitatAE, degradedhabitatAE, eraAE, worldclimAE

include("Movement.jl")
export GaussianKernel, LongTailKernel, BirthOnlyMovement, AlwaysMovement, NoMovement, getkernel, Torus, Cylinder, NoBoundary

include("Traits.jl")
export GaussTrait, DiscreteTrait, TempBin,RainBin,
TraitCollection2, TraitCollection3,DiscreteEvolve,
ContinuousEvolve

include("Demographics.jl")
export PopGrowth, EqualPop, NoGrowth

include("SpeciesList.jl")
export SpeciesList

include("Landscape.jl")
export GridLandscape

include("Ecosystem.jl")
export Ecosystem, getsize, gethabitat,
gettraitrel, getgridsize, getdispersaldist, getdispersalvar,
resetrate!,resettime!, getbudget, addspecies!

include("Traitfuns.jl")
export TraitFun, getpref, gettraitrel, gethabitat

include("HabitatUpdate.jl")
export getchangefun, TempChange, RainChange, TempFluct, eraChange, worldclimChange

include("Scenarios.jl")
export SimpleScenario, TempIncrease, RandHabitatLoss!, ClustHabitatLoss!, DisturbanceScenario,
 HabitatDisturbance!, UniformDecline, ProportionalDecline, LargeDecline, RareDecline,
 CommonDecline, HabitatReplacement, Invasive, SusceptibleDecline

include("Generate.jl")
export populate!, repopulate!, traitpopulate!, traitrepopulate!, emptypopulate!, reenergise!, randomniches, update!, update_birth_move!, convert_coords, get_neighbours

include("Helper.jl")
export simulate!, simulate_record!,simulate_record_diversity!, expected_counts, generate_storage

include("AdditionalDiversity.jl")
export meta_simpson, meta_shannon, meta_speciesrichness, mean_abun, geom_mean_abun, sorenson, pd, makeunique

end
