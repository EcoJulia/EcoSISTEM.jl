# SPDX-License-Identifier: LGPL-3.0-or-later

# ---------------------------------------------------------------------------
# Niche axes
# ---------------------------------------------------------------------------
# A `NicheAxis` names *what* a layer measures — an axis of the (Hutchinsonian) niche — in
# plain English, independent of its physical unit, and orthogonal to `Role` (which says how
# a layer is *used*). A layer and the trait matched to it share the same `NicheAxis`, so
# "same axis" is the matching rule rather than "same unit". Users add their own with
# `struct MyAxis <: NicheAxis end` plus, if needed, the interface overrides below.

"""
    NicheAxis

Abstract supertype of the niche-axis markers naming what a layer measures (e.g.
[`MeanTemperature`](@ref), [`Precipitation`](@ref)), independent of its unit and orthogonal
to [`Role`](@ref). Extend it with your own `struct MyAxis <: NicheAxis end`. Related axes are
grouped under an abstract intermediate ([`AbstractTemperature`](@ref),
[`AbstractPrecipitation`](@ref)) that carries their shared interface.
"""
abstract type NicheAxis end

# --- Temperature axes (all measured in K) ----------------------------------
"""    AbstractTemperature <: NicheAxis — temperature niche axes (all in K) """
abstract type AbstractTemperature <: NicheAxis end
struct MeanTemperature <: AbstractTemperature end            # Climate tavg, BioClim 1
struct MinTemperature <: AbstractTemperature end             # Climate tmin, BioClim 6
struct MaxTemperature <: AbstractTemperature end             # Climate tmax, BioClim 5
struct DiurnalTemperatureRange <: AbstractTemperature end    # BioClim 2
struct TemperatureSeasonality <: AbstractTemperature end     # BioClim 4
struct TemperatureAnnualRange <: AbstractTemperature end     # BioClim 7
struct WettestQuarterTemperature <: AbstractTemperature end  # BioClim 8
struct DriestQuarterTemperature <: AbstractTemperature end   # BioClim 9
struct WarmestQuarterTemperature <: AbstractTemperature end  # BioClim 10
struct ColdestQuarterTemperature <: AbstractTemperature end  # BioClim 11

# --- Precipitation axes (all mm; provide a water budget/requirement) --------
"""    AbstractPrecipitation <: NicheAxis — precipitation niche axes (all in mm) """
abstract type AbstractPrecipitation <: NicheAxis end
struct Precipitation <: AbstractPrecipitation end                # Climate prec, BioClim 12
struct WettestMonthPrecipitation <: AbstractPrecipitation end    # BioClim 13
struct DriestMonthPrecipitation <: AbstractPrecipitation end     # BioClim 14
struct WettestQuarterPrecipitation <: AbstractPrecipitation end  # BioClim 16
struct DriestQuarterPrecipitation <: AbstractPrecipitation end   # BioClim 17
struct WarmestQuarterPrecipitation <: AbstractPrecipitation end  # BioClim 18
struct ColdestQuarterPrecipitation <: AbstractPrecipitation end  # BioClim 19

# --- Standalone axes -------------------------------------------------------
"""    SolarRadiation <: NicheAxis """
struct SolarRadiation <: NicheAxis end
"""    WindSpeed <: NicheAxis """
struct WindSpeed <: NicheAxis end
"""    VaporPressure <: NicheAxis """
struct VaporPressure <: NicheAxis end
"""    Isothermality <: NicheAxis (dimensionless) """
struct Isothermality <: NicheAxis end
"""    PrecipitationSeasonality <: NicheAxis (dimensionless) """
struct PrecipitationSeasonality <: NicheAxis end
"""    Heterogeneity <: NicheAxis (spatial habitat-heterogeneity metrics; dimensionless) """
struct Heterogeneity <: NicheAxis end
"""    LandType <: NicheAxis (categorical land-cover / land-use classes) """
struct LandType <: NicheAxis end
"""    Altitude <: NicheAxis (elevation above sea level) """
struct Altitude <: NicheAxis end

# Retained pending the water-axis design — water spans several dimensions (depth / mass per
# area / flux / volume) that still need a unified treatment (see the plan).
"""    VolumetricWater <: NicheAxis """
struct VolumetricWater <: NicheAxis end

# --- The NicheAxis interface -----------------------------------------------
# The small set of hooks that replace unit-based inference of a layer's meaning. Each has a
# safe default; a group (via its abstract intermediate) or an individual axis overrides only
# what applies. The budget/requirement/change types they name live in later-included files
# and are resolved at call time (nothing here runs at load).

"""
    canonicalunit(::NicheAxis)

The unit a layer on this axis is normalised to when materialised (e.g. temperature → `K`),
or `nothing` to leave values as-is.
"""
canonicalunit(::NicheAxis) = nothing
canonicalunit(::AbstractTemperature) = K
canonicalunit(::AbstractPrecipitation) = mm
canonicalunit(::Altitude) = m

"""
    budgettype(::NicheAxis)

The `AbstractBudget` concrete type for this axis when used as a `Budget` resource, or
`nothing` if the axis is not a consumable resource (so it errors clearly rather than
silently guessing one if a budget is nonetheless requested).
"""
budgettype(::NicheAxis) = nothing
budgettype(::AbstractPrecipitation) = WaterBudget
budgettype(::SolarRadiation) = SolarBudget
budgettype(::VolumetricWater) = VolWaterBudget

"""
    requirementtype(::NicheAxis)

The `AbstractRequirement` concrete type a species uses to consume this axis' resource, or
`nothing` if the axis is not a resource.
"""
requirementtype(::NicheAxis) = nothing
requirementtype(::AbstractPrecipitation) = WaterRequirement
requirementtype(::SolarRadiation) = SolarRequirement
requirementtype(::VolumetricWater) = VolWaterRequirement

"""
    dynamics(::NicheAxis)

The default per-timestep change function for a layer on this axis (default
[`NoChange`](@ref)); overridden per layer at build time.
"""
dynamics(::NicheAxis) = NoChange
