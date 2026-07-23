# SPDX-License-Identifier: LGPL-3.0-or-later

using Unitful
using Unitful.DefaultSymbols

# ---------------------------------------------------------------------------
# Layer roles
# ---------------------------------------------------------------------------
# A materialised layer plays one of two roles in a `GridHabitat`: a `Condition`
# (an environmental condition matched to species tolerances) or a `Resource` (a resource
# consumed by species). The role is a phantom marker used to keep the two
# type-distinguishable while sharing one storage implementation.

"""
    Role

Abstract supertype of the layer-role markers [`Condition`](@ref) and [`Resource`](@ref).
"""
abstract type Role end

"""
    Condition <: Role

Marker for a layer used as an environmental condition (matched to species tolerances).
"""
struct Condition <: Role end

"""
    Resource <: Role

Marker for a layer used as a consumable resource.
"""
struct Resource <: Role end

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
"""    MeanTemperature <: AbstractTemperature — annual mean temperature, the mean of the 12 monthly mean temperatures (Climate `tavg`, BioClim 1; K). """
struct MeanTemperature <: AbstractTemperature end
"""    MinTemperature <: AbstractTemperature — minimum temperature of the coldest month (Climate `tmin`, BioClim 6; K). """
struct MinTemperature <: AbstractTemperature end
"""    MaxTemperature <: AbstractTemperature — maximum temperature of the warmest month (Climate `tmax`, BioClim 5; K). """
struct MaxTemperature <: AbstractTemperature end
"""    DiurnalTemperatureRange <: AbstractTemperature — mean diurnal range, the mean of the monthly (max − min) temperatures (BioClim 2; K). """
struct DiurnalTemperatureRange <: AbstractTemperature end
"""    TemperatureSeasonality <: AbstractTemperature — temperature seasonality, the standard deviation of the monthly mean temperatures (BioClim 4). """
struct TemperatureSeasonality <: AbstractTemperature end
"""    TemperatureAnnualRange <: AbstractTemperature — temperature annual range, the warmest month's max minus the coldest month's min (BioClim 5 − BioClim 6; BioClim 7; K). """
struct TemperatureAnnualRange <: AbstractTemperature end
"""    WettestQuarterTemperature <: AbstractTemperature — mean temperature of the wettest 3 consecutive months (BioClim 8; K). """
struct WettestQuarterTemperature <: AbstractTemperature end
"""    DriestQuarterTemperature <: AbstractTemperature — mean temperature of the driest 3 consecutive months (BioClim 9; K). """
struct DriestQuarterTemperature <: AbstractTemperature end
"""    WarmestQuarterTemperature <: AbstractTemperature — mean temperature of the warmest 3 consecutive months (BioClim 10; K). """
struct WarmestQuarterTemperature <: AbstractTemperature end
"""    ColdestQuarterTemperature <: AbstractTemperature — mean temperature of the coldest 3 consecutive months (BioClim 11; K). """
struct ColdestQuarterTemperature <: AbstractTemperature end

# --- Precipitation axes (all mm; provide a water supply/demand) --------
"""    AbstractPrecipitation <: NicheAxis — precipitation niche axes (all in mm) """
abstract type AbstractPrecipitation <: NicheAxis end
"""    Precipitation <: AbstractPrecipitation — annual precipitation, the sum of the 12 monthly precipitation totals (Climate `prec`, BioClim 12; mm). """
struct Precipitation <: AbstractPrecipitation end
"""    WettestMonthPrecipitation <: AbstractPrecipitation — precipitation of the wettest month (BioClim 13; mm). """
struct WettestMonthPrecipitation <: AbstractPrecipitation end
"""    DriestMonthPrecipitation <: AbstractPrecipitation — precipitation of the driest month (BioClim 14; mm). """
struct DriestMonthPrecipitation <: AbstractPrecipitation end
"""    WettestQuarterPrecipitation <: AbstractPrecipitation — total precipitation of the wettest 3 consecutive months (BioClim 16; mm). """
struct WettestQuarterPrecipitation <: AbstractPrecipitation end
"""    DriestQuarterPrecipitation <: AbstractPrecipitation — total precipitation of the driest 3 consecutive months (BioClim 17; mm). """
struct DriestQuarterPrecipitation <: AbstractPrecipitation end
"""    WarmestQuarterPrecipitation <: AbstractPrecipitation — total precipitation of the warmest 3 consecutive months (BioClim 18; mm). """
struct WarmestQuarterPrecipitation <: AbstractPrecipitation end
"""    ColdestQuarterPrecipitation <: AbstractPrecipitation — total precipitation of the coldest 3 consecutive months (BioClim 19; mm). """
struct ColdestQuarterPrecipitation <: AbstractPrecipitation end

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
"""    Heterogeneity <: NicheAxis (spatial regime-heterogeneity metrics; dimensionless) """
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
# what applies. The supply/demand/change types they name live in later-included files
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
    supplytype(::NicheAxis)

The `AbstractSupply` concrete type for this axis when used as a `Resource` resource, or
`nothing` if the axis is not a consumable resource (so it errors clearly rather than
silently guessing one if a supply is nonetheless requested).
"""
supplytype(::NicheAxis) = nothing
supplytype(::AbstractPrecipitation) = WaterSupply
supplytype(::SolarRadiation) = SolarSupply
supplytype(::VolumetricWater) = VolWaterSupply

"""
    demandtype(::NicheAxis)

The `AbstractDemand` concrete type a species uses to consume this axis' resource, or
`nothing` if the axis is not a resource.
"""
demandtype(::NicheAxis) = nothing
demandtype(::AbstractPrecipitation) = WaterDemand
demandtype(::SolarRadiation) = SolarDemand
demandtype(::VolumetricWater) = VolWaterDemand

"""
    dynamics(::NicheAxis)

The default per-timestep change function for a layer on this axis (default
[`NoChange`](@ref)); overridden per layer at build time.
"""
dynamics(::NicheAxis) = NoChange
