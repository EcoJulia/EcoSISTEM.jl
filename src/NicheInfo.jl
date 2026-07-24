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
[`Temperature`](@ref), [`Precipitation`](@ref)), independent of its unit and orthogonal
to [`Role`](@ref). Extend it with your own `struct MyAxis <: NicheAxis end`. Related axes are
grouped under an abstract intermediate named `…Axis` ([`TemperatureAxis`](@ref),
[`WaterAxis`](@ref)) that carries their shared interface.
"""
abstract type NicheAxis end

# --- Temperature axes ------------------------------------------------------
# Group default is K; the two dimensionless leaves and the degree-day sum override it (see the
# `canonicalunit` overrides below). Level, range and seasonality follow the aggregation policy:
# max/min/mean collapse to one level axis, a range gets its own axis, seasonality (an SD) stays apart.
"""    TemperatureAxis <: NicheAxis — temperature-related niche axes (K unless a leaf overrides). """
abstract type TemperatureAxis <: NicheAxis end
"""    Temperature <: TemperatureAxis — an absolute temperature (annual/monthly/quarter mean, min, max, or growing-season mean; BioClim 1/5/6/8-11, Climate tavg/tmin/tmax/tas/tasmin/tasmax, BioClimPlus gst; K). """
struct Temperature <: TemperatureAxis end
"""    TemperatureRange <: TemperatureAxis — a temperature range (max − min): mean diurnal range or annual range (BioClim 2 and 7; K interval). """
struct TemperatureRange <: TemperatureAxis end
"""    TemperatureSeasonality <: TemperatureAxis — temperature seasonality, the standard deviation of the monthly mean temperatures (BioClim 4; K). """
struct TemperatureSeasonality <: TemperatureAxis end
"""    CumulativeHeat <: TemperatureAxis — accumulated growing degree-days above a threshold (BioClimPlus gdd0/gdd5/gdd10; K·day). """
struct CumulativeHeat <: TemperatureAxis end
"""    Isothermality <: TemperatureAxis — isothermality, 100 × (mean diurnal range / annual range) (BioClim 3; dimensionless). """
struct Isothermality <: TemperatureAxis end
"""    FrostChangeFrequency <: TemperatureAxis — frost change frequency, the number of days the temperature crosses 0 °C (BioClimPlus fcf; dimensionless). """
struct FrostChangeFrequency <: TemperatureAxis end

# --- Water axes ------------------------------------------------------------
# `WaterAxis` is the umbrella for everything that measures water in any form. It is a pure grouping
# node — it carries NO unit or resource of its own, so it can span mm precipitation, kPa/dimensionless
# humidity and kg·m⁻²[·month⁻¹] fluxes/stocks without conflict; each leaf (or sub-group) declares its own.
"""    WaterAxis <: NicheAxis — umbrella for all water-measuring niche axes (precipitation, humidity, evapotranspiration, moisture, water stocks). A pure grouping node with no unit or resource of its own. """
abstract type WaterAxis <: NicheAxis end
"""    PrecipitationAxis <: WaterAxis — precipitation niche axes. A topical group; the `mm` unit and the water supply/demand live on the [`Precipitation`](@ref) leaf. """
abstract type PrecipitationAxis <: WaterAxis end
"""    Precipitation <: PrecipitationAxis — a precipitation amount over a month, quarter or year (BioClim 12-19, Climate prec/pr; mm). Carries the `mm` unit and provides a water supply/demand. """
struct Precipitation <: PrecipitationAxis end
"""    PrecipitationSeasonality <: PrecipitationAxis — precipitation seasonality, the coefficient of variation of the monthly precipitation totals (BioClim 15; dimensionless). """
struct PrecipitationSeasonality <: PrecipitationAxis end
"""    HumidityAxis <: WaterAxis — atmospheric-water niche axes (vapour pressure, its deficit, relative humidity). """
abstract type HumidityAxis <: WaterAxis end
"""    VaporPressure <: HumidityAxis — near-surface water-vapour partial pressure (Climate vapr; kPa). """
struct VaporPressure <: HumidityAxis end
"""    VaporPressureDeficitAxis <: HumidityAxis — vapour-pressure-deficit niche axes (level and range). """
abstract type VaporPressureDeficitAxis <: HumidityAxis end
"""    VaporPressureDeficit <: VaporPressureDeficitAxis — near-surface vapour pressure deficit (Climate vpd, BioClimPlus vpd_max/mean/min; Pa). """
struct VaporPressureDeficit <: VaporPressureDeficitAxis end
"""    VaporPressureDeficitRange <: VaporPressureDeficitAxis — annual range of the monthly vapour pressure deficit (BioClimPlus vpd_range; Pa). """
struct VaporPressureDeficitRange <: VaporPressureDeficitAxis end
"""    RelativeHumidityAxis <: HumidityAxis — relative-humidity niche axes (level and range). """
abstract type RelativeHumidityAxis <: HumidityAxis end
"""    RelativeHumidity <: RelativeHumidityAxis — near-surface relative humidity (Climate hurs, BioClimPlus hurs_max/mean/min; dimensionless). """
struct RelativeHumidity <: RelativeHumidityAxis end
"""    RelativeHumidityRange <: RelativeHumidityAxis — annual range of the monthly relative humidity (BioClimPlus hurs_range; dimensionless). """
struct RelativeHumidityRange <: RelativeHumidityAxis end
"""    EvapotranspirationAxis <: WaterAxis — evapotranspiration niche axes (level and range). """
abstract type EvapotranspirationAxis <: WaterAxis end
"""    Evapotranspiration <: EvapotranspirationAxis — potential evapotranspiration (Climate pet, BioClimPlus pet_penman_max/mean/min; kg·m⁻²·month⁻¹). """
struct Evapotranspiration <: EvapotranspirationAxis end
"""    EvapotranspirationRange <: EvapotranspirationAxis — annual range of the monthly potential evapotranspiration (BioClimPlus pet_penman_range; kg·m⁻²·month⁻¹). """
struct EvapotranspirationRange <: EvapotranspirationAxis end
"""    ClimateMoistureAxis <: WaterAxis — climate-moisture-index niche axes (level and range). """
abstract type ClimateMoistureAxis <: WaterAxis end
"""    ClimateMoisture <: ClimateMoistureAxis — climate moisture index, precipitation minus potential evapotranspiration (Climate cmi, BioClimPlus cmi_max/mean/min; kg·m⁻²·month⁻¹). """
struct ClimateMoisture <: ClimateMoistureAxis end
"""    ClimateMoistureRange <: ClimateMoistureAxis — annual range of the monthly climate moisture index (BioClimPlus cmi_range; kg·m⁻²·month⁻¹). """
struct ClimateMoistureRange <: ClimateMoistureAxis end
"""    SnowWaterEquivalent <: WaterAxis — snow water equivalent, the liquid-water amount of the snowpack (BioClimPlus swe; kg·m⁻²). """
struct SnowWaterEquivalent <: WaterAxis end
"""    SiteWaterBalance <: WaterAxis — cumulative site water balance over the year, capped by soil water-holding capacity (BioClimPlus swb; kg·m⁻²). """
struct SiteWaterBalance <: WaterAxis end
"""    GrowingSeasonPrecipitation <: WaterAxis — precipitation accumulated over the growing season (BioClimPlus gsp; kg·m⁻²). """
struct GrowingSeasonPrecipitation <: WaterAxis end

# --- Solar radiation / wind speed / cloud cover ----------------------------
"""    SolarRadiationAxis <: NicheAxis — solar-radiation niche axes (level and range). """
abstract type SolarRadiationAxis <: NicheAxis end
"""    SolarRadiation <: SolarRadiationAxis — surface downwelling shortwave (solar) radiation flux (Climate srad/rsds, BioClimPlus rsds_max/mean/min). Provides a solar supply/demand. """
struct SolarRadiation <: SolarRadiationAxis end
"""    SolarRadiationRange <: SolarRadiationAxis — annual range of the monthly solar radiation flux (BioClimPlus rsds_range). """
struct SolarRadiationRange <: SolarRadiationAxis end
"""    WindSpeedAxis <: NicheAxis — wind-speed niche axes (level and range). """
abstract type WindSpeedAxis <: NicheAxis end
"""    WindSpeed <: WindSpeedAxis — near-surface (10 m) wind speed (Climate wind/sfcWind, BioClimPlus sfcWind_max/mean/min; m·s⁻¹). """
struct WindSpeed <: WindSpeedAxis end
"""    WindSpeedRange <: WindSpeedAxis — annual range of the monthly near-surface wind speed (BioClimPlus sfcWind_range; m·s⁻¹). """
struct WindSpeedRange <: WindSpeedAxis end
"""    CloudCoverAxis <: NicheAxis — cloud-cover niche axes (level and range). """
abstract type CloudCoverAxis <: NicheAxis end
"""    CloudCover <: CloudCoverAxis — total cloud cover (Climate clt, BioClimPlus clt_max/mean/min; dimensionless). """
struct CloudCover <: CloudCoverAxis end
"""    CloudCoverRange <: CloudCoverAxis — annual range of the monthly total cloud cover (BioClimPlus clt_range; dimensionless). """
struct CloudCoverRange <: CloudCoverAxis end

# --- Day-based, carbon and categorical axes --------------------------------
"""    DayAxis <: NicheAxis — day-based phenology niche axes (day-of-year positions and day counts; day). """
abstract type DayAxis <: NicheAxis end
"""    DayOfYear <: DayAxis — an ordinal day-of-year: first/last growing day or first/last day above a degree-day threshold (BioClimPlus fgd/lgd/gdgfgd*/gddlgd*; day). """
struct DayOfYear <: DayAxis end
"""    DayRange <: DayAxis — a count or span of days: growing-season length, snow-cover days or number of days above a threshold (BioClimPlus gsl/scd/ngd*; day). """
struct DayRange <: DayAxis end
"""    CarbonAxis <: NicheAxis — carbon-cycle niche axes. """
abstract type CarbonAxis <: NicheAxis end
"""    CarbonFlux <: CarbonAxis — net primary productivity, a carbon flux (BioClimPlus npp; g·m⁻²·year⁻¹). """
struct CarbonFlux <: CarbonAxis end
"""    TypologyAxis <: NicheAxis — categorical classification niche axes (discrete class labels; dimensionless, used with a `DiscreteRegime`). """
abstract type TypologyAxis <: NicheAxis end
"""    LandCoverTypology <: TypologyAxis — categorical land-cover / land-use classes (EarthEnv LandCover 1-12). """
struct LandCoverTypology <: TypologyAxis end
"""    ClimateTypology <: TypologyAxis — categorical climate-classification classes (BioClimPlus kg0-5: Köppen-Geiger, Wissmann, Thornthwaite, Troll-Pfaffen). """
struct ClimateTypology <: TypologyAxis end

# --- Other standalone axes -------------------------------------------------
"""    Heterogeneity <: NicheAxis — spatial habitat-heterogeneity metrics of EVI (EarthEnv HabitatHeterogeneity; dimensionless). """
struct Heterogeneity <: NicheAxis end
"""    Altitude <: NicheAxis — elevation above sea level (WorldClim Elevation; m). """
struct Altitude <: NicheAxis end

# --- The NicheAxis interface -----------------------------------------------
# The small set of hooks that replace unit-based inference of a layer's meaning. Each has a
# safe default; a group (via its abstract intermediate) or an individual axis overrides only
# what applies. The supply/demand/change types they name live in later-included files
# and are resolved at call time (nothing here runs at load).

"""
    canonicalunit(::NicheAxis)
    canonicalunit(::Type{<:Role}, ::NicheAxis)

The unit a layer on this axis is normalised to when materialised (e.g. temperature → `K`),
or `nothing` to leave values as-is. The 2-argument form is role-aware: a `Condition` (a
niche tolerance — a descriptive climatological normal) and a `Resource` (a literal
consumable pool, replenished over time) can legitimately want a different canonical unit for
the same axis — e.g. `Precipitation` stays a bare depth (`mm`) as a `Condition`, but is a
genuine volumetric flow (`L/day`) as a `Resource`. It defaults to the 1-argument form for
any role/axis combination without a specific override, so existing (implicitly
`Condition`-role) call sites are unaffected.
"""
canonicalunit(::NicheAxis) = nothing
canonicalunit(::Type{<:Role}, A::NicheAxis) = canonicalunit(A)
# TemperatureAxis defaults to K, covering the three real Kelvin leaves; the two dimensionless
# leaves and the degree-day sum override it. mm/Water live on the `Precipitation` leaf (not the
# topical `PrecipitationAxis`), so `PrecipitationSeasonality` inherits nothing.
canonicalunit(::TemperatureAxis) = K
canonicalunit(::CumulativeHeat) = K * Unitful.d
canonicalunit(::Isothermality) = NoUnits
canonicalunit(::FrostChangeFrequency) = NoUnits
canonicalunit(::Precipitation) = mm
canonicalunit(::Altitude) = m
canonicalunit(::DayAxis) = Unitful.d
canonicalunit(::TypologyAxis) = NoUnits
# WorldClim's srad is kJ·m⁻²·day⁻¹, CHELSA's rsds*/BioClimPlus's rsds_* are MJ·m⁻²·day⁻¹ — without
# this override the two sources' regimes are never reconciled to one scale (see _tocanon). kJ matches
# the existing Resource-role choice below and WorldClim's native scale. SolarRadiationRange shares the
# same physical unit as the level (an annual range of the same flux), so this lives on the group.
canonicalunit(::SolarRadiationAxis) = kJ / (m^2 * Unitful.d)

# Resource-role canonical units: only axes with a dedicated Supply type need one (currently
# Precipitation, SolarRadiation) — a genuine consumption rate, not the Condition-role's bare
# amount. `day` is the shared canonical timebase (see NEWS/docs): it's the unit solar
# irradiance is natively defined in, and water is matched to it for consistency, not because
# water has its own intrinsic daily period. (A third Resource-role override, for the free/simple
# per-individual case with no physical axis, lives in Layer.jl next to the `Unclassified` axis.)
canonicalunit(::Type{Resource}, ::Precipitation) = Unitful.L / Unitful.d
canonicalunit(::Type{Resource}, ::SolarRadiation) = kJ / Unitful.d

"""
    supplytype(::NicheAxis)

The `AbstractSupply` concrete type for this axis when used as a `Resource` resource, or
`nothing` if the axis is not a consumable resource (so it errors clearly rather than
silently guessing one if a supply is nonetheless requested).
"""
supplytype(::NicheAxis) = nothing
supplytype(::Precipitation) = WaterSupply
supplytype(::SolarRadiation) = SolarSupply

"""
    demandtype(::NicheAxis)

The `AbstractDemand` concrete type a species uses to consume this axis' resource, or
`nothing` if the axis is not a resource.
"""
demandtype(::NicheAxis) = nothing
demandtype(::Precipitation) = WaterDemand
demandtype(::SolarRadiation) = SolarDemand

# A genuine, dispatchable dimension (like the built-in `Unitful.Power`) for a volumetric flow
# rate (L³/T) — Unitful has no named alias for this, unlike energy/time → power.
Unitful.@derived_dimension VolumeFlow (Unitful.𝐋^3 * Unitful.𝐓^-1)

# `supplytype`/`demandtype` overloaded on a *quantity type* rather than a `NicheAxis`: given a
# bare Resource-role value (no axis attached — e.g. a user-supplied `supply =`/`resource =`
# keyword), pick the concrete Supply/Demand type from its physical dimension via ordinary
# dispatch — no `Dict`, no `if`/`elseif` chain. `Unitful.Power` covers `kJ/day` (energy per
# time); `VolumeFlow` (above) covers `L/day` (volume per time); anything real and
# dimensionless is the free/simple case. The fallback method is the "no match" case,
# resolved the same way `Base.show`'s fallback is — the least specific method in the table,
# not a separate test-then-error branch.
supplytype(::Type{<:Unitful.Power}) = SolarSupply
supplytype(::Type{<:VolumeFlow}) = WaterSupply
supplytype(::Type{<:Real}) = SimpleSupply
function supplytype(::Type{T}) where {T}
    return error("cannot pick a supply type for values of type $T — expected a Power (kJ/day), a VolumeFlow (L/day), or a bare Real (dimensionless).")
end

demandtype(::Type{<:Unitful.Power}) = SolarDemand
demandtype(::Type{<:VolumeFlow}) = WaterDemand
demandtype(::Type{<:Real}) = SimpleDemand
function demandtype(::Type{T}) where {T}
    return error("cannot pick a demand type for values of type $T — expected a Power (kJ/day), a VolumeFlow (L/day), or a bare Real (dimensionless).")
end

"""
    _basedimension(::Type)

The physical "substance" dimension of a *stored* Resource element type, with its rate
component removed — e.g. `dimension(kJ)` for `typeof(1.0*kJ/day)`, `dimension(L)` for
`typeof(1.0*L/day)`, `NoDims` for the free/simple case's bare `typeof(1.0/day)`
(`Unitful.Frequency`, the built-in T⁻¹ alias — the stored type already carries the rate, not
a bare `Real`, unlike `supplytype`/`demandtype`'s constructor-input dispatch above). Used
only by the rate-dimension guard test (`test/test_Energy.jl`): every concrete Demand/Supply
element type must satisfy `dimension(T) / _basedimension(T) == 𝐓^-1` exactly — a genuine
rate, not merely "any negative time exponent" (energy already embeds 𝐓^-2 in its own
definition, so a bare exponent check would be wrong; verified directly: `dimension(kJ/d) /
dimension(kJ) == 𝐓^-1`).
"""
_basedimension(::Type{<:Unitful.Power}) = dimension(kJ)
_basedimension(::Type{<:VolumeFlow}) = dimension(Unitful.L)
_basedimension(::Type{<:Unitful.Frequency}) = NoDims

"""
    dynamics(::NicheAxis)

The default per-timestep change function for a layer on this axis (default
[`NoChange`](@ref)); overridden per layer at build time.
"""
dynamics(::NicheAxis) = NoChange
