# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Which niche axes exist - the shipped catalogue, grouped by family, each declared with
# [`@nicheaxis`](@ref). Every one is abstract, so a leaf can be subdivided later without changing
# the hierarchy.
#
# Needs `Ecology.jl` and `nicheaxis_macro.jl` before it: the macro emits methods naming
# `Resource`, and the interface hooks they define must already exist.

using Unitful
using Unitful.DefaultSymbols

# --- Temperature axes ------------------------------------------------------
# **`TemperatureAxis` is a pure grouping node and declares no unit** - it groups *things to do with
# temperature*, not *things measured in K*. Its six leaves do not share a unit: three are in kelvin,
# `CumulativeHeat` is in `K*d` and two are dimensionless.
#
# That is a rule rather than a tidiness preference. Promotion can move a layer from a group onto a
# leaf, and it **rewraps without converting**, so a group declaring a unit its leaves override would
# let kelvin data be relabelled as degree-days. There is no meaningful conversion between those: they
# are different *quantities*, not one quantity in two units. The macro refuses such a declaration, so
# the trap is unrepresentable rather than merely remembered.
#
# Level, range and seasonality follow the aggregation policy: max, min and mean collapse to one level
# axis, a range gets its own axis, and seasonality - a standard deviation - stays apart.
"""    TemperatureAxis <: NicheAxis - temperature-related niche axes. A pure grouping node: each leaf states its own unit, because they do not share one. """
@nicheaxis(TemperatureAxis<:NicheAxis)
"""    Temperature <: TemperatureAxis - an absolute temperature (annual/monthly/quarter mean, min, max, or growing-season mean; BioClim 1/5/6/8-11, Climate tavg/tmin/tmax/tas/tasmin/tasmax, BioClimPlus gst; K). """
@nicheaxis(Temperature<:TemperatureAxis, condition=K, bounds=(0.0K, nothing),
           densitywidth=1.0K)
"""    TemperatureRange <: TemperatureAxis - a temperature range (max - min): mean diurnal range or annual range (BioClim 2 and 7; K interval). """
@nicheaxis(TemperatureRange<:TemperatureAxis, condition=K,
           bounds=(0.0K, nothing), densitywidth=1.0K)
"""    TemperatureSeasonality <: TemperatureAxis - temperature seasonality, the standard deviation of the monthly mean temperatures (BioClim 4; K). """
@nicheaxis(TemperatureSeasonality<:TemperatureAxis, condition=K,
           bounds=(0.0K, nothing), densitywidth=1.0K)
"""    CumulativeHeat <: TemperatureAxis - accumulated growing degree-days above a threshold (BioClimPlus gdd0/gdd5/gdd10; K*day). """
@nicheaxis(CumulativeHeat<:TemperatureAxis, condition=K*Unitful.d,
           bounds=(0.0K*Unitful.d, nothing), densitywidth=1.0K*Unitful.d)
"""    Isothermality <: TemperatureAxis - isothermality, 100 * (mean diurnal range / annual range) (BioClim 3; dimensionless). """
@nicheaxis(Isothermality<:TemperatureAxis, condition=NoUnits,
           bounds=(0.0, nothing), densitywidth=1.0NoUnits)
"""    FrostChangeFrequency <: TemperatureAxis - frost change frequency, the number of days the temperature crosses 0 °C (BioClimPlus fcf; dimensionless). """
@nicheaxis(FrostChangeFrequency<:TemperatureAxis, condition=NoUnits,
           bounds=(0.0, nothing), densitywidth=1.0NoUnits)

# --- Water axes ------------------------------------------------------------
# `WaterAxis` is the umbrella for everything that measures water in any form. It is a pure grouping
# node - it carries NO unit or resource of its own, so it can span mm precipitation, kPa/dimensionless
# humidity and L*m^-2 fluxes/stocks without conflict; each leaf (or sub-group) declares its own.
"""    WaterAxis <: NicheAxis - umbrella for all water-measuring niche axes (precipitation, humidity, evapotranspiration, moisture, water stocks). A pure grouping node with no unit or resource of its own. """
@nicheaxis(WaterAxis <: NicheAxis)
"""    PrecipitationAxis <: WaterAxis - precipitation niche axes. A topical group; the unit and the water supply/demand live on the [`Precipitation`](@ref) leaf. """
@nicheaxis(PrecipitationAxis <: WaterAxis)
# A **rate**, not a depth: every shipped precipitation layer is read per unit time. The unit belongs on
# the *leaf* rather than on `WaterAxis` or `PrecipitationAxis`, because the water **stocks** -
# `SnowWaterEquivalent`, `SiteWaterBalance`, `GrowingSeasonPrecipitation` - genuinely are `L*m^-2`, a
# length, which is what makes a depth look plausible here.
#
# Both roles at once, which is the case that shows `resource` and `condition` are independent: `mm/day`
# as a condition, a climatological normal, and `L/day` as a resource, a pool species draw on.
"""    Precipitation <: PrecipitationAxis - a precipitation amount over a month, quarter or year (BioClim 12-19, Climate prec/pr; mm/day). Carries the unit and provides a water supply/demand. """
@nicheaxis(Precipitation<:PrecipitationAxis, condition=mm/Unitful.d,
           resource=Unitful.L/Unitful.d, supply=Supply{Precipitation},
           demand=Demand{Precipitation},
           bounds=(0.0mm/Unitful.d, nothing), densitywidth=1.0mm/Unitful.d)
"""    PrecipitationSeasonality <: PrecipitationAxis - precipitation seasonality, the coefficient of variation of the monthly precipitation totals (BioClim 15; dimensionless). """
@nicheaxis(PrecipitationSeasonality<:PrecipitationAxis, condition=NoUnits,
           densitywidth=1.0NoUnits)
"""    HumidityAxis <: WaterAxis - atmospheric-water niche axes (vapour pressure, its deficit, relative humidity). """
@nicheaxis(HumidityAxis <: WaterAxis)
"""    VapourPressure <: HumidityAxis - near-surface water-vapour partial pressure (Climate vapr; kPa). """
@nicheaxis(VapourPressure<:HumidityAxis, condition=kPa, densitywidth=1.0kPa)
"""    VapourPressureDeficitAxis <: HumidityAxis - vapour-pressure-deficit niche axes (level and range). """
@nicheaxis(VapourPressureDeficitAxis<:HumidityAxis, condition=Unitful.Pa,
           densitywidth=1.0Unitful.Pa)
"""    VapourPressureDeficit <: VapourPressureDeficitAxis - near-surface vapour pressure deficit (Climate vpd, BioClimPlus vpd_max/mean/min; Pa). """
@nicheaxis(VapourPressureDeficit <: VapourPressureDeficitAxis)
"""    VapourPressureDeficitRange <: VapourPressureDeficitAxis - annual range of the monthly vapour pressure deficit (BioClimPlus vpd_range; Pa). """
@nicheaxis(VapourPressureDeficitRange <: VapourPressureDeficitAxis)
# A **fraction**, 0 to 1. The shipped layers are published on the 0 to 100 scale, which their
# `Units = percent` cell records; `_tocanon` divides it out once at read, so nothing downstream
# carries a `%` unit around.
"""    RelativeHumidityAxis <: HumidityAxis - relative-humidity niche axes (level and range). """
@nicheaxis(RelativeHumidityAxis<:HumidityAxis, condition=NoUnits,
           bounds=(0.0, 1.0), densitywidth=1.0NoUnits)
"""    RelativeHumidity <: RelativeHumidityAxis - near-surface relative humidity (Climate hurs, BioClimPlus hurs_max/mean/min; dimensionless). """
@nicheaxis(RelativeHumidity <: RelativeHumidityAxis)
"""    RelativeHumidityRange <: RelativeHumidityAxis - annual range of the monthly relative humidity (BioClimPlus hurs_range; dimensionless). """
@nicheaxis(RelativeHumidityRange <: RelativeHumidityAxis)
# **Read as a flow.** Evaporative demand happens continuously, so the accumulated total divided by
# its period *is* a meaningful mean daily rate - which is the test for dividing, not the `Category`
# column (see `docs/src/units.md`). Declaring the rate is the whole of it: `layerrate` then reports
# `L*m^-2*day^-1`, and `_readdivisors` - which asks `layerrate` whether the unit changed - divides each
# slice by *its own* month. Unit and values move together, or not at all.
# Stated in `mm/day` rather than `L/(m^2*day)` to match [`Precipitation`](@ref): `cmi` is literally
# `prec - pet`, so the three should not print in two different spellings of the same unit. They are
# exactly equal (1 L*m^-2 ≡ 1 mm), so this is a display choice, not a conversion.
# Non-negative, and the floor lives on the group so the range inherits it: a max - min of a
# non-negative flow cannot be negative either. Unlike `ClimateMoistureAxis` below - no row on this
# axis is a `balance`.
"""    EvapotranspirationAxis <: WaterAxis - evapotranspiration niche axes (level and range). """
@nicheaxis(EvapotranspirationAxis<:WaterAxis, condition=mm/Unitful.d,
           bounds=(0.0mm/Unitful.d, nothing), densitywidth=1.0mm/Unitful.d)
"""    Evapotranspiration <: EvapotranspirationAxis - potential evapotranspiration (Climate pet, BioClimPlus pet_penman_max/mean/min; mm/day). """
@nicheaxis(Evapotranspiration <: EvapotranspirationAxis)
"""    EvapotranspirationRange <: EvapotranspirationAxis - annual range of the monthly potential evapotranspiration (BioClimPlus pet_penman_range; mm/day). """
@nicheaxis(EvapotranspirationRange <: EvapotranspirationAxis)
# **Also read as a flow, and it is a `balance`** - which is precisely why "is it a rate?" cannot be
# read off the `Category` column. `cmi` is `prec - pet`, both per month, so it is a *balance whose
# dimension is a rate*: a mean daily surplus or deficit is as real a quantity as a monthly one.
# Dividing it does not stop it being signed, and it is still never a supply.
# **No bounds, deliberately**: sign-indefinite by construction, so it must declare no floor - the
# same fact the catalogue records as `Category = balance`, and asserted against it in
# `test_nicheaxis_macro.jl`.
"""    ClimateMoistureAxis <: WaterAxis - climate-moisture-index niche axes (level and range). """
@nicheaxis(ClimateMoistureAxis<:WaterAxis, condition=mm/Unitful.d,
           densitywidth=1.0mm/Unitful.d)
"""    ClimateMoisture <: ClimateMoistureAxis - climate moisture index, precipitation minus potential evapotranspiration (Climate cmi, BioClimPlus cmi_max/mean/min; mm/day). """
@nicheaxis(ClimateMoisture <: ClimateMoistureAxis)
"""    ClimateMoistureRange <: ClimateMoistureAxis - annual range of the monthly climate moisture index (BioClimPlus cmi_range; mm/day). """
@nicheaxis(ClimateMoistureRange <: ClimateMoistureAxis)
"""    SnowWaterEquivalent <: WaterAxis - snow water equivalent, the liquid-water amount of the snowpack (BioClimPlus swe; L*m^-2). """
@nicheaxis(SnowWaterEquivalent<:WaterAxis, condition=Unitful.L/m^2,
           densitywidth=1.0Unitful.L/m^2)
# A **balance** (precipitation - PET, capped by soil capacity), so it is sign-indefinite and
# declares no floor - see the `bounds` docstring, and the testset asserting this against the
# catalogue's `Category = balance` marker.
"""    SiteWaterBalance <: WaterAxis - cumulative site water balance over the year, capped by soil water-holding capacity (BioClimPlus swb; L*m^-2). """
@nicheaxis(SiteWaterBalance<:WaterAxis, condition=Unitful.L/m^2,
           densitywidth=1.0Unitful.L/m^2)
# **A water *stock*, so an amount per area with no time in it** - unlike `Precipitation`. `gsp` is a
# total precisely because its period (`gsl`) varies per cell, so it has no fixed denominator to divide
# by. Declared because an axis with no canonical unit cannot carry a tolerance at all:
# `NicheTolerance` compares `dimension(support)` against `dimension(canonicalunit(axis))`.
"""    GrowingSeasonPrecipitation <: WaterAxis - precipitation accumulated over the growing season (BioClimPlus gsp; L*m^-2). """
@nicheaxis(GrowingSeasonPrecipitation<:WaterAxis, condition=Unitful.L/m^2,
           bounds=(0.0Unitful.L/m^2, nothing), densitywidth=1.0Unitful.L/m^2)

# --- Solar radiation / wind speed / cloud cover ----------------------------
# WorldClim's srad is kJ*m^-2*day^-1, CHELSA's rsds*/BioClimPlus's rsds_* are MJ*m^-2*day^-1 - without
# this the two sources' regimes are never reconciled to one scale (see `_tocanon`). kJ matches the
# Resource-role choice below and WorldClim's native scale. The range shares the level's physical unit
# (an annual range of the same flux), so it lives on the group.
"""    SolarRadiationAxis <: NicheAxis - solar-radiation niche axes (level and range). """
@nicheaxis(SolarRadiationAxis<:NicheAxis, condition=kJ/(m^2*Unitful.d),
           bounds=(0.0kJ/(m^2*Unitful.d), nothing),
           densitywidth=1.0kJ/(m^2*Unitful.d))
# Deliberately no `condition` of its own - it **inherits** the group's, which is the whole point of the
# group existing. Declaring `resource` says nothing about the condition reading.
"""    SolarRadiation <: SolarRadiationAxis - surface downwelling shortwave (solar) radiation flux (Climate srad/rsds, BioClimPlus rsds_max/mean/min). Provides a solar supply/demand. """
@nicheaxis(SolarRadiation<:SolarRadiationAxis, resource=kJ/Unitful.d,
           supply=Supply{SolarRadiation}, demand=Demand{SolarRadiation})
"""    SolarRadiationRange <: SolarRadiationAxis - annual range of the monthly solar radiation flux (BioClimPlus rsds_range). """
@nicheaxis(SolarRadiationRange <: SolarRadiationAxis)
"""    WindSpeedAxis <: NicheAxis - wind-speed niche axes (level and range). """
@nicheaxis(WindSpeedAxis<:NicheAxis, condition=m/Unitful.s,
           densitywidth=1.0m/Unitful.s)
"""    WindSpeed <: WindSpeedAxis - near-surface (10 m) wind speed (Climate wind/sfcWind, BioClimPlus sfcWind_max/mean/min; m*s^-1). """
@nicheaxis(WindSpeed <: WindSpeedAxis)
"""    WindSpeedRange <: WindSpeedAxis - annual range of the monthly near-surface wind speed (BioClimPlus sfcWind_range; m*s^-1). """
@nicheaxis(WindSpeedRange <: WindSpeedAxis)
# A **fraction**, 0 to 1. The shipped layers are published on the 0 to 100 scale, which their
# `Units = percent` cell records; `_tocanon` divides it out once at read, so nothing downstream
# carries a `%` unit around.
"""    CloudCoverAxis <: NicheAxis - cloud-cover niche axes (level and range). """
@nicheaxis(CloudCoverAxis<:NicheAxis, condition=NoUnits, bounds=(0.0, 1.0),
           densitywidth=1.0NoUnits)
"""    CloudCover <: CloudCoverAxis - total cloud cover (Climate clt, BioClimPlus clt_max/mean/min; dimensionless). """
@nicheaxis(CloudCover <: CloudCoverAxis)
"""    CloudCoverRange <: CloudCoverAxis - annual range of the monthly total cloud cover (BioClimPlus clt_range; dimensionless). """
@nicheaxis(CloudCoverRange <: CloudCoverAxis)

# --- Day-based, carbon and categorical axes --------------------------------
"""    DayAxis <: NicheAxis - day-based phenology niche axes (day-of-year positions and day counts; day). """
@nicheaxis(DayAxis<:NicheAxis, condition=Unitful.d,
           bounds=(0.0Unitful.d, nothing), densitywidth=1.0Unitful.d)
"""    DayOfYear <: DayAxis - an ordinal day-of-year: first/last growing day or first/last day above a degree-day threshold (BioClimPlus fgd/lgd/gdgfgd*/gddlgd*; day). """
@nicheaxis(DayOfYear <: DayAxis)
"""    DayCount <: DayAxis - a number of days: growing-season length, snow-cover days or days above a threshold (BioClimPlus gsl/scd/ngd*; day). Counts days rather than measuring a spread, so deliberately not named `*Range` as every genuine max-min axis is. """
@nicheaxis(DayCount <: DayAxis)
# **The bound is stated in the *Resource* unit, because this group has no Condition one** (see
# `CarbonFlux`). A supply cannot be negative by definition; declared here so the fact does not depend
# on which half of the interface is consulted, and so that giving the axis a Condition unit later
# cannot silently leave it unbounded. The shipped `npp` is *potential* NPP (Miami model),
# non-negative by construction - a net carbon flux, which genuinely goes both ways, would be a
# `Category = balance` layer and thus a different axis.
"""    CarbonAxis <: NicheAxis - carbon-cycle niche axes. """
@nicheaxis(CarbonAxis <: NicheAxis, bounds = (0.0g / Unitful.d, nothing))
# **`condition = nothing` is a decision, not an omission: `CarbonFlux` is supply-only.** `npp` is
# potential productivity - a resource species compete for, not a condition they are matched to a
# tolerance against. Written explicitly so that the axis-wide canonical-unit audit cannot "fix" it by
# reflex, and so it stays correct when the fallback becomes an error.
# Consequence worth knowing: `layerrate` falls back to the Resource unit when an axis declares no
# Condition one, which is what makes `npp` read as `g*m^-2*day^-1` rather than its annual total.
"""    CarbonFlux <: CarbonAxis - net primary productivity, a carbon flux (BioClimPlus npp; g*m^-2*day^-1). """
@nicheaxis(CarbonFlux<:CarbonAxis, condition=nothing, resource=g/Unitful.d,
           supply=Supply{CarbonFlux}, demand=Demand{CarbonFlux},
           densitywidth=1.0nothing)
# The one categorical group the package ships: its values are class codes, so a resampler must take
# the nearest class rather than interpolate. `categorical = true` is inherited by every descendant,
# so `LandCoverTypology` and `ClimateTypology` below declare nothing and answer `true` too.
"""    TypologyAxis <: NicheAxis - categorical classification niche axes (class labels; dimensionless, used with a `CategoricalRegime`). """
@nicheaxis(TypologyAxis<:NicheAxis, condition=NoUnits, densitywidth=1.0NoUnits,
           categorical=true)
"""    LandCoverTypology <: TypologyAxis - categorical land-cover / land-use classes (EarthEnv LandCover 1-12). """
@nicheaxis(LandCoverTypology <: TypologyAxis)
# **`SurfaceArea` is GONE; these layers are `SurfaceArea`** (see the space family below). Each
# EarthEnv band is the proportion of a cell covered by one class - an **area per unit area** - which
# is the same quantity a space supply measures, not a kind of class label. It never appeared in
# v0.4.0, so it is renamed outright with no shim. What it leaves behind is a tidier `TypologyAxis`,
# holding only the two axes that really are class labels.
"""    ClimateTypology <: TypologyAxis - categorical climate-classification classes (BioClimPlus kg0-5: Köppen-Geiger, Wissmann, Thornthwaite, Troll-Pfaffen). """
@nicheaxis(ClimateTypology <: TypologyAxis)

# --- Space ----------------------------------------------------------------
# **The most physically direct resource the model has.** Every other supply mediates crowding
# through an abstract pool; space *is* the rivalry - two plants cannot occupy the same square metre.
#
# **Areal, not per-cell**, exactly like solar and water: a layer holds a **fraction** (m^2/m^2, so
# dimensionless), and `cancel` multiplies by the cell's own area to give the per-cell supply in m^2.
# That is what makes "all ground" (1.0) and "available ground" (0.4) the *same kind of layer*,
# differing only in values - and it is why the twelve land-cover bands belong here.
#
# **"Areal" is invisible in the unit.** A cover fraction is `NoDims`, so no unit-based test can
# tell a space layer from any other dimensionless one. Only the axis knows - which is this
# subproject's whole thesis, and a mistake I made once while writing it.
#
# **One axis per stratum, not one axis with a stratum parameter**, because the strata do not share
# a dimension: surface and canopy are areas (m^2/m^2 -> m^2), soil is a volume (m^3/m^2 -> m^3). A single
# axis cannot carry two canonical units. `CanopyArea`/`SoilVolume` are not declared yet; the family
# exists so they can be added without disturbing anything.
"""    SpaceAxis <: NicheAxis - the physical space niche axes: ground surface now, canopy and soil later. """
@nicheaxis(SpaceAxis <: NicheAxis)
# `bounds = (0.0, 1.0)` is safe **only because axis bounds are consulted for the `Condition` role
# alone** - measured: `_enforcebounds!`'s `Resource` method clamps at zero and never reads them. So a
# derived available-ground *supply* that sums to more than 1 is never checked against this, while an
# individual land-cover band used as a *regime* maxes at exactly 1.0 and passes.
"""    SurfaceArea <: SpaceAxis - ground surface, as a fraction of a cell (EarthEnv LandCover 1-12, or a `SurfaceSpec`). """
@nicheaxis(SurfaceArea<:SpaceAxis, condition=NoUnits, resource=m^2,
           supply=Supply{SurfaceArea}, demand=Demand{SurfaceArea},
           bounds=(0.0, 1.0), densitywidth=1.0NoUnits)

# --- Other standalone axes -------------------------------------------------
"""    Heterogeneity <: NicheAxis - spatial habitat-heterogeneity metrics of EVI (EarthEnv HabitatHeterogeneity; dimensionless). """
@nicheaxis(Heterogeneity<:NicheAxis, condition=NoUnits, densitywidth=1.0NoUnits)
# Deliberately unbounded: below sea level in places. It is *not* a special case - bounds are
# declared, never derived, so there is nothing for it to be an exception to.
"""    Altitude <: NicheAxis - elevation above sea level (WorldClim Elevation; m). """
@nicheaxis(Altitude<:NicheAxis, condition=m, densitywidth=1.0m)

# == Functions ==================================================================================

# ---------------------------------------------------------------------------
# Layer roles
# ---------------------------------------------------------------------------
# A materialised layer plays one of two roles in a `GridHabitat`: a `Condition`, an environmental
# condition matched against species tolerances, or a `Resource`, something species draw on. The role
# is a phantom marker, which keeps the two type-distinguishable while they share one storage
# implementation.

# --- The NicheAxis interface -----------------------------------------------
# The small set of hooks that say what a layer on an axis means. Each has a safe default, and a group
# - through its abstract intermediate - or an individual axis overrides only what applies. The
# supply, demand and change types they name live in later-included files and are resolved at call
# time, so nothing here runs at load.

"""
    bounds(axis::Type{<:NicheAxis})

Return the physical range `axis`' values must stay within, as a `(floor, ceiling)` tuple in the
axis' **canonical unit**. Either end is `nothing` where no bound applies, and
`(nothing, nothing)` - the default - means the axis is unbounded in both directions.

**Stated in canonical units, and that is load-bearing**: a temperature is `≥ 0` in kelvin but not in
degrees Celsius, so a bound compared against a layer still holding °C would be nonsense. Every caller
converts first.

**The floor is the same fact the catalogue records as `Category = balance`.** An axis with a
`balance` row is sign-indefinite - a difference of two quantities, which may fall either side
of zero - and is *therefore* not supply-eligible either. The two are one fact, and a testset asserts that the
catalogue and these methods agree. Declared as methods rather than looked up per call, so that no
bound check costs a catalogue search.

Ceilings are declared here and **cannot** be read from the catalogue: its `OfficialUnit` column says
`°C` for `Isothermality`, `TemperatureSeasonality` and `PrecipitationSeasonality`, which are a ratio,
a standard deviation and a coefficient of variation respectively. That column is wrong for them, so
it is no evidence either way.
"""
bounds(::Type{<:NicheAxis}) = (nothing, nothing)

# **Deliberately not `public`.** `bounds` is exported by both `DimensionalData` and `Rasters`, so
# `using EcoSISTEM, Rasters` would be a collision - measured, not assumed. Declaring an axis's bounds
# goes through `@nicheaxis`, which is the public surface; reading them is internal.

"""
    supplytype(::Type{<:NicheAxis})

The `AbstractSupply` concrete type for this axis when used as a `Resource` resource, or
`nothing` if the axis is not a consumable resource (so it errors clearly rather than
silently guessing one if a supply is nonetheless requested).
"""
supplytype(::Type{<:NicheAxis}) = nothing

"""
    demandtype(::Type{<:NicheAxis})

The `AbstractDemand` concrete type a species uses to consume this axis' resource, or
`nothing` if the axis is not a resource.
"""
demandtype(::Type{<:NicheAxis}) = nothing

# The one wording for "this axis is not a consumable resource", raised from every route that can
# reach it: `_supplytype` and `_demandtype` when a supply or demand **type** is wanted, and
# `_canonicalresource` when a supply **value** is being converted.
#
# One function rather than one message per site, because a user who meets this has made a single
# mistake - declaring a condition as a supply - and should not be told about it in two voices
# depending on which conversion happened to run first. A claim that a wording is raised everywhere is
# worth checking against the call sites rather than assuming.
function _notaresource(axis, noun::AbstractString = "supply")
    return error("axis $axis is not a consumable resource, so it cannot be built as a " *
                 "$noun - it declares no resource unit. An axis becomes a resource by " *
                 "declaring one: `@nicheaxis($axis <: ..., resource = ..., supply = ..., " *
                 "demand = ...)`. (A condition species are *matched against* is used as " *
                 "`regime =`, not `supply =`.)")
end

# A genuine, dispatchable dimension (like the built-in `Unitful.Power`) for a volumetric flow
# rate (L^3/T) - Unitful has no named alias for this, unlike energy/time -> power.
Unitful.@derived_dimension VolumeFlow (Unitful.𝐋^3 * Unitful.𝐓^-1)

# **A layer's meaning comes from its axis, never from the dimension of its values**, and this is the
# defect the whole interface above exists to remove: `m/s` and `mm/day` are both `𝐋𝐓^-1`, so choosing a
# supply type by dimension would silently build a `WindSpeed` layer as a water supply - a wrong
# answer rather than an error. Both the supply and the demand path dispatch on the axis alone,
# `supplytype(axis)` and `demandtype(axis)`, which is the direction that cannot be wrong.
#
# The three shipped resources happen not to collide by dimension today, so nothing would fail
# immediately were this reversed. A fourth resource sharing a dimension with an existing one would,
# and silently - which is why the rule is stated rather than left to hold by luck.

# The physical "substance" dimension of a *stored* Resource element type, with its rate
# component removed - e.g. `dimension(kJ)` for `typeof(1.0*kJ/day)`, `dimension(L)` for
# `typeof(1.0*L/day)`, `dimension(g)` for `typeof(1.0*g/day)`, `NoDims` for the free/simple
# case's bare `typeof(1.0/day)`
# (`Unitful.Frequency`, the built-in T^-1 alias - the stored type already carries the rate, not
# a bare `Real`, unlike `supplytype`/`demandtype`'s constructor-input dispatch above). Used
# only by the rate-dimension guard test (`test/test_Demand.jl`): every concrete Demand/Supply
# element type must satisfy `dimension(T) / _basedimension(T) == 𝐓^-1` exactly - a genuine
# rate, not merely "any negative time exponent" (energy already embeds 𝐓^-2 in its own
# definition, so a bare exponent check would be wrong; verified directly: `dimension(kJ/d) /
# dimension(kJ) == 𝐓^-1`).
_basedimension(::Type{<:Unitful.Power}) = dimension(kJ)
_basedimension(::Type{<:VolumeFlow}) = dimension(Unitful.L)
_basedimension(::Type{<:Unitful.MassFlow}) = dimension(g)
_basedimension(::Type{<:Unitful.Frequency}) = NoDims
