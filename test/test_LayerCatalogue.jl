# SPDX-License-Identifier: LGPL-3.0-or-later

module TestLayerCatalogue

using EcoSISTEM
# Abstract types are `public`, not exported, so `using EcoSISTEM` does not bring them in.
using EcoSISTEM: NicheAxis, TemperatureAxis, WaterAxis, PrecipitationAxis
using EcoSISTEM.Units
using RasterDataSources
using Rasters
using Unitful
using Unitful.DefaultSymbols
using Test

@testset "layerunit" begin
    @testset "temperature and precipitation (BioClim)" begin
        @test layerunit(WorldClim{BioClim}, 1) == °C        # annual mean temperature
        @test layerunit(WorldClim{BioClim}, 2) == K            # mean diurnal range
        # `layerunit` is the **table's** unit - the accumulated amount, with the period split out
        # into `AccumulationPeriod`. What a *read* yields is the rate, and that is `SourceSpec.unit`;
        # the two are asserted side by side below so the distinction cannot quietly collapse again.
        @test layerunit(WorldClim{BioClim}, 12) == Unitful.L * m^-2          # annual precipitation
        # A dimensionless index is not necessarily *unscaled*: isothermality is published as a
        # ratio ×100, so its unit is `percent` and `_tocanon` divides it out to a 0-1 fraction when
        # a layer is built. Dimension `NoDims` either way - only the scale differs.
        @test layerunit(WorldClim{BioClim}, 3) == Unitful.percent      # isothermality
        @test dimension(layerunit(WorldClim{BioClim}, 3)) == NoDims
    end

    @testset "compound units and Symbol/String keys (Climate)" begin
        @test layerunit(WorldClim{Climate}, :srad) == kJ * m^-2
        @test layerunit(WorldClim{Climate}, "wind") == m * s^-1
        @test dimension(layerunit(WorldClim{Climate}, :vapr)) ==
              dimension(u"kPa")
    end

    @testset "land cover is a dimensionless percentage" begin
        # Published per-class as 0-100, which `Units = percent` records; the canonical form is the
        # **fraction**, so a built layer is 0-1 and no `%` unit travels downstream.
        @test layerunit(EarthEnv{LandCover}, 1) == Unitful.percent
        @test layerunit(EarthEnv{LandCover}, 12) == Unitful.percent
        @test dimension(layerunit(EarthEnv{LandCover}, 1)) == NoDims
    end

    @testset "other tables resolve (BOM, month, elevation)" begin
        @test layerunit(WorldClim{Elevation}, :elev) == m
        @test layerunit(WorldClim{BioClimPlus}, :cmi_max) == Unitful.L * m^-2
    end

    @testset "_datasettype extracts the wrapped dataset" begin
        @test EcoSISTEM._datasettype(WorldClim{BioClim}) === BioClim
        @test EcoSISTEM._datasettype(EarthEnv{LandCover}) === LandCover
        @test EcoSISTEM._datasettype(CHELSA{BioClim}) === BioClim
    end

    @testset "dataset type resolves through any wrapper" begin
        # CHELSA{BioClim} shares the BioClim layer table
        @test layerunit(CHELSA{BioClim}, 1) == °C
    end

    @testset "code aliases + reconciled codes resolve" begin
        # BioClim layers take both an int and a `bioN` symbol (`;`-separated `Code` aliases)
        @test layerunit(WorldClim{BioClim}, :bio1) ==
              layerunit(WorldClim{BioClim}, 1) == °C
        @test layerunit(WorldClim{BioClim}, :bio12) == Unitful.L * m^-2
        # LandCover takes an int or a descriptive symbol
        @test layerunit(EarthEnv{LandCover}, :open_water) ==
              layerunit(EarthEnv{LandCover}, 12)
        # BioClimPlus documents bio1-19
        @test layerunit(CHELSA{BioClimPlus}, :bio1) == °C
        # CHELSA{Climate} monthly layers live in Climate.csv; water-mass rates canonicalise to
        # L m^-2 month^-1 (litres, not kg, per the water-density identity: 1 kg water ≈ 1 L)
        @test layerunit(CHELSA{Climate}, :cmi) == Unitful.L * m^-2
        @test layerunit(CHELSA{Climate}, :pr) == Unitful.L * m^-2
        @test layerunit(CHELSA{Climate}, :tas) == °C
    end

    @testset "unknown layer / source errors" begin
        @test_throws ErrorException layerunit(WorldClim{BioClim}, 99)
        @test_throws ErrorException layerunit(RasterDataSources.AWAP, 1)
    end
end

# The shipped `ValueType` column is what says whether a layer may be *averaged*, and the split
# that matters is categorical against the other two - not the everyday sense of the words. A
# `discrete` layer here is a day-of-year or a count of days: an ordinary number. Only `categorical`
# marks a class code, whose mean is meaningless, and there are exactly six of those. Until now the
# column was shipped but parsed away entirely - nothing in `src/` read it - which is why
# `_resamplemethod` still interpolates everything.
@testset "ValueType reaches the catalogue" begin
    CP = EcoSISTEM
    vt(code, dataset) = only(filter(r -> r.dataset === dataset,
                                    CP.layerinfo(code))).valuetype

    # The only categorical layers in the shipped tables: the climate typologies.
    @test vt("kg0", :BioClimPlus) === :categorical
    @test vt("kg5", :BioClimPlus) === :categorical
    # `discrete` is numeric and averageable - the case the old `_isdiscrete` name got backwards.
    @test vt("gsl", :BioClimPlus) === :discrete
    @test vt("scd", :BioClimPlus) === :discrete
    @test vt("bio1", :BioClim) === :continuous

    # Exactly six categorical layers, all of them `kg*`; everything else is numeric.
    cats = filter(r -> r.valuetype === :categorical, CP._catalogue())
    @test length(cats) == 6
    @test all(r -> all(startswith("kg"), r.aliases), cats)
    @test all(r -> r.valuetype in (:categorical, :discrete, :continuous),
              CP._catalogue())

    # An unrecognised value must fail loudly rather than default to "safe to interpolate", which
    # would silently mangle a class-code layer.
    @test CP._parsevaluetype("Categorical") === :categorical   # case-insensitive
    @test CP._parsevaluetype("") === :continuous               # blank = not stated
    @test_throws ErrorException CP._parsevaluetype("ordinal")
end

# `Category` is the only machine-readable record of what *kind* of quantity a layer is, and
# nothing in `src/` reads it yet - which is exactly why it is worth pinning here. A miscategorised
# layer costs nothing today and stops being free the moment anything dispatches on it (§4 of the
# layer-units plan makes `Category` decide supply-eligibility).
@testset "Category says what kind of quantity a layer is" begin
    CP = EcoSISTEM
    cat(code, dataset) = only(filter(r -> r.dataset === dataset,
                                     CP.layerinfo(code))).category

    # A climate moisture index is precipitation *minus* PET, so it can be negative - a `balance`,
    # never a `rate`. The distinction is what stops a sign-indefinite layer being read as a supply:
    # a negative water supply is meaningless.
    @test cat("cmi", :Climate) === :balance
    @test cat("cmi_max", :BioClimPlus) === :balance
    @test cat("cmi_mean", :BioClimPlus) === :balance
    @test cat("cmi_min", :BioClimPlus) === :balance

    # ...but a *range* of a balance is a spread, not a balance: `max - min` cannot be negative. Every
    # `*Range`-axis layer is a `range`, whatever the quantity it is a range of - a spread is not a
    # flow, so none of them is supply-eligible either.
    @test cat("cmi_range", :BioClimPlus) === :range
    @test cat("pet_penman_range", :BioClimPlus) === :range
    @test cat("rsds_range", :BioClimPlus) === :range
    spread = filter(r -> !isnothing(r.axis) &&
                         endswith(string(r.axis), "Range"), CP._catalogue())
    @test !isempty(spread)
    @test all(r -> r.category === :range, spread)

    # `DayCount` counts days (`gsl`, `scd`, `ngd*`) rather than measuring a spread between two
    # values, so `count` is honest for it and the sweep above must not catch it.
    #
    # `gsl` is catalogued `range`, not `stock`, and the reason is worth stating because it accrues
    # qualifying days and so invites the wrong one. So does `ngd0`, which is a `count`; and `gsl` is
    # `lgd - fgd`, a *span* between two `DayOfYear` positions rather than an accumulation.
    # Decisively, `gsl` is the layer other layers measure their
    # accumulation period *by* (`gsp` is precipitation per growing season), and an interval cannot
    # itself be measured over an interval. It was named
    # `DayRange` until 2026-08-04, which put it in that sweep by its name alone while meaning
    # something else entirely - so the name is pinned here too, and renaming it back would fail this
    # test loudly rather than silently start demanding `:range` of a day count.
    @test !endswith(string(DayCount), "Range")
    @test all(r -> r.category === :count,
              filter(r -> r.axis === DayCount, CP._catalogue()))
    @test cat("gsl", :BioClimPlus) === :count

    @test all(r -> r.category in CP._CATEGORIES, CP._catalogue())
    # `Category` is cross-checked against the **accumulation period**, not against a separate
    # rate flag: such a flag would say exactly `category === :rate` in every row and so could never
    # disagree. See the `AccumulationPeriod` testset below.
end

# The accumulation period: the interval a value accumulated over, which is what turns a total into
# an honest rate. It is a *different question* from `temporal` (how often the layer is sampled), and
# they must not share a column - that is invisible for `prec`, where both are months, and wrong for
# `srad`, which is sampled monthly but accumulates per day.
@testset "AccumulationPeriod" begin
    CP = EcoSISTEM
    rec(code, dataset) = only(filter(r -> r.dataset === dataset &&
                                          first(r.aliases) == code,
                                     CP._catalogue()))

    # A constant period is a real, parseable unit.
    @test rec("gdd0", :BioClimPlus).period ==
          CP.ConstantAccumulationPeriod(year)
    @test rec("bio13", :BioClimPlus).period ==
          CP.ConstantAccumulationPeriod(month_mean_duration)
    @test rec("bio16", :BioClimPlus).period ==
          CP.ConstantAccumulationPeriod(quarter_mean_duration)

    # The four rows a single unit cannot describe: each slice accumulated over its own real month.
    for c in ("prec", "pr", "cmi", "pet")
        @test rec(c, :Climate).period isa CP.PerSliceAccumulationPeriod
    end
    # ...and the one whose period is another *layer*, varying by cell.
    @test rec("gsp", :BioClimPlus).period == CP.PerCellAccumulationPeriod(:gsl)

    # `srad`/`rsds` are the rows that prove sampling and accumulation are different questions:
    # twelve monthly slices, but the value is already a per-*day* flux, so no month divides it.
    @test rec("srad", :Climate).period == CP.ConstantAccumulationPeriod(day)
    @test rec("srad", :Climate).temporal == month_mean_duration        # sampled monthly
    @test rec("srad", :Climate).numslices == 12

    # `Temporal Resolution` now means sampling only, so a single-slice layer has none at all - even
    # though `bio13` is emphatically a monthly quantity. That is the whole point of the split.
    @test isnothing(rec("bio13", :BioClimPlus).temporal)
    @test rec("bio13", :BioClimPlus).numslices == 1

    # The cross-column invariant: a period is declared exactly where something accumulates. This is
    # what turns ~100 blank cells from "not filled in" into "asserted not to apply" - and it is what
    # exposed `gsl` being catalogued a `stock` when it is a span between two dates.
    # One assertion, not 139: a per-row `@test` in a loop would bury the suite count in near
    # duplicates. Collecting the offenders first keeps a failure just as informative - it prints the
    # codes - while staying a single test.
    # `:range` is excluded because it is neither required nor forbidden: a range **inherits**
    # whether it has a period from the quantity it ranges over, so no per-row rule can decide it.
    # That case is checked against its siblings instead - see the range testset below.
    offenders = [first(r.aliases)
                 for r in CP._catalogue()
                 if r.category !== :range &&
        (r.category in (:rate, :stock, :balance)) !=
        !isnothing(r.period)]
    @test isempty(offenders)
    @test isnothing(rec("gsl", :BioClimPlus).period)     # a count, so no period
    @test isnothing(rec("tmin", :Climate).period)        # instantaneous, but still sampled monthly
    @test rec("tmin", :Climate).temporal == month_mean_duration

    # The grammar itself. `=` is a safe discriminator because a unit can never contain one.
    @test isnothing(CP._parseperiod(""))
    @test CP._parseperiod("year") == CP.ConstantAccumulationPeriod(year)
    @test CP._parseperiod("perslice=calendar_month") isa
          CP.PerSliceAccumulationPeriod
    @test CP._parseperiod("percell=gsl") == CP.PerCellAccumulationPeriod(:gsl)
    @test CP._parseperiod("percell=7") == CP.PerCellAccumulationPeriod(7)
    # A period must be an interval, and an unknown form must fail loudly rather than default.
    @test_throws ErrorException CP._parseperiod("m")                    # a length, not a time
    @test_throws ErrorException CP._parseperiod("perslice=calendar_week")
    @test_throws ErrorException CP._parseperiod("persplice=calendar_month")
    @test_throws ErrorException CP._parseperiod("percell=")
    # Both directions of the invariant are enforced, not just the missing one.
    @test_throws ErrorException CP._checkperiod(nothing, :rate, "x")
    @test_throws ErrorException CP._checkperiod(CP.ConstantAccumulationPeriod(year),
                                                :instantaneous, "x")

    # `IsRate` is gone: it said exactly `category === :rate` in every row, and two columns holding
    # one fact is how they drift apart.
    @test !(:IsRate in CP._REQUIRED_COLUMNS)

    # Each kind renders as what it means, not as a struct.
    detail(r) = sprint(show, MIME"text/plain"(), r)
    @test occursin("accumulated: over 1 d", detail(rec("srad", :Climate)))
    @test occursin("own calendar month", detail(rec("prec", :Climate)))
    @test occursin("`gsl` layer, which varies by cell",
                   detail(rec("gsp", :BioClimPlus)))
    @test occursin("sampled    : month_mean_duration (12 slices)",
                   detail(rec("prec", :Climate)))
end

# A range is the difference of two values of what it ranges over, and subtraction preserves both
# unit and accumulation period. So a `*Range` row must carry `absoluteunit` of its siblings' unit and
# exactly their period - no judgement about what a range "should" mean, just that subtraction is
# unit-preserving. The catalogue's own definitions say so outright ("Difference between maximum and
# minimum monthly X").
#
# This is the rule three rows silently broke. `cmi_range`, `pet_penman_range` and `rsds_range` kept
# a rate's unit with a blank period, having fallen between two steps: the categories were corrected
# `rate` -> `range` while their `Units` cells were deliberately left for later, and the later pass
# stripped periods from "every row *with* a period" - which, by then, they no longer had.
@testset "a range carries its siblings' unit, made absolute, and their period" begin
    CP = EcoSISTEM
    rec(code, dataset) = only(filter(r -> r.dataset === dataset &&
                                          first(r.aliases) == code,
                                     CP._catalogue()))
    axisname(r) = isnothing(r.axis) ? "" : String(nameof(r.axis))
    cat = CP._catalogue()
    ranges = [r for r in cat if endswith(axisname(r), "Range")]
    # The shipped tables really do have range rows to check; an empty sweep would pass vacuously.
    @test length(ranges) == 11

    # One assertion over every range row, for the same reason the invariant above is one: a loop of
    # `@test`s would inflate the count with near-duplicates. Offenders are named on failure.
    offenders = String[]
    for r in ranges
        base = axisname(r)[1:(end - length("Range"))]
        sibs = [s for s in cat if s.dataset == r.dataset && axisname(s) == base]
        isempty(sibs) && continue
        units = unique(s.unit for s in sibs)
        length(units) == 1 || continue
        (r.unit == Unitful.absoluteunit(only(units)) &&
         CP._periodsagree(r.period, first(sibs).period)) ||
            push!(offenders, first(r.aliases))
    end
    @test isempty(offenders)

    # The three that were wrong, pinned individually so a regression names itself.
    @test rec("cmi_range", :BioClimPlus).unit == u"L*m^-2"
    @test rec("cmi_range", :BioClimPlus).period ==
          CP.ConstantAccumulationPeriod(month_mean_duration)
    @test rec("pet_penman_range", :BioClimPlus).unit == u"L*m^-2"
    @test rec("rsds_range", :BioClimPlus).unit == u"MJ*m^-2"
    @test rec("rsds_range", :BioClimPlus).period ==
          CP.ConstantAccumulationPeriod(day)
    # The affine case, and the reason the rule says *absolute*: a difference of temperatures is an
    # interval, so `bio2` is K where the `Temperature` it ranges over is °C.
    @test rec("bio2", :BioClimPlus).unit == u"K"
    @test rec("bio1", :BioClimPlus).unit == u"°C"
    @test Unitful.absoluteunit(u"°C") == u"K"

    # Periods compare by what they mean, so a `Constant` is never quietly equal to a `PerSlice`.
    @test CP._periodsagree(nothing, nothing)
    @test CP._periodsagree(CP.ConstantAccumulationPeriod(day),
                           CP.ConstantAccumulationPeriod(day))
    @test !CP._periodsagree(CP.ConstantAccumulationPeriod(day), nothing)
    @test !CP._periodsagree(CP.ConstantAccumulationPeriod(month_mean_duration),
                            CP.PerSliceAccumulationPeriod())
    @test !CP._periodsagree(CP.PerCellAccumulationPeriod(:gsl),
                            CP.PerCellAccumulationPeriod(:scd))

    # And the check bites. Rebuilding one record with one field changed is the whole negative test:
    # if `_checkrangerows` accepted these, the shipped tables passing it would mean nothing.
    with(r, field, value) = CP.LayerRecord((f === field ? value :
                                            getfield(r, f)
                                            for f in fieldnames(CP.LayerRecord))...)
    rng = rec("cmi_range", :BioClimPlus)
    fam = [s
           for s in cat
           if s.dataset == :BioClimPlus && axisname(s) == "ClimateMoisture"]
    # The unit this row actually had was `L*m^-2*month^-1`, which can no longer even be *written*:
    # `Units.month` is gone, and `u"...month..."` fails to parse. `d^-1` reproduces the shape of the bug -
    # one spurious time denominator - in a unit that still exists.
    @test_throws ErrorException CP._checkrangerows([fam;
                                                    with(rng, :unit,
                                                         u"L*m^-2*d^-1")])
    @test_throws ErrorException CP._checkrangerows([fam;
                                                    with(rng, :period, nothing)])
    @test_throws ErrorException CP._checkrangerows([fam;
                                                    with(rng, :period,
                                                         CP.PerSliceAccumulationPeriod())])
    # The unmodified family passes, so the failures above are the perturbation and not the setup.
    @test isnothing(CP._checkrangerows([fam; rng]))
end

# Option C, stated directly: the two questions and their two answers. `layerunit` says what the
# shipped table declares; `SourceSpec.unit` says what materialising it actually yields. Before the
# accumulation period was split out, one field meant both - which is how a monthly total came to be
# divided by a fixed 30.4375-day month for every month of the year.
@testset "layerunit is the table's amount; a read yields the rate" begin
    CP = EcoSISTEM

    # A layer whose canonical reading is a rate: the table holds the amount, the read divides.
    @test layerunit(WorldClim{Climate}, :prec) == Unitful.L * m^-2
    @test SourceSpec(WorldClim{Climate}, :prec).unit ==
          Unitful.L * m^-2 * day^-1
    @test layerunit(WorldClim{BioClim}, :bio13) == Unitful.L * m^-2
    @test SourceSpec(WorldClim{BioClim}, :bio13).unit ==
          Unitful.L * m^-2 * day^-1

    # Having a period does NOT mean a layer is read as a rate - the **axis** decides that.
    # `CumulativeHeat`'s canonical unit is the heat sum itself, and §4b established heat is a
    # condition rather than a consumable resource, so `gdd0` is never divided into a temperature.
    @test layerunit(CHELSA{BioClimPlus}, :gdd0) == day * K
    @test SourceSpec(CHELSA{BioClimPlus}, :gdd0).unit == day * K

    # Nor is a per-cell period divided here: `gsp`'s stock-and-rate readings are role-dependent
    # and the role is unknown when a spec is built, so it keeps the declared amount.
    @test SourceSpec(CHELSA{BioClimPlus}, :gsp).unit == Unitful.L * m^-2

    # A **resource-only** axis has no `Condition` unit to consult, so its `Resource`-role one
    # stands in: `CarbonFlux` declares `g/day`, so `npp`'s annual total is read as a daily rate
    # like any other rate layer. Without the fallback it would stay `g*m^-2` and no carbon supply
    # could be built from it at all (`cancel` would be handed a mass per area, not a mass flux).
    @test layerunit(CHELSA{BioClimPlus}, :npp) == g * m^-2
    @test SourceSpec(CHELSA{BioClimPlus}, :npp).unit == g * m^-2 * day^-1

    # A layer with no period at all is untouched by any of this.
    @test SourceSpec(WorldClim{BioClim}, :bio1).unit ==
          layerunit(WorldClim{BioClim}, :bio1)

    # The derivation itself, independent of any table row.
    @test CP.layerrate(Unitful.L * m^-2, nothing, Precipitation) ==
          Unitful.L * m^-2
    @test CP.layerrate(Unitful.L * m^-2,
                       CP.PerSliceAccumulationPeriod(), Precipitation) ==
          Unitful.L * m^-2 * day^-1
    @test CP.layerrate(Unitful.L * m^-2,
                       CP.PerCellAccumulationPeriod(:gsl), Precipitation) ==
          Unitful.L * m^-2
end

# **The divisors are the whole point of the subproject**, so they are pinned here directly rather
# than only through a read: dividing a monthly total by a fixed 30.4375-day month is 7.7% wrong for
# February, and these numbers are what make it right. No data is downloaded - the divisors come from
# the catalogue and the calendar, not from any raster.
@testset "accumulation-period divisors" begin
    CP = EcoSISTEM
    rec(ds, c) = only(filter(r -> r.dataset === ds && string(c) in r.aliases,
                             CP._catalogue()))

    # Per-slice: one divisor per month, in days, February at its 28.25-day mean.
    d = CP._readdivisors(rec(:Climate, "prec"), 1:12)
    @test d ≈ [31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    @test sum(d) ≈ 365.25                     # the year is conserved exactly

    # A **partial** read must use the months actually asked for, not 1:n. This is the case the
    # `Ti` axis cannot answer (it would say 1,2,3 for a `month = 2:4` read), and getting it wrong
    # would divide February's data by January's length - silently, and 9% out.
    @test CP._readdivisors(rec(:Climate, "prec"), 2:4) ≈ [28.25, 31, 30]
    @test CP._readdivisors(rec(:Climate, "prec"), [2]) ≈ [28.25]

    # A constant period is one scalar for the whole layer.
    @test CP._readdivisors(rec(:BioClimPlus, "bio13"), 1:1) ≈ 30.4375
    @test CP._readdivisors(rec(:BioClimPlus, "bio16"), 1:1) ≈ 91.3125
    @test CP._readdivisors(rec(:BioClim, "12"), 1:1) ≈ 365.25
    @test CP._readdivisors(rec(:BioClimPlus, "npp"), 1:1) ≈ 365.25

    # And the two kinds that must NOT be divided: a heat sum (the axis says the accumulated
    # amount is the canonical reading) and a per-cell period (role-dependent, resolved later).
    @test isnothing(CP._readdivisors(rec(:BioClimPlus, "gdd0"), 1:1))
    @test isnothing(CP._readdivisors(rec(:BioClimPlus, "gsp"), 1:1))
    @test isnothing(CP._readdivisors(rec(:BioClim, "1"), 1:1))   # no period at all
end

@testset "layeraxis / _resolveaxis" begin
    # populated Axis columns resolve to the right axis (incl. the hierarchy leaves)
    @test layeraxis(WorldClim{BioClim}, 1) === Temperature
    @test layeraxis(WorldClim{BioClim}, 2) === TemperatureRange
    @test layeraxis(WorldClim{BioClim}, 12) === Precipitation
    @test layeraxis(WorldClim{Climate}, :srad) === SolarRadiation
    @test layeraxis(EarthEnv{LandCover}, 1) === SurfaceArea
    @test layeraxis(WorldClim{Elevation}, :elev) === Altitude
    @test layeraxis(EarthEnv{HabitatHeterogeneity}, :shannon) === Heterogeneity
    # the max/min/mean collapse: quarter temps all fold into Temperature
    @test layeraxis(CHELSA{BioClimPlus}, :bio8) === Temperature
    # the aggregate families: level vs range; and the new axes
    @test layeraxis(CHELSA{BioClimPlus}, :cmi_max) === ClimateMoisture
    @test layeraxis(CHELSA{BioClimPlus}, :cmi_range) === ClimateMoistureRange
    @test layeraxis(CHELSA{Climate}, :vpd) === VapourPressureDeficit
    @test layeraxis(CHELSA{Climate}, :hurs) === RelativeHumidity
    @test layeraxis(CHELSA{Climate}, :pet) === Evapotranspiration
    @test layeraxis(CHELSA{Climate}, :clt) === CloudCover
    @test layeraxis(CHELSA{BioClimPlus}, :gdd0) === CumulativeHeat
    @test layeraxis(CHELSA{BioClimPlus}, :gst) === Temperature
    @test layeraxis(CHELSA{BioClimPlus}, :swe) === SnowWaterEquivalent
    @test layeraxis(CHELSA{BioClimPlus}, :gsp) === GrowingSeasonPrecipitation
    @test layeraxis(CHELSA{BioClimPlus}, :fgd) === DayOfYear
    @test layeraxis(CHELSA{BioClimPlus}, :gsl) === DayCount
    @test layeraxis(CHELSA{BioClimPlus}, :npp) === CarbonFlux
    @test layeraxis(CHELSA{BioClimPlus}, :fcf) === FrostChangeFrequency
    @test layeraxis(CHELSA{BioClimPlus}, :kg0) === ClimateTypology
    # a name resolves to its concrete NicheAxis by autodiscovery (through abstract groups)
    @test EcoSISTEM._resolveaxis("Temperature") === Temperature
    @test EcoSISTEM._resolveaxis("Heterogeneity") === Heterogeneity
    @test EcoSISTEM._resolveaxis("Altitude") === Altitude
    # an abstract grouping supertype is NOT a leaf, so it does not resolve as a layer's axis
    @test_throws ErrorException EcoSISTEM._resolveaxis("TemperatureAxis")
    @test_throws ErrorException EcoSISTEM._resolveaxis("WaterAxis")
    @test_throws ErrorException EcoSISTEM._resolveaxis("DayAxis")
    # a name that isn't a loaded NicheAxis errors clearly
    @test_throws ErrorException EcoSISTEM._resolveaxis("NotAnAxis")
end

@testset "guard: shipped data/ tables are well-formed" begin
    datadir = pkgdir(EcoSISTEM, "data", "RasterDataSources")
    csvs = filter(endswith(".csv"), readdir(datadir))
    @test !isempty(csvs)
    for f in csvs
        base = first(splitext(f))
        # every table's basename names a real RasterDataSources dataset type
        @test isdefined(RasterDataSources, Symbol(base))
        rows = EcoSISTEM._layertable(joinpath(datadir, f))
        for (_, cell) in rows
            # every non-blank Units parses ...
            isempty(cell.units) ||
                @test uparse(cell.units, unit_context = [Unitful, Units]) isa
                      Unitful.Units
            # ... and every non-blank Axis resolves to a loaded NicheAxis
            isempty(cell.axis) ||
                @test EcoSISTEM._resolveaxis(cell.axis) <: NicheAxis
        end
    end
end

@testset "guard: tables reconcile to RasterDataSources layers(T)" begin
    RDS = RasterDataSources
    # each shipped table ↔ the RDS source type(s) that map to it (via `_datasettype` -> CSV name)
    tablesources = Dict("BioClim" => Any[WorldClim{BioClim}, CHELSA{BioClim}],
                        "BioClimPlus" => Any[CHELSA{BioClimPlus}],
                        "Climate" => Any[WorldClim{Climate}, CHELSA{Climate}],
                        "Elevation" => Any[WorldClim{Elevation}],
                        "HabitatHeterogeneity" =>
                            Any[EarthEnv{HabitatHeterogeneity}],
                        "LandCover" => Any[EarthEnv{LandCover}])
    # `ncdf` is listed by `layers(CHELSA{Climate})` but is a spurious RDS entry (no such CHELSA
    # variable / file), so it is deliberately absent from the table - carve it out here.
    known_spurious = Set(["ncdf"])
    datadir = pkgdir(EcoSISTEM, "data", "RasterDataSources")
    # every code any source accepts: `layers(T)` plus its `layerkeys(T)` aliases
    accepted(T) = union(Set(string(c) for c in RDS.layers(T)),
                        Set(string(c)
                            for c in (try
                                          RDS.layerkeys(T)
                                      catch
                                          RDS.layers(T)
                                      end)))
    for (base, Ts) in tablesources
        acc = setdiff(union((accepted(T) for T in Ts)...), known_spurious)
        csvcodes = Set(keys(EcoSISTEM._layertable(joinpath(datadir,
                                                           "$base.csv"))))
        # exactly reconciled: every fetchable code is documented, and no phantom rows exist
        @test acc == csvcodes
    end
end

@testset "catalogue query helpers (layerinfo / layersbyaxis / layeraxes)" begin
    CP = EcoSISTEM
    # layerinfo(T, code): one detailed record
    r = CP.layerinfo(WorldClim{BioClim}, 1)
    @test r isa CP.LayerRecord
    @test r.name == "Annual Mean Temperature"
    @test r.unit == °C
    @test r.axis === Temperature
    @test "bio1" in r.aliases && "1" in r.aliases
    @test r.sources == ["WorldClim", "CHELSA"]     # full list of supporting sources
    @test_throws ErrorException CP.layerinfo(WorldClim{BioClim}, 99)

    # layerinfo(code): every dataset that has the code
    @test Set(x.dataset for x in CP.layerinfo(:bio1)) ==
          Set([:BioClim, :BioClimPlus])
    @test only(CP.layerinfo(:tas)).sources == ["CHELSA"]   # a CHELSA-only Climate code

    # temporal metadata surfaced from the NumSlices / Temporal Resolution columns
    @test CP.layerinfo(WorldClim{Climate}, :tmin).temporal ==
          month_mean_duration
    @test CP.layerinfo(WorldClim{Climate}, :tmin).numslices == 12
    @test CP.layerinfo(WorldClim{BioClim}, 1).numslices == 1       # a static, single-slice layer
    @test isnothing(CP.layerinfo(WorldClim{BioClim}, 1).temporal)

    # layersbyaxis: concrete leaf, and abstract group spanning several axes
    @test Set(x.dataset for x in CP.layersbyaxis(Temperature)) ⊇
          Set([:BioClim, :Climate])
    temp = CP.layersbyaxis(TemperatureAxis)
    @test all(x.axis <: TemperatureAxis for x in temp)
    @test length(temp) > length(CP.layersbyaxis(Temperature))
    # the water umbrella spans precipitation + humidity + fluxes/stocks
    water = CP.layersbyaxis(WaterAxis)
    @test all(x.axis <: WaterAxis for x in water)
    @test Set(x.axis for x in water) ⊇
          Set([Precipitation, RelativeHumidity, Evapotranspiration,
                  ClimateMoisture])

    # The no-argument form is the whole catalogue, and the reason it exists is the trap on the next
    # line: `layersbyaxis(NicheAxis)` spans every *axis*, but a layer with a blank `Axis` cell has none,
    # so it looks complete and silently is not. Today nothing is unclassified, which is precisely what
    # makes that easy to miss - so the relationship is asserted rather than the counts.
    every = CP.layersbyaxis()
    @test length(every) ==
          length(CP.layersbyaxis(NicheAxis)) + length(CP.layersbyaxis(nothing))
    @test all(isnothing(x.axis) for x in CP.layersbyaxis(nothing))
    @test all(!isnothing(x.axis) for x in CP.layersbyaxis(NicheAxis))
    # No shipped layer is unclassified today - pinned so that shipping one is a deliberate act that
    # shows up here, rather than a quiet gap in every axis-based sweep.
    @test isempty(CP.layersbyaxis(nothing))
    @test length(every) == 139
    # A copy, not the cached vector: sorting the result in place must not corrupt later lookups.
    sort!(every, by = x -> first(x.aliases))
    @test length(CP.layersbyaxis()) == 139

    # layeraxes: nested tree of NicheAxis subtypes down to the concrete leaves, each leaf
    # carrying the names of the shipped layers that use it
    tree = CP.layeraxes()
    @test tree.axis === NicheAxis
    @test isempty(tree.names)                       # NicheAxis itself is never a layer's axis
    water = only(c for c in tree.children if c.axis === WaterAxis)
    precip = only(c for c in water.children if c.axis === PrecipitationAxis)
    precipleaf = only(c for c in precip.children if c.axis === Precipitation)
    @test isempty(precip.names)                     # abstract group: no layers of its own
    @test "Annual Mean Temperature" ∉ precipleaf.names # sanity: wrong branch has no bio1
    @test "Annual Precipitation" in precipleaf.names
    @test isempty(precipleaf.children)               # a concrete leaf has no subtypes

    # scoping to a subtree with the optional argument
    tempaxis = CP.layeraxes(TemperatureAxis)
    @test tempaxis.axis === TemperatureAxis
    @test any(c -> c.axis === Temperature, tempaxis.children)
    @test "Annual Mean Temperature" in
          only(c for c in tempaxis.children if c.axis === Temperature).names

    # public but NOT exported
    @test :layerinfo in names(EcoSISTEM, all = false)
    @test !Base.isexported(EcoSISTEM, :layerinfo)
end

# The check that needs no declaration, and so covers precisely the axes `_checkaxisunit` gives up
# on. Tested by *doctoring* a copy of the real catalogue rather than by hand-building records: the
# point is that it runs over the shipped rows, and a synthetic pair would not prove that.
@testset "cross-row dimension homogeneity" begin
    CP = EcoSISTEM
    catalogue = CP._catalogue()

    # The shipped tables must already satisfy it - this is what the load-time call asserts.
    @test isnothing(CP._checkaxishomogeneity(catalogue))

    # Rebuild one record on a well-populated axis with a dimensionally wrong unit, and check it is
    # caught and that the message names the axis and both layers.
    victim = first(CP.layersbyaxis(Temperature))
    fields = [getfield(victim, f) for f in fieldnames(CP.LayerRecord)]
    fields[findfirst(==(:unit), collect(fieldnames(CP.LayerRecord)))] = u"kg"
    broken = CP.LayerRecord(fields...)
    err = try
        CP._checkaxishomogeneity(vcat(catalogue, [broken]))
        nothing
    catch e
        sprint(showerror, e)
    end
    @test !isnothing(err)
    @test occursin("Temperature", err)
    @test occursin(first(victim.aliases), err)

    # Differing *scale* on one axis is legal - it is what `canonicalunit` exists to reconcile -
    # so `SolarRadiation`'s kJ*m^-2 and MJ*m^-2 rows must NOT be reported.
    @test isnothing(CP._checkaxishomogeneity(CP.layersbyaxis(SolarRadiation)))
end

# Step 12: `pet` and `cmi` are read as flows, so declaring the rate must move the **unit and the
# values together** - `_readdivisors` asks `layerrate` whether the unit changed, so a declaration that
# changed only one of them would be a silent scaling error.
@testset "pet and cmi read as flows; swb and gsp do not" begin
    CP = EcoSISTEM
    rate_axes = (Evapotranspiration, EvapotranspirationRange, ClimateMoisture,
                 ClimateMoistureRange)
    for A in rate_axes, r in CP.layersbyaxis(A)
        @test dimension(CP.layerrate(r.unit, r.period, r.axis)) ==
              dimension(u"mm" / day)
        @test !isnothing(CP._readdivisors(r, 1:12))     # values divided too
    end
    # **`Category` does not decide this, and these rows prove it**: `cmi_mean` is a `balance` and
    # *does* divide, while `swb` is also a `balance` and does not; `pet_penman_mean` is a `rate` and
    # divides, while `gdd0` is a `stock` and does not. Two `balance` rows with opposite answers in one
    # testset is the whole content of the discriminator table in `docs/src/units.md`.
    #
    # The three that must NOT divide, each for a different reason: `gdd0` is a heat **sum** (the
    # accumulated total *is* the quantity), `swb` is a *capped* cumulative (the mean is not a rate of
    # anything), and `gsp`'s period varies per cell so there is no fixed divisor.
    for A in (CumulativeHeat, SiteWaterBalance, GrowingSeasonPrecipitation),
        r in CP.layersbyaxis(A)
        @test CP.layerrate(r.unit, r.period, r.axis) === r.unit
        @test isnothing(CP._readdivisors(r, 1:12))
    end
end

# Axis-name resolution is memoised, and the cache must not change what resolves.
#
# Why it exists: `_resolveaxis` runs once per catalogue row, and walking `_leafaxes()` each time
# costs `subtypes()` all the way down - ~0.6 s a call, 83 s a build. Every axis being abstract is
# what makes the walk expensive, since it cannot stop at concrete leaves. Nothing fails without the
# cache, so the only symptom is slower gates, which is exactly why it needs a test of its own.
@testset "axis names resolve through the cache" begin
    CP = EcoSISTEM

    @test CP._resolveaxis("Temperature") === Temperature
    @test CP._resolveaxis("Temperature") === Temperature      # ...and again, from the cache
    @test CP._resolveaxis("Precipitation") === Precipitation

    # An unknown name still errors, and says how many it found.
    @test_throws ErrorException CP._resolveaxis("NoSuchAxisAnywhere")

    # The path that makes the cache safe for downstream packages: an axis declared *after* the map
    # was built is still found, because a miss re-walks before giving up.
    @nicheaxis(CachedLateAxis <: EcoSISTEM.NicheAxis, condition = K)
    @test CP._resolveaxis("CachedLateAxis") === CachedLateAxis
end

# **A failure in this testset is GOOD NEWS**, and the only one in the suite of which that is true.
#
# Two providers publish values that are not in the unit they document, and nothing in the file says
# so - the discrepancy is in the provider's own metadata. `PublishedScaleFactor` in the shipped
# tables records how many times too large the published values are, so that a **fix upstream** is
# noticed instead of silently turning any compensation we apply into a second error the other way.
#
# So these assertions deliberately pin the **broken** state. If one fails, the provider has
# probably corrected its data: clear that layer's `PublishedScaleFactor` cell, remove any
# compensation keyed on it, and re-bless anything downstream - do **not** adjust the threshold here.
#
# The catalogue half needs no download and runs everywhere. Only the one *read* costs anything,
# and it is deliberately the cheapest anomaly available (EarthEnv heterogeneity at 25 km, ~0.6 MB)
# rather than the most interesting - the runner's time budget is the binding constraint, and the
# runtime check in `read` covers every other dataset for anyone who downloads one.
@testset "known upstream scaling defects are still present" begin
    CP = EcoSISTEM
    @testset "recorded in the catalogue (no download)" begin
        # WorldClim documents bio4 as "standard deviation ×100"; CHELSA gives its unit as `°C/100`.
        @test CP.layerinfo(WorldClim{BioClim}, 4).publishedscale == 100
        # EarthEnv documents `Scale = 0.0001` for every heterogeneity measure and applies none.
        @test CP.layerinfo(EarthEnv{HabitatHeterogeneity}, :evenness).publishedscale ==
              10_000
        @test CP.layerinfo(EarthEnv{HabitatHeterogeneity}, :contrast).publishedscale ==
              10_000
        # The column must stay **empty** everywhere else: a filled cell is a claim that the
        # provider is wrong, and one made by accident would silently rescale real data.
        @test isnothing(CP.layerinfo(WorldClim{BioClim}, 1).publishedscale)
        @test isnothing(CP.layerinfo(EarthEnv{LandCover}, 1).publishedscale)
        @test isnothing(CP.layerinfo(CHELSA{BioClim}, 1).publishedscale)
    end

    @testset "the published file still needs the factor (one small read)" begin
        # **The pin, taken from the FILE rather than from our own arithmetic.** Dividing a read by
        # 10 000 and then multiplying back would only test our own multiplication; this reads what
        # EarthEnv actually published. `evenness` is documented `0 to 1` and the file holds 0-9937,
        # because nothing applies the `Scale = 0.0001` they document - not EarthEnv, not GDAL (their
        # GeoTIFFs carry no scale tag, unlike CHELSA's), and until now not us.
        path = RasterDataSources.getraster(EarthEnv{HabitatHeterogeneity},
                                           :evenness)
        rawmax = maximum(skipmissing(Raster(path isa AbstractString ? path :
                                            first(path))))
        # The same **geometric midpoint** the runtime check uses, `√factor × ceiling` = 100 here.
        # That tolerates the data moving by up to 100× - any ordinary revision - while failing on a
        # rescale. If this fails, EarthEnv has started applying its own scale: clear that layer's
        # `PublishedScaleFactor` cell, and the read-time division disappears with it.
        @test rawmax > sqrt(10_000) * 1.0

        # ...and with the factor applied, a read now lands inside the documented range.
        vals = filter(isfinite,
                      vec(read(EarthEnv{HabitatHeterogeneity}, :evenness).array))
        @test maximum(vals) <= 1.0
        @test maximum(vals) > 0.5          # corrected, not flattened

        # The alarm itself is exercised, not merely trusted: values that look *already fixed*
        # must make the runtime check speak up.
        @test_logs (:warn, r"no longer looks") CP._checkupstreamscale(EarthEnv{HabitatHeterogeneity},
                                                                      :evenness,
                                                                      vals)
    end
end

end
