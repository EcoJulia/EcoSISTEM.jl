# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Report every shipped layer `Code` (across all `data/RasterDataSources/*.csv` tables) with the one
# quantity the catalogue does **not** store: how much a cell actually receives per unit area if the
# timestep is one month. Writes `rate_table.csv` alongside this script.
#
# **This script reads the catalogue; it does not re-derive it.** Everything except
# `AmountPerAreaPerMonth` and its dimension is read straight off the `LayerRecord` - unit, axis,
# category, accumulation period, value type, official unit, notes. That is the whole design: an earlier
# version carried its own hand-written axis->category map, substance table and notes dictionary, and
# those drifted out of agreement with the shipped tables the moment the catalogue was corrected. A
# second source of truth for data that already has one is a bug waiting for someone to re-run it.
#
# It also used to **rewrite every shipped CSV in place** from those maps. That block was removed on
# 2026-08-04; re-running it would have silently reverted the catalogue corrections. This script is now
# read-only apart from its own `rate_table.csv`.
#
# ## The columns
#
#   - OfficialUnit: the unit as given by the layer's primary reference documentation (CHELSA V2.1 /
#     CHELSA-BIOCLIM+ technical specifications, or WorldClim's own docs). NOT necessarily the same as
#     the shipped `Units` cell, and NOT necessarily valid Unitful syntax - e.g. `kg*m^-2*gsl^-1` for
#     `gsp` and `g*C*m^-2*year^-1` for `npp`, since carbon is not a Unitful dimension and `gsl` is not
#     a time unit. Read from the catalogue's own `OfficialUnit` column.
# - Unit / UnitDimension: the shipped `Units` cell and its dimension. Since the accumulation
#     period was split into its own column, this is the accumulated **amount** (`L*m^-2`), not a rate
#     - the period is in AccumulationPeriod, and the rate is the two together.
#   - AccumulationPeriod: the interval that amount accumulated over, or blank where none applies.
# - AmountPerAreaPerMonth: the derived figure, and the only reason this script exists. For a
#     `:rate` it is `amount / its accumulation period × one month` - a genuine `uconvert`, so an annual
#     rate picks up a real ≈0.083 coefficient and a daily flux ≈30.44. Everything else repeats its own
#     unit unchanged. `uparse`-compatible: a leading numeric coefficient parses fine.
#   - Dimension: the dimension of AmountPerAreaPerMonth - for `bio19`, with the time dimension resolved
#     to a concrete month, just `L` (a depth), against UnitDimension's own `L`.
#
# ## Why `:stock` and `:balance` are excluded from the per-month conversion
#
# They are dimensionally accumulations too, but their reporting window is part of what they MEAN, not a
# resolution choice the way a month or year is for precipitation/PET/CMI/NPP/solar radiation (those are
# quantities you could equally well have measured at finer resolution, so "annual / 12" is a genuine if
# coarse monthly estimate). `:stock` - CumulativeHeat `gdd0/5/10`, GrowingSeasonPrecipitation `gsp` -
# is a flow integrated up to the point of reading, heat or water accumulated toward a phenological
# threshold or over a bounded season, analogous to `swe` (a literal physical stock of frozen water).
# Dividing by 12 would produce a number matching no real month, since the accumulation concentrates in
# particular parts of the year. `:balance` (SiteWaterBalance `swb`) is different again: a
# two-directional running balance, **capped by soil water-holding capacity** and able to go negative (a
# deficit), not a monotonic accumulation - so a mean rate would not be meaningful even if asked for.
# `:range` is excluded for the same reason a spread is not a flow.
#
# ## Two conventions the shipped tables already encode
#
# The **water-density identity** - 1 kg water / 1 m^2 ≈ 1 L / m^2 ≈ 1 mm, at ≈1000 kg/m^3 - is why every
# water-substance layer ships as a volume per area (`L*m^-2`) whatever its documentation says
# (`kg*m^-2` for CHELSA, `mm` for WorldClim). It is the same identity behind `Precipitation`'s own `mm`
# canonical unit (`src/NicheInfo.jl`). For evapotranspiration, the moisture index built from it, and
# snow water equivalent this is a *liquid-water-equivalent* depth, not the substance's own volume -
# which is exactly what "snow water equivalent" conventionally means. The conversion is already done in
# the shipped `Units` cells, so this script performs none; the tables' own `Notes` say where it applies.
#
# **`gsp` is treated as a plain total.** Its documentation gives `kg*m^-2*gsl^-1`, normalised by
# growing-season length, but `gsl` is not a fixed physical time unit - it varies by cell - so no
# conversion to it is attempted here. The catalogue records that honestly as
# `AccumulationPeriod = percell=gsl`, and deriving a rate from it needs the `gsl` layer read alongside.
#
# ## Running it
#
# A manual diagnostic, NOT part of the test suite - its name does not start with `test_`, so
# `runtests.jl` does not auto-run it, and nothing else runs `examples/scripts/` either.
#
#     julia --project=examples examples/scripts/rate_table.jl

using EcoSISTEM
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using CSV

const CP = EcoSISTEM

const SUPERSCRIPT = Dict('⁰' => '0', '¹' => '1', '²' => '2', '³' => '3',
                         '⁴' => '4', '⁵' => '5',
                         '⁶' => '6', '⁷' => '7', '⁸' => '8', '⁹' => '9',
                         '⁻' => '-')

# Convert Unitful's pretty-printed unit/dimension string (space-separated factors, unicode superscript
# exponents - e.g. "L m⁻² yr⁻¹") into ASCII syntax (e.g. "L*m^-2*yr^-1") parseable by `Unitful.uparse`
# - the same call `_catalogue()` uses to read the shipped `Units` column, so these cells are usable as
# direct replacements for it. Unitful's short forms (`d`, `yr`) are fine as they stand: `uparse`
# resolves them through the `Unitful` module in its `unit_context`, which is how the shipped
# `MJ*d^-1*m^-2` parses today.
function asciiform(s::AbstractString; sep = "*")
    parts = String[]
    for f in split(s)
        chars = collect(f)
        i = findfirst(c -> haskey(SUPERSCRIPT, c), chars)
        if isnothing(i)
            push!(parts, f)
        else
            push!(parts,
                  String(chars[1:(i - 1)]) * "^" *
                  String(map(c -> SUPERSCRIPT[c], chars[i:end])))
        end
    end
    return isempty(parts) ? "" : join(parts, sep)
end

# Format a Unitful units/dimensions object as `uparse`-compatible ASCII text, or "" for a dimensionless one.
unitstr(u) = u == Unitful.NoUnits ? "" : asciiform(string(u))
function dimstr(u)
    return (d = dimension(u);
            d == NoDims ? "" : asciiform(string(d), sep = " * "))
end

# The duration an accumulation period resolves to for a single-figure report. One method per kind.
_perioddur(p::CP.ConstantAccumulationPeriod) = p.duration
# A per-slice month is twelve different divisors; one figure has to pick one, and the mean month is
# the honest choice - it is the declared approximation `month_mean_duration` exists to be.
_perioddur(::CP.PerSliceAccumulationPeriod) = month_mean_duration
# `percell=gsl` varies by cell, so there is no figure to give without reading that layer.
_perioddur(::CP.PerCellAccumulationPeriod) = nothing
_perioddur(::Nothing) = nothing

# How the AccumulationPeriod cell reads in the report. Blank where none applies.
_periodstr(p::CP.ConstantAccumulationPeriod) = string(p.duration)
_periodstr(::CP.PerSliceAccumulationPeriod) = "perslice=calendar_month"
_periodstr(p::CP.PerCellAccumulationPeriod) = "percell=" * string(p.code)
_periodstr(::Nothing) = ""

# The derived figure. Only a `:rate` converts; everything else repeats its own unit, since there is
# no meaningful area/timestep decomposition (see the header for why stocks and balances are excluded).
function amount_per_month(r, code)
    # A non-rate returns its *unit*, not `1.0 × unit`: there is no conversion to report, and a
    # spurious coefficient would imply one happened. The return type is what the formatter dispatches
    # on, so the two cases cannot be confused.
    r.category === :rate || return r.unit
    # The one genuine invariant left in this script. It is now *also* enforced at catalogue parse
    # time by `_checkperiod`, so this should be unreachable - kept because a report that silently
    # invents a period would be worse than one that stops.
    isnothing(r.period) &&
        error("`$code` is Category = :rate but declares no AccumulationPeriod; a rate without an " *
              "interval cannot be converted to a per-month amount.")
    d = _perioddur(r.period)
    isnothing(d) && return nothing        # per-cell: not derivable without reading that layer
    return uconvert(r.unit,
                    (1.0 * r.unit) / (1 * d) * (1 * month_mean_duration))
end

# `nothing` (a per-cell period) reports as blank rather than as a wrong number; a bare unit is a
# non-rate repeating itself; a quantity is a genuine conversion and carries its coefficient.
_amountstr(::Nothing) = ""
_amountstr(u::Unitful.Units) = unitstr(u)
function _amountstr(q::Unitful.Quantity)
    return "$(round(ustrip(q), sigdigits = 4))*" *
           unitstr(unit(q))
end
_dimof(::Nothing) = ""
_dimof(u::Unitful.Units) = dimstr(u)
_dimof(q::Unitful.Quantity) = dimstr(unit(q))

# The public route to the whole catalogue. **Not `layersbyaxis(EcoSISTEM.NicheAxis)`**, which looks complete
# and is not: it spans every *axis*, but a layer with a blank `Axis` cell has none and would be silently
# omitted. No shipped layer is unclassified today, which is exactly what makes that trap easy to fall
# into. The no-argument form means the whole catalogue by construction.
rows = NamedTuple[]
for r in CP.layersbyaxis()
    code = join(r.aliases, ";")
    q = amount_per_month(r, code)
    push!(rows,
          (Dataset = string(r.dataset), Code = code,
           Axis = isnothing(r.axis) ? "" : string(nameof(r.axis)),
           OfficialUnit = r.officialunit, Category = string(r.category),
           Unit = unitstr(r.unit), UnitDimension = dimstr(r.unit),
           AccumulationPeriod = _periodstr(r.period),
           AmountPerAreaPerMonth = _amountstr(q), Dimension = _dimof(q),
           ValueType = string(r.valuetype), Notes = r.notes))
end

CSV.write(joinpath(@__DIR__, "rate_table.csv"), rows)
println("wrote ", length(rows), " rows to ",
        joinpath(@__DIR__, "rate_table.csv"))
println("\n=== RATE rows only (the ones that actually convert) ===")
for r in rows
    r.Category == "rate" &&
        println(rpad(r.Dataset, 14), rpad(r.Code, 18), rpad(r.Axis, 24),
                rpad(r.Unit, 14), rpad(r.AccumulationPeriod, 22), " -> ",
                rpad(r.AmountPerAreaPerMonth, 22), r.Dimension)
end
