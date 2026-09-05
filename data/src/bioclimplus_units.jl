# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Establish the physical units of two "accumulated" CHELSA-BIOCLIM+ layers by downloading them and
# comparing magnitudes against well-defined references over contrasting climates:
#
#   * gsp ("Growing season precipitation") - a seasonal TOTAL (kg*m^-2 = mm) or a per-day rate
#     (kg*m^-2*day^-1)?  Compared to bio12 (annual precipitation, mm) and gsl (growing-season length, days).
#   * gdd0/gdd5/gdd10 ("Growing degree days above 0/5/10 C") - labelled `K`, but a degree-day SUM is `K*day`.
#     Compared to gst (growing-season mean temperature, K) and gsl.
#
# It also shows that the reader now masks CHELSA's raw integer fill sentinels (typemax of the band, which the
# file's declared nodata misses) - before the fix gsp's `0xffffffff` fill leaked as ~4.29e8.
#
# **Both questions are SETTLED, and this script is the evidence behind the answers rather than an
# open investigation.** `data/RasterDataSources/BioClimPlus.csv` records them:
#
#   * `gsp` ships as `L*m^-2` with `AccumulationPeriod = percell=gsl` - a seasonal total, as Part A
# concludes. Confirmed empirically 2026-07-23.
#   * `gdd0`/`gdd5`/`gdd10` ship as `d*K` over `year`, which is Part B's recommendation **already
#     applied** - they were `K` when this was written.
#
# So a re-run is a *regression check* on the catalogue, not a decision: both verdicts below should
# still come out the way the shipped table says. One that does not is a finding.
#
# It downloads seven global CHELSA layers, each ~100-600 MB, cached after the first run.
#
# This is a manual diagnostic, NOT part of the test suite (its name does not start with `test_`, so
# runtests.jl does not auto-run it).
#
#     julia --project=data/src data/src/bioclimplus_units.jl

using EcoSISTEM
using RasterDataSources
using Unitful, Unitful.DefaultSymbols
using EcoSISTEM.Units
using DimensionalData: DimensionalData, X, Y
using Statistics
using Printf

const SCALE = 20            # mean-aggregate the global 30-arcsec layers (magnitude-preserving)

# (name, (lat_lo, lat_hi), (long_lo, long_hi)) - small ° boxes spanning contrasting growing seasons
const REGIONS = [
    ("Amazon (tropical wet, gsl≈365)", (-4.0, 0.0), (-66.0, -62.0)),
    ("England (temperate)", (51.0, 53.0), (-2.0, 0.0)),
    ("C. Siberia (boreal, short season)", (60.0, 63.0), (90.0, 93.0))
]

# Read one layer once, globally, mean-aggregated (downloads via getraster on first call, then cached).
readlayer(T, layer) = read(T, layer, scale = SCALE, fn = mean).array

# Plain finite Float values (strip any unit / Missing / NaN).
function finitevals(a)
    return filter(isfinite,
                  [Float64(ustrip(float(x)))
                   for x in vec(Array(a)) if !ismissing(x)])
end

# Mean of finite values inside a lat/long box.
function boxmean(arr, (latlo, lathi), (longlo, longhi))
    # The lookup's own values, asked for by dimension name rather than by axis position.
    coords(d) = [Float64(ustrip(float(x)))
                 for x in parent(DimensionalData.lookup(arr, d))]
    lat, long = coords(Y), coords(X)
    sub = Array(arr)[(lat .>= latlo) .& (lat .<= lathi),
                     (long .>= longlo) .& (long .<= longhi)]
    vals = filter(isfinite,
                  [Float64(ustrip(float(x))) for x in vec(sub) if !ismissing(x)])
    return isempty(vals) ? NaN : mean(vals)
end

@info "downloading/reading CHELSA layers (cached after first run; ~106 MiB each) ..."
gsp = readlayer(CHELSA{BioClimPlus}, :gsp)
gsl = readlayer(CHELSA{BioClimPlus}, :gsl)
bio12 = readlayer(CHELSA{BioClim}, 12)
gst = readlayer(CHELSA{BioClimPlus}, :gst)                       # growing-season mean temp (K)
gdd = Dict(x => readlayer(CHELSA{BioClimPlus}, Symbol("gdd$x"))
           for x in (0, 5, 10))

println("\n================ FILL-SENTINEL CHECK (reader masks raw integer typemax) ================")
gspvals = finitevals(gsp)
@printf("gsp global (scale=%d): max = %.1f mm, median = %.1f mm\n",
        SCALE, maximum(gspvals), median(gspvals))
println("-> physical (hundreds-thousands of mm). Before the reader masked CHELSA's raw `0xffffffff` fill,")
println("  the unmasked UInt32 sentinel scaled to 4.29e8 (= (2^32-1)*0.1) and polluted the data.")

println("\n================ PART A - gsp: seasonal TOTAL vs per-day RATE ================")
for (nm, latb, longb) in REGIONS
    A = boxmean(gsp, latb, longb)
    G = boxmean(gsl, latb, longb)
    P = boxmean(bio12, latb, longb)
    @printf("%-34s gsp=%.0f  gsl=%.0f d  bio12=%.0f mm | TOTAL: gsp/gsl=%.2f mm/d, gsp/bio12=%.3f | RATE: gsp*gsl=%.0f\n",
            nm, A, G, P, A / G, A / P, A * G)
end
gsp_ratios = [boxmean(gsp, r[2], r[3]) / boxmean(bio12, r[2], r[3])
              for r in REGIONS]
println(all(x -> 0.1 <= x <= 1.2, filter(isfinite, gsp_ratios)) ?
        "VERDICT: gsp is a SEASONAL TOTAL -> `kg*m^-2` (= mm) is correct. ✅ This is what BioClimPlus.csv " *
        "already ships (`L*m^-2`, via the water-density identity, with `percell=gsl`)." :
        all(x -> x < 0.05, filter(isfinite, gsp_ratios)) ?
        "🔴 VERDICT CHANGED: gsp is a PER-DAY RATE -> `kg*m^-2*day^-1`. BioClimPlus.csv ships it as a total." :
        "VERDICT: inconclusive - inspect the numbers.")

println("\n================ PART B - gdd: `K` vs `K*day` ================")
# A degree-day SUM over the season is `K*day` (= °C*day). The clincher: gdd_x ≈ (gst - x)*gsl - the mean daily
# excess above x°C times the number of days - so gdd_x ≫ any single temperature. If instead gdd_x were O(10s)
# it would be a genuine `K`.
# NB the reader returns CHELSA BioClimPlus temperatures in °C-magnitude (it applies no °C->K shift for this
# source: `_layerunit` is `NoUnits`), so `gst` prints as e.g. ~26 for the Amazon. That °C-vs-`K` labelling of
# the CHELSA temperature layers is a SEPARATE reader/CSV issue, out of scope here.
for (nm, latb, longb) in REGIONS
    G = boxmean(gsl, latb, longb)
    Tgst = boxmean(gst, latb, longb)                   # growing-season mean temp (°C-magnitude - see note)
    println(nm)
    @printf("  gsl=%.0f d   gst=%.1f (°C)\n", G, Tgst)
    for x in (0, 5, 10)
        D = boxmean(gdd[x], latb, longb)
        @printf("  gdd%-2d = %.0f    gdd%d/gsl = %.2f /day   |  (gst-%d)*gsl = %.0f  (≈ gdd%d if `K*day`)\n",
                x, D, x, D / G, x, max(Tgst - x, 0.0) * G, x)
    end
end
gdd0_over_gst = [boxmean(gdd[0], r[2], r[3]) / boxmean(gst, r[2], r[3])
                 for r in REGIONS]
println(all(x -> x > 3, filter(isfinite, gdd0_over_gst)) ?
        "VERDICT: gdd0/5/10 are ACCUMULATED degree-day sums (gdd_x ≈ (gst-x)*gsl, gdd0 ≫ gst) -> `K*day` " *
        "(°C*day), NOT `K`. ✅ This is what BioClimPlus.csv already ships (`d*K` over `year`); gst is a mean and stays `K`." :
        "🔴 VERDICT CHANGED: gdd values are O(gst) -> genuine `K`. BioClimPlus.csv ships `d*K`, so one of the two is now wrong.")
