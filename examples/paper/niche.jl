# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Figure 3, `NicheWidth2.pdf`: the analytic curves behind the demographics - how a species'
# resource demand scales its rates, and how its niche width sets its match to the environment. No
# simulation, so the figure is identical at both scales.
#
# Included by `examples/paper.jl`, never run directly.

using Distributions: Normal
using Plots

# The suitability the demographics see for species of optimum `opt` and the given niche widths
# across `temperatures`: the package's own tolerance and nichefit, so the normalisation is exactly
# the model's rather than a hand-written Gaussian. One column per width.
function suitability_curves(opt, widths, temperatures)
    tolerance = NicheTolerance(Temperature, Normal, fill(opt, length(widths)),
                               widths)
    fit = NicheSuitability(tolerance)
    return [fit(tolerance.dists[i], T)
            for T in temperatures, i in eachindex(widths)]
end

function niche_figure()
    # A: the demand term. A species' birth and death rates are both scaled by
    # `(demand / mean demand)^-longevity`, so a species demanding less than the community average
    # turns over faster and one demanding more turns over more slowly, with the birth/death ratio
    # unchanged. Longevity is 1 here, as in the paper.
    ratio = range(0.01, 1.0, length = 200)
    panela = plot(ratio, ratio .^ -1.0,
                  xlab = "Resource demand, relative to the community mean",
                  ylab = "Multiplier on birth and death rates",
                  xlim = (0, 1), label = "", title = "A", titleloc = :left)

    # B: match to the environment across 20-30 C for niche widths of 1-5 K about an optimum of
    # 25 C. Light to dark is narrow to wide, on the colour bar.
    opt = uconvert(K, 25.0°C)
    # Collected first: a range of Celsius values is affine and cannot be indexed.
    temperatures = uconvert.(K, collect(20.0:0.05:30.0) .* °C)
    widths = (1.0:0.5:5.0) .* K
    curves = suitability_curves(opt, widths, temperatures)
    shifts = 0.5:0.5:2.5
    panelb = plot(ustrip.(°C, temperatures), curves,
                  line_z = permutedims(ustrip.(widths)), c = :blues,
                  clim = (1, 5), colorbar_title = "Niche width (C)",
                  colorbar_ticks = ([1, 5],
                                    ["1 (specialist)", "5 (generalist)"]),
                  xlab = "Temperature (C)", ylab = "Match to the environment",
                  label = "", title = "B", titleloc = :left)
    vline!(panelb, 25.0 .+ shifts, ls = :dash, c = :grey, label = "")

    # C: match against niche width, for an environment 0.5-2.5 C away from the optimum - the
    # dashed lines in B. The widths start at 0.1 K rather than 0, as figure 7 does.
    finewidths = range(0.1, 5.0, length = 200) .* K
    panelc = plot(xlab = "Niche width (C)", ylab = "Match to the environment",
                  legendtitle = "Shift in temperature (C)", title = "C",
                  titleloc = :left)
    for shift in shifts
        curve = vec(suitability_curves(opt, finewidths, [opt + shift * K]))
        plot!(panelc, ustrip.(finewidths), curve, label = string(shift))
    end

    # Recorded so the log carries the check: a 1 K width peaks at 0.399, a 5 K one near 0.08.
    println("    peak match at width 1 K: ",
            round(maximum(curves[:, 1]), digits = 3),
            ", at 5 K: ", round(maximum(curves[:, end]), digits = 3))

    return plot(panela, panelb, panelc, layout = (1, 3), size = (1800, 550),
                margin = 8Plots.mm, guidefontsize = 12, tickfontsize = 10,
                titlefontsize = 16)
end

timed("NicheWidth2") do
    return savefigure(niche_figure(), "NicheWidth2")
end
