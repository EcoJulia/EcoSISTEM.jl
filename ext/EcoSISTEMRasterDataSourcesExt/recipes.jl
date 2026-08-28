# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Plot recipes for monthly climate rasters -------------------------------------------------------
#
# Here rather than in `EcoSISTEMPlotsExt` because they dispatch on `WorldClim{Climate}`/
# `CHELSA{Climate}`, and an extension can only see its own trigger packages and the parent's
# dependencies, so a Plots extension cannot name a dataset at all. A `@recipe` needs only
# `RecipesBase`, which the parent depends on, so defining them here costs nothing: without a plotting
# backend loaded a recipe is simply inert.

@recipe function f(wc::ClimateRaster{<:Union{WorldClim{Climate},
                                             CHELSA{Climate}}},
                   time::Union{Integer, Unitful.Time})
    ind = _monthindex(time)
    mnth = Dates.monthabbr(ind)
    A = transpose(ustrip.(wc.array[Ti(At(ind * month_mean_duration))]))
    x = ustrip.(parent(DimensionalData.lookup(wc.array, 1)))
    y = ustrip.(parent(DimensionalData.lookup(wc.array, 2)))
    seriestype := :heatmap
    grid --> false
    background_color_inside --> :grey
    title --> "$mnth"
    return x, y, A
end

@recipe function f(wc::ClimateRaster{<:Union{WorldClim{Climate},
                                             CHELSA{Climate}}},
                   time::Union{Integer, Unitful.Time}, xrange, yrange)
    ind = _monthindex(time)
    mnth = Dates.monthabbr(ind)
    A = transpose(ustrip.(wc.array[Y(xrange), X(yrange),
                                   Ti(At(ind * month_mean_duration))]))
    step1 = ustrip(parent(DimensionalData.lookup(wc.array, 1))[2] -
                   parent(DimensionalData.lookup(wc.array, 1))[1])
    step2 = ustrip(parent(DimensionalData.lookup(wc.array, 2))[2] -
                   parent(DimensionalData.lookup(wc.array, 2))[1])
    x = ustrip(xrange.left):step1:ustrip(xrange.right)
    y = ustrip(yrange.left):step2:ustrip(yrange.right)
    seriestype := :heatmap
    grid --> false
    background_color_inside --> :grey
    title --> "$mnth"
    return x, y, A
end
