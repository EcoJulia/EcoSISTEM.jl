# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The inputs of the data pipeline example, shared by its two runs, `fetch.jl` and `model.jl`: one
# of each kind of file a tracked spec reads. WorldClim's two layers are downloads of a source cited
# by its paper, Natural Earth's outline of Scotland a download of a source identified by its web
# page, and a month of ERA5 over Scotland an extract the Climate Data Store makes to a request.
#
# Each spec is tracked by `handle`, the run's handle, under a root of its own; a `handle` of
# `nothing`, for a script started any other way, tracks nothing, so the same specs read the asset
# cache as an untracked script's do.

using EcoSISTEM
using EcoSISTEM: Track
using RasterDataSources

function inputs(handle)
    return (temperature = SourceSpec(WorldClim{BioClim}, :bio1,
                                     track = Track(handle, "worldclim")),
            rainfall = SourceSpec(WorldClim{BioClim}, :bio12,
                                  track = Track(handle, "worldclim")),
            summer = SourceSpec(EcoSISTEM.ERA, "t2m", years = 2015:2015,
                                months = 7:7, area = "Scotland",
                                atend = HoldAtEnd(),
                                track = Track(handle, "era5")),
            region = NaturalEarthSpec("Scotland", coverage = LargestLandmass(),
                                      track = Track(handle, "naturalearth")))
end
