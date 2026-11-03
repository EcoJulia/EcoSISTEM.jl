# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The inputs of the data pipeline example, shared by its two runs, `fetch.jl` and `model.jl`: one
# of each kind of file a tracked spec reads, and a file of your own. WorldClim's two layers are
# downloads of a source cited by its paper, Natural Earth's outline of Scotland a download of a
# source identified by its web page, a month of ERA5 over Scotland an extract the Climate Data
# Store makes to a request, and a box around Scotland a GeoJSON file of your own, cut from Natural
# Earth's table of regions and recorded as such. The model run is built from WorldClim's two
# layers within the outline and the box, so needs no Climate Data Store key; the ERA5 month's
# fetch run shows how an extract is registered.
#
# Each catalogued spec is tracked by `handle`, the run's handle, under a root of its own; a
# `handle` of `nothing`, for a script started any other way, tracks nothing, so the same specs read
# the asset cache as an untracked script's do. The box, a file of your own, enters its spec through
# its record, which `record_provenance` writes and, under a run, registers or reads from the store,
# so the spec itself carries no track.

using EcoSISTEM
using EcoSISTEM: Track
using EcoSISTEM.Provenancer: SourceFacts, OtherOrigin, record_provenance
using RasterDataSources
using Unitful: ustrip

# The box, a file of your own standing for one you made yourself: written under the asset cache
# from the shipped table's bounding box of Scotland where it is not there, and recorded with
# Natural Earth's facts as a file cut from its data, which is what the table was made from.
const BOX = joinpath(EcoSISTEM.assetdir(), "pipeline-example",
                     "scotland_box.geojson")

function writebox(path = BOX)
    isfile(path) && return path
    box = EcoSISTEM.boundingbox("Scotland")
    (y1, y2), (x1, x2) = ustrip.(box.Y), ustrip.(box.X)
    ring = join(("[$x, $y]"
                 for (x, y) in ((x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1))),
                ", ")
    mkpath(dirname(path))
    write(path,
          """{"type": "FeatureCollection", "features": [{"type": "Feature", """ *
          """"properties": {"name": "Scotland box"}, "geometry": {"type": "Polygon", """ *
          """"coordinates": [[$ring]]}}]}\n""")
    return path
end

# The box, written where it is not there, recorded once as a file cut from Natural Earth's data,
# and its record returned. Under a run's handle the fetch run registers it and a model run reads it
# from the run's store; with none the record is written beside the file and left there.
function recordbox(handle)
    writebox()
    facts = SourceFacts(EcoSISTEM.datasetinfo(EcoSISTEM.NaturalEarthLevel),
                        title = "A box around Scotland (scotland_box)")
    return record_provenance(BOX, facts, role = :region,
                             route = OtherOrigin(),
                             track = Track(handle, "box"))
end

# The box as a spec, built from its record, which says where the file is: in the asset cache, or,
# under a run, the run's copy put there.
boxspec(handle) = ShapeSpec(recordbox(handle))

# The specs, with the box's built from its record unless another is given: `nothing` for a fetch
# run of another source, whose config neither writes nor reads the box.
function inputs(handle; box = boxspec(handle))
    return (temperature = SourceSpec(WorldClim{BioClim}, :bio1,
                                     track = Track(handle, "worldclim")),
            rainfall = SourceSpec(WorldClim{BioClim}, :bio12,
                                  track = Track(handle, "worldclim")),
            summer = SourceSpec(EcoSISTEM.ERA, "t2m", years = 2015:2015,
                                months = 7:7, area = "Scotland",
                                atend = HoldAtEnd(),
                                track = Track(handle, "era5")),
            region = NaturalEarthSpec("Scotland", coverage = LargestLandmass(),
                                      track = Track(handle, "naturalearth")),
            box = box)
end
