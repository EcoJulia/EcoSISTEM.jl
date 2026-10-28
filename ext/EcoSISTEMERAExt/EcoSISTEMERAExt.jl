# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEMERAExt

using EcoSISTEM
using CDSAPI

# Fetch a request's file from the Climate Data Store, through CDSAPI's client, which submits the
# job, polls it and downloads the result. The reply lands under a temporary name beside the target
# and is renamed in, so an interrupted transfer never looks complete; `assetpath` is what decides a
# file is present, and it checks nothing but presence. A netCDF reply has its `history` replaced
# before the rename, since it records when and where the CDS converted the file and would make two
# fetches of the same data differ. Returns the job's final status document, whose `jobID` the
# file's provenance record keeps, with the client's version and the replaced `history` added.
function EcoSISTEM._fetchcds(request::EcoSISTEM.CDSRequest)
    mkpath(dirname(abspath(request.path)))
    part = request.path * ".part"
    status, history = try
        s = CDSAPI.retrieve(request.dataset, request.request, part, wait = 20.0)
        h = get(request.request, "data_format", nothing) == "netcdf" ?
            EcoSISTEM._normalisehistory!(part) : nothing
        mv(part, request.path, force = true)
        s, h
    finally
        rm(part, force = true)
    end
    status isa AbstractDict || return status
    added = Dict{String, Any}("client" => "CDSAPI $(pkgversion(CDSAPI))")
    isnothing(history) || (added["history"] = history)
    return merge(Dict{String, Any}(status), added)
end

end
