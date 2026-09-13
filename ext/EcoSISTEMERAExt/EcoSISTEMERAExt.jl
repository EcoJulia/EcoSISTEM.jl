# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEMERAExt

using EcoSISTEM
using CDSAPI

# Fetch a request's file from the Climate Data Store, through CDSAPI's client, which submits the
# job, polls it and downloads the result. The reply lands under a temporary name beside the target
# and is renamed in, so an interrupted transfer never looks complete; `assetpath` is what decides a
# file is present, and it checks nothing but presence. Returns the job's final status document,
# whose `jobID` the file's provenance record keeps, with the client's version added.
function EcoSISTEM._fetchcds(request::EcoSISTEM.CDSRequest)
    mkpath(dirname(abspath(request.path)))
    part = request.path * ".part"
    status = try
        s = CDSAPI.retrieve(request.dataset, request.request, part, wait = 20.0)
        mv(part, request.path, force = true)
        s
    finally
        rm(part, force = true)
    end
    status isa AbstractDict || return status
    return merge(Dict{String, Any}(status),
                 Dict{String, Any}("client" => "CDSAPI $(pkgversion(CDSAPI))"))
end

end
