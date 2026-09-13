# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEMERAExt

using EcoSISTEM
using CDSAPI

# Fetch a request's file from the Climate Data Store, through CDSAPI's client, which submits the
# job, polls it and downloads the result. The reply lands under a temporary name beside the target
# and is renamed in, so an interrupted transfer never looks complete; `assetpath` is what decides a
# file is present, and it checks nothing but presence.
function EcoSISTEM._fetchcds(request::EcoSISTEM.CDSRequest)
    mkpath(dirname(abspath(request.path)))
    part = request.path * ".part"
    try
        CDSAPI.retrieve(request.dataset, request.request, part, wait = 20.0)
        mv(part, request.path, force = true)
    finally
        rm(part, force = true)
    end
    return request.path
end

end
