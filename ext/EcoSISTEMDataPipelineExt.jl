# SPDX-License-Identifier: LGPL-3.0-or-later

module EcoSISTEMDataPipelineExt

using EcoSISTEM
import p7zip_jll

@info "Creating data pipeline interface for EcoSISTEM..."

function EcoSISTEM.unziptemp(path::String)
    # unzip into a fresh dir under EcoSISTEM's DataPipeline scratch subdirectory
    newpath = mktempdir(EcoSISTEM.assetdir(DataPipeline))
    run(`$(p7zip_jll.p7zip()) x -tzip -y -o$(newpath) $(path)`)
    @debug "Unzipped to $newpath"
    return newpath
end

end
