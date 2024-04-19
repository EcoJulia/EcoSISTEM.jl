module EcoSISTEMDataPipelineExt

using EcoSISTEM
import p7zip_jll

@info "Creating data pipeline interface for EcoSISTEM..."

function EcoSISTEM.unziptemp(path::String)
    newpath = mktempdir()
    run(`$(p7zip_jll.p7zip()) x -tzip -y -o$(newpath) $(path)`)
    @debug "Unzipped to $newpath"
    return newpath
end

end
