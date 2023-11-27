module EcoSISTEMDataPipelineExt

using DataPipeline
using EcoSISTEM
using EcoSISTEM.ClimatePref
using Tar

@info "Creating data pipeline interface ..."

function EcoSISTEM.ClimatePref.unzip(path::String)
    newpath = joinpath(splitpath(path)[1:end-1])
    Tar.extract($path, $newpath)
    @info "Unzipped to $newpath"
    return newpath
end

end
