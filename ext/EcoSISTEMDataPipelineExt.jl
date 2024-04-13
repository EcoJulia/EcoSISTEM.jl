module EcoSISTEMDataPipelineExt

using EcoSISTEM
using Tar

@info "Creating data pipeline interface ..."

function EcoSISTEM.unzip(path::String)
    newpath = joinpath(splitpath(path)[1:(end - 1)])
    Tar.extract(path, newpath)
    @info "Unzipped to $newpath"
    return newpath
end

end
