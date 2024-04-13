module EcoSISTEMDataPipelineExt

using EcoSISTEM
using ZipArchives

@info "Creating data pipeline interface for EcoSISTEM..."

function EcoSISTEM.unziptemp(path::String)
    zr = ZipReader(read(open(path)))
    newpath = mktempdir()
    files = zip_names(zr)
    for file in files
        if zip_isdir(zr, file)
            mkdir(joinpath(newpath, file))
        end
    end
    for file in files
        if !zip_isdir(zr, file)
            write(joinpath(newpath, file), zip_readentry(zr, file, String))
        end
    end
    @debug "Unzipped to $newpath"
    return newpath
end

end
