using Requires
function __init__()
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" begin
        println("Creating RCall interface ...")
        include("ProcessData.jl")
        export processMet, writeMet
        include("DownloadClimate.jl")
        export MetOfficeDownload, getMetparams, getMetdata
    end
end

@warn "This functionality remains under development!"

include("ClimateTypes.jl")
export Worldclim, Bioclim, ERA, CERA, CRUTS, CHELSA, Reference

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA, readCERA, readfile, readCHELSA, readMet

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, upresolution, downresolution

include("Plotting.jl")
export getprofile
