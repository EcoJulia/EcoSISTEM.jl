using Requires
function __init__()
    # Only try loading this when we're not on travis, and we've specified SIMULATION_ECMWF
    if !env_bool("TRAVIS") && env_bool("SIMULATION_ECMWF")
        @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
            println("Creating ECMWF interface ...")
            include("ERA_interim_tools.jl")
            export retrieve_era_interim
            include("ECMWF_tools.jl")
            export retrieve_ECMWF
        end
    end
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" begin
        println("Creating RCall interface ...")
        include("ProcessData.jl")
        export processMet, writeMet
        include("DownloadClimate.jl")
        export MetOfficeDownload, getMetparams, getMetdata
    end
end

include("ClimateTypes.jl")
export Worldclim, Bioclim, ERA, CERA, Reference

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA, readCERA, readfile, readCHELSA, readMet

include("ReadGBIF.jl")
export ReadGBIF

include("ReadTPL.jl")
export ReadTPL

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, gardenmask, genus_worldclim_average,
    genus_worldclim_monthly, upresolution, downresolution

include("Conversion.jl")
export worldclim_to_DB, era_to_DB, CHELSA_to_DB

include("Plotting.jl")
export getprofile

include("PhyloModels.jl")
export Brownian, Lambda, varcovar
