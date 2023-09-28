if !isdefined(Base, :get_extension)
    using Requires
end

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
            println("Creating ECMWF interface ...")
            include("ERA5_tools.jl")
            export retrieve_era5
        end
        @require DataPipeline = "9ced6f0a-eb77-43a8-bbd1-bbf3031b0d12" begin
            println("Creating data pipeline interface ...")
            include("Pipeline.jl")
            export unzip 
        end
    end
end

@warn "This functionality remains under development!"

include("Extensions.jl")
export unzip, retrieve_era5 

include("ClimateTypes.jl")
export Worldclim_monthly, Worldclim_bioclim, ERA, CERA, CRUTS, CHELSA_bioclim, CHELSA_monthly, Reference, Landcover

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA, 
readCERA, readfile, readCHELSA_monthly, readCHELSA_bioclim, readCRUTS, readlc

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, upresolution, downresolution, downresolution!, compressLC

include("Plotting.jl")
export getprofile

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

