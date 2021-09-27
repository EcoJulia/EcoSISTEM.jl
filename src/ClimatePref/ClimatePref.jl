using Requires
function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" begin
        println("Creating ECMWF interface ...")
        include("ERA_interim_tools.jl")
        export retrieve_era_interim
        include("ECMWF_tools.jl")
        export retrieve_ECMWF
    end
end

@warn "This functionality remains under development!"

include("ClimateTypes.jl")
export Worldclim_monthly, Worldclim_bioclim, ERA, CERA, CRUTS, CHELSA_bioclim, CHELSA_monthly, Reference

include("ReadData.jl")
export read, searchdir, readworldclim, readbioclim, readERA, 
readCERA, readfile, readCHELSA_monthly, readCHELSA_bioclim, readCRUTS

include("ExtractClimate.jl")
export extractvalues

include("DataCleaning.jl")
export create_reference, upresolution, downresolution, downresolution!

include("Plotting.jl")
export getprofile

include("PhyloModels.jl")
export Brownian, Lambda, fitBrownian, fitLambda, varcovar

include("SimpleSDMInterface.jl")
export Worldclim_monthly, Worldclim_bioclim, CHELSA_bioclim, Landcover
