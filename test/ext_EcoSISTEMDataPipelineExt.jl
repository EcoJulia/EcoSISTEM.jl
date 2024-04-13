module ExtEcoSISTEMDataPipelineExt

using EcoSISTEM
using DataPipeline
using Test

if !Sys.iswindows()
    if isdir("assets")
        @testset "DataPipeline unziptemp() test" begin
            @test isdir(EcoSISTEM.unziptemp("assets/WorldClim/BioClim/zips/wc2.1_10m_bio.zip"))
        end
    end
end

end
