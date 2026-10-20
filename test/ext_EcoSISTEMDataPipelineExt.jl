# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `ext/EcoSISTEMDataPipelineExt/`, the FAIR data pipeline's unzip helper. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["ext_EcoSISTEMDataPipelineExt.jl"])'
# or as part of `core_ext.jl`.
#
# The pipeline calls themselves - `initialise`, `link_read!`, `link_write!`, `finalise` - need a
# local registry and a `fair run` around the process, so they are exercised by
# `examples/pipeline/AfricaRun.jl` rather than here. What this extension supplies is `unziptemp`,
# which needs neither, and is tested in full.

module TestEcoSISTEMDataPipelineExt

using Test
using EcoSISTEM
using DataPipeline
using p7zip_jll

# Build a zip holding a file at the top level and one in a subdirectory, and return its path. The
# subdirectory is what says the extraction preserves an archive's shape rather than flattening it,
# which is what a claimed dataset's layout depends on.
function _samplezip()
    src = mktempdir()
    mkpath(joinpath(src, "sub"))
    write(joinpath(src, "top.txt"), "top")
    write(joinpath(src, "sub", "nested.txt"), "nested")
    zipfile = joinpath(mktempdir(), "bundle.zip")
    cd(src) do
        return run(pipeline(`$(p7zip_jll.p7zip()) a -tzip $(zipfile) .`,
                            stdout = devnull))
    end
    return zipfile
end

@testset "the extension loads on DataPipeline and supplies the unzip hook" begin
    @test !isnothing(Base.get_extension(EcoSISTEM,
                                        :EcoSISTEMDataPipelineExt))
    # The parent declares the hook method-less, so the extension supplies the sole method - the
    # pattern that lets an extension add rather than overwrite.
    @test length(methods(EcoSISTEM.unziptemp)) == 1
    @test hasmethod(EcoSISTEM.unziptemp, Tuple{String})
end

@testset "an archive is unpacked whole, under the pipeline's own asset directory" begin
    out = EcoSISTEM.unziptemp(_samplezip())
    @test isdir(out)
    @test read(joinpath(out, "top.txt"), String) == "top"
    @test read(joinpath(out, "sub", "nested.txt"), String) == "nested"
    # Under `assetdir(DataPipeline)` rather than the system temporary directory, so a claimed
    # archive unpacks beside the rest of this package's cached data.
    @test startswith(out, EcoSISTEM.assetdir(DataPipeline))

    # A second claim gets its own directory: two data products unpacked in one run must not
    # overwrite each other's files, which a shared destination would let them do silently.
    @test EcoSISTEM.unziptemp(_samplezip()) != out
end

@testset "an archive that cannot be read fails loudly" begin
    notazip = joinpath(mktempdir(), "notazip.zip")
    write(notazip, "plain text, whatever the extension says")
    # p7zip reports a refusal on its own streams, which is twenty lines of alarming output in a
    # passing run, so the two calls that are meant to fail are silenced here rather than in
    # `unziptemp`, where the same text is what tells a user what went wrong.
    redirect_stdio(stdout = devnull, stderr = devnull) do
        # A non-zero exit must surface: returning an empty directory would leave the caller
        # reading nothing where it asked for a dataset.
        @test_throws ProcessFailedException EcoSISTEM.unziptemp(joinpath(mktempdir(),
                                                                         "absent.zip"))
        return @test_throws ProcessFailedException EcoSISTEM.unziptemp(notazip)
    end
end

end
