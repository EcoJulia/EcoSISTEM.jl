# SPDX-License-Identifier: LGPL-3.0-or-later

module TestAsset

using EcoSISTEM
using Test

# The download helper is exercised against `file://` URLs, which curl serves without a network and
# for which it reports no HTTP status. What that does not cover, and cannot without a server: a
# transfer cut off part way and resumed from a byte range, and a server answering a range with the
# whole file. Both branches are read, not run.
@testset "a download resumes from what it holds and fails after its attempts" begin
    E = EcoSISTEM
    @test E._resumeheaders(0) == Pair{String, String}[]
    @test E._resumeheaders(1024) == ["Range" => "bytes=1024-"]
    dir = mktempdir()
    src = joinpath(dir, "source.bin")
    write(src, rand(UInt8, 10_000))
    dest = joinpath(dir, "dest.bin")
    @test E._download("file://" * src, dest) == dest
    @test read(dest) == read(src)
    # A URL that cannot be fetched is retried as many times as asked and then rethrown.
    missing = joinpath(dir, "missing.bin")
    @test_throws Exception E._download("file://" * missing,
                                       joinpath(dir, "never.bin"),
                                       attempts = 2)
    # The downloader is held to HTTP/1.1.
    @test E._http11downloader() isa E.Downloads.Downloader
end

@testset "an asset resolves under its owner and is fetched once" begin
    E = EcoSISTEM
    dir = mktempdir()
    src = joinpath(dir, "asset_fixture.bin")
    write(src, rand(UInt8, 100))
    asset = E.CachedAsset(TwentyCR, "file://" * src)
    path = E.assetpath(asset)
    @test dirname(path) == E.assetdir(owner = TwentyCR)
    @test basename(path) == "asset_fixture.bin" && read(path) == read(src)
    # Already there: nothing is fetched, and no part file is left beside it.
    rm(src)
    @test E.assetpath(asset) == path
    @test !any(f -> occursin(".part-", f), readdir(dirname(path)))
    rm(path)
end

end
