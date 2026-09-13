# SPDX-License-Identifier: LGPL-3.0-or-later

module TestAsset

using EcoSISTEM
using EcoSISTEM: CachedAsset, assetpath, provenance
using Sockets
using Test

# A one-file HTTP/1.1 server on the loopback, serving `data` and honouring byte ranges, with the
# two faults the resumable download must survive available on demand: cut the first transfer off
# after `dropafter` bytes, or ignore ranges and answer every request with the whole file. Closes
# the connection after each reply, as a static file host does.
function servefile(data::Vector{UInt8}; honourranges = true,
                   dropafter = nothing)
    server = listen(Sockets.localhost, 0)
    port = getsockname(server)[2]
    dropped = Ref(false)
    function handle(sock)
        lines = String[]
        while true
            line = readline(sock)
            isempty(line) && break
            push!(lines, line)
        end
        method = first(split(first(lines)))
        range = nothing
        for line in lines
            m = match(r"^Range:\s*bytes=(\d+)-"i, line)
            isnothing(m) || (range = parse(Int, m[1]))
        end
        start = (honourranges && !isnothing(range)) ? range : 0
        body = data[(start + 1):end]
        status = start > 0 ? "206 Partial Content" : "200 OK"
        head = "HTTP/1.1 $status\r\nContent-Length: $(length(body))\r\n" *
               "ETag: \"fixture-etag\"\r\nLast-Modified: Sat, 13 Sep 2026 00:00:00 GMT\r\n" *
               (start > 0 ?
                "Content-Range: bytes $start-$(length(data) - 1)/$(length(data))\r\n" :
                "") * "Connection: close\r\n\r\n"
        write(sock, head)
        if method == "HEAD"
            # nothing more
        elseif !isnothing(dropafter) && !dropped[]
            dropped[] = true
            write(sock, body[1:dropafter])
        else
            write(sock, body)
        end
        return close(sock)
    end
    @async while isopen(server)
        sock = try
            accept(server)
        catch
            break
        end
        @async handle(sock)
    end
    return (url = "http://127.0.0.1:$port/fixture.bin", server = server)
end

@testset "a download resumes from what it holds and survives a cut-off transfer" begin
    E = EcoSISTEM
    @test E._resumeheaders(0) == Pair{String, String}[]
    @test E._resumeheaders(1024) == ["Range" => "bytes=1024-"]
    data = rand(UInt8, 20_000)
    dir = mktempdir()
    # The plain path: one transfer, the response's final URL and validators come back.
    whole = servefile(data)
    part = joinpath(dir, "whole.bin.part")
    response = E._download(whole.url, part)
    @test read(part) == data && response.status == 200
    @test any(h -> lowercase(first(h)) == "etag", response.headers)
    close(whole.server)
    # Cut off after 7 000 bytes: the retry asks for the rest by range and gets a 206.
    cut = servefile(data, dropafter = 7_000)
    part = joinpath(dir, "cut.bin.part")
    response = @test_logs (:warn, r"download attempt 1 of 5 failed") E._download(cut.url,
                                                                                 part)
    @test read(part) == data && response.status == 206
    close(cut.server)
    # A server that ignores the range answers 200 to a resume: the part file is restarted, not
    # appended to.
    deaf = servefile(data, honourranges = false)
    part = joinpath(dir, "deaf.bin.part")
    write(part, rand(UInt8, 5_000))
    response = @test_logs (:warn, r"download attempt 1 of 5 failed") E._download(deaf.url,
                                                                                 part)
    @test read(part) == data && response.status == 200
    close(deaf.server)
    # A URL that cannot be fetched is retried as many times as asked and then rethrown, and the
    # part file is left for a later run to resume from.
    missing = joinpath(dir, "missing.bin")
    @test_throws Exception E._download("file://" * missing,
                                       joinpath(dir, "never.bin.part"),
                                       attempts = 2)
    @test E._http11downloader() isa E.Downloads.Downloader
end

@testset "an asset lands where it is told, once, with its provenance beside it" begin
    E = EcoSISTEM
    data = rand(UInt8, 10_000)
    site = servefile(data)
    dir = mktempdir()
    dest = joinpath(dir, "data", "fixture.bin")
    asset = CachedAsset(TwentyCR, site.url, path = dest)
    @test assetpath(asset) == dest
    @test read(dest) == data
    @test !isfile(dest * ".part") && !isfile(dest * ".lock")
    # The sidecar: the base fields, the URL and validators, the source's row (the URL is none of
    # its layers' files, so no code), and nothing that names this machine.
    sidecar = E._sidecarpath(dest)
    @test isfile(sidecar)
    text = read(sidecar, String)
    @test !occursin(homedir(), text) && !occursin(dir, text) &&
          !occursin(gethostname(), text)
    r = provenance(dest)
    @test r.role === :habitat && r.dataset == "TwentyCR" && isnothing(r.code)
    @test r.path == "fixture.bin" && r.url == site.url
    @test r.bytes == 10_000 && r.sha256 == E._sha256(dest)
    @test !isnothing(r.fetched) && r.doi == "10.1002/qj.3598" &&
          !isempty(r.citation)
    @test occursin("etag = \"\\\"fixture-etag\\\"\"", text) ||
          occursin("fixture-etag", text)
    # Already there: nothing is fetched (the server is gone), and verification passes.
    close(site.server)
    @test assetpath(asset) == dest
    @test assetpath(asset, verify = true) == dest
    # A replaced file fails verification by its checksum.
    write(dest, rand(UInt8, 10_000))
    @test_throws "does not match the checksum" assetpath(asset, verify = true)
    # An asset with no path lands in its owner's cache directory under the URL's name.
    @test E._localpath(CachedAsset(TwentyCR, "https://example.org/x/y.nc")) ==
          joinpath(E.assetdir(owner = TwentyCR), "y.nc")
    # A URL that is one of the source's layers' files is recorded as that layer.
    air = E.layerinfo(TwentyCR, "air").file
    rec = E._assetrecord(CachedAsset(TwentyCR, air, path = dest), dest,
                         (url = air, headers = Pair{String, String}[],
                          status = 200))
    @test rec["code"] == "air" && rec["source"] == "TwentyCR"
    @test rec["licence"] == E.datasetinfo(TwentyCR).licence
    # A region's zip records a region, and no dataset row.
    nerec = E._assetrecord(CachedAsset(E.NaturalEarthLevel,
                                       "https://example.org/ne.zip",
                                       path = dest), dest,
                           (url = "https://example.org/ne.zip",
                            headers = Pair{String, String}[], status = 200))
    @test nerec["role"] == "region" && !haskey(nerec, "source")
    # So does a shape file's, whose owner is a parametric type rather than a source.
    shrec = E._assetrecord(CachedAsset(ShapeSpec,
                                       "https://example.org/site.zip",
                                       path = dest), dest,
                           (url = "https://example.org/site.zip",
                            headers = Pair{String, String}[], status = 200))
    @test shrec["role"] == "region" && !haskey(shrec, "doi")
end

@testset "a fetch waits on another's lock and then finds the file" begin
    E = EcoSISTEM
    dir = mktempdir()
    dest = joinpath(dir, "shared.bin")
    calls = Ref(0)
    held = E.mkpidlock(dest * ".lock")
    waiter = @async E._fetchlocked(dest) do
        calls[] += 1
        return write(dest, "fetched by the waiter")
    end
    sleep(0.5)
    @test !istaskdone(waiter)
    # The holder writes the file and lets go: the waiter takes the lock, sees the file, fetches
    # nothing.
    write(dest, "fetched by the holder")
    close(held)
    wait(waiter)
    @test calls[] == 0 && read(dest, String) == "fetched by the holder"
    @test !isfile(dest * ".lock")
    # With no file to find, the fetch runs.
    rm(dest)
    E._fetchlocked(dest) do
        calls[] += 1
        return write(dest, "x")
    end
    @test calls[] == 1
end

end
