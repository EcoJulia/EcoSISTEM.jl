# SPDX-License-Identifier: LGPL-3.0-or-later

module TestAsset

using EcoSISTEM
using EcoSISTEM: CachedAsset, assetpath, provenance
using Sockets
using Test

include("tomlstrings.jl")

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
    @test !any(s -> occursin(homedir(), s) || occursin(dir, s) ||
                    occursin(gethostname(), s), tomlstrings(sidecar))
    r = provenance(dest)
    @test r.role === :habitat && r.dataset == "TwentyCR" && isnothing(r.code)
    @test r.path == "fixture.bin" && r.url == site.url
    @test r.bytes == 10_000 && r.sha256 == E._sha256(dest)
    @test !isnothing(r.fetched) && r.doi == "10.1002/qj.3598" &&
          !isempty(r.citation)
    text = read(sidecar, String)
    @test occursin("etag = \"\\\"fixture-etag\\\"\"", text) ||
          occursin("fixture-etag", text)
    # Already there: nothing is fetched (the server is gone), and verification passes.
    close(site.server)
    @test assetpath(asset) == dest
    @test assetpath(asset, verify = true) == dest
    # A replaced file fails verification by its checksum, and verification says what it checked.
    @test E._verifyfile(dest)
    write(dest, rand(UInt8, 10_000))
    @test_throws "does not match the checksum" assetpath(asset, verify = true)
    # A record with no checksum, or none of ours, is nothing to check against.
    bare = joinpath(dir, "bare.bin")
    write(bare, "x")
    @test !E._verifyfile(bare)
    write(E._sidecarpath(bare),
          "writer = \"EcoSISTEM 0.8.0\"\nfile = \"bare.bin\"\n")
    @test !E._verifyfile(bare)
    write(E._sidecarpath(bare), "writer = \"someone else\"\nsha256 = \"00\"\n")
    @test !E._verifyfile(bare)
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
    # A region's zip records a region, with Natural Earth's row.
    nerec = E._assetrecord(CachedAsset(E.NaturalEarthLevel,
                                       "https://example.org/ne.zip",
                                       path = dest), dest,
                           (url = "https://example.org/ne.zip",
                            headers = Pair{String, String}[], status = 200))
    @test nerec["role"] == "region" && nerec["source"] == "NaturalEarth"
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

@testset "a present file fetched for another request is refused, never served" begin
    E = EcoSISTEM
    dir = mktempdir()
    path = joinpath(dir, "era5_t2m_1990s.nc")
    write(path, "not a netCDF file, and never opened")
    asked(; kw...) = E.CDSRequest(E.ERA, "t2m"; years = 1990:1999, path = path,
                                  kw...)
    # No record: presence is all there is, so the file is served and recorded as this request.
    @test assetpath(asked(area = (40, -25, -36, 56))) == path
    @test provenance(path).request["area"] == [40, -25, -36, 56]
    # The same request again, and one asking for the same data in another order, are served.
    @test assetpath(asked(area = [40.0, -25.0, -36.0, 56.0])) == path
    @test assetpath(E.CDSRequest(E.ERA, "t2m", years = reverse(1990:1999),
                                 area = (40, -25, -36, 56), path = path)) ==
          path
    # Another box, other months, or the globe are refused, naming what differs, and the file and
    # its record are left as they were.
    sidecar = read(E._sidecarpath(path), String)
    @test_throws "area: recorded" assetpath(asked(area = (41, -25, -36, 56)))
    @test_throws "month: recorded" assetpath(asked(area = (40, -25, -36, 56),
                                                   months = 6:8))
    @test_throws "different request" assetpath(asked())
    @test read(E._sidecarpath(path), String) == sidecar
    @test read(path, String) == "not a netCDF file, and never opened"
    # A record of ours written in another order or shape - years out of order, a single value not
    # in a list, as an older release might write it - asks for the same data, so is served.
    older = joinpath(dir, "older.nc")
    write(older, "x")
    rec = E._requestrecord(E.CDSRequest(E.ERA, "t2m", years = 1990:1991,
                                        path = older),
                           older, nothing)
    rec["request"]["year"] = ["1991", "1990"]
    rec["request"]["product_type"] = "monthly_averaged_reanalysis"
    E._writesidecar(older, rec)
    @test assetpath(E.CDSRequest(E.ERA, "t2m", years = 1990:1991, path = older)) ==
          older
    @test_throws "year: recorded" assetpath(E.CDSRequest(E.ERA, "t2m",
                                                         years = 1990:1992,
                                                         path = older))
    # A record another program wrote says nothing about our request, so it is not compared.
    other = joinpath(dir, "theirs.nc")
    write(other, "x")
    write(E._sidecarpath(other),
          "writer = \"someone else\"\nrequest = {area = [1, 2, 3, 4]}\n")
    @test assetpath(E.CDSRequest(E.ERA, "t2m", years = 1990:1990, path = other)) ==
          other
end

@testset "the asset cache moves to a directory of the user's, and back" begin
    E = EcoSISTEM
    # The preference is written into the first environment on the load path that has EcoSISTEM,
    # so one made here, put first, keeps it out of every real environment.
    env = mktempdir()
    write(joinpath(env, "Project.toml"),
          "[deps]\nEcoSISTEM = \"$(Base.PkgId(EcoSISTEM).uuid)\"\n")
    scratch = E.assetdir()
    pushfirst!(LOAD_PATH, env)
    try
        cache = joinpath(mktempdir(), "cache")
        @test E.set_assetdir!(cache) == cache && isdir(cache)
        @test E.assetdir() == joinpath(cache, "EcoSISTEM")
        @test E.assetdir(owner = E.ERA) == joinpath(cache, "EcoSISTEM", "ERA")
        # Parsed, since TOML escapes the backslashes of a Windows path in the file's text.
        prefs = E.TOML.parsefile(joinpath(env, "LocalPreferences.toml"))
        @test prefs["EcoSISTEM"]["assetdir"] == cache
        # A spec built now fetches into it.
        spec = SourceSpec(E.ERA, "t2m", years = 1990:1999)
        @test startswith(only(spec.files).path, cache)
        # `nothing` returns the cache to the scratch space.
        @test E.set_assetdir!(nothing) == dirname(scratch)
        @test E.assetdir() == scratch
    finally
        E.set_assetdir!(nothing)
        filter!(!=(env), LOAD_PATH)
    end
    @test E.assetdir() == scratch
end

end
