# SPDX-License-Identifier: LGPL-3.0-or-later

module TestTrack

using EcoSISTEM
using EcoSISTEM: Track, ReadKey, LayerCache, ERA, fetchfiles
using Extents
using RasterDataSources
using Test

# A handle no pipeline knows, standing in for a run: its route methods record what they are asked
# and then stop the call, so the routing is tested with no registry and nothing read or fetched.
struct Recording
    calls::Vector{Any}
end

struct Stopped <: Exception end

function EcoSISTEM._trackedprovide!(track::Track{Recording}, paths)
    push!(track.handle.calls, (:provide, track.root, copy(paths)))
    return throw(Stopped())
end

function EcoSISTEM._trackedfetch(track::Track{Recording}, spec)
    push!(track.handle.calls, (:fetch, track.root, spec))
    return throw(Stopped())
end

# A run that answers and puts nothing in place, which a read must notice rather than then fetch.
struct Empty end

EcoSISTEM._trackedprovide!(::Track{Empty}, paths) = nothing

@testset "a track names a run and the root its files are filed under" begin
    @test Track(nothing, "era5").root == "era5"
    @test Track(nothing, "inputs/era5.v2_1").root == "inputs/era5.v2_1"
    # The root is checked whatever the handle, so a script fails alike with the pipeline or without.
    for bad in ("", "/era5", "era5/", "a//b", ".era5", "era5./x", "era 5", "\u00e9r")
        @test_throws "cannot begin a data product name" Track(nothing, bad)
    end
    # No run means no tracking: the spec stores none, so there is one untracked state.
    @test isnothing(SourceSpec(WorldClim{BioClim}, :bio1,
                               track = Track(nothing, "wc")).track)
    @test isnothing(NaturalEarthSpec("Scotland",
                                     track = Track(nothing, "ne")).track)
    # A track among the read keywords is a misspelt `track`, refused rather than read with.
    @test_throws "`track = Track(handle, root)`" SourceSpec(WorldClim{BioClim},
                                                            :bio1,
                                                            trak = Track(Recording([]),
                                                                         "wc"))
end

@testset "a tracked spec's files take its run's route, an untracked spec's do not" begin
    run = Recording([])
    track = Track(run, "worldclim")
    spec = SourceSpec(WorldClim{BioClim}, :bio1, track = track)
    # A read: every file the source would fetch is put in place by the run, named where it would be.
    @test_throws Stopped read(spec)
    kind, root, paths = only(run.calls)
    @test kind === :provide && root == "worldclim"
    @test paths ==
          EcoSISTEM._sourcepaths(WorldClim{BioClim}, :bio1;
                                 EcoSISTEM._readkw(spec)...)
    # Deciding a grid opens the files before any read, so it takes the same route.
    empty!(run.calls)
    @test_throws Stopped StudyArea(regime = spec, verbosity = :silent)
    @test first(only(run.calls)) === :provide
    # The fetch run's call takes the run's own fetch, which registers what it fetches.
    empty!(run.calls)
    @test_throws Stopped fetchfiles(spec)
    @test only(run.calls) == (:fetch, "worldclim", spec)
    # A spec naming its files asks for them by where they would be, before anything opens them.
    empty!(run.calls)
    file = tempname() * ".nc"
    @test_throws Stopped read(SourceSpec(ERA, "t2m", file = file,
                                         track = track))
    @test last(only(run.calls)) == [file]
    # A named region asks for its level's zip.
    empty!(run.calls)
    @test_throws Stopped read(NaturalEarthSpec("Scotland", track = track))
    @test endswith(only(last(only(run.calls))), ".zip")
    # Untracked, nothing is asked of a run.
    @test isnothing(EcoSISTEM._provide!(nothing,
                                        SourceSpec(WorldClim{BioClim}, :bio1)))
end

@testset "a file the run did not put in place is never fetched instead" begin
    file = tempname() * ".nc"
    @test_throws "did not provide $(basename(file))" read(SourceSpec(ERA, "t2m",
                                                                     file = file,
                                                                     track = Track(Empty(),
                                                                                   "era5")))
    # And a handle that is not a pipeline run is refused, naming its type.
    @test_throws "a `Int64` is not one" read(SourceSpec(ERA, "t2m", file = file,
                                                        track = Track(42,
                                                                      "era5")))
end

@testset "a spec rebuilt inside the package keeps its track" begin
    run = Recording([])
    track = Track(run, "chelsa")
    # The divisor a per-cell supply is desugared into reads through the same run.
    gsp = EcoSISTEM._desugarsupply(SourceSpec(CHELSA{BioClimPlus}, :gsp,
                                              track = track))
    @test all(layer.track === track for layer in gsp.layers)
    # So does each layer of a multi-layer spec read one at a time.
    wc = SourceSpec(WorldClim{BioClim}, [1, 12], track = track)
    @test_throws Stopped EcoSISTEM._stackcached(wc, LayerCache(), nothing, 1)
    @test basename(only(last(only(run.calls)))) == "wc2.1_10m_bio_1.tif"
end

@testset "a tracked read is not an untracked one" begin
    track = Track(Recording([]), "wc")
    spec = SourceSpec(WorldClim{BioClim}, :bio1)
    tracked = SourceSpec(WorldClim{BioClim}, :bio1, track = track)
    @test ReadKey(tracked) != ReadKey(spec)
    @test ReadKey(tracked) ==
          ReadKey(SourceSpec(WorldClim{BioClim}, :bio1, track = track))
    @test ReadKey(NaturalEarthSpec("Scotland", track = track)) !=
          ReadKey(NaturalEarthSpec("Scotland"))
    # The two hash alike, the track being left out of the hash, and are still two entries.
    @test length(Dict(ReadKey(spec) => 1, ReadKey(tracked) => 2)) == 2
    # A whole untracked read in the cache serves an untracked window of it, never a tracked one.
    cache = LayerCache()
    cache.reads[ReadKey(spec)] = :whole
    box = Extents.Extent(X = (-5.0, 0.0), Y = (55.0, 58.0))
    @test EcoSISTEM._servedread(cache, ReadKey(spec, cut = box)) === :whole
    @test isnothing(EcoSISTEM._servedread(cache, ReadKey(tracked, cut = box)))
end

@testset "the spelling shows a track, a request's region does not" begin
    track = Track(Recording([]), "era5")
    @test occursin("track = Track(..., \"era5\")",
                   string(SourceSpec(WorldClim{BioClim}, :bio1, track = track)))
    region = NaturalEarthSpec("Scotland", track = track)
    @test occursin("track = Track(..., \"era5\")", string(region))
    # A tracked region as an ERA5 area is accepted, its track ignored: no file of it is read, and the
    # request records the ground it covers.
    request = only(SourceSpec(ERA, "t2m", years = 2000:2000, area = region).files)
    @test request.region == string(NaturalEarthSpec("Scotland"))
    @test string(Track(nothing, "x")) == "Track(nothing, \"x\")"
end

end
