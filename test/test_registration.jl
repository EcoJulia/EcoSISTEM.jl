# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/registration.jl`: what a tracked spec's files are registered as with a data pipeline
# run, decided from plain data, so with no registry - and `fetchfiles` and `missingfiles` over the
# inputs a habitat is built from. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_registration.jl"])'

module TestRegistration

using EcoSISTEM
using EcoSISTEM: Track, InputRecord, ERA, fetchfiles, missingfiles
using EcoSISTEM.Units
using Dates: Dates
using RasterDataSources
using TOML: TOML
using Test
using Unitful: K

const E = EcoSISTEM

# A record of a file as a fetch writes it, read back, for a source's row.
function _record(owner, layer, path; route, extra = Dict{String, Any}())
    record = E._baserecord(path, :habitat)
    merge!(record, E._cataloguerecord(owner, layer, path, route = route), extra)
    E._writesidecar(path, record)
    return provenance(path)
end

@testset "a tracked file's data products are named by its stem under the root" begin
    track = Track(nothing, "inputs/chelsa")
    path = "/somewhere/CHELSA_bio1_1981-2010_V.2.1.tif"
    @test E._filestem(path) == "CHELSA_bio1_1981-2010_V.2.1"
    @test E._dataname(track, path) ==
          "inputs/chelsa/data/CHELSA_bio1_1981-2010_V.2.1"
    @test E._sidecarname(track, path) ==
          "inputs/chelsa/provenance/CHELSA_bio1_1981-2010_V.2.1"
    @test E._depositname(track) == "inputs/chelsa/deposit"
    @test E._filestem("air.2m.mon.mean.nc") == "air.2m.mon.mean"
    @test isnothing(E._checkstems([path, "x/ne_10m_admin_0_countries.zip"]))
    # A stem the pipeline cannot name is refused before anything is fetched, naming the file.
    @test_throws "cru ts 4.08" E._checkstems(["data/cru ts 4.08.tif"])
    @test_throws "Rename the file" E._checkstems(["data/.hidden.tif"])
end

@testset "each file of a spec has the route it comes by" begin
    published = only(E._specroutes(SourceSpec(WorldClim{BioClim}, :bio1)))
    @test published.route === :published &&
          endswith(published.path, "wc2.1_10m_bio_1.tif")
    @test only(E._specroutes(SourceSpec(TwentyCR, "air"))).route === :published
    @test only(E._specroutes(SourceSpec(ERA, "t2m", years = 1990:1999))).route ===
          :extract
    file = tempname() * ".tif"
    @test only(E._specroutes(SourceSpec(CRUTS, "tavg", file = file))) ==
          (path = file, route = :given, code = "tavg", url = nothing)
    @test endswith(only(E._specroutes(NaturalEarthSpec("Scotland"))).path,
                   ".zip")
end

@testset "a file's version is its source's, its bytes' or a stamped one" begin
    day = Dates.Date(2023, 1, 4)
    # Not yet registered, or registered with these very bytes: the source's own, in full.
    @test E._registeredversion("2.1", nothing, false, day) == "2.1.0"
    @test E._registeredversion("", nothing, false, day) == "1.0.0"
    @test E._registeredversion("3", nothing, false, day) == "3.0.0"
    # Its source's version held by other bytes: a provider's silent change, dated by its day.
    stamped = E._registeredversion("2.1", nothing, true, day)
    @test stamped == "2.1.0+20230104" && E._isstamped(stamped)
    @test !E._isstamped("2.1.0")
    # These bytes registered under a stamp already: that version, so a third fetch adds nothing.
    @test E._registeredversion("2.1", "2.1.0+20230104", true,
                               Dates.Date(2026)) ==
          "2.1.0+20230104"
    # A registry takes a stamp as it stands.
    @test E._isregistryversion(stamped)
end

@testset "a file's registration keywords come from its record and its row" begin
    dir = mktempdir()
    # A CHELSA layer as published: named by its deposit's DOI, its source's own.
    chelsa = joinpath(dir, "CHELSA_bio1_1981-2010_V.2.1.tif")
    write(chelsa, "x")
    rec = E.datasetinfo(CHELSA{BioClim})
    r = _record(CHELSA{BioClim}, E.layerinfo(CHELSA{BioClim}, :bio1), chelsa,
                route = :published,
                extra = Dict{String, Any}("modified" => "2023-01-04T00:00:00Z"))
    kw = E._sourcekeywords(r, rec, :published, version = "2.1.0",
                           date = rec.releasedate)
    @test kw.doi == "10.16904/envidat.228.v2.1" && !haskey(kw, :bibtex)
    @test kw.primary === true && kw.licence == "CC0-1.0"
    @test kw.title == r.title && kw.release_version == "2.1.0"
    @test kw.release_date == "2021-06-01" &&
          kw.description == E._requestpage(rec)
    # Its stamp is dated by when the provider changed it, else when it was fetched.
    @test E._stampday(r, Dates.Date(2026, 10, 27)) == Dates.Date(2023, 1, 4)
    @test E._registereddate(r, rec, Dates.Date(2023, 1, 4)) == "2023-01-04"
    @test E._registereddate(r, rec, nothing) == "2021-06-01"
    # A paper's file is not its source's own, whatever route it came by.
    tcr = joinpath(dir, "air.2m.mon.mean.nc")
    write(tcr, "x")
    p = _record(TwentyCR, E.layerinfo(TwentyCR, "air"), tcr, route = :published)
    @test p.primary === false
    @test E._sourcekeywords(p, E.datasetinfo(TwentyCR), :published,
                            version = "3.0.0", date = nothing).primary === false
    # A row identified by its web page names its source by its citation's key, and a file stating
    # its own version is dated by when its provider changed it.
    ne = E.datasetinfo(E.NaturalEarthLevel)
    zip = joinpath(dir, "ne_10m_admin_0_countries.zip")
    write(zip, "x")
    z = _record(E.NaturalEarthLevel, nothing, zip, route = :published)
    @test E._needsdate(z, ne)
    zk = E._sourcekeywords(z, ne, :published, version = "5.1.2", date = nothing)
    @test zk.bibtex == ne.request && !haskey(zk, :doi) && zk.primary === true
    dated = InputRecord(role = :region, dataset = "NaturalEarth",
                        modified = "2022-05-13T05:20:49Z")
    @test !E._needsdate(dated, ne)
    @test E._registereddate(dated, ne, nothing) == "2022-05-13"
    # Failing that, when it was fetched; with neither, the citation's date stands.
    fetched = InputRecord(role = :region, dataset = "NaturalEarth",
                          fetched = "2026-09-03T00:49:00Z")
    @test E._registereddate(fetched, ne, nothing) == "2026-09-03"
    @test isnothing(E._registereddate(InputRecord(role = :region,
                                                  dataset = "NaturalEarth"), ne,
                                      nothing))
    # An extract says what it was made to, and its deposit is the dataset's own.
    request = only(E.era5requests("t2m", 1990:1999, dir = dir))
    write(request.path, "x")
    E._writesidecar(request.path,
                    E._requestrecord(request, request.path, nothing))
    x = provenance(request.path)
    era = E.datasetinfo(ERA)
    xk = E._sourcekeywords(x, era, :extract, version = "1.0.0", date = nothing)
    @test xk.primary === false && occursin("2m_temperature", xk.description) &&
          occursin(E._requestpage(era), xk.description)
    dk = E._depositkeywords(era)
    @test dk.primary === true && dk.title == era.title && dk.doi == era.doi &&
          isnothing(dk.release_version) && dk.release_date == era.releasedate
    @test E._depositkeywords(E.datasetinfo(CERA)).primary === false
    # A row giving no release leaves its deposit's to its citation, whichever file comes first.
    nk = E._depositkeywords(ne)
    @test isnothing(nk.release_version) && isnothing(nk.release_date) &&
          nk.bibtex == "naturalearth"
end

@testset "a file the user gives is recorded with the release its track names" begin
    dir = mktempdir()
    file = joinpath(dir, "cru_ts4.08.1901.2023.tmp.dat.tif")
    write(file, "x")
    spec = SourceSpec(CRUTS, "tavg", file = file,
                      track = Track(:run, "cru", release_version = "4.8",
                                    release_date = "2024-05"))
    entry = only(E._specroutes(spec))
    @test entry.route === :given && entry.code == "tavg" && isnothing(entry.url)
    @test E._recordgiven(spec, entry) == file
    r = provenance(file)
    @test r.primary === false && r.version == "4.8" &&
          r.releasedate == "2024-05"
    @test isnothing(r.fetched) && r.dataset == "CRUTS"
    # A record there already is left alone, and a file that is not there is refused.
    @test E._recordgiven(spec, entry) == file
    @test_throws "is not there" E._recordgiven(spec,
                                               (; entry...,
                                                path = joinpath(dir,
                                                                "absent.tif")))
end

@testset "a record from an earlier version is brought up to date before it is registered" begin
    dir = mktempdir()
    path = joinpath(dir, "air.2m.mon.mean.nc")
    write(path, "x")
    # What a fetch by an earlier version recorded: the file, its checksum, no title, no URL.
    old = Dict{String, Any}("writer" => "EcoSISTEM 0.8.0", "role" => "habitat",
                            "file" => basename(path), "source" => "TwentyCR",
                            "code" => "air", "licence" => "", "sha256" => "abc",
                            "bytes" => 1)
    E._writesidecar(path, old)
    spec = SourceSpec(TwentyCR, "air", file = path)
    entry = (; only(E._specroutes(spec))..., route = :published,
             url = E.layerinfo(TwentyCR, "air").file)
    facts = E._currentfacts(spec, entry)
    # Read with the catalogue's facts, the record on disk is untouched.
    seen = E._recordwith(path, facts)
    @test !isempty(seen.title) && seen.url == entry.url &&
          seen.primary === false
    @test isempty(provenance(path).title)
    # Brought up to date, what the fetch recorded is kept and the writer is this version.
    r = E._refreshrecord!(path, facts)
    @test r.title == seen.title && r.licence == "CC0-1.0" && r.sha256 == "abc"
    @test startswith(TOML.parsefile(E._sidecarpath(path))["writer"],
                     "EcoSISTEM $(pkgversion(EcoSISTEM))")
    # Up to date already, it is not written again.
    before = mtime(E._sidecarpath(path))
    sleep(1.1)
    E._refreshrecord!(path, facts)
    @test mtime(E._sidecarpath(path)) == before
    @test isnothing(E._refreshrecord!(joinpath(dir, "none.nc"), facts))
end

@testset "a config's missing lines and a file's reasons are written out" begin
    track = Track(nothing, "era5")
    lines = E._writelines(track, "ECMWF", E._filetype("x/era5_t2m_1990s.nc"))
    # The source's namespace for the run's inputs and outputs, since it registers that source alone.
    @test occursin("default_input_namespace: ECMWF\n  default_output_namespace: ECMWF",
                   lines)
    # The files and their records, each entry with the description DataPipeline needs of one.
    @test occursin("- data_product: era5/data/*\n  description: ", lines)
    @test occursin("  file_type: nc\n- data_product: era5/provenance/*\n  description: ",
                   lines)
    @test endswith(lines, "  file_type: toml") &&
          count("description: ", lines) == 2
    @test E._filetype("ne_10m_admin_0_map_units.zip") == "zip"
    # A model reads both with one entry, a `*` matching one segment, which leaves the deposit out.
    @test E._readlines(track, "ECMWF") ==
          "read:\n- data_product: era5/*/*\n  use:\n    namespace: ECMWF"
    reason(; present = true, held = "2.1.0", same = true, recorded = true,
           named = true) = E._fetchreason(present = present, held = held,
                                          same = same, recorded = recorded,
                                          named = named)
    @test isnothing(reason())
    @test startswith(reason(present = false), "not held")
    @test startswith(reason(held = nothing, named = false), "not registered")
    @test occursin("stamped version", reason(held = nothing))
    @test occursin("other facts", reason(same = false))
    @test occursin("provenance record", reason(recorded = false))
    @test occursin("no `read:` entry", E._modelreason(nothing, track, "ECMWF"))
    @test occursin("fair pull", E._modelreason(true, track, "ECMWF"))
    @test occursin("run the fetch run first",
                   E._modelreason(false, track, "ECMWF"))
end

@testset "fetchfiles and missingfiles take whatever a habitat is built from" begin
    dir = mktempdir()
    held = joinpath(dir, "held.tif")
    write(held, "x")
    absent = joinpath(dir, "absent.tif")
    given = SourceSpec(CRUTS, "tavg", file = held)
    notheld = SourceSpec(CRUTS, "tavg", file = absent)
    # Wrappers and combinations are walked to the specs reading files, synthetic ones adding none.
    combined = ConstructedRasterSpec((a, b) -> a, given,
                                     UniformSpec(1.0K, axis = Temperature),
                                     axis = Temperature)
    inputs = (regime = Varying(combined,
                               IncrementBy(0.01K / month_mean_duration)),
              supply = [notheld], within = trues(2, 2))
    @test [only(s.files) for s in E._filespecs(inputs)] == [held, absent]
    @test fetchfiles((given, given)) == [held, held]
    @test missingfiles(inputs) == [(file = absent, reason = "not held")]
    @test isempty(missingfiles(given))
    # A dry run reports each file: its entry, whether it is present, and its bytes.
    dry = fetchfiles(inputs, dryrun = true)
    @test [d.entry for d in dry] == [held, absent]
    @test [d.present for d in dry] == [true, false]
    @test all(keys(d) == (:entry, :present, :bytes) for d in dry)
    # A shape's file and a named region's zip are found through the mask and coverage they sit in.
    shape = ShapeSpec(held)
    region = NaturalEarthSpec("Scotland")
    @test E._filespecs(ShapeMaskSpec(ConstructedShapeSpec(ShapeUnion(), shape,
                                                          region))) ==
          Any[shape, region]
    @test fetchfiles(ShapeCoverage(shape, axis = SurfaceArea)) == [held]
end

end
