# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/registration.jl`: what a tracked spec's files are registered as with a data pipeline
# run, decided from plain data, so with no registry - and `fetchfiles` and `missingfiles` over the
# inputs a habitat is built from. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_registration.jl"])'

module TestRegistration

using EcoSISTEM
using EcoSISTEM: Track, ProvenanceRecord, ERA, fetchfiles, missingfiles
using EcoSISTEM.Units
using Dates: Dates
using RasterDataSources
using TOML: TOML
using Test
using Unitful: K

const E = EcoSISTEM

# A record of a file as a fetch writes it, read back, for a source's row.
function _record(owner, layer, path; route, extra = Dict{String, Any}())
    record = E.Provenancer.baserecord(path, :habitat)
    merge!(record, E._cataloguerecord(owner, layer, path, route = route), extra)
    E.Provenancer.writesidecar(path, record)
    return provenance(path)
end

@testset "a tracked file's data products are named by its stem under the root" begin
    track = Track(nothing, "inputs/chelsa")
    path = "/somewhere/CHELSA_bio1_1981-2010_V.2.1.tif"
    @test E.Provenancer._filestem(path) == "CHELSA_bio1_1981-2010_V.2.1"
    @test E.Provenancer._dataname(track, path) ==
          "inputs/chelsa/data/CHELSA_bio1_1981-2010_V.2.1"
    @test E.Provenancer._sidecarname(track, path) ==
          "inputs/chelsa/provenance/CHELSA_bio1_1981-2010_V.2.1"
    @test E.Provenancer._depositname(track) == "inputs/chelsa/deposit"
    @test E.Provenancer._filestem("air.2m.mon.mean.nc") == "air.2m.mon.mean"
    @test isnothing(E.Provenancer._checkstems([
                                                  path,
                                                  "x/ne_10m_admin_0_countries.zip"
                                              ]))
    # A stem the pipeline cannot name is refused before anything is fetched, naming the file.
    @test_throws "cru ts 4.08" E.Provenancer._checkstems(["data/cru ts 4.08.tif"])
    @test_throws "Rename the file" E.Provenancer._checkstems(["data/.hidden.tif"])
end

@testset "each file of a spec has the route it comes by" begin
    published = only(E._specroutes(SourceSpec(WorldClim{BioClim}, :bio1)))
    @test published.route isa E.VersionOfRecord &&
          endswith(published.path, "wc2.1_10m_bio_1.tif")
    @test only(E._specroutes(SourceSpec(TwentyCR, "air"))).route isa
          E.VersionOfRecord
    @test only(E._specroutes(SourceSpec(ERA, "t2m", years = 1990:1999))).route isa
          E.GeneratedOnRequest
    file = tempname() * ".tif"
    @test only(E._specroutes(SourceSpec(CRUTS, "tavg", file = file))) ==
          (path = file, route = E.OtherOrigin(), code = "tavg",
           url = nothing)
    @test endswith(only(E._specroutes(NaturalEarthSpec("Scotland"))).path,
                   ".zip")
    # A source resolving its own files answers for them as a type and a shape spec for its own;
    # anything saying nothing of its route is refused by name rather than failing in dispatch.
    @test E.Provenancer.route(WorldClim{BioClim}) == E.VersionOfRecord()
    @test E.Provenancer.route(NaturalEarthSpec("Scotland")) ==
          E.VersionOfRecord()
    @test E._fileroute(WorldClim{BioClim}) == E.VersionOfRecord()
    @test_throws "IsFetchable" E._fileroute(1)
    @test_throws "IsFetchable" E._fileroute(Int)
end

@testset "a name authority's release is registered as a catalogued source" begin
    wfo = WFOPlantList(release = "2026-06")
    entry = only(E._specroutes(wfo))
    @test entry.route == E.VersionOfRecord() && isnothing(entry.code)
    @test entry.url ==
          "https://zenodo.org/records/20782718/files/_DwC_backbone_R.zip"
    @test entry.path == only(E._specpaths(wfo)) == only(E._wouldbepaths(wfo))
    @test E._filespecs(wfo) == Any[wfo]
    # Registered with the catalogue row's facts: its version DOI, a deposit of several files.
    src = E._trackedsource(wfo)
    @test src.facts.doi == "10.5281/zenodo.20782718" && !src.fileversioned
    facts = E._currentfacts(wfo, entry)
    @test facts["url"] == entry.url
    @test facts["version"] == "2026-06" && facts["licence"] == "CC0-1.0"
    @test E._catalogued(wfo).namespace == "WFO"
    # A track saying which release a fetched file is from is refused, as for any fetched source.
    @test_throws "Leave them out" WFOPlantList(release = "2026-06",
                                               track = Track(nothing, "wfo",
                                                             release_version = "1"))
end

@testset "a file's version is its source's, its bytes' or a stamped one" begin
    day = Dates.Date(2023, 1, 4)
    # Not yet registered, or registered with these very bytes: the source's own, in full.
    @test E.Provenancer._registeredversion("2.1", nothing, false, day) ==
          "2.1.0"
    @test E.Provenancer._registeredversion(nothing, nothing, false, day) ==
          "1.0.0"
    @test E.Provenancer._registeredversion("3", nothing, false, day) == "3.0.0"
    # Its source's version held by other bytes: a provider's silent change, dated by its day.
    stamped = E.Provenancer._registeredversion("2.1", nothing, true, day)
    @test stamped == "2.1.0+20230104" && E.Provenancer._isstamped(stamped)
    @test !E.Provenancer._isstamped("2.1.0")
    # These bytes registered under a stamp already: that version, so a third fetch adds nothing.
    @test E.Provenancer._registeredversion("2.1", "2.1.0+20230104", true,
                                           Dates.Date(2026)) ==
          "2.1.0+20230104"
    # A registry takes a stamp as it stands.
    @test E.Provenancer._isregistryversion(stamped)
end

@testset "a file's registration keywords come from its record and its row" begin
    dir = mktempdir()
    # A CHELSA layer as published: named by its deposit's DOI, its source's own.
    chelsa = joinpath(dir, "CHELSA_bio1_1981-2010_V.2.1.tif")
    write(chelsa, "x")
    rec = E.datasetinfo(CHELSA{BioClim})
    r = _record(CHELSA{BioClim}, E.layerinfo(CHELSA{BioClim}, :bio1), chelsa,
                route = E.VersionOfRecord(),
                extra = Dict{String, Any}("modified" => "2023-01-04T00:00:00Z"))
    kw = E._sourcekeywords(r, rec, E.VersionOfRecord(), version = "2.1.0",
                           date = rec.releasedate)
    @test kw.doi == "10.16904/envidat.228.v2.1" && !haskey(kw, :bibtex)
    @test kw.primary === true && kw.licence == "CC0-1.0"
    @test kw.title == r.title && kw.release_version == "2.1.0"
    @test kw.release_date == "2021-06-01" &&
          kw.description == E._requestpage(rec)
    # Its stamp is dated by when the provider changed it, else when it was fetched.
    @test E.Provenancer._stampday(r, Dates.Date(2026, 10, 27)) ==
          Dates.Date(2023, 1, 4)
    @test E._registereddate(r, rec, Dates.Date(2023, 1, 4)) == "2023-01-04"
    @test E._registereddate(r, rec, nothing) == "2021-06-01"
    # A paper's file is not its source's own, whatever route it came by.
    tcr = joinpath(dir, "air.2m.mon.mean.nc")
    write(tcr, "x")
    p = _record(TwentyCR, E.layerinfo(TwentyCR, "air"), tcr,
                route = E.VersionOfRecord())
    @test p.primary === false
    @test E._sourcekeywords(p, E.datasetinfo(TwentyCR), E.VersionOfRecord(),
                            version = "3.0.0", date = nothing).primary === false
    # A row identified by its web page names its source by its citation's key, and a file stating
    # its own version is dated by when its provider changed it.
    ne = E.datasetinfo(E.NaturalEarthLevel)
    zip = joinpath(dir, "ne_10m_admin_0_countries.zip")
    write(zip, "x")
    z = _record(E.NaturalEarthLevel, nothing, zip,
                route = E.VersionOfRecord())
    @test E._needsdate(z, ne)
    zk = E._sourcekeywords(z, ne, E.VersionOfRecord(), version = "5.1.2",
                           date = nothing)
    @test zk.bibtex == ne.request && !haskey(zk, :doi) && zk.primary === true
    dated = ProvenanceRecord(role = :region, source = "NaturalEarth",
                             modified = "2022-05-13T05:20:49Z")
    @test !E._needsdate(dated, ne)
    @test E._registereddate(dated, ne, nothing) == "2022-05-13"
    # Failing that, when it was fetched; with neither, the citation's date stands.
    fetched = ProvenanceRecord(role = :region, source = "NaturalEarth",
                               fetched = "2026-09-03T00:49:00Z")
    @test E._registereddate(fetched, ne, nothing) == "2026-09-03"
    @test isnothing(E._registereddate(ProvenanceRecord(role = :region,
                                                       source = "NaturalEarth"),
                                      ne,
                                      nothing))
    # An extract says what it was made to, and its deposit is the dataset's own.
    request = only(E.era5requests("t2m", 1990:1999, dir = dir))
    write(request.path, "x")
    E.Provenancer.writesidecar(request.path,
                               E._requestrecord(request, request.path, nothing))
    x = provenance(request.path)
    era = E.datasetinfo(ERA)
    xk = E._sourcekeywords(x, era, E.GeneratedOnRequest(), version = "1.0.0",
                           date = nothing)
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
    @test entry.route isa E.OtherOrigin && entry.code == "tavg" &&
          isnothing(entry.url)
    @test E._recordgiven(spec, entry) == file
    r = provenance(file)
    @test r.primary === false && r.version == "4.8" &&
          r.releasedate == "2024-05"
    @test isnothing(r.fetched) && r.source == "CRUTS"
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
    E.Provenancer.writesidecar(path, old)
    spec = SourceSpec(TwentyCR, "air", file = path)
    entry = (; only(E._specroutes(spec))..., route = E.VersionOfRecord(),
             url = E.layerinfo(TwentyCR, "air").file)
    facts = E._currentfacts(spec, entry)
    # Read with the catalogue's facts, the record on disk is untouched.
    seen = E.Provenancer._recordwith(path, facts)
    @test !isnothing(seen.title) && seen.url == entry.url &&
          seen.primary === false
    @test isnothing(provenance(path).title)
    # Brought up to date, what the fetch recorded is kept and the writer is this version.
    r = E.Provenancer._refreshrecord!(path, facts)
    @test r.title == seen.title && r.licence == "CC0-1.0" && r.sha256 == "abc"
    @test startswith(TOML.parsefile(E.Provenancer.sidecarpath(path))["writer"],
                     "EcoSISTEM.Provenancer $(pkgversion(EcoSISTEM))")
    # Up to date already, it is not written again.
    before = mtime(E.Provenancer.sidecarpath(path))
    sleep(1.1)
    E.Provenancer._refreshrecord!(path, facts)
    @test mtime(E.Provenancer.sidecarpath(path)) == before
    @test isnothing(E.Provenancer._refreshrecord!(joinpath(dir, "none.nc"),
                                                  facts))
end

@testset "a config's missing lines and a file's reasons are written out" begin
    track = Track(nothing, "era5")
    lines = E.Provenancer._writelines(track, "ECMWF",
                                      E.Provenancer._filetype("x/era5_t2m_1990s.nc"))
    # The source's namespace for the run's inputs and outputs, since it registers that source alone.
    @test occursin("default_input_namespace: ECMWF\n  default_output_namespace: ECMWF",
                   lines)
    # The files and their records, each entry with the description DataPipeline needs of one.
    @test occursin("- data_product: era5/data/*\n  description: ", lines)
    @test occursin("  file_type: nc\n- data_product: era5/provenance/*\n  description: ",
                   lines)
    @test endswith(lines, "  file_type: toml") &&
          count("description: ", lines) == 2
    @test E.Provenancer._filetype("ne_10m_admin_0_map_units.zip") == "zip"
    # A model reads both with one entry, a `*` matching one segment, which leaves the deposit out.
    @test E.Provenancer._readentry(track, "ECMWF") ==
          "read:\n- data_product: era5/*/*\n  use:\n    namespace: ECMWF"
    reason(; present = true, held = "2.1.0", same = true, recorded = true,
           named = true) = E.Provenancer._fetchreason(present = present,
                                                      held = held,
                                                      same = same,
                                                      recorded = recorded,
                                                      named = named)
    @test isnothing(reason())
    @test startswith(reason(present = false), "not held")
    @test startswith(reason(held = nothing, named = false), "not registered")
    @test occursin("stamped version", reason(held = nothing))
    @test occursin("other facts", reason(same = false))
    @test occursin("provenance record", reason(recorded = false))
    @test occursin("no `read:` entry",
                   E.Provenancer._modelreason(nothing, track, "ECMWF"))
    @test occursin("fair pull",
                   E.Provenancer._modelreason(true, track, "ECMWF"))
    @test occursin("run the fetch run first",
                   E.Provenancer._modelreason(false, track, "ECMWF"))
    # A file registered as its source itself is the run's input, which no `write:` entry names:
    # its config writes the record alone, and its reasons say how it is registered, and whether the
    # run fetches it or it is given by hand.
    given = E.Provenancer._writelines(Track(nothing, "given"), "Mine", "csv",
                                      E.DepositIsFile())
    @test !occursin("given/data/*", given) &&
          occursin("- data_product: given/provenance/*", given) &&
          count("description: ", given) == 1
    @test E.Provenancer._writelines(track, "ECMWF", "nc",
                                    E.DepositSharesIdentifier()) == lines
    asitself(; present, held = nothing,
             fetched = true) = E.Provenancer._fetchreason(present = present,
                                                          held = held,
                                                          same = true,
                                                          recorded = true,
                                                          named = false,
                                                          deposit = E.DepositIsFile(),
                                                          fetched = fetched)
    @test asitself(present = false) ==
          "not held: the fetch run fetches and registers it as its source itself"
    @test asitself(present = false, fetched = false) ==
          "not there: a file given by hand is not fetched, so put it in place; the fetch run " *
          "registers it as its source itself"
    @test asitself(present = true) ==
          "not registered: the fetch run registers it as its source itself"
    @test reason(present = false) ==
          "not held: the fetch run fetches and registers it"
    # Primary is the route's and the identifier's, as for any file: the version of record of such
    # a source is the cited item, a copy of other origin is not.
    file = E.SourceFacts(source = "Given", key = "given", identifies = :data,
                         deposit = E.DepositIsFile())
    @test E.Provenancer.isprimary(file, E.VersionOfRecord()) &&
          !E.Provenancer.isprimary(file, E.OtherOrigin())
    # A deposit is registered from the source's facts, or from its own where it is cited apart.
    paper = E.SourceFacts(source = "Paper", doi = "10.1/p", identifies = :paper)
    apart = E.SourceFacts(file, deposit = E.DepositCitedApart(paper))
    @test E.Provenancer._depositfacts(file.deposit, file) === file &&
          E.Provenancer._depositfacts(apart.deposit, apart) === paper
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

@testset "a shape or raster file of your own enters a spec through its record" begin
    dir = mktempdir()
    geojson = joinpath(dir, "box.geojson")
    write(geojson, "{}")
    spec = ShapeSpec(geojson)
    @test sprint(show, spec) == "ShapeSpec($(repr(geojson)))"
    entry = only(E._specroutes(spec))
    @test entry == (path = geojson, route = E.OtherOrigin(), code = nothing,
           url = nothing)
    @test E._specpaths(spec) == [geojson] && E._fetchspec(spec) == [geojson]
    @test fetchfiles(spec) == [geojson] && isempty(missingfiles(spec))
    # A spec of a file no catalogue row describes has no source to register the file with: the
    # file enters a spec through its record, which names one.
    @test_throws "no catalogued source" E._trackedsource(spec)
    facts = E.Provenancer.SourceFacts(E.datasetinfo(E.NaturalEarthLevel),
                                      title = "A box")
    record = E.Provenancer.record_provenance(geojson, facts, role = :region,
                                             route = E.OtherOrigin())
    recorded = ShapeSpec(record)
    @test recorded.path === record && E._specpaths(recorded) == [geojson] &&
          isnothing(E._spectrack(recorded))
    @test sprint(show, recorded) == "ShapeSpec($(repr(record)))"
    @test only(E._specroutes(recorded)) ==
          (path = geojson, route = E.OtherOrigin(), code = nothing,
           url = nothing)
    @test only(provenance(recorded)) == record &&
          provenance(recorded) == provenance(spec) &&
          E._specinputs(recorded) == [record]
    # The record carries the route it was written with and where its sidecar is, wherever that is,
    # so a spec of it takes both from the record.
    other = joinpath(dir, "other.geojson")
    write(other, "{}")
    elsewhere = joinpath(dir, "records", "other.toml")
    moved = E.Provenancer.record_provenance(other, facts, role = :region,
                                            sidecar = elsewhere)
    @test moved.route == E.VersionOfRecord() && moved.sidecar == elsewhere
    ospec = ShapeSpec(moved)
    @test only(E._specroutes(ospec)).route == E.VersionOfRecord() &&
          only(provenance(ospec)) == moved &&
          isnothing(only(provenance(ShapeSpec(other))))
    # A record naming no file is refused.
    @test_throws "names no file" ShapeSpec(ProvenanceRecord(role = :region,
                                                            source = "x"))
    # A raster file likewise, as a layer.
    tif = joinpath(dir, "layer.tif")
    write(tif, "x")
    trecord = E.Provenancer.record_provenance(tif, facts, role = :habitat,
                                              route = E.OtherOrigin())
    rspec = RasterFileSpec(trecord, axis = E.NicheAxis)
    @test only(rspec.files) === trecord && E._specpaths(rspec) == [tif] &&
          isnothing(rspec.track)
    @test sprint(show, rspec) ==
          "RasterFileSpec($(repr(trecord)), axis = NicheAxis)"
    @test only(E._specroutes(rspec)).route == E.OtherOrigin() &&
          only(provenance(rspec)) == trecord &&
          E._readinputs(rspec) == [trecord] && E.verifyassets(rspec) == 1
    @test_throws "no catalogued source" E._trackedsource(rspec)
    @test_throws "names no file" RasterFileSpec(ProvenanceRecord(role = :habitat,
                                                                 source = "x"),
                                                axis = E.NicheAxis)
    # A file the user gives by path is recorded by a fetch run with a layer's role where it has no
    # record.
    given = joinpath(dir, "given.tif")
    write(given, "x")
    gspec = RasterFileSpec(given, axis = E.NicheAxis)
    gentry = only(E._specroutes(gspec))
    @test gentry.route == E.OtherOrigin() && isnothing(E._catalogued(gspec))
    @test E._recordgiven(gspec, gentry) == given &&
          provenance(given).role === :habitat
    @test E._catalogued(SourceSpec(WorldClim{BioClim}, :bio1)).namespace ==
          "WorldClim" &&
          E._catalogued(NaturalEarthSpec("Scotland")).dataset == "NaturalEarth"
    # The config lines for a source with no namespace of its own leave the run's defaults alone.
    lines = E.Provenancer._writelines(Track(nothing, "box"), nothing, "geojson")
    @test startswith(lines, "write:") && occursin("box/data/*", lines) &&
          occursin("file_type: geojson", lines)
    @test startswith(E.Provenancer._writelines(Track(nothing, "box"),
                                               "NaturalEarth", "geojson"),
                     "run_metadata:")
    @test occursin("the fetch run's default namespace",
                   E.Provenancer._readentry(Track(nothing, "box"), nothing))
    # Whether a file states its own version is asked of the source directly too.
    r = ProvenanceRecord(role = :region, source = "x",
                         modified = "2022-05-13T05:20:49Z")
    @test E._registereddate(r, true, nothing) == "2022-05-13" &&
          isnothing(E._registereddate(r, false, nothing))
    @test E._needsdate(ProvenanceRecord(role = :region, source = "x"), true) &&
          !E._needsdate(r, false)
    # An extract recorded with no request describes itself without one.
    @test E._extractdescription(ProvenanceRecord(role = :abundance,
                                                 source = "Occurrence Download")) ==
          "An extract of Occurrence Download."
end
end
