# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/NameAuthority.jl`: the WFO Plant List as a name authority - which release it takes, the
# rules a name resolves to a species by, the classification above a species, and the release's file
# fetched, recorded and cached - on a made-up release written for the test, and on the real 2026-06
# release where the network allows. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_NameAuthority.jl"])'

module TestNameAuthority

using EcoSISTEM
using EcoSISTEM: resolve_name, lineage, fetchfiles, missingfiles
using Test

include("wfofixtures.jl")

# Whether the test reading the real release runs: never on a CI runner, and not when switched off.
# The release is 122 MB, fetched into the asset cache once.
function wforelease()
    haskey(ENV, "RUNNER_OS") && return false
    return get(ENV, "ECOSISTEM_WFO_FETCH", "true") != "false"
end

@testset "a WFO Plant List names a catalogued release" begin
    @test_throws "catalogues release `2026-06`" WFOPlantList()
    @test_throws "`2025-12` is not catalogued" WFOPlantList(release = "2025-12")
    wfo = WFOPlantList(release = "2026-06")
    @test (wfo.release, wfo.allow_unchecked, wfo.prefer_accepted) ==
          ("2026-06", false, false)
    @test string(wfo) == "WFOPlantList(release = \"2026-06\")"
    path = joinpath("data", "wfo.zip")
    @test string(WFOPlantList(release = "2026-06", allow_unchecked = true,
                              prefer_accepted = true, path = path)) ==
          "WFOPlantList(release = \"2026-06\", allow_unchecked = true, " *
          "prefer_accepted = true, path = $(repr(path)))"
    @test wfo isa EcoSISTEM.AbstractNameAuthority
    # Its file is the release's Darwin Core file, from the Zenodo record the catalogued DOI names,
    # under the release's own directory of the asset cache.
    asset = EcoSISTEM._wfoasset(wfo)
    @test asset.url ==
          "https://zenodo.org/records/20782718/files/_DwC_backbone_R.zip"
    @test EcoSISTEM._localpath(asset) ==
          joinpath(EcoSISTEM.assetdir(owner = WFOPlantList), "2026-06",
                   "_DwC_backbone_R.zip")
end

@testset "names resolve to species by the release's rules" begin
    withfakewfo() do zip
        for au in (false, true), pa in (false, true)
            wfo = WFOPlantList(release = "2026-06", path = zip,
                               allow_unchecked = au, prefer_accepted = pa)
            for (name, expected) in FAKEWFO_EXPECTED[(au, pa)]
                @test resolve_name(wfo, name) == expected
            end
        end
    end
end

@testset "an authorship picks among a homonym's rows" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        # Mill.'s "Acacia beta" is a synonym of Acacia gamma; L.'s is accepted.
        @test resolve_name(wfo, "Acacia beta", authorship = "Mill.") ==
              (name = "Acacia gamma", id = "wfo-0000000005", match = :synonym)
        @test resolve_name(wfo, "Acacia beta", authorship = "L.") ==
              (name = "Acacia beta", id = "wfo-0000000003", match = :accepted)
        # Spaces are trimmed and collapsed before comparing.
        @test resolve_name(wfo, "Acacia delta", authorship = "  B ").id ==
              "wfo-0000000007"
        # An authorship no row has leaves the homonym as it is, and one is not read for a name
        # that is no homonym.
        @test resolve_name(wfo, "Acacia delta", authorship = "C").id ==
              "wfo-0000000006+wfo-0000000007"
        @test resolve_name(wfo, "Acacia alpha", authorship = "C").id ==
              "wfo-0000000001"
    end
end

@testset "a lineage walks up to the kingdom" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        path = lineage(wfo, "wfo-0000000001")
        @test [l.rank for l in path] ==
              ["species", "genus", "family", "order", "kingdom"]
        @test [l.name for l in path] ==
              ["Acacia alpha", "Acacia", "Fabaceae", "Fabales", "Plantae"]
        @test first(path).id == "wfo-0000000001"
        # A synonym's identifier walks from its accepted species.
        @test lineage(wfo, "wfo-0000000002") == path
        # Joined identifiers share what their species share: the genus where both are Acacia,
        # the family where one is Mimosa.
        @test first(lineage(wfo, "wfo-0000000001+wfo-0000000005")).name ==
              "Acacia"
        @test [l.name for l in lineage(wfo, "wfo-0000000006+wfo-0000000007")] ==
              ["Fabaceae", "Fabales", "Plantae"]
        @test_throws "no identifier" lineage(wfo, "wfo-9999999999")
    end
end

@testset "the release's file is recorded, its names cached and read once" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        @test isempty(missingfiles(wfo))
        @test only(fetchfiles(wfo)) == zip
        resolve_name(wfo, "Acacia alpha")
        cache = EcoSISTEM._wfocachepath(zip)
        @test isfile(cache)
        # The table cached holds the rows the reader keeps, Unchecked ones too, so one table serves
        # every setting.
        @test countlines(cache) == length(FAKEWFO_ROWS) + 1
        @test EcoSISTEM._wfonames(wfo) === EcoSISTEM._wfonames(wfo)
        # A file changed since its table was cached is read afresh, and the older table removed, so a
        # file keeps one however often the reader or the file changes.
        touch(zip)
        resolve_name(wfo, "Acacia alpha")
        @test !isfile(cache) && isfile(EcoSISTEM._wfocachepath(zip))
        record = only(provenance(wfo))
        @test record.role === :species
        @test record.source == "WFOPlantList"
        @test record.doi == "10.5281/zenodo.20782718"
        @test record.version == "2026-06"
        # The record of a file found present, not fetched, has no fetch time.
        @test isnothing(record.fetched)
    end
    # A release not yet held is a file to do, and a dry run says where it would land.
    mktempdir() do dir
        zip = joinpath(dir, "_DwC_backbone_R.zip")
        wfo = WFOPlantList(release = "2026-06", path = zip)
        @test only(missingfiles(wfo)).file == zip
        @test isnothing(only(provenance(wfo)))
    end
end

if wforelease()
    @testset "the 2026-06 release" begin
        wfo = WFOPlantList(release = "2026-06")
        t = EcoSISTEM._wfonames(wfo)
        # The counts read when the release was catalogued.
        @test length(t.id) == 1_664_014
        @test count(==(:Accepted), t.status) == 454_635
        @test count(==(:Synonym), t.status) == 1_026_257
        @test count(==(:Unchecked), t.status) == 183_122
        @test resolve_name(wfo, "Acacia binervata") ==
              (name = "Acacia binervata", id = "wfo-0000185290+wfo-0000202067",
               match = :homonym)
        @test resolve_name(WFOPlantList(release = "2026-06",
                                        prefer_accepted = true),
                           "Acacia binervata").id == "wfo-0000185290"
        @test [(l.rank, l.name) for l in lineage(wfo, "wfo-0000185290")] ==
              [("species", "Acacia binervata"), ("genus", "Acacia"),
            ("subfamily", "Caesalpinioideae"), ("family", "Fabaceae"),
            ("order", "Fabales"), ("phylum", "Angiosperms"),
            ("subkingdom", "Pteridobiotina"), ("kingdom", "Plantae")]
        record = only(provenance(wfo))
        @test record.doi == "10.5281/zenodo.20782718"
        @test !isnothing(record.sha256)
    end
else
    @info "WFO Plant List release not read: on a CI runner, or ECOSISTEM_WFO_FETCH=false"
end

end
