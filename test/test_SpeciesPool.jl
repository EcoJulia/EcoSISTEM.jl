# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/SpeciesPool.jl`: a pool of species made from a count, names or a table, standardised
# through a name authority (a made-up WFO release), merged by its rule, changed by its verbs, given
# values by `add_trait!`, read as a table, and saying how each name was matched and what it was
# read from; and the error severity its refusals are raised at. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_SpeciesPool.jl"])'

module TestSpeciesPool

using CSV
using DataFrames
using EcoSISTEM
using EcoSISTEM: AbstractMergeRule, ProblemError
using EcoSISTEM.Provenancer
using Tables
using Test

include("wfofixtures.jl")

# A pool's species' names and identifiers, in order.
names_of(pool) = DataFrame(pool).name

ids_of(pool) = DataFrame(pool).id

@testset "a pool from a count, names or a table" begin
    pool = SpeciesPool(3)
    @test length(pool) == 3
    @test names_of(pool) == ["1", "2", "3"]
    @test all(ismissing, ids_of(pool))
    @test all(DataFrame(pool).active)
    @test Tables.istable(pool)
    @test propertynames(DataFrame(pool)) ==
          [:name, :id, :active, :copyof, :copyofid]
    @test string(pool) == "SpeciesPool(3 species, 3 active)"
    @test occursin("traits    none", sprint(show, MIME"text/plain"(), pool))
    # A name given twice is one species, with nothing to report.
    @test_logs SpeciesPool(["b", "a", "b"])
    @test names_of(SpeciesPool(["b", "a", "b"])) == ["b", "a"]
    # A table's other columns are the species' values, `missing` kept where a value is unknown.
    table = DataFrame(species = ["x", "y", "z"], tmin = [1.0, 2.0, missing])
    pool = SpeciesPool(table)
    @test names_of(pool) == ["x", "y", "z"]
    @test isequal(DataFrame(pool).tmin, [1.0, 2.0, missing])
    @test_throws "no column `name`" SpeciesPool(table, column = :name)
    @test_throws "would hide the pool's own" SpeciesPool(DataFrame(species = ["x"],
                                                                   active = [true]))
    @test_throws "names no species" SpeciesPool(DataFrame(species = [
                                                              "x",
                                                              missing
                                                          ]))
    @test_throws "a Tables.jl table" SpeciesPool(1.5)
    @test_throws "`verbosity`" SpeciesPool(table, verbosity = :loud)
end

@testset "rows reaching one species are merged by the rule" begin
    table = DataFrame(species = ["x", "y", "x"], v = [1, 2, 3])
    @test RefuseMerge() isa AbstractMergeRule
    @test_throws "reached by several rows: `x` (rows 1, 3)" SpeciesPool(table)
    @test_logs (:info, r"merged: `x` \(rows 1, 3\)") SpeciesPool(table,
                                                                 merge = MergeWith(first))
    @test DataFrame(SpeciesPool(table, merge = MergeWith(first),
                                verbosity = :silent)).v == [1, 2]
    # A merge may build its row from those it is given.
    total(rows) = (v = sum(rows.v),)
    @test DataFrame(SpeciesPool(table, merge = MergeWith(total),
                                verbosity = :silent)).v == [4, 2]
    @test string(MergeWith(first)) == "MergeWith(first)"
    @test_throws "has no `v`" SpeciesPool(table,
                                          merge = MergeWith(_ -> (w = 1,)),
                                          verbosity = :silent)
    @test_throws "returned `nothing`" SpeciesPool(table,
                                                  merge = MergeWith(_ -> nothing),
                                                  verbosity = :silent)
end

@testset "a pool's table may be a file, with its record" begin
    mktempdir() do dir
        path = joinpath(dir, "envelopes.csv")
        CSV.write(path, DataFrame(species = ["x", "y"], tmin = [270.0, 280.0]))
        pool = SpeciesPool(path)
        @test names_of(pool) == ["x", "y"]
        record = only(provenance(pool))
        @test (record.role, record.path) == (:species, "envelopes.csv")
        @test unique(name_table(pool).source) == ["envelopes.csv"]
        # A record of the file is kept as it is.
        facts = SourceFacts(source = "Envelopes", doi = "10.1234/envelopes",
                            title = "Envelopes", licence = "CC-BY-4.0",
                            identifies = :data)
        given = record_provenance(path, facts, role = :species)
        @test only(provenance(SpeciesPool(given))).doi == "10.1234/envelopes"
        # A table in memory has none, unless one is given.
        @test isempty(provenance(SpeciesPool(DataFrame(species = ["x"]))))
        @test only(provenance(SpeciesPool(DataFrame(species = ["x"]),
                                          provenance = given))) == given
        # A file's names are read as text, so names made of digits keep their spelling, in whichever
        # column holds them.
        numbered = joinpath(dir, "numbered.csv")
        CSV.write(numbered,
                  DataFrame(code = ["007", "1e3"], tmin = [270.0, 280.0]))
        @test names_of(SpeciesPool(numbered, column = :code,
                                   verbosity = :silent)) == ["007", "1e3"]
    end
end

@testset "an authority standardises the pool's names" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        table = DataFrame(species = ["Acacia alpha", "Mimosa alpha",
                              "Acacia delta",
                              "Mimosa delta", "Acacia eta", "Algae x",
                              "Acacia zeta"],
                          v = 1:7)
        pool = @test_logs((:warn,
                           r"2 names discarded.*1 only unchecked, 1 not in it.*`Acacia eta`, `Algae x`"),
                          (:info, r"2 synonyms translated"),
                          (:info, r"1 homonym kept"),
                          (:info,
                           r"2 species were reached by several rows, merged"),
                          SpeciesPool(table, authority = wfo,
                                      merge = MergeWith(first)))
        @test names_of(pool) == ["Acacia alpha", "Acacia delta", "Acacia zeta"]
        @test ids_of(pool) ==
              [
            "wfo-0000000001",
            "wfo-0000000006+wfo-0000000007",
            "wfo-0000000010"
        ]
        @test DataFrame(pool).v == [1, 3, 7]
        # Every name the table used, with what it reached and how.
        nt = name_table(pool)
        @test nt.name == sort(table.species)
        row(n) = only(eachrow(nt[nt.name .== n, :]))
        @test (row("Mimosa alpha").species, row("Mimosa alpha").match,
               row("Mimosa alpha").merged) == ("Acacia alpha", :synonym, true)
        @test ismissing(row("Algae x").species) &&
              row("Algae x").match === :not_found
        @test row("Acacia eta").match === :unchecked_only
        @test !row("Acacia zeta").merged
        # The release's record is among the pool's inputs.
        @test any(r -> r.doi == "10.5281/zenodo.20782718", provenance(pool))
        # Refused by default, the merge names the species and its rows.
        @test_throws "`Acacia alpha` (rows 1, 2)" SpeciesPool(table,
                                                              authority = wfo,
                                                              verbosity = :silent)
        # Unchecked names kept where they are allowed.
        lax = WFOPlantList(release = "2026-06", path = zip,
                           allow_unchecked = true)
        @test "Acacia eta" in
              names_of(SpeciesPool(["Acacia eta"], authority = lax,
                                   verbosity = :silent))
        # Rows a homonym's authorship splits, reaching one name, are one species: merged, its
        # identifier joining both.
        split = DataFrame(species = ["Acacia delta", "Acacia delta"],
                          authors = ["A", "B"], v = [1, 2])
        @test_throws "`Acacia delta` (rows 1, 2)" SpeciesPool(split,
                                                              authority = wfo,
                                                              authorship = :authors,
                                                              verbosity = :silent)
        merged = SpeciesPool(split, authority = wfo, authorship = :authors,
                             merge = MergeWith(first), verbosity = :silent)
        @test ids_of(merged) == ["wfo-0000000006+wfo-0000000007"]
        @test name_table(merged).authorship == ["A", "B"]
    end
end

@testset "the verbs change a pool before the build" begin
    pool = SpeciesPool(["a", "b", "c", "d"])
    @test_logs (:info, r"2 species deactivated; 2 of 4 active") deactivate_species!(pool,
                                                                                    "b",
                                                                                    ["c"])
    @test DataFrame(pool).active == [true, false, false, true]
    @test_logs (:info, r"1 species activated") activate_species!(pool, "c")
    @test DataFrame(pool).active == [true, false, true, true]
    # No names means every species, either way.
    deactivate_species!(pool, verbosity = :silent)
    @test !any(DataFrame(pool).active)
    @test_logs (:info, r"4 species activated; 4 of 4 active") activate_species!(pool)
    @test all(DataFrame(pool).active)
    @test_logs (:warn, r"could not deactivate `q`") match_mode = :any deactivate_species!(pool,
                                                                                          "q")
    @test_throws "needs the pool's authority" activate_species!(pool, "a",
                                                                rank = "genus")
    # A species chosen by a function of its row, which must answer `true` or `false`.
    deactivate_species!(pool, verbosity = :silent) do row
        return row.name in ("a", "d")
    end
    @test DataFrame(pool).active == [false, true, true, false]
    @test_throws "returned missing for `a`" activate_species!(_ -> missing,
                                                              pool)
    activate_species!(pool, verbosity = :silent)
    # Named species are removed, whether active or not; with no names nothing is, unless only the
    # inactive are asked for. Every other species keeps its key.
    @test_throws "would remove every species" remove_species!(pool)
    @test_logs (:info, r"1 species removed; 3 remain") remove_species!(pool,
                                                                       "b")
    @test names_of(pool) == ["a", "c", "d"]
    @test pool.species.key == [1, 3, 4]
    deactivate_species!(pool, "d", verbosity = :silent)
    @test_logs (:info, r"kept 1 active species: `c`") (:info,
                                                       r"1 species removed") remove_species!(pool,
                                                                                             "c",
                                                                                             "d",
                                                                                             keep_active = true)
    @test names_of(pool) == ["a", "c"]
    deactivate_species!(pool, "c", verbosity = :silent)
    remove_species!(pool, keep_active = true, verbosity = :silent)
    @test names_of(pool) == ["a"] && pool.species.key == [1]
    # A pool is never left empty, however the species to remove are given.
    @test_throws "leave the pool with no species" remove_species!(pool, "a")
    deactivate_species!(pool, verbosity = :silent)
    @test_throws "leave the pool with no species" remove_species!(pool,
                                                                  keep_active = true)
    @test_throws "leave the pool with no species" remove_species!(_ -> true,
                                                                  pool)
    @test length(pool) == 1
    # A copy has a key after the last row, its own label, and its original's values and name.
    pool = SpeciesPool(DataFrame(species = ["a", "b"], v = [1, 2]))
    copy_species!(pool, "b", as = "b copy")
    copy_species!(pool, "b copy", as = "b copy again")
    table = DataFrame(pool)
    @test table.name == ["a", "b", "b copy", "b copy again"]
    @test table.v == [1, 2, 2, 2]
    @test isequal(table.copyof, [missing, missing, "b", "b"])
    @test pool.species.key == [1, 2, 3, 4]
    @test isequal(pool.species.copykey, [missing, missing, 2, 2])
    @test_throws "already names a species" copy_species!(pool, "a", as = "b")
    @test_throws "needs the copy's label" copy_species!(pool, "a")
    @test_throws "no species `q`" copy_species!(pool, "q", as = "q2")
    # A copy is matched by its label.
    deactivate_species!(pool, "b copy again", verbosity = :silent)
    @test DataFrame(pool).active == [true, true, true, false]
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        pool = SpeciesPool(["Acacia alpha", "Acacia zeta"], authority = wfo,
                           verbosity = :silent)
        # A synonym finds its species; a recognised species not in the pool is a notice; a name
        # the authority does not know, a warning.
        deactivate_species!(pool, "Mimosa alpha", verbosity = :silent)
        @test DataFrame(pool).active == [false, true]
        @test_logs (:info, r"`Acacia gamma`: recognised") match_mode = :any deactivate_species!(pool,
                                                                                                "Acacia gamma")
        @test_logs (:warn, r"`Algae x`.*does not recognise") match_mode = :any deactivate_species!(pool,
                                                                                                   "Algae x")
        # A copy's label may not be a name the authority recognises.
        @test_throws "recognises" copy_species!(pool, "Acacia alpha",
                                                as = "Acacia gamma")
        copy_species!(pool, "Acacia alpha", as = "alpha copy")
        table = DataFrame(pool)
        @test isequal(table.id, ["wfo-0000000001", "wfo-0000000010", missing])
        @test table.copyofid[3] == "wfo-0000000001"
        activate_species!(pool, "alpha copy", verbosity = :silent)
        @test last(DataFrame(pool).active)
    end
end

@testset "species are chosen by the group a rank names" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        pool = SpeciesPool([
                               "Acacia alpha",
                               "Acacia zeta",
                               "Acacia delta",
                               "Mimosa beta"
                           ],
                           authority = wfo, verbosity = :silent)
        copy_species!(pool, "Acacia alpha", as = "alpha copy")
        # The genus Acacia, named by a synonym: its species and the copy of one; not the homonym
        # whose two species are in different genera, nor the species of the other genus. "Acacia"
        # is a section as well as a genus, and the rank asked for picks the genus.
        deactivate_species!(pool, verbosity = :silent)
        activate_species!(pool, "Racosperma", rank = "genus",
                          verbosity = :silent)
        @test DataFrame(pool).active == [true, true, false, false, true]
        deactivate_species!(pool, "Acacia", rank = "genus", verbosity = :silent)
        @test !any(DataFrame(pool).active)
        # The family holds every species, the homonym's two included.
        activate_species!(pool, "Fabaceae", rank = "family",
                          verbosity = :silent)
        @test all(DataFrame(pool).active)
        # A name found only at another rank, or not at all, is warned of and chooses nothing.
        @test_logs (:warn, r"`Fabaceae` \(a family\): not a genus") match_mode = :any deactivate_species!(pool,
                                                                                                          "Fabaceae",
                                                                                                          rank = "genus")
        @test_logs (:warn, r"`Algae`: .*does not recognise") match_mode = :any deactivate_species!(pool,
                                                                                                   "Algae",
                                                                                                   rank = "genus")
        @test all(DataFrame(pool).active)
        # Keep one genus and remove the rest.
        deactivate_species!(pool, verbosity = :silent)
        activate_species!(pool, "Mimosa", rank = "genus", verbosity = :silent)
        remove_species!(pool, keep_active = true, verbosity = :silent)
        @test names_of(pool) == ["Mimosa beta"]
        # A rank removes a group outright.
        pool = SpeciesPool(["Acacia alpha", "Mimosa beta"], authority = wfo,
                           verbosity = :silent)
        remove_species!(pool, "Acacia", rank = "genus", verbosity = :silent)
        @test names_of(pool) == ["Mimosa beta"]
    end
end

@testset "values are added by species, never by position" begin
    pool = SpeciesPool(["x", "y", "z"])
    @test add_trait!(pool, DataFrame(species = ["z", "x", "q"], w = [3, 1, 9]),
                     verbosity = :silent) === pool
    @test isequal(DataFrame(pool).w, [1, missing, 3])
    @test_throws "already has `w`" add_trait!(pool,
                                              DataFrame(species = ["x"],
                                                        w = [0]))
    add_trait!(pool, DataFrame(species = ["x"], w = [0]), replace = true,
               verbosity = :silent)
    @test isequal(DataFrame(pool).w, [0, missing, missing])
    @test_throws "column of the pool's own" add_trait!(pool,
                                                       DataFrame(species = ["x"],
                                                                 name = ["n"]))
    # A refused merge leaves the pool as it was.
    @test_throws "reached by several rows" add_trait!(pool,
                                                      DataFrame(species = [
                                                                    "x",
                                                                    "x"
                                                                ],
                                                                u = [1, 2]))
    @test !(:u in propertynames(DataFrame(pool)))
    @test_logs (:warn, r"`q`") add_trait!(pool,
                                          DataFrame(species = ["q"], s = [1]))
    # A value computed from each species' row, with a do-block.
    add_trait!(pool, :double) do row
        return coalesce(row.w, 0) * 2
    end
    @test DataFrame(pool).double == [0, 0, 0]
    @test_throws "already has `double`" add_trait!(row -> 1, pool, :double)
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        pool = SpeciesPool(["Acacia alpha", "Acacia zeta"], authority = wfo,
                           verbosity = :silent)
        # A table naming a synonym reaches the species, and its names join the name table.
        widths = DataFrame(name = ["Mimosa alpha", "Acacia gamma"],
                           width = [2.0, 3.0])
        @test_logs (:info, r"1 name of the table recognised, but of no species") add_trait!(pool,
                                                                                            widths,
                                                                                            column = :name,
                                                                                            role = :dispersal)
        @test isequal(DataFrame(pool).width, [2.0, missing])
        nt = name_table(pool)
        @test nrow(nt) == 4
        @test only(nt[nt.name .== "Acacia gamma", :match]) === :not_in_pool
        @test only(nt[nt.name .== "Mimosa alpha", :species]) == "Acacia alpha"
    end
end

@testset "a problem raised at error level stops" begin
    @test EcoSISTEM._severitytag(ProblemError()) == "error"
    @test_throws "cannot go on" EcoSISTEM._report(ProblemError(),
                                                  "cannot go on")
end

end
