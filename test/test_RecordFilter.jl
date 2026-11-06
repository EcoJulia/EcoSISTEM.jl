# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests `src/RecordFilter.jl`: `IsActiveSpecies` keeping the rows of a table that name an active species
# of a pool or a species list, in place or into a new table, matching names through a name authority (a
# made-up WFO release), and its report - each kind of name dropped, its level and its detail. Run with
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_RecordFilter.jl"])'

module TestRecordFilter

using CSV
using DataFrames
using EcoSISTEM
using EcoSISTEM: AbstractRecordFilter, ProblemWarning, ProblemError
using EcoSISTEM.Units
using Test
using Unitful
using Unitful.DefaultSymbols

include("wfofixtures.jl")

# The species of `pool`, every one alike.
function _species(pool::SpeciesPool)
    return build_species(pool,
                         tolerance = ToleranceSpec(Temperature, μ = 298.0K,
                                                   σ = 2.0K),
                         demand = DemandSpec(SolarRadiation, 1.0e-4kJ / day),
                         seed = 1, verbosity = :silent)
end

# A pool of three species, the third not active, and a copy of the first.
function _pool()
    pool = SpeciesPool(["a", "b", "c"], verbosity = :silent)
    deactivate_species!(pool, "c", verbosity = :silent)
    return copy_species!(pool, "a", as = "a copy")
end

@testset "a filter keeps the rows naming an active species, and reports the rest" begin
    pool = _pool()
    records = DataFrame(species = ["a", "b", "c", "zz", missing, "a copy", "zz",
                            "a"],
                        n = 1:8)
    kept = @test_logs((:info,
                       r"^2 records of 1 name of no species in the pool were removed: `zz` \(2\)\.$"),
                      (:info,
                       r"^1 record of 1 species not active was removed: `c` \(1\)\. `activate_species!` activates them\.$"),
                      (:info, r"^1 record naming no species was removed\.$"),
                      filter!(IsActiveSpecies(pool), records))
    @test kept === records
    @test records.n == [1, 2, 6, 8]
    @test IsActiveSpecies(pool) isa AbstractRecordFilter
    # A species list filters as the pool it was built from does, naming the way to activate a species.
    species = _species(pool)
    records = DataFrame(species = ["a", "c", "a copy", "zz"])
    @test_logs((:info, r"of no species in the species list was removed"),
               (:info,
                r"Activate them in their pool with `activate_species!` and build the species again"),
               filter!(IsActiveSpecies(species), records))
    @test records.species == ["a", "a copy"]
end

@testset "filter returns a new table, of the kind it was given" begin
    pool = _pool()
    quiet = IsActiveSpecies(pool, verbosity = :silent)
    records = DataFrame(species = ["a", "zz"], n = [1, 2])
    new = @test_logs filter(quiet, records)
    @test new.n == [1] && records.n == [1, 2]
    @test filter(quiet, (species = ["zz", "b"], n = [1, 2])) ==
          (species = ["b"], n = [2])
    mktempdir() do dir
        path = joinpath(dir, "records.csv")
        CSV.write(path, records)
        file = CSV.File(path)
        filtered = filter(quiet, file)
        @test filtered isa DataFrame
        @test filtered.species == ["a"] && filtered.n == [1]
        @test_throws "cannot be changed in place" filter!(quiet, file)
    end
end

@testset "a report is raised at its level, in more detail at `:verbose`" begin
    pool = _pool()
    records = DataFrame(species = vcat(["x$i" for i in 1:12], "a"))
    # The first ten names at `:normal`, the heaviest first and then by name.
    @test_logs((:warn,
                r"^12 records of 12 names .*: `x1` \(1\), `x10` \(1\), .*`x7` \(1\) and 2 more\.$"),
               filter(IsActiveSpecies(pool, level = ProblemWarning()), records))
    @test_logs((:info, r"`x8` \(1\), `x9` \(1\)\.$"),
               filter(IsActiveSpecies(pool, verbosity = :verbose), records))
    heavy = DataFrame(species = ["y", "x", "y"])
    @test_logs((:info, r": `y` \(2\), `x` \(1\)\.$"),
               filter(IsActiveSpecies(pool), heavy))
    # An error stops the filter before the table changes.
    @test_throws "12 records of 12 names" filter!(IsActiveSpecies(pool,
                                                                  level = ProblemError()),
                                                  records)
    @test nrow(records) == 13
    # The column naming the species, and what is refused.
    @test filter(IsActiveSpecies(pool, column = :taxon, verbosity = :silent),
                 DataFrame(taxon = ["a", "zz"])).taxon == ["a"]
    @test_throws "names no species in a `species` column; it has `taxon`" filter(IsActiveSpecies(pool),
                                                                                 DataFrame(taxon = ["a"]))
    @test_throws "row 1 of a table names its species as `3`" filter(IsActiveSpecies(pool),
                                                                    DataFrame(species = [3]))
    @test_throws "`verbosity` must be" IsActiveSpecies(pool, verbosity = :loud)
end

@testset "names are found through the pool's name authority" begin
    withfakewfo() do zip
        wfo = WFOPlantList(release = "2026-06", path = zip)
        pool = SpeciesPool(["Acacia alpha", "Acacia zeta"], authority = wfo,
                           verbosity = :silent)
        # A synonym reaches its accepted species, an identifier its species; a name the authority
        # recognises but the pool does not hold is reported apart from one it does not recognise.
        given = ["Mimosa alpha", "Acacia alpha", "wfo-0000000010",
            "Acacia theta", "Acacia theta", "Algae x"]
        for holder in (pool, _species(pool))
            records = DataFrame(species = given)
            @test_logs((:info,
                        r"^1 record of 1 name not recognised by the WFO Plant List 2026-06 was removed: `Algae x` \(1\)\.$"),
                       (:info,
                        r"^2 records of 1 name recognised by the WFO Plant List 2026-06 but of no species in the (pool|species list) were removed: `Acacia theta` \(2\)\.$"),
                       filter!(IsActiveSpecies(holder), records))
            @test records.species == given[1:3]
        end
        # Names all found directly - an own name, an identifier - never read the authority's names.
        empty!(EcoSISTEM._WFOLOADED)
        @test filter(IsActiveSpecies(pool),
                     DataFrame(species = ["Acacia alpha", "wfo-0000000010"])).species ==
              ["Acacia alpha", "wfo-0000000010"]
        @test isempty(EcoSISTEM._WFOLOADED)
        filter(IsActiveSpecies(pool), DataFrame(species = ["Mimosa alpha"]))
        @test !isempty(EcoSISTEM._WFOLOADED)
    end
end

@testset "a filter shows as it is written" begin
    pool = _pool()
    shown = sprint(show, pool)
    @test sprint(show, IsActiveSpecies(pool)) == "IsActiveSpecies($shown)"
    @test sprint(show,
                 IsActiveSpecies(pool, column = :taxon,
                                 level = ProblemWarning(),
                                 verbosity = :verbose)) ==
          "IsActiveSpecies($shown, column = :taxon, level = ProblemWarning(), verbosity = :verbose)"
    plain = sprint(show, MIME"text/plain"(), IsActiveSpecies(pool))
    @test startswith(plain, "IsActiveSpecies\n") &&
          occursin("level     ProblemNotice()", plain)
end

end
