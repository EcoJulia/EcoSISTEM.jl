# SPDX-License-Identifier: LGPL-3.0-or-later

module TestDeprecations

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Test
using Unitful
using Unitful.DefaultSymbols

# Coverage for `src/deprecations.jl` and the dataset-typed shims in the RasterDataSources
# extension, which holds none at present. A shim is checked on both halves: it warns and its
# result matches the current API it forwards to. A method's warning is a log record
# (`@test_deprecated`); a deprecated binding's is printed by the runtime to stderr, so it is caught
# with `@test_warn`, and the access is evaluated inside it rather than written in the test's own
# code, which on Julia 1.11 would warn when that code is compiled, before anything is captured.

# A kernel's positional form builds the very kernel its keyword form does, and so the same lookup.
# That both match the released code to the bit is for the canonical tests, which still build their
# kernels positionally.
@testset "kernels' positional forms forward to their keyword forms" begin
    regime = build_habitat(DefaultEcosystem(), verbosity = :silent).regime
    function samelookup(a, b)
        la, lb = EcoSISTEM.genlookups(regime, a),
                 EcoSISTEM.genlookups(regime, b)
        return (la.y, la.x, la.p) == (lb.y, lb.x, lb.p)
    end
    old = @test_deprecated GaussianKernel(2.4u"km", 1.0e-9)
    new = GaussianKernel(meandist = 2.4u"km", thresh = 1.0e-9)
    @test old == new && samelookup(old, new)
    old = @test_deprecated LongTailKernel(2.4u"km", 2.0, 1.0e-6)
    new = LongTailKernel(scale = (2 * 2.4u"km") / sqrt(pi), shape = 2.0,
                         thresh = 1.0e-6)
    @test old == new && samelookup(old, new)
end

# A list-wide movement is the vector of its species' movements, each in the mode the old type
# named, with its own kernel, the kernel's type as the family and its own flag as the dead-cell
# fate.
@testset "list-wide movements forward to each species' movement" begin
    kernels = [GaussianKernel(meandist = 1.0u"km"),
        GaussianKernel(meandist = 2.0u"km")]
    @test (@test_deprecated BirthOnlyMovement(kernels)) ==
          SeedDispersal(GaussianKernel, disperse_safely = true).(kernels)
    @test (@test_deprecated BirthOnlyMovement(kernels, [true, false])) ==
          [SeedDispersal(GaussianKernel, disperse_safely = true)(kernels[1]),
        SeedDispersal(GaussianKernel, disperse_safely = false)(kernels[2])]
    @test (@test_deprecated AlwaysMovement(kernels)) ==
          AlwaysDispersal().(kernels)
    @test (@test_deprecated NoMovement(kernels)) == fill(NoDispersal()(), 2)
    long = [LongTailKernel(scale = 1.0u"km", shape = 2.0)]
    @test (@test_deprecated BirthOnlyMovement(long)) ==
          SeedDispersal(LongTailKernel, disperse_safely = true).(long)
    @test_throws "one entry per species" EcoSISTEM._listmovement(SeedDispersal,
                                                                 kernels,
                                                                 [true])
end

# A form that took a species list's starting abundances, or read them from one, raises an error saying
# how to rewrite the call; `build_species(numspecies; ...)`'s is written from the keywords it was given.
@testset "forms reading starting abundances are removed, each saying how to rewrite it" begin
    rewrite(f) =
        try
            f()
            ""
        catch err
            sprint(showerror, err)
        end
    message = rewrite() do
        return build_species(5, tolerance = (290.0u"K", 2.0u"K"),
                             toleranceaxis = Temperature,
                             demand = 1.0e9u"kJ/d", demandaxis = SolarRadiation,
                             dispersal = 5.0u"km", abundance = 500, seed = 1)
    end
    @test occursin("`build_species(5; ...)` was removed in v0.11.0", message)
    @test occursin("build_species(SpeciesPool(5); tolerance = ToleranceSpec(Temperature, " *
                   "μ = 290.0 K, σ = 2.0 K)", message)
    @test occursin("movement = MovementSpec(InsectDispersal(), meandist = 5.0 km), seed = 1)",
                   message)
    @test occursin("build_ecosystem(BySupply(total = 500), species, habitat)",
                   message)
    message = rewrite() do
        return build_species(2, tolerance = ([290.0, 291.0]u"K", 2.0u"K"),
                             toleranceaxis = Temperature, demand = 1.0u"kJ/d",
                             demandaxis = SolarRadiation, abundance = [3, 4],
                             native = [true, false])
    end
    @test occursin("μ = TraitColumn(:μ)", message)
    @test occursin("BySupply(each = row -> abundance[row.key])", message)
    @test occursin("deactivate_species!", message)
    # The `DefaultEcosystem()` shorthand builds from a pool, refusing what only the removed form took.
    @test_throws "builds species from a `SpeciesPool` since v0.11.0" build_species(DefaultEcosystem(),
                                                                                   abundance = 10,
                                                                                   verbosity = :silent)
    @test_throws "builds species from a `SpeciesPool` since v0.11.0" build_species(DefaultEcosystem(),
                                                                                   movement = BirthOnlyMovement,
                                                                                   verbosity = :silent)
    species = build_species(DefaultEcosystem(), numspecies = 3,
                            verbosity = :silent)
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    eco = build_ecosystem(NoAbundance(), species, habitat, verbosity = :silent)
    @test_throws "`abun` and `native` arguments were removed" SpeciesList(3,
                                                                          species.tolerance,
                                                                          [
                                                                              1,
                                                                              2,
                                                                              3
                                                                          ],
                                                                          species.demand,
                                                                          species.movement,
                                                                          species.params,
                                                                          fill(true,
                                                                               3))
    @test_throws "`Ecosystem(spplist, habitat, nichefit)` was removed" Ecosystem(species,
                                                                                 habitat,
                                                                                 eco.nichefit)
    @test_throws "`Ecosystem(popfun, spplist, habitat, nichefit)` was removed" Ecosystem(identity,
                                                                                         species,
                                                                                         habitat,
                                                                                         eco.nichefit)
    @test_throws "`MPIEcosystem(spplist, habitat, nichefit)` was removed" MPIEcosystem(species,
                                                                                       habitat,
                                                                                       eco.nichefit)
    for (populate, placement) in ((populate!, "BySupply()"),
        (populate_by_tolerance!, "BySuitability()"),
        (emptypopulate!, "NoAbundance()"))
        @test_throws "`$placement` places as" populate(eco.abundances, species,
                                                       habitat, eco.nichefit,
                                                       eco.rngs)
    end
end

# Placing again by suitability is `repopulate!` with the suitability placement.
@testset "repopulate_by_tolerance! places again by suitability" begin
    species = build_species(DefaultEcosystem(), numspecies = 3,
                            verbosity = :silent)
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    eco = build_ecosystem(BySupply(total = 300), species, habitat, seed = 1,
                          verbosity = :silent)
    @test (@test_deprecated repopulate_by_tolerance!(eco)) === eco
    @test sum(eco.abundances.matrix) == 3000
    @test eco.placement == BySuitability()
end

@testset "InputRecord is ProvenanceRecord" begin
    @test Base.isdeprecated(EcoSISTEM, :InputRecord)
    @test (@test_warn "deprecated" Core.eval(@__MODULE__,
                                             :(EcoSISTEM.InputRecord))) ===
          EcoSISTEM.ProvenanceRecord
    @test Base.ispublic(EcoSISTEM, :InputRecord) &&
          !Base.isexported(EcoSISTEM, :InputRecord)
end

# Shims removed at a breaking release must leave nothing behind: a stale `export` or `public` of a
# deleted binding is only a precompile warning, and a method left on a function that stays is
# silent.
# A small ecosystem on a 10 x 10 grid, every cell holding 100 of any species at capacity.
function _eco(; abundance = 3000)
    area = StudyArea(extent = (10.0km, 10.0km), cellsize = 1.0km,
                     verbosity = :silent)
    env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                      supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                           axis = SolarRadiation),
                      area = area)
    spp = build_species(SpeciesPool(3),
                        tolerance = ToleranceSpec(Temperature, μ = 285.0K,
                                                  σ = 5.0K),
                        demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                        seed = 1, verbosity = :silent)
    return build_ecosystem(BySupply(total = abundance), spp, env, seed = 1,
                           verbosity = :silent)
end

# A deprecated call's result, its warning kept out of the test output.
_quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

# The table operation is applied by hand here, at chosen elapsed times, so each step's additions can
# be read off the abundance matrix exactly: the dynamics never run between them.
const MONTH = 1.0month_mean_duration

function _applytable!(eco, op, elapsed; timestep = MONTH, region = AllCells())
    return EcoSISTEM.applyinterventions!(eco,
                                         Intervention(EveryStep(), region, op),
                                         elapsed, timestep, 0)
end

# Rows out of time order, naming species by index and by name and cells by linear index and by
# `CartesianIndex`, with two fractional rows of one species and cell in the same step.
function _abundancerows(eco)
    return (species = Any[2, 1, 1, eco.spplist.names[3], 1],
            cell = Any[5, CartesianIndex(3, 4), 5, 7, 5],
            count = [2, 3, 0.3, 0.4, 0.3],
            time = [2.0, 0.0, 1.0, 1.0, 1.0] .* month_mean_duration)
end

# What each of the elapsed times adds, as the change in the abundance matrix after each.
function _tableadditions(eco, op, elapsed; kw...)
    before = copy(eco.abundances.matrix)
    return map(elapsed) do t
        _applytable!(eco, op, t; kw...)
        return eco.abundances.matrix - before
    end
end

# The table operation forwards to a dated placement, converting what only it takes - a species by its
# position, a cell by its linear index, fractional counts - so the released fixture's additions are
# unchanged: its two fractional rows of one species and cell share a time, and sum to one.
@testset "AddAbundanceTable forwards to AddAbundance of a dated placement" begin
    eco = _eco()
    op = @test_deprecated AddAbundanceTable(_abundancerows(eco))
    @test sprint(show, op) == "AddAbundance(the rows of an `AddAbundanceTable`)"
    start, one, two = _tableadditions(eco, op, (0.0, 1.0, 2.0) .* MONTH)
    # The start takes every row at or before it; `CartesianIndex(3, 4)` is column-major cell 33.
    @test start[1, 33] == 3 && sum(start) == 3
    @test one[1, 5] == 1 && one[3, 7] == 0 && sum(one) == 4
    @test two[2, 5] == 2 && sum(two) == 6
    # Another ecosystem, another kind of table, and its region, as released.
    again = _tableadditions(_eco(), op, (0.0, 1.0, 2.0) .* MONTH)
    @test last(again) == two
    mask = falses(10, 10)
    mask[5] = true
    masked = last(_tableadditions(_eco(), op, (0.0, 1.0, 2.0) .* MONTH,
                                  region = CellMask(mask)))
    @test masked[1, 5] == 1 && masked[2, 5] == 2 && sum(masked) == 3
    # What it refuses, as it resolves.
    table(; kw...) = _quiet() do
        return AddAbundanceTable((species = [1], cell = [1],
                                  time = [0.0 * MONTH], kw...))
    end
    @test_throws "species list has 3" _applytable!(_eco(), table(species = [4]),
                                                   0.0 * MONTH)
    @test_throws "grid has 100 cells" _applytable!(_eco(), table(cell = [101]),
                                                   0.0 * MONTH)
    @test_throws "off the 10 by 10 grid" _applytable!(_eco(),
                                                      table(cell = [CartesianIndex(11,
                                                                                   1)]),
                                                      0.0 * MONTH)
    @test_throws "a count must be" _applytable!(_eco(), table(count = [-1]),
                                                0.0 * MONTH)
    @test_throws "of its own" Intervention(AtTime(MONTH), AllCells(), table())
end

@testset "build_abundance_table writes the table the released operation reads" begin
    # Non-square, so a `y` and an `x` exchanged give a different cell rather than the same one.
    area = StudyArea(extent = (7.0km, 5.0km), cellsize = 1.0km,
                     verbosity = :silent)
    env = GridHabitat(regime = UniformSpec(285.0K, axis = Temperature),
                      supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                           axis = SolarRadiation),
                      area = area)
    spp = build_species(SpeciesPool(3),
                        tolerance = ToleranceSpec(Temperature, μ = 285.0K,
                                                  σ = 5.0K),
                        demand = DemandSpec(SolarRadiation, 1.0e9kJ / day),
                        seed = 1, verbosity = :silent)
    eco = build_ecosystem(BySupply(total = 3000), spp, env, seed = 1,
                          verbosity = :silent)
    height, width = Base.size(parent(eco.habitat.active))
    @test height != width
    cells = LinearIndices((height, width))
    names = eco.spplist.names
    records = (species = [names[3], names[1], names[3], names[3]],
               y = [height, 1, height, 1], x = [1, width, 1, width],
               time = [1.0, 0.0, 1.0, 2.0] .* MONTH)
    table = @test_deprecated build_abundance_table(eco, records)
    @test table.species == [1, 3, 3]
    @test table.cell == [cells[1, width], cells[height, 1], cells[1, width]]
    @test table.count == [1, 2, 1] && table.count isa Vector{Int}
    @test table.time == uconvert.(s, [0.0, 1.0, 2.0] .* MONTH)

    # Its output is the operation's input.
    added = _tableadditions(eco, _quiet(() -> AddAbundanceTable(table)),
                            (0.0, 1.0) .* MONTH)
    @test last(added)[3, cells[height, 1]] == 2

    fractional = _quiet() do
        return build_abundance_table(eco,
                                     (species = [1, 1], y = [1, 1],
                                      x = [1, 1],
                                      count = [0.25, 0.5],
                                      time = [MONTH, MONTH]))
    end
    @test fractional.count == [0.75]

    @test_throws "not in the species list" _quiet(() -> build_abundance_table(eco,
                                                                              (species = ["nosuch"],
                                                                               y = [1],
                                                                               x = [1],
                                                                               time = [MONTH])))
    @test_throws "off the" _quiet(() -> build_abundance_table(eco,
                                                              (species = [1],
                                                               y = [height],
                                                               x = [width + 1],
                                                               time = [MONTH])))
    @test_throws "no rows" _quiet(() -> build_abundance_table(eco,
                                                              (species = Int[],
                                                               y = Int[],
                                                               x = Int[],
                                                               time = typeof(MONTH)[])))
end

@testset "removed shims leave nothing behind" begin
    for n in (:simulate_action!, :simulate_record!, :simulate_record_diversity!,
        :retrieve_era5, :readfile, :ConstructedSpec,
        :emptyMPIgridlandscape,
        :empty_mpi_gridlandscape)
        @test !isdefined(EcoSISTEM, n)
        @test !Base.isexported(EcoSISTEM, n) && !Base.ispublic(EcoSISTEM, n)
    end
    # Methods taken from functions that stay: the six-argument `simulate!` that cached to files,
    # the five-argument demographic constructors, `build_species`'s `boost`, and `read` on the
    # dataset types.
    @test !hasmethod(simulate!,
                     Tuple{Ecosystem, Unitful.Time, Unitful.Time, Unitful.Time,
                           String, String})
    @test !hasmethod(EqualPop, NTuple{5, Any})
    @test_throws MethodError build_species(DefaultEcosystem(), numspecies = 2,
                                           boost = 1.0, verbosity = :silent)
    # Another package's generic method may still accept these (FlatBuffers has `read(::Type{T},
    # bytes)`), so the question is whether any method is the package's own.
    ext = Base.get_extension(EcoSISTEM, :EcoSISTEMRasterDataSourcesExt)
    ours(sig) = [m for m in methods(read, sig) if m.module in (EcoSISTEM, ext)]
    @test isempty(ours(Tuple{Type{WorldClim{BioClim}}, Any}))
    @test isempty(ours(Tuple{Type{CHELSA{Climate}}, String, String}))
    @test isempty(ours(Tuple{Type{EcoSISTEM.ERA}, String, String}))
    @test isempty(ours(Tuple{Type{EcoSISTEM.CRUTS}, String, String}))
end

end
