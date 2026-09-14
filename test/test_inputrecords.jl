# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Tests for `src/inputrecords.jl` - what a study area, a habitat and an ecosystem say they were built
# from.
#
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["test_inputrecords.jl"])'

module TestInputRecords

using Test
using EcoSISTEM
using EcoSISTEM: Provenance
using EcoSISTEM.Units
using Unitful, Unitful.DefaultSymbols
using ArchGDAL
using Dates: Dates
using TOML

# A 5 x 7 WGS84 GeoTIFF of one degree cells whose top-left corner is at 10 degrees east, 55 north.
function _geotiff(path)
    ArchGDAL.create(path, driver = ArchGDAL.getdriver("GTiff"), width = 7,
                    height = 5, nbands = 1, dtype = Float32) do ds
        ArchGDAL.write!(ds, Float32[280 + i + 10j for i in 1:7, j in 1:5], 1)
        ArchGDAL.setgeotransform!(ds, [10.0, 1.0, 0.0, 55.0, 0.0, -1.0])
        return ArchGDAL.setproj!(ds,
                                 ArchGDAL.toWKT(ArchGDAL.importEPSG(4326)))
    end
    return path
end

@testset "an area and its habitat record what they read, what cut the grid, the grid and the software" begin
    dir = mktempdir()
    spec = RasterFileSpec(_geotiff(joinpath(dir, "field.tif")),
                          axis = Temperature, unit = K)
    mask = ConstructedRasterSpec(r -> .!isnan.(r),
                                 RasterFileSpec(_geotiff(joinpath(dir,
                                                                  "mask.tif")),
                                                axis = EcoSISTEM.NicheAxis),
                                 axis = EcoSISTEM.NicheAxis)
    area = StudyArea(regime = spec, within = mask, verbosity = :silent)
    p = provenance(area)
    @test p isa Provenance && isnothing(p.run)
    field = only(r for r in p.inputs if r.path == "field.tif")
    @test field.role === :habitat && field.dataset == "file"
    @test any(r -> r.path == "mask.tif" && r.role === :region, p.inputs)
    @test p.grid.cellsize == 1.0° && p.grid.cells == 35 && p.grid.active == 35
    @test p.grid.extent == EcoSISTEM.Extents.Extent(Y = (50.0°, 55.0°),
                                   X = (10.0°, 17.0°))
    @test p.software.package == "EcoSISTEM"
    @test p.software.version == string(pkgversion(EcoSISTEM))
    @test startswith(p.software.doi, "10.5281/zenodo.")
    @test p.software.julia == string(VERSION)
    # A habitat built on the area gives the same answer, its reads long discarded.
    habitat = GridHabitat(regime = spec,
                          supply = UniformSpec(1.0e5kJ / (m^2 * day),
                                               axis = SolarRadiation),
                          area = area)
    q = provenance(habitat)
    @test isnothing(habitat.area.report.cache)
    @test [r.path for r in q.inputs] == [r.path for r in p.inputs]
    @test q.grid == p.grid && isnothing(q.run)
    # A combination of layers lists its data members' records and none for a synthetic one.
    both = ConstructedRasterSpec((a, b) -> a, spec,
                                 UniformSpec(290.0K, axis = Temperature),
                                 axis = Temperature)
    @test provenance(both) == [nothing]
end

@testset "an ecosystem adds its run, and a synthetic one reads nothing" begin
    eco = build_ecosystem(DefaultEcosystem(), seed = 7)
    p = provenance(eco)
    @test isempty(p.inputs) && isnothing(p.grid.crs)
    @test p.run.seed == eco.seed && p.run.epoch === eco.epoch &&
          p.run.elapsed == eco.elapsed
    @test p.run.calendar == ExactDates()
    @test p.grid == provenance(eco.habitat).grid
end

@testset "records a caller attaches reach the ecosystem, and an intervention's once it has acted" begin
    E = EcoSISTEM
    traits = E.InputRecord(role = :species, dataset = "trait table",
                           path = "traits.csv")
    seeding = E.InputRecord(role = :abundance,
                            dataset = "GBIF occurrence download",
                            doi = "10.15468/dl.abc123")
    species = build_species(DefaultEcosystem(), numspecies = 3,
                            verbosity = :silent, provenance = traits)
    @test only(species.inputs).dataset == "trait table"
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    eco = build_ecosystem(species, habitat, seed = 1, provenance = [seeding])
    datasets() = [r.dataset for r in provenance(eco).inputs]
    @test issetequal(datasets(), ["trait table", "GBIF occurrence download"])
    # An intervention's records join once it acts, and never for one that does not.
    cull = E.InputRecord(role = :intervention, dataset = "cull plan")
    idle = E.InputRecord(role = :intervention, dataset = "never applied")
    plan = InterventionSet(Intervention(AtTime(0s), AllCells(),
                                        RemoveAbundance(1, 1),
                                        provenance = cull),
                           Intervention(NeverScheduled(), AllCells(),
                                        RemoveAbundance(1, 1),
                                        provenance = idle))
    @test !("cull plan" in datasets())
    simulate!(eco, 1month_mean_duration, 1month_mean_duration,
              intervention = plan)
    @test "cull plan" in datasets() && !("never applied" in datasets())
    # Stripping the similarity keeps both the species' and the run's records.
    uniqueeco = EcoSISTEM.makeunique(eco)
    @test [r.dataset for r in uniqueeco.spplist.inputs] == ["trait table"]
    @test [r.dataset for r in uniqueeco.inputs] ==
          [r.dataset for r in eco.inputs]
    # Anything but a record, or a vector of them, is refused naming the keyword.
    @test_throws "`provenance` takes" build_ecosystem(species, habitat,
                                                      provenance = "GBIF")
end

@testset "write_provenance writes what an area and an ecosystem were built from, as TOML" begin
    dir = mktempdir()
    spec = RasterFileSpec(_geotiff(joinpath(dir, "field.tif")),
                          axis = Temperature, unit = K)
    area = StudyArea(regime = spec, verbosity = :silent)
    path = write_provenance(joinpath(dir, "out", "provenance.toml"), area)
    @test isfile(path)
    toml = TOML.parsefile(path)
    @test toml["software"]["package"] == "EcoSISTEM"
    @test startswith(toml["software"]["doi"], "10.5281/zenodo.")
    @test toml["grid"]["cells"] == 35 && toml["grid"]["active"] == 35
    @test toml["grid"]["cellsize"] == string(1.0°)
    @test toml["grid"]["extent"]["Y"] == [string(50.0°), string(55.0°)]
    @test !haskey(toml, "run")
    input = only(toml["inputs"])
    @test input["path"] == "field.tif" && input["role"] == "habitat"
    @test !haskey(input, "url") && !haskey(input, "doi")
    # Nothing in the file belongs to the machine that wrote it.
    text = read(path, String)
    @test !occursin(homedir(), text) && !occursin(gethostname(), text) &&
          !occursin(dir, text)
    # An ecosystem adds its run, and writing it gives the same file as writing its provenance.
    seeding = EcoSISTEM.InputRecord(role = :abundance,
                                    dataset = "GBIF occurrence download",
                                    doi = "10.15468/dl.abc123",
                                    fetched = "2026-09-13T07:21:07Z")
    eco = build_ecosystem(DefaultEcosystem(), seed = 7, provenance = seeding)
    whole = write_provenance(joinpath(dir, "eco.toml"), eco)
    parts = write_provenance(joinpath(dir, "parts.toml"), provenance(eco))
    @test read(whole, String) == read(parts, String)
    run = TOML.parsefile(whole)
    @test run["run"]["seed"] == 7
    @test only(run["inputs"])["fetched"] ==
          Dates.DateTime(2026, 9, 13, 7, 21, 7)
    # A CRS is written as its name, not as the type's own display.
    @test EcoSISTEM._tomlvalue(EcoSISTEM.Rasters.EPSG(27700)) == "EPSG:27700"
end

@testset "a provenance lists its inputs under their datasets, each citation once" begin
    E = EcoSISTEM
    cite = "Smith, A. (2020). A trait table. Journal of Traits, 1, 1-2."
    traits = [E.InputRecord(role = :species, dataset = "trait table",
                            path = file, citation = cite)
              for file in ("traits.csv", "more.csv")]
    seeding = E.InputRecord(role = :abundance,
                            dataset = "GBIF occurrence download",
                            doi = "10.15468/dl.abc123")
    species = build_species(DefaultEcosystem(), numspecies = 3,
                            verbosity = :silent, provenance = traits)
    habitat = build_habitat(DefaultEcosystem(), verbosity = :silent)
    eco = build_ecosystem(species, habitat, seed = 1, provenance = seeding)
    shown() = sprint(show, MIME("text/plain"), provenance(eco))
    before = shown()
    @test length(findall(cite, before)) == 1
    # The dataset heads its group, its citation next, then its files.
    at(text) = first(findfirst(text, before))
    @test at("    trait table") < at("cite  " * cite) < at("more.csv")
    @test at("cite  " * cite) < at("traits.csv")
    @test occursin("    GBIF occurrence download (doi 10.15468/dl.abc123)",
                   before)
    @test occursin("0.0 s elapsed", before)
    # Elapsed time reads in days once the run has gone on for some.
    simulate!(eco, 1month_mean_duration, 1month_mean_duration)
    @test occursin(r"seed 1, [0-9.]+ d elapsed", shown())
end

end
