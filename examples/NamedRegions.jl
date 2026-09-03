# SPDX-License-Identifier: LGPL-3.0-or-later
#
# Finding a study area by name, and building one from several.
#
# **The question this answers.** You have a coordinate, or a raster, or a study area, and you want to
# know what it is *in* - or you want to simulate "the British Isles" and have no shapefile for it.
# Both are lookups against Natural Earth's own polygons, shipped with the package as a table of
# 2 444 named regions.
#
# **Three things it shows, in order:**
#
# 1. **Discovery.** `investigate_regions` asks which named regions relate to something you have, as
#    `investigate_study_area` reports on a grid before one is built. The answer is a report you can
#    read, and whose rows convert straight into specs.
# 2. **Composition.** Regions combine as *geometry* through `ConstructedShapeSpec`, which is exact
#    and carries no resolution - so the grid is still free to be decided afterwards.
# 3. **A run.** The composed region becomes a study area and a short simulation, which is what all
#    of it is for.
#
# **What it is careful about.** The discovery query compares **bounding boxes**, because that is what
# costs no download - and a box can be far larger than the ground it names. The run below prints a
# case of that, and the difference is the reason `NaturalEarthSpec` exists beside `boundingbox`.
#
# **A module, deliberately** - see `examples/README.md`.
#
#     julia --project=examples examples/NamedRegions.jl

module NamedRegions

using EcoSISTEM
using EcoSISTEM.Units
using Rasters: EPSG
using Unitful, Unitful.DefaultSymbols
using Diversity
using Random

const SEED = 20260903

# --- 1. What is this coordinate in? -----------------------------------------------------------
#
# Edinburgh. The report is ordered smallest-enclosing-first, so the tightest region containing the
# point comes first and the coarsest last.
const EDINBURGH = LatLong(55.95°, -3.19°)

println("Which named regions enclose Edinburgh?\n")
const ENCLOSING = EcoSISTEM.investigate_regions(EDINBURGH, kind = :political,
                                                limit = 6)
display(ENCLOSING)
println()

# The box tier being honest rather than wrong: **Norway** is on that list, because its territory runs
# west to Jan Mayen and north to Svalbard, and the box around all of it contains Edinburgh. Nothing
# is broken - the query said which *boxes* enclose the point, and Norway's does. It is why `only`
# is the safe way to take a single answer, and why an actual shape is a separate step.
println("Norway on the list? ",
        any(m -> m.name == "Norway", ENCLOSING),
        " - its box runs from Jan Mayen to Svalbard.\n")

# Narrowing at source is the first tool: one level, one answer.
const SCOTLAND = only(EcoSISTEM.investigate_regions(EDINBURGH,
                                                    level = "SUBUNIT"))
println("At SUBUNIT level there is exactly one: ", SCOTLAND.name, "\n")

# --- 2. Building a region from several ---------------------------------------------------------
#
# The British Isles. Natural Earth *has* a polygon of that name - and it stops at 59.80 degrees
# north, which cuts off Shetland. The union of the three countries does not, so the union is what
# this uses. That is the general lesson about these outlines: they are drawn to be printed at a
# scale, not to be authoritative boundaries.
const UK = NaturalEarthSpec("United Kingdom", level = "ADMIN")
const IRELAND = NaturalEarthSpec("Ireland", level = "ADMIN")
const MAN = NaturalEarthSpec("Isle of Man", level = "ADMIN")

# `LandmassesAbove` drops Rockall, which is 0.031 km2 of uninhabitable granite 300 km out in the
# Atlantic and would otherwise stretch the grid west by three degrees for one cell's worth of rock.
const ISLES = ConstructedShapeSpec(ShapeUnion(), UK, IRELAND, MAN,
                                   coverage = LandmassesAbove(1km^2))

println("The shipped BRITISH ISLES polygon vs the union of the three countries:")
println("  polygon: ",
        EcoSISTEM.boundingbox("BRITISH ISLES", level = "Physical Island group"))
println("  union  : reaches Shetland, and drops Rockall\n")

# --- 3. Simulating on it ------------------------------------------------------------------------
#
# `within` takes the composed spec directly: it both restricts which cells are simulated and sets the
# grid's extent, so no separate bounding box is needed. The grid is **projected**, which a simulation
# requires - dispersal is expressed against one cell size, which only a projected grid has.
const TEMPERATURE = UniformSpec(285.0K, axis = Temperature)
const SUNLIGHT = UniformSpec(1.0e4kJ / (km^2 * day), axis = SolarRadiation)

const AREA = StudyArea(regime = TEMPERATURE, supply = SUNLIGHT, within = ISLES,
                       crs = EPSG(27700), cellsize = 20.0km,
                       verbosity = :silent)

const ENVIRONMENT = GridHabitat(regime = TEMPERATURE, supply = SUNLIGHT,
                                area = AREA, topology = Island())

println("Grid built on the composed region:")
println("  cells: ", size(ENVIRONMENT.active), "  active: ",
        count(ENVIRONMENT.active), " of ", length(ENVIRONMENT.active))
println("  cell size: ", getcellsize(ENVIRONMENT), "\n")

const SPECIES = build_species(8,
                              tolerance = NicheTolerance(Temperature, Normal,
                                                         fill(285.0K, 8),
                                                         fill(3.0K, 8)),
                              demand = SimpleDemand(SolarRadiation,
                                                    fill(1.0kJ / day, 8)))

const ECO = build_ecosystem(SPECIES, ENVIRONMENT, seed = SEED)
simulate!(ECO, 2years, 1month)

println("After two simulated years on the British Isles:")
println("  total abundance : ", sum(ECO.abundances.matrix))
println("  occupied cells  : ",
        count(>(0), sum(ECO.abundances.matrix, dims = 1)))
println("  species present : ",
        count(>(0), sum(ECO.abundances.matrix, dims = 2)))

end
