# SPDX-License-Identifier: LGPL-3.0-or-later
#
# **A categorical environment**, where a cell *is* a land-cover class rather than holding a number.
# Species are matched to a **set of acceptable classes**, not to an optimum and a width - so this is
# the whole other half of the niche machinery from the continuous examples.
#
# **The axis is what makes a layer categorical, not the values.** `compress_landcover` collapses
# EarthEnv's twelve per-class cover fractions into the single winning class *code* per cell - an
# integer - and it is `axis = LandCoverTypology` that says those integers are **labels**. Declare a
# continuous axis instead and the same integers are canonicalised to `Float64`, after which the
# builder refuses the pairing outright: *"layer `tolerance` is Int64 in the species tolerance but
# Float64 in the environment regime"*. That refusal is the axes rule working - meaning comes from the
# declared axis, never from the values.
#
# **The three pieces that differ from a continuous run**, and they follow from each other:
#
# | | continuous | categorical |
# |---|---|---|
# | regime | `ContinuousRegime` | `CategoricalRegime` |
# | tolerance | `NicheTolerance` (optimum, width) | `SimpleCategoricalTolerance` (a set of classes) |
# | fit | `NicheSuitability` | `CategoricalSuitability` |
#
# Only the tolerance is written here. `build_ecosystem` **infers** the fit from the tolerance type,
# which is why nothing below names `CategoricalSuitability` - and why getting the axis wrong is caught
# as a type mismatch rather than as a wrong answer.
#
# **Classes by name, never by number.** The numeric codes are an EarthEnv implementation detail;
# the shipped catalogue (`data/RasterDataSources/LandCover.csv`) maps names to them and
# `EcoSISTEM.landcoverclass` is the lookup.
#
#     julia --project=examples examples/CategoricalLandCover.jl

# **A module, deliberately** - see `examples/README.md`.
module CategoricalLandCover

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Extents: Extent
using Unitful
using Unitful.DefaultSymbols

const SMALL = get(ENV, "ECOSISTEM_SCALE", "large") == "small"
const CELLSIZE = SMALL ? 20.0km : 5.0km
const YEARS = SMALL ? 2year : 10year

# A bare dataset with no code hands the combine the whole multi-band raster, which is what
# `compress_landcover` needs - it must see all twelve classes to pick a winner.
const LANDCOVER = ConstructedSpec(compress_landcover, EarthEnv{LandCover},
                                  axis = LandCoverTypology)
const SUNLIGHT = UniformSpec(1.0e4kJ / (km^2 * day), axis = SolarRadiation)

# A small projected box over Britain. Projected because dispersal assumes one uniform cell size.
const AREA = StudyArea(regime = LANDCOVER, supply = SUNLIGHT,
                       within = Extent(Y = (54.0°, 58.5°), X = (-6.5°, -1.0°)),
                       crs = EPSG(27700), cellsize = CELLSIZE,
                       verbosity = :silent)

const ENVIRONMENT = GridHabitat(regime = LANDCOVER, supply = SUNLIGHT,
                                area = AREA)

# Two guilds that partition the landscape between them: one on wooded ground, one on open ground.
# A `SimpleCategoricalTolerance` takes **one vector of acceptable classes per species**, so this is a
# vector of two.
const WOODLAND = (:needleleaf_trees, :evergreen_broadleaf_trees,
                  :deciduous_broadleaf_trees, :other_trees, :shrubs)
const OPEN = (:herbaceous, :cultivated_and_managed, :regularly_flooded)

_classes(names) = collect(EcoSISTEM.landcoverclass.(names))

# A **pre-built** tolerance, which `build_species` takes as it stands: a set of classes is not a
# mean and a width, so there is no `(mean, width)` form for it and no `toleranceaxis` to give - the
# tolerance already carries its own meaning. It does still declare its `axis`, which is what pairs
# it with the `LandCoverTypology` regime above.
const TOLERANCE = SimpleCategoricalTolerance([
                                                 _classes(WOODLAND),
                                                 _classes(OPEN)
                                             ],
                                             axis = LandCoverTypology)

const SPECIES = build_species(2, tolerance = TOLERANCE,
                              demand = 5.0kJ / day,
                              demandaxis = SolarRadiation,
                              dispersal = 2CELLSIZE, abundance = 200_000,
                              seed = 1)

const ECO = build_ecosystem(SPECIES, ENVIRONMENT, seed = 1)

# `println`, not `@info`: logging goes to stderr, which `@test_nowarn` fails on.
ny, nx = size(ENVIRONMENT.active)
println("Britain on a $(ny) × $(nx) grid of $(CELLSIZE) cells ",
        "($(count(ENVIRONMENT.active)) active).")

# The regime really is categorical, and its element type is what says so.
println("regime is a ", nameof(typeof(ENVIRONMENT.regime)), " of ",
        eltype(ENVIRONMENT.regime), " class codes; the fit inferred from the ",
        "tolerance is a ", nameof(typeof(ECO.nichefit)), ".")

simulate!(ECO, YEARS, 1month_mean_duration)

const FINAL = ECO.abundances.matrix
const WOOD, GRASS = sum(FINAL[1, :]), sum(FINAL[2, :])
println("After $(YEARS): woodland guild $(WOOD), open-ground guild $(GRASS).")

# --- what this must show -------------------------------------------------------------------

# The point of a categorical niche: each guild does **better** where the cell's class is one it
# accepts. Not "occupies only those cells" - dispersal puts individuals into unsuitable cells every
# step, where they simply fare badly. A suitability of zero makes a cell a poor place to be, not an
# unreachable one, so "disjoint occupancy" is a property the model does not claim and does not have.
#
# Mean abundance per cell, split by whether that cell's class is in the guild's set. Read through the
# `grid` view so the `(y, x)` of an abundance and of a regime cell are the same `(y, x)`.
function _meanbyacceptance(guild, accepted)
    classes = ENVIRONMENT.regime.matrix
    inside, outside = Int[], Int[]
    for y in axes(classes, 1), x in axes(classes, 2)
        ENVIRONMENT.active[y, x] || continue
        push!(classes[y, x] in accepted ? inside : outside,
              ECO.abundances.grid[guild, y, x])
    end
    return (inside = isempty(inside) ? 0.0 : sum(inside) / length(inside),
            outside = isempty(outside) ? 0.0 : sum(outside) / length(outside))
end

for (guild, name, names) in ((1, "woodland", WOODLAND), (2, "open", OPEN))
    m = _meanbyacceptance(guild, _classes(names))
    println("  ", rpad(name, 9), "mean ", round(m.inside, digits = 1),
            " per accepted cell vs ", round(m.outside, digits = 1),
            " per rejected one")
    m.inside > m.outside ||
        error("the $(name) guild is no better off in the classes it accepts, so the " *
              "categorical niche is not doing anything")
end

println("✓ each guild is commonest in the land-cover classes its tolerance names.")

end
