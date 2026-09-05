# SPDX-License-Identifier: LGPL-3.0-or-later
#
#### SINGLE SPECIES ####
# One species across Africa on WorldClim data, with the run **bracketed by the FAIR data pipeline**:
# inputs are claimed through `link_read!` and the figure is registered with `link_write!`, so the
# result is traceable to the exact data that produced it. That bracketing is the point of this file -
# the ecology is deliberately the simplest thing that exercises it.
#
# **A `StudyArea` decides a projected grid first, and `GridHabitat` samples onto it.** Building the
# layers by hand and handing them to `Ecosystem` decides no grid at all, leaving the run on
# WorldClim's own **geographic** (°) coordinates, which cannot be simulated: a degree cell's real
# extent changes with latitude while dispersal assumes one uniform cell size.
#
# **The pipeline read does real work rather than being decorative.** `link_read!` returned a
# path that nothing used - the script went on to call `read(WorldClim{BioClim}, ...)` and download the
# data itself. The claimed path is now what the layers are actually read from.
#
# **Not runnable without a configured pipeline.** `DataPipeline.initialise()` needs a
# `config.yaml` and a local registry, so this cannot be part of any automatic run and is not in
# `examples/` proper. With one:
#
#     julia --project=examples examples/pipeline/AfricaRun.jl

using EcoSISTEM
using EcoSISTEM.Units
using RasterDataSources
using Rasters: EPSG
using Unitful
using Unitful.DefaultSymbols
using Plots
using DataPipeline

# --- claim the inputs -------------------------------------------------------------------

handle = DataPipeline.initialise()

# The claimed directory is handed to `RasterDataSources` as its download path, so the layers named
# below resolve to the *pipeline's* copy rather than to one fetched independently. That is what makes
# the run reproducible from the registry rather than merely repeatable.
const CLAIMED = EcoSISTEM.unziptemp(link_read!(handle, "AfricaModel/WorldClim"))
ENV["RASTERDATASOURCES_PATH"] = CLAIMED

# --- the layers -------------------------------------------------------------------------

# Named, not read: bio1 is annual mean temperature and bio12 annual precipitation, and both the
# unit and the axis come from the shipped catalogue (`data/RasterDataSources/BioClim.csv`) rather
# than being asserted here - which is what the old `uconvert.(K, africa_temp .* °C)` was doing by
# hand, and could get wrong.
const TEMPERATURE = SourceSpec(WorldClim{BioClim}, :bio1)
const RAINFALL = SourceSpec(WorldClim{BioClim}, :bio12)

# --- decide the grid, before anything is built on it ------------------------------------

# A projected CRS is required to simulate. EPSG:10592 (WGS 84 / GLANCE Africa) is the package's own
# advice for this extent. `within` positions the area - WorldClim is global - and the land mask now
# comes from the layers' own coverage rather than from a hand-built `.!isnan.(...)` matrix.
const AREA = StudyArea(regime = TEMPERATURE, supply = RAINFALL,
                       within = EcoSISTEM.boundingbox("Africa",
                                                      level = "CONTINENT",
                                                      coverage = LargestLandmass()),
                       crs = EPSG(10592), cellsize = 50.0km)

const ENVIRONMENT = GridHabitat(regime = TEMPERATURE, supply = RAINFALL,
                                area = AREA)

# --- one species, tolerant of warmth and limited by water --------------------------------

# The tolerance is a `Temperature` and the demand a `Precipitation`, and each says so: meaning
# comes from the declared axis, never from the values or their units.
const SPECIES = build_species(1,
                              tolerance = (280.0K, 10.0K),
                              toleranceaxis = Temperature,
                              demand = 0.1Unitful.L / day,
                              demandaxis = Precipitation,
                              dispersal = 15.0km,
                              movement = AlwaysMovement,
                              abundance = 20_000,
                              seed = 1)

const ECO = build_ecosystem(SPECIES, ENVIRONMENT, seed = 1)

# --- run, recording monthly --------------------------------------------------------------

const YEARS = 10year
const INTERVAL = 1month_mean_duration
const LENSIM = length((0year):INTERVAL:YEARS)

abuns = zeros(Int64, 1, length(ENVIRONMENT.active), LENSIM)
@time simulate_record!(abuns, ECO, YEARS, INTERVAL, INTERVAL)

# --- what happened ------------------------------------------------------------------------

ny, nx = size(ENVIRONMENT.active)
frames = reshape(abuns[1, :, :], ny, nx, LENSIM)
masked(i) = replace(Float64.(frames[:, :, i]), 0.0 => NaN)

println("Africa on a $(ny) × $(nx) grid of 50 km cells ",
        "($(count(ENVIRONMENT.active)) active) over $(YEARS).")
println("Total abundance: ", sum(frames[:, :, 1]), " -> ",
        sum(frames[:, :, end]), ".")

heatmap(masked(1), clim = (0, maximum(frames)),
        background_color = :lightblue, background_color_outside = :white,
        grid = false, color = cgrad(:algae, scale = :exp), title = "start",
        layout = (@layout [a b; c d]))
heatmap!(masked(LENSIM), clim = (0, maximum(frames)),
         background_color = :lightblue, background_color_outside = :white,
         grid = false, color = cgrad(:algae, scale = :exp), title = "end",
         subplot = 2)
# The layers as they were actually sampled onto the grid - read off the built habitat rather than
# re-read, so the figure cannot disagree with what the simulation saw.
heatmap!(Array(ENVIRONMENT.regime.matrix), grid = false,
         title = "temperature", subplot = 3)
heatmap!(Array(ENVIRONMENT.supply.matrix), grid = false, title = "rainfall",
         subplot = 4)

# --- register the output -------------------------------------------------------------------

Plots.pdf(link_write!(handle, "Africa-plot"))
DataPipeline.finalise(handle)
