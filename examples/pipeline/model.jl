# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The model run of the data pipeline example: one species over Scotland, built from WorldClim's
# temperature and rainfall within Natural Earth's outline, as their fetch runs registered them.
# Launched by `fair run` from `model.yaml`, which runs
#
#     julia --threads=${JULIA_NUM_THREADS:-auto} --project=examples/pipeline examples/pipeline/model.jl
#
# (`JULIA_NUM_THREADS` threads where it is set, else every core): every file is read from the
# run's store, never fetched, and the figure is registered as the run's output. Started any other
# way it reads the asset cache and writes the figure to a temporary directory.

using DataPipeline
using Plots
using Rasters: EPSG
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units

include("inputs.jl")
using EcoSISTEM: getcellareas

launched = haskey(ENV, "FDP_CONFIG_DIR")
handle = launched ? DataPipeline.initialise() : nothing
const INPUTS = inputs(handle)

# --- decide the grid, before anything is built on it ------------------------------------

# British National Grid, which simulating needs a projected CRS for, within Scotland's outline.
const AREA = StudyArea(regime = INPUTS.temperature, supply = INPUTS.rainfall,
                       within = INPUTS.region, crs = EPSG(27700),
                       cellsize = 10.0km)
const HABITAT = GridHabitat(regime = INPUTS.temperature,
                            supply = INPUTS.rainfall, area = AREA)

# --- one species, tolerant of a cool climate and limited by water ------------------------

const SPECIES = build_species(1, tolerance = (281.0K, 3.0K),
                              toleranceaxis = Temperature,
                              demand = 0.1Unitful.L / day,
                              demandaxis = Precipitation, dispersal = 10.0km,
                              movement = AlwaysMovement, abundance = 10_000,
                              seed = 1)
const ECO = build_ecosystem(SPECIES, HABITAT, seed = 1)

# --- run, recording yearly ------------------------------------------------------------------

const YEARS = 5year
const INTERVAL = 1month_mean_duration
const STEPS = Int(round(uconvert(NoUnits, YEARS / 1year))) + 1
abuns = zeros(Int64, 1, length(HABITAT.active), STEPS)
simulate!(RecordAbundance(abuns), ECO, YEARS, INTERVAL,
          every = EveryInterval(1year))

# --- the figure -------------------------------------------------------------------------------

ny, nx = size(HABITAT.active)
frames = reshape(abuns[1, :, :], ny, nx, STEPS)
masked(i) = replace(Float64.(frames[:, :, i]), 0.0 => NaN)
println("Scotland on a $(ny) x $(nx) grid of 10 km cells ",
        "($(count(HABITAT.active)) active): total abundance ",
        sum(frames[:, :, 1]), " -> ", sum(frames[:, :, end]), ".")
heatmap(masked(1), clim = (0, maximum(frames)), title = "start",
        layout = (@layout [a b; c d]), size = (800, 600))
heatmap!(masked(STEPS), clim = (0, maximum(frames)), title = "end",
         subplot = 2)
heatmap!(Array(HABITAT.regime.matrix), title = "annual temperature",
         subplot = 3)
# The supply is the volume of rain a cell receives a day, shown as a depth a year over the land.
rainfall = uconvert.(u"mm/yr",
                     Array(HABITAT.supply.matrix) ./ getcellareas(HABITAT))
heatmap!(ifelse.(Array(HABITAT.active), rainfall, NaN * u"mm/yr"),
         title = "annual rainfall", subplot = 4)

# --- register the figure ----------------------------------------------------------------------

figure = isnothing(handle) ? joinpath(mktempdir(), "Scotland-plot.png") :
         link_write!(handle, "Scotland-plot")
Plots.png(figure)
isnothing(handle) || DataPipeline.finalise(handle)
