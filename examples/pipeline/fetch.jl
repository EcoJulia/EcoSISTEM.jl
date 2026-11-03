# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The fetch runs of the data pipeline example, one for each source `inputs.jl` reads from: the
# files of the source the script's argument names are fetched where they are not held and
# registered with the run, each an output derived from the source's deposit, with its provenance
# record. The registry derives every output of a run from all of the run's inputs, so each source
# is fetched in a run of its own, its deposit the run's one input. The box, a file of your own, is
# not fetched: its run writes it where it is not there, records it with the facts of the source it
# was cut from, and registers it the same way, through its record rather than a spec, so its run
# asks about it as a file. Launched by `fair run` from that source's config,
# which runs, for WorldClim,
#
#     julia --project=examples/pipeline examples/pipeline/fetch.jl worldclim
#
# Run that way again, with everything registered, it stops before a code run exists: the rehearsal
# asks what is missing and registers nothing. Started any other way it fetches into the asset cache
# and registers nothing, every source with no argument. Fetching the ERA5 month needs a Climate
# Data Store key in `~/.cdsapirc`; a file already in the asset cache needs none.

using DataPipeline
using CDSAPI

include("inputs.jl")

# Which of the inputs each source's run fetches, by the root each is tracked under.
const SOURCES = (worldclim = (:temperature, :rainfall),
                 naturalearth = (:region,),
                 era5 = (:summer,),
                 box = (:box,))

launched = haskey(ENV, "FDP_CONFIG_DIR")
wanted = if isempty(ARGS)
    launched &&
        error("name the source this run fetches: one of $(join(keys(SOURCES), ", ")).")
    reduce((a, b) -> (a..., b...), values(SOURCES))
else
    SOURCES[Symbol(only(ARGS))]
end
boxwanted = :box in wanted
# The specs this run fetches, the box left out: it enters the model through its record, so its
# run records and registers it with `recordbox`.
chosen(handle) = inputs(handle, box = nothing)[filter(!=(:box), wanted)]

# What this run has still to do: its specs' files, and the box, asked about as a file.
function todo(handle)
    files = EcoSISTEM.missingfiles(chosen(handle))
    boxwanted &&
        append!(files,
                EcoSISTEM.missingfiles(BOX, track = Track(handle, "box"),
                                       route = OtherOrigin()))
    return files
end

if launched
    nothing_to_fetch = DataPipeline.rehearse() do handle
        return isempty(todo(handle))
    end
    nothing_to_fetch && exit()
end
handle = launched ? DataPipeline.initialise() : nothing
boxwanted && recordbox(handle)
EcoSISTEM.fetchfiles(chosen(handle))
isnothing(handle) || DataPipeline.finalise(handle)
