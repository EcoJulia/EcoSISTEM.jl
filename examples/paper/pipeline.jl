# SPDX-License-Identifier: LGPL-3.0-or-later
#
# The paper's figures as a FAIR data pipeline run. Launched by `fair run` from `paper.yaml`, which
# runs
#
#     julia --threads=${JULIA_NUM_THREADS:-auto} --project=examples/paper examples/paper/pipeline.jl
#
# on `JULIA_NUM_THREADS` threads where it is set, else every core - each simulation is threaded, and a
# `julia` started without `--threads` has one - in this directory's environment, after `fair pull` has registered `parameters.toml` as
# `Examples/Parameters`. It reads the parameters from the registry, includes `examples/paper.jl`,
# which draws everything from them exactly as a plain run does, into `outputs/` here, and registers
# each output in namespace `EcoSISTEM`, derived from the parameters: every figure as
# `Examples/Figures/<name>`, Tables S1 and S2 as `Examples/Tables/<name>` and `figures.toml` as
# `Examples/Config`.

using DataPipeline

const HANDLE = DataPipeline.initialise()

# Read by `common.jl` in place of the repository's file.
const PAPER_PARAMETERS = link_read!(HANDLE, "Examples/Parameters")

include(normpath(joinpath(@__DIR__, "..", "paper.jl")))

# The data product an output of `paper.jl` is registered as, by its kind.
function productname(output)
    output.kind === :figure && return "Examples/Figures/" * output.name
    output.kind === :table && return "Examples/Tables/" * output.name
    output.kind === :config && return "Examples/Config"
    return error("no data product for an output of kind `$(output.kind)`")
end

for output in Paper.PAPER_OUTPUTS
    link_write!(HANDLE, productname(output), path = output.path)
    println("    registered ", productname(output))
end
DataPipeline.finalise(HANDLE)
