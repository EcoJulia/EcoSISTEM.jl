function _compartment_idx(compartment, names)
    idx = findfirst(names .== compartment)
    isnothing(idx) && throw(ArgumentError("Compartment $compartment not in $names"))
    return idx
end

function _default_steps(abuns)
    N_steps = size(abuns, 3)
    return Int.(floor.(range(1, N_steps; length=4)))
end

function _default_clim(abuns, idx, steps)
    max_val = maximum(maximum(abuns[idx, :, step]) for step in steps)
    return (0, max_val)
end

"""
    epiheatmaps(
        abuns::AbstractArray{<:Integer, 3},
        epilist::EpiList,
        epienv::AbstractEpiEnv;
        compartment="Infected",
        steps=[],
    )

Plot heatmaps of `abuns` for `compartment` at `steps`.
"""
epiheatmaps
@userplot EpiHeatmaps
function _check_args(h::EpiHeatmaps)
    correct_args = (
        length(h.args) == 3 &&
        isa(h.args[1], AbstractArray{<:Integer, 3}) &&
        isa(h.args[2], EpiList) && isa(h.args[3], AbstractEpiEnv)
    )
    if !correct_args
        throw(ArgumentError(
            "epiheatmaps requires (abuns, EpiList, EpiEnv); got: $(typeof(h.args))"
        ))
    end
end
@recipe function f(
    h::EpiHeatmaps;
    compartment="Infected",
    steps=[],
)
    _check_args(h)
    abuns, epilist, epienv = h.args
    idx = _compartment_idx(compartment, epilist.names)
    if isempty(steps)
        steps = _default_steps(abuns)
    end

    layout := length(steps)

    subplot = 1
    gridsize = size(epienv.habitat.matrix)
    for step in steps
        data = reshape(abuns[idx, :, step], gridsize...)
        @series begin
            seriestype := :heatmap
            title := "Day $step ($compartment)"
            subplot := subplot
            data
        end
        subplot += 1
    end
end
