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
    plot_epiheatmaps(
        abuns::AbstractArray{<:Integer, 3},
        epilist::EpiList,
        epienv::AbstractEpiEnv;
        compartment="Infected",
        steps=[],
    )

Plot heatmaps of `abuns` for `compartment` at `steps`.
"""
plot_epiheatmaps
@userplot Plot_EpiHeatmaps
function _check_args(h)
    correct_args = (
        length(h.args) == 3 &&
        isa(h.args[1], AbstractArray{<:Integer, 3}) &&
        isa(h.args[2], EpiList) && isa(h.args[3], AbstractEpiEnv)
    )
    if !correct_args
        throw(ArgumentError(
            "$(typeof(h)) requires (abuns, EpiList, EpiEnv); got: $(typeof(h.args))"
        ))
    end
end
@recipe function f(
    h::Plot_EpiHeatmaps;
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
        data = Float64.(reshape(abuns[idx, :, step], gridsize...))
        data[.!epienv.active] .= NaN
        @series begin
            seriestype := :heatmap
            title := "Day $step ($compartment)"
            subplot := subplot
            background_color --> :lightblue
            background_color_outside --> :white
            aspect_ratio --> 1
            grid --> false
            data
        end
        subplot += 1
    end
end

"""
    plot_epidynamics(
        abuns::AbstractArray{<:Integer, 3},
        epilist::EpiList,
        epienv::AbstractEpiEnv,
    )

Plot the dynamics over time of `abuns`.
"""
plot_epidynamics
@userplot Plot_EpiDynamics
@recipe function f(h::Plot_EpiDynamics)
    _check_args(h)
    abuns, epilist, epienv = h.args

    for (idx, name) in enumerate(epilist.names)
        data = vec(mapslices(sum, abuns[idx, :, :], dims = 1))
        title --> "Infection dynamics"
        xguide --> "Step"
        yguide --> "Totals"
        @series begin
            label := name
            data
        end
    end
end
