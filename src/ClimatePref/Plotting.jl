using RecipesBase
import Plots: px

# Recipe for plotting ERA data from a particular time period.
@recipe function f(era::ERA, time::Unitful.Time)
    tm = ustrip.(uconvert(year, time))
    yr = floor(Int64, tm)
    ind = round(Int64, (tm - yr)/(1/12))
    typeof(ind) <: Int64 || error("NO")
    mnth = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][ind+1]
    A = transpose(ustrip.(era.array[:, :, time]))
    x = ustrip.(era.array.axes[1].val)
    y = ustrip.(era.array.axes[2].val)
    seriestype  :=  :heatmap
    grid --> false
    title --> "$yr $mnth"
    x, y, A
end

@recipe function f(era::ERA, time::Unitful.Time, xmin::typeof(1.0°),
    xmax::typeof(1.0°), ymin::typeof(1.0°), ymax::typeof(1.0°))
    tm = ustrip.(uconvert(year, time))
    yr = floor(Int64, tm)
    ind = round(Int64, (tm - yr)/(1/12))
    typeof(ind) <: Int64 || error("NO")
    mnth = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][ind+1]
    A = transpose(ustrip.(era.array[xmin .. xmax, ymin .. ymax, time]))
    step1 = ustrip(era.array.axes[1].val[2] - era.array.axes[1].val[1])
    step2 = ustrip(era.array.axes[2].val[2] - era.array.axes[2].val[1])
    x = (ustrip.(xmin) + step1):step1:(ustrip.(xmax) - step1)
    y = (ustrip.(ymin)+ step1):step1:ustrip.(ymax)
    seriestype  :=  :heatmap
    grid --> false
    title --> "$yr $mnth"
    x, y, A
end

@recipe function f(cera::CERA, time::Unitful.Time)
    tm = ustrip.(uconvert(year, time))
    yr = floor(Int64, tm)
    ind = round(Int64, (tm - yr)/(1/12))
    typeof(ind) <: Int64 || error("NO")
    mnth = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][ind + 1]
    A = transpose(ustrip.(cera.array[:, :, time + 1month]))
    x = ustrip.(cera.array.axes[1].val)
    y = ustrip.(cera.array.axes[2].val)
    seriestype  :=  :heatmap
    grid --> false
    title --> "$yr $mnth"
    x, y, A
end

@recipe function f(era::CERA, time::Unitful.Time, xmin::typeof(1.0°),
    xmax::typeof(1.0°), ymin::typeof(1.0°), ymax::typeof(1.0°))
    tm = ustrip.(uconvert(year, time))
    yr = floor(Int64, tm)
    ind = round(Int64, (tm - yr)/(1/12))
    typeof(ind) <: Int64 || error("NO")
    mnth = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][ind+1]
    A = transpose(ustrip.(era.array[xmin .. xmax, ymin .. ymax, time + 1month]))
    step1 = ustrip(era.array.axes[1].val[2] - era.array.axes[1].val[1])
    step2 = ustrip(era.array.axes[2].val[2] - era.array.axes[2].val[1])
    x = (ustrip.(xmin) + step1):step1:ustrip.(xmax)
    y = (ustrip.(ymin)+ step1):step1:ustrip.(ymax)
    seriestype  :=  :heatmap
    grid --> false
    title --> "$yr $mnth"
    x, y, A
end

@recipe function f(wc::Worldclim_monthly, time::Unitful.Time)
    ind = (time + 1month) / month
    typeof(ind) <: Int64 || error("NO")
    mnth = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][ind]
    A = transpose(ustrip.(wc.array[:, :, time + 1month]))
    x = ustrip.(AxisArrays.axes(wc.array, 1).val)
    y = ustrip.(AxisArrays.axes(wc.array, 2).val)
    seriestype  :=  :heatmap
    grid --> false
    background_color_inside --> :grey
    title --> "$mnth"
    x, y, A
end

@recipe function f(wc::Worldclim_monthly, time::Unitful.Time, xmin::typeof(1.0°),
    xmax::typeof(1.0°), ymin::typeof(1.0°), ymax::typeof(1.0°))
    ind = (time + 1month) / month
    typeof(ind) <: Int64 || error("NO")
    mnth = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][ind]
    A = transpose(ustrip.(wc.array[xmin .. xmax, ymin .. ymax, time + 1month]))
    step1 = ustrip(wc.array.axes[1].val[2] - wc.array.axes[1].val[1])
    step2 = ustrip(wc.array.axes[2].val[2] - wc.array.axes[2].val[1])
    x = (ustrip.(xmin)+ step1):step1:ustrip.(xmax)
    y = (ustrip.(ymin) + step1):step1:ustrip.(ymax)
    seriestype  :=  :heatmap
    grid --> false
    background_color_inside --> :grey
    title --> "$mnth"
    x, y, A
end


"""
    getprofile(spp_names::Array{String, 1}, data::IndexedTable, variable_name::String, dims::Tuple{Int64, Int64} = (1,1))

Function to plot climate profiles for a vector of species names, `spp_names`, using a JuliaDB table of GBIF records, `data`, column containing climate variable of interest, `var`, and dimensions over which it should be plotted, `dims`.
"""
# function getprofile(spp_names::Array{String, 1}, data::IndexedTable, var::Symbol, label::String, dims::Tuple{Int64, Int64} = (1,1))

#     # Check for dimensions to be greater or the same as the length of species, or for all to be plotted in one plot pane.
#     (dims[1] * dims[2] >= length(spp_names) || dims == (1,1)) || error("Dimensions not big enough for number of species")
#     first_plot = true
#     # Loop through species names, adding profiles to plot
#     hist = histogram(layout = dims)
#     for i in eachindex(spp_names)
#         # Filter data for species in question
#         spp = filter(p-> p[:species] == spp_names[i], data)
#         # Select climate variable of interest
#         vals = select(spp, var)
#         # Remove NaNs
#         res = vcat(vals...)
#         res = res[.!isnan.(res)]

#         # If dimensions are (1,1)...
#         # ... all plotted to the same subplot window, else individual for each species
#         sp = ifelse(dims == (1,1), 1, i)
#         # ... put legend on the left, otherwise none
#         lg = ifelse(dims == (1,1), :left, :none)
#         # ... title as species name, else blank
#         title = ifelse(dims == (1,1), "", spp_names[i])
#         # For the first species plot histogram, otherwise add to previous plot
#         histogram!(res, bins = -20:2:30, grid = false, xlabel = label, label=spp_names[i], subplot = sp, legend= lg, top_margin = 20px,
#         bottom_margin = 20px, title = title)
#     end
#     return hist
# end
