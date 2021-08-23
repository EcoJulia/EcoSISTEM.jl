using HDF5
using BritishNationalGrid
using AxisArrays
using TOML
"""
    parse_hdf5(path; grid="10k", component="scotland_2018")

Parse HDF5-format file at `path`, containing `component`. The data is
read on a grid with cells of size `grid` meters x `grid` meters and ages are assumed to be
binned into 5 years intervals, starting at zero and going up to 90+.
"""
function parse_hdf5(path; component="grid area/age/persons", aggregate = true)
    data = h5read(path, component)
    gridrefs = data["Dimension_1_values"]
    east,north = get_en(gridrefs)
    vals = data["array"]
    
    if aggregate
        new_ages = ceil(Int, size(vals, 2)/10)
        vals_new = zeros(Float64, size(vals, 1), new_ages)
        vals_new[:, end] = vals[:, end]
        for i in 0:(new_ages -2)
            vals_new[:, i+1] = sum(vals[:, (i*10 + 1):((i+1)*10)], dims = 2)
        end
        @assert sum(vals_new) == sum(vals)
    else
        vals_new = vals
    end

    age_cats = collect(1:size(vals_new, 2))
    output = create_BNG_grid(east, north, vals_new, age_cats)
    return output
end

function get_en(bng_name::String)
    bng = BritishNationalGrid.BNGPoint(parse(Int, bng_name[3:4]) * 1000, parse(Int, bng_name[5:6]) * 1000, uppercase(bng_name[1:2]))
    return bng.e, bng.n
end
function get_en(bng_names::Vector{String})
    ens = get_en.(bng_names)
    easts = [ens[i][1] for i in eachindex(ens)]
    norths = [ens[i][2] for i in eachindex(ens)]
    return easts, norths
end
function get_bng(east::Real, north::Real, ref::Int64 = 4)
    pnt = BNGPoint(east, north)
    return gridref(pnt, ref, true)
end

function create_BNG_grid(east::Vector{Int64}, north::Vector{Int64}, vals::Matrix{T}, ages::Vector{Int64}) where T
    easts = collect(minimum(east):1_000:maximum(east)) .* m
    norths = collect(minimum(north):1_000:maximum(north)) .* m
    grid_a = AxisArray(zeros(typeof(vals[1]), length(norths), length(easts), length(unique(ages))), Axis{:northing}(norths), Axis{:easting}(easts), Axis{:age}(unique(ages)))
    for i in eachindex(east)
        for j in eachindex(ages)
            grid_a[atvalue(north[i] * m), atvalue(east[i] * m), atvalue(ages[j])] = vals[i, j]
        end
    end
    return grid_a
end
function create_BNG_grid(east::Vector{Int64}, north::Vector{Int64}, vals::Array{T, 1}) where T
    easts = collect(minimum(east):1_000:maximum(east)) .* m
    norths = collect(minimum(north):1_000:maximum(north)) .* m
    grid_a = AxisArray(zeros(typeof(vals[1]), length(norths), length(easts)), Axis{:northing}(norths), Axis{:easting}(easts))
    for i in eachindex(east)
        grid_a[atvalue(north[i] * m), atvalue(east[i] * m)] = vals[i]
    end
    return grid_a
end

function parse_toml(path, component)
    file = TOML.parsefile(path)
    return file[component]["value"]
end

function parse_table(path, component)
    tab = h5read(path, component)["table"]
    return DataFrame(tab)
end

function parse_data_toml(path)
    file = TOML.parsefile(path)
    file_keys = keys(file)
    for f in file_keys
        isfile(file[f]["target"]) ||
            download(file[f]["url"], file[f]["target"])
    end
    return 
end