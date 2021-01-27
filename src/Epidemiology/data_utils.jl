using HDF5
using SQLite
"""
    parse_hdf5(path; grid="10k", component="scotland_2018")

Parse HDF5-format file at `path`, containing `component`. The data is
read on a grid with cells of size `grid` meters x `grid` meters and ages are assumed to be
binned into 5 years intervals, starting at zero and going up to 90+.
"""
function parse_hdf5(path; grid="10km", component="grid10km/10year/persons")
    # TODO write a method that automatically downloads the file from Github if no path is
    # provided.
    findfirst(grid, component) === nothing && throw(
        ArgumentError("Invalid argument value grid=$grid. Allowed values are 1km or 10km.")
    )

    data = h5read(path, component)

    # Get rows and columns
    distance = data["Dimension_1_units"][1]
    cell_ids = map(x -> parse.(Int, split(x, "-")), data["Dimension_1_names"])
    max_x = maximum(getindex.(cell_ids, 1))
    max_y = maximum(getindex.(cell_ids, 2))

    # Create empty grid
    cell_size = parse(Int, split(distance, "km")[1])km
    n_ages = length(data["Dimension_2_names"])
    n_years = parse(Int, split(split(component, "/")[2], "year")[1])year
    population = AxisArray(
        zeros(Int, max_x, max_y, n_ages),
        # Placing the coordinates of the cell as the distance between its
        # origin and origin of the grid.
        grid_x=[cell_size * (i - 1) for i in 1:max_x],
        grid_y=[cell_size * (i - 1) for i in 1:max_y],
        age=[i * n_years for i in 0:(n_ages-1)],
    )

    # Populate grid
    for (id, pop) in zip(cell_ids, eachrow(data["array"]))
        population[grid_x = id[1], grid_y = id[2]] = Int.(pop)
    end

    return population
end

# """
#     parse_scottish_population(api::DataPipelineAPI; product="human/demographics/population/scotland", component="grid1km/age/persons", aggregate_age=10)
#
# Parse HDF5-format file of Scottish Population sizes at `product`, containing `component`. The data can be aggregated by age, `aggregate_age`, but otherwise is binned in yearly intervals, starting at zero and going up to 90+.
# """
# function parse_scottish_population(api::DataPipelineAPI; product="human/demographics/population/scotland", component="grid1km/age/persons", aggregate_age=10)
#
#     data = read_array(api, product, component)
#     data.data[data.data .< 0] .= 0
#
#     # Get rows and columns
#     distance = data.dimensions[1].units
#     cell_ids = map(x -> parse.(Int, split(x, "-")), data.dimensions[1].names)
#     max_x = maximum(getindex.(cell_ids, 1))
#     max_y = maximum(getindex.(cell_ids, 2))
#
#     # Create empty grid
#     cell_size = parse(Int, split(distance, "km")[1])km
#     n_ages = length(data.dimensions[2].names)
#     agg_years = collect(1:aggregate_age:n_ages)
#     push!(agg_years, n_ages + 1)
#     expanded_agg_years = [agg_years[i]:agg_years[i+1]-1 for i in 1:length(agg_years) - 1]
#     population = AxisArray(
#         zeros(Int, max_x, max_y, length(agg_years)-1),
#         # Placing the coordinates of the cell as the distance between its
#         # origin and origin of the grid.
#         grid_x=[cell_size * (i - 1) for i in 1:max_x],
#         grid_y=[cell_size * (i - 1) for i in 1:max_y],
#         age=[i for i in 0:aggregate_age:(n_ages-1)],
#     )
#
#     # Populate grid
#     data_age = hcat([sum(data.data[i, :], dims = 1)[1, :] for i in expanded_agg_years]...)
#     for (id, pop) in zip(cell_ids, eachrow(data_age))
#         population[grid_x = id[1], grid_y = id[2]] = Int.(pop)
#     end
#
#     return population
# end
#
#
function get_3d_km_grid_axis_array(cn::SQLite.DB, dims::Array{String,1}, msr::String, tbl::String)
    sel_sql = ""
    dim_ax = []
    for i in eachindex(dims)
        sel_sql = string(sel_sql, dims[i], ",")
        dim_st = SQLite.Stmt(cn, string("SELECT DISTINCT ", dims[i], " AS val FROM ", tbl, " ORDER BY ", dims[i]))
        dim_vals = SQLite.DBInterface.execute(dim_st) |> DataFrames.DataFrame
        av = i < 3 ? [(v)km for v in dim_vals.val] : dim_vals.val   # unit conversion
        push!(dim_ax, AxisArrays.Axis{Symbol(dims[i])}(av))
    end
    sel_sql = string("SELECT ", sel_sql, " SUM(", msr, ") AS val\nFROM ", tbl, "\nGROUP BY ", rstrip(sel_sql, ','))
    stmt = SQLite.Stmt(cn, sel_sql)
    df = SQLite.DBInterface.execute(stmt) |> DataFrames.DataFrame
    ## scottish population AxisArray
    axis_size = Tuple(Int64[length(d) for d in dim_ax])
    # data = zeros(typeof(df.val[1]), axis_size)
    output = AxisArrays.AxisArray(zeros(typeof(df.val[1]), axis_size), Tuple(dim_ax))
    for row in eachrow(df)
        output[AxisArrays.atvalue(row[Symbol(dims[1])]km), AxisArrays.atvalue(row[Symbol(dims[2])]km), AxisArrays.atvalue(row[Symbol(dims[3])])] = row.val
    end
    return output
end
