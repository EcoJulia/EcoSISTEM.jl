using HDF5
using SQLite
using BritishNationalGrid
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
function get_bng(east_north::Tuple{Int64, Int64}, ref::Int64)
    pnt = BNGPoint(east_north[1], east_north[2])
    return gridref(pnt, ref, true)
end
function create_BNG_grid(east::Vector{Int64}, north::Vector{Int64}, vals::Array{Float64, 1}, ages::Vector{Int64})
    easts = collect(minimum(east):1_000:maximum(east)) .* m
    norths = collect(minimum(north):1_000:maximum(north)) .* m
    grid_a = AxisArray(zeros(Float64, length(norths), length(easts), length(unique(ages))), Axis{:northing}(norths), Axis{:easting}(easts), Axis{:age}(unique(ages)))
    for i in eachindex(east)
        grid_a[atvalue(north[i] * m), atvalue(east[i] * m), atvalue(ages[i])] = vals[i]
    end
    return grid_a
end
"""
    get_3d_km_grid_axis_array(cn::SQLite.DB, dims::Array{String,1}, msr::String, tbl::String)

Function to take Scottish Population data from an SQLite database and convert to an axis array.

"""
function get_3d_km_grid_axis_array(cn::SQLite.DB, dims::Array{String,1}, msr::String, tbl::String)
    sel_sql = ""
    dim_ax = []
    for i in eachindex(dims)
        sel_sql = string(sel_sql, dims[i], ",")
        dim_st = SQLite.Stmt(cn, string("SELECT DISTINCT ", dims[i], " AS val FROM ", tbl, " ORDER BY ", dims[i]))
        dim_vals = SQLite.DBInterface.execute(dim_st) |> DataFrames.DataFrame
        # av = i < 3 ? [(v)km for v in dim_vals.val] : dim_vals.val   # unit conversion
        push!(dim_ax, AxisArrays.Axis{Symbol(dims[i])}(dim_vals.val))
    end
    sel_sql = string("SELECT ", sel_sql, " SUM(", msr, ") AS val\nFROM ", tbl, "\nGROUP BY ", rstrip(sel_sql, ','))
    stmt = SQLite.Stmt(cn, sel_sql)
    df = SQLite.DBInterface.execute(stmt) |> DataFrames.DataFrame
    grid_area = df[!, :grid_area]
    east,north = get_en(grid_area)
    vals = df[!, :val]
    ages = df[!, :age_aggr]
    output = create_BNG_grid(east, north, vals, ages)
    return output
end
