using HDF5

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
