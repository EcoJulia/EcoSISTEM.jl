using HDF5

"""
    parse_hdf5(path; grid="10k", component="scotland_2018")

Parse HDF5-format file at `path`, containing `component`. The data is
read on a grid with cells of size `grid` meters x `grid` meters and ages are assumed to be
binned into 5 years intervals, starting at zero and going up to 90+.
"""
function parse_hdf5(path; grid="10k", component="scotland_2018")
    # TODO write a method that automatically downloads the file from Github if no path is
    # provided.
    grid in ["1k", "10k"] || throw(
        ArgumentError("Invalid argument value grid=$grid. Allowed values are 1k or 10k.")
    )

    data = h5read(path, component)

    # Get rows and columns
    field = "grid" * grid
    cell_ids = map(x -> parse.(Int, split(x, "-")), data[field*"_id"])
    max_row = maximum(getindex.(cell_ids, 1))
    max_column = maximum(getindex.(cell_ids, 2))

    # Create empty grid
    cell_size = parse(Int, chop(grid))km
    population = AxisArray(
        zeros(max_row, max_column, 19),
        # Placing the coordinates of the cell as the distance between its origin and the
        # origin of the grid.
        grid_x=[cell_size * (i - 1) for i in 1:max_row],
        grid_y=[cell_size * (i - 1) for i in 1:max_column],
        age=[5 * i * year for i in 0:18],
    )

    # Populate grid
    field = "grid" * grid
    for (id, pop) in zip(cell_ids, eachrow(data[field*"_binned"]))
        population[grid_x = id[1], grid_y = id[2]] = pop
    end

    return population
end
