using HDF5

"""
    parse_scotpop(path; grid="10k")

Parse HDF5-format file at `path`, containing the Scottish population by age. The data is
read on a grid with cells of size `grid` meters x `grid` meters.
"""
function parse_scotpop(path; grid="10k")
    # TODO write a method that automatically downloads the file from Github if no path is
    # provided.
    grid in ["1k", "10k"] || throw(
        ArgumentError("Invalid argument value grid=$grid. Allowed values are 1k or 10k.")
    )

    data = h5read(path, "scotland_2018")

    # Get age brackets
    ages = [string(5*i)*"-"*string(5*i+4) for i in 0:17]
    push!(ages, "90+")

    # Get rows and columns
    field = "grid" * grid
    cell_ids = map(x -> parse.(Int, split(x, "-")), data[field*"_id"])
    max_row = maximum(getindex.(cell_ids, 1))
    max_column = maximum(getindex.(cell_ids, 2))

    # Create empty grid
    cell_size = parse(Int, chop(grid))km
    population = AxisArray(
        zeros(max_row, max_column, length(ages)),
        # Placing the coordinates of the cell as the distance between its origin and the
        # origin of the grid.
        grid_x=[cell_size * (i - 1) for i in 1:max_row],
        grid_y=[cell_size * (i - 1) for i in 1:max_column],
        age=ages,
    )

    # Populate grid
    field = "grid" * grid
    for (id, pop) in zip(cell_ids, eachrow(data[field*"_binned"]))
        population[grid_x = id[1], grid_y = id[2]] = pop
    end

    return population
end
