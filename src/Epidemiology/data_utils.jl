using HDF5

"""
    parse_scotpop(path; grid="10k")

Parse HDF5-format file at `path`, containing the Scottish population by age. The data is
read on a `grid`x`grid` square grid.
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

    # Create empty grid
    grid_size = parse(Int, chop(grid)) * 1000
    population = AxisArray(
        zeros(grid_size, grid_size, length(ages)),
        grid_x=collect(1:grid_size),
        grid_y=collect(1:grid_size),
        age=ages,
    )

    # Populate grid
    field = "grid" * grid
    for (id, pop) in zip(data[field*"_id"], eachrow(data[field*"_binned"]))
        x, y = parse.(Int, split(id, "-"))
        population[grid_x = x, grid_y = y] = pop
    end

    return population
end
