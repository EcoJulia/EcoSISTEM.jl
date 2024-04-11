function unzip(path::String)
    newpath = joinpath(splitpath(path)[1:(end - 1)])
    run(`tar xvf $path -C $newpath`)
    println("Unzipped to $newpath")
    return newpath
end
