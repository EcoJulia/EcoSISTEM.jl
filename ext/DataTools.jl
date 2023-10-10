module DataTools

if isdefined(Base, :get_extension)
    using DataPipeline
else
    using ..DataPipeline
end

println("Creating data pipeline interface ...")

using EcoSISTEM
using EcoSISTEM.ClimatePref

function EcoSISTEM.ClimatePref.unzip(path::String)
    newpath = joinpath(splitpath(path)[1:end-1])
    run(`tar xvf $path -C $newpath`)
    println("Unzipped to $newpath")
    return newpath
end

end