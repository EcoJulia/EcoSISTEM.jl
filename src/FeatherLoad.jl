using Feather
using AxisArrays
using DataFrames

function Feathersave(A::AxisArray{Int64,2,Array{Int64,2},
    Tuple{AxisArrays.Axis{:Species,Array{String,1}},
    AxisArrays.Axis{:Param,Array{String,1}}}}, file::String)
    DF = DataFrame(Array(A))
    names!(DF, Symbol.(axes(A)[2].val))
    DF[:Species] = axes(A)[1].val
    Feather.write(file, DF)
end

function Featherload(file::String)
    DF = Feather.read(file)
    A = AxisArray(convert(Array, DF[1:end-1]), Axis{:Species}(DF[:Species]),
    Axis{:Param}(string.(names(DF)[1:end-1])))
    return A
end
