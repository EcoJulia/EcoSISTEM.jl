using AxisArrays
using Unitful
abstract type AbstractClimate end

mutable struct Worldclim <: AbstractClimate
    array::AxisArray
    function Worldclim(array::AxisArray)
        size(array, 3) == 12 ||
            error("There should be 12 months of data for worldclim")
        new(array)
    end
end

mutable struct Bioclim <: AbstractClimate
    array::AxisArray
    function Bioclim(array::AxisArray)
        size(array, 3) == 19 ||
            error("There should 19 climate variables for bioclim")
        new(array)
    end
end

mutable struct ERA <: AbstractClimate
    array::AxisArray
    function ERA(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end
mutable struct CERA <: AbstractClimate
    array::AxisArray
    function CERA(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end

mutable struct Reference <: AbstractClimate
    array::AxisArray
end
