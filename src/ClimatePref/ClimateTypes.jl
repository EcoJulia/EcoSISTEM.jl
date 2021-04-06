using AxisArrays
using Unitful
using EcoSISTEM.Units
using RecipesBase

import AxisArrays.axes

"""
    AbstractClimate

Abstract supertype of all climate data.
"""
abstract type AbstractClimate end

"""
    Worldclim <: AbstractClimate

Type that houses data extracted from Worldclim raster files.
"""
mutable struct Worldclim <: AbstractClimate
    array::AxisArray
    function Worldclim(array::AxisArray)
        size(array, 3) == 12 ||
            error("There should be 12 months of data for worldclim")
        new(array)
    end
end

"""
    Bioclim <: AbstractClimate

Type that houses data extracted from Bioclim raster files.
"""
mutable struct Bioclim <: AbstractClimate
    array::AxisArray
    function Bioclim(array::AxisArray)
        size(array, 3) == 19 ||
            error("There should 19 climate variables for bioclim")
        new(array)
    end
end

"""
    ERA <: AbstractClimate

Type that houses data extracted from ERA raster files.
"""
mutable struct ERA <: AbstractClimate
    array::AxisArray
    function ERA(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end

"""
    CERA <: AbstractClimate

Type that houses data extracted from CERA-20C raster files.
"""
mutable struct CERA <: AbstractClimate
    array::AxisArray
    function CERA(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end

"""
    Reference <: AbstractClimate

Type that houses a reference data array.
"""
mutable struct Reference <: AbstractClimate
    array::AxisArray
end

"""
    CRUTS <: AbstractClimate

Type that houses data extracted from CRUTS raster files.
"""
mutable struct CRUTS <: AbstractClimate
    array::AxisArray
    function CRUTS(array::AxisArray)
        typeof(collect(axes(array, 3).val)[1])<: Unitful.Time ||
            error("Third dimension of array must be time")
        new(array)
    end
end

"""
    CHELSA <: AbstractClimate

Type that houses data extracted from CHELSA raster files.
"""
mutable struct CHELSA <: AbstractClimate
    array::AxisArray
    function CHELSA(array::AxisArray)
        size(array, 3) == 12 ||
            error("There should be 12 months of data for CHELSA")
        new(array)
    end
end
