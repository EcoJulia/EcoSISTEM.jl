# SPDX-License-Identifier: LGPL-3.0-or-later

using AxisArrays
using Unitful
using EcoSISTEM.Units
using RecipesBase

"""
    AbstractClimate

Abstract supertype of all climate data.
"""
abstract type AbstractClimate end

"""
    Worldclim_monthly <: AbstractClimate

Type that houses data extracted from Worldclim raster files.
"""
struct Worldclim_monthly{A <: AxisArray} <: AbstractClimate
    array::A
    function Worldclim_monthly(array::A) where {A <: AxisArray}
        size(array, 3) == 12 ||
            error("There should be 12 months of data for worldclim")
        return new{A}(array)
    end
end

"""
    Worldclim_bioclim <: AbstractClimate

Type that houses data extracted from Bioclim raster files.
"""
struct Worldclim_bioclim{A <: AxisArray} <: AbstractClimate
    array::A
end

"""
    Landcover <: AbstractClimate

Type that houses data extracted from EarthEnv Landcover raster files.
"""
struct Landcover{A <: AxisArray} <: AbstractClimate
    array::A
end

"""
    ERA <: AbstractClimate

Type that houses data extracted from ERA raster files.
"""
struct ERA{A <: AxisArray} <: AbstractClimate
    array::A
    function ERA(array::A) where {A <: AxisArray}
        collect(AxisArrays.axes(array, 3).val)[1] isa Unitful.Time ||
            error("Third dimension of array must be time")
        return new{A}(array)
    end
end

"""
    CERA <: AbstractClimate

Type that houses data extracted from CERA-20C raster files.
"""
mutable struct CERA{A <: AxisArray} <: AbstractClimate
    array::A
    function CERA(array::A) where {A <: AxisArray}
        collect(AxisArrays.axes(array, 3).val)[1] isa Unitful.Time ||
            error("Third dimension of array must be time")
        return new{A}(array)
    end
end

"""
    Reference <: AbstractClimate

Type that houses a reference data array.
"""
struct Reference{A <: AxisArray} <: AbstractClimate
    array::A
end

"""
    CRUTS <: AbstractClimate

Type that houses data extracted from CRUTS raster files.
"""
struct CRUTS{A <: AxisArray} <: AbstractClimate
    array::A
    function CRUTS(array::A) where {A <: AxisArray}
        collect(AxisArrays.axes(array, 3).val)[1] isa Unitful.Time ||
            error("Third dimension of array must be time")
        return new{A}(array)
    end
end

"""
    CHELSA_monthly <: AbstractClimate

Type that houses data extracted from CHELSA raster files.
"""
struct CHELSA_monthly{A <: AxisArray} <: AbstractClimate
    array::A
    function CHELSA_monthly(array::A) where {A <: AxisArray}
        size(array, 3) == 12 ||
            error("There should be 12 months of data for CHELSA")
        return new{A}(array)
    end
end

"""
    CHELSA_bioclim <: AbstractClimate

Type that houses data extracted from CHELSA raster files.
"""
struct CHELSA_bioclim{A <: AxisArray} <: AbstractClimate
    array::A
end
