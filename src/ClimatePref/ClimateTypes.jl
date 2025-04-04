# SPDX-License-Identifier: LGPL-3.0-or-later

using AxisArrays
using Unitful
using EcoSISTEM.Units
using RecipesBase
using RasterDataSources
const RDS = RasterDataSources

import Base: size, length, eltype

"""
    AbstractClimate

Abstract supertype of all climate data.
"""
abstract type AbstractClimate end

"""
    ClimateRaster{<:RDS.RasterDataSource, <: AxisArray} <: AbstractClimate

Type for climate data derived from `RasterDataSource`s.
"""
mutable struct ClimateRaster{R <: RDS.RasterDataSource, A <: AxisArray} <:
               AbstractClimate
    array::A
    function ClimateRaster(T::Type{<:RDS.RasterDataSource},
                           a::A) where {A <: AxisArray}
        return new{T, A}(a)
    end
end

Base.size(cr::ClimateRaster) = size(cr.array)
Base.length(cr::ClimateRaster) = length(cr.array)
Base.eltype(::ClimateRaster{RDS, A}) where {RDS, A} = eltype(A)

"""
    Worldclim_monthly <: AbstractClimate

Type that houses data extracted from Worldclim raster files.
"""
mutable struct Worldclim_monthly{A <: AxisArray} <: AbstractClimate
    array::A
    function Worldclim_monthly(array::A) where {A <: AxisArray}
        size(array, 3) == 12 ||
            error("There should be 12 months of data for worldclim")
        return new{A}(array)
    end
end

"""
    ERA <: AbstractClimate

Type that houses data extracted from ERA raster files.
"""
mutable struct ERA{A <: AxisArray} <: AbstractClimate
    array::A
    function ERA(array::A) where {A <: AxisArray}
        typeof(collect(AxisArrays.axes(array, 3).val)[1]) <: Unitful.Time ||
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
        typeof(collect(AxisArrays.axes(array, 3).val)[1]) <: Unitful.Time ||
            error("Third dimension of array must be time")
        return new{A}(array)
    end
end

"""
    Reference <: AbstractClimate

Type that houses a reference data array.
"""
mutable struct Reference{A <: AxisArray} <: AbstractClimate
    array::A
end

"""
    CRUTS <: AbstractClimate

Type that houses data extracted from CRUTS raster files.
"""
mutable struct CRUTS{A <: AxisArray} <: AbstractClimate
    array::A
    function CRUTS(array::A) where {A <: AxisArray}
        typeof(collect(AxisArrays.axes(array, 3).val)[1]) <: Unitful.Time ||
            error("Third dimension of array must be time")
        return new{A}(array)
    end
end

"""
    CHELSA_monthly <: AbstractClimate

Type that houses data extracted from CHELSA raster files.
"""
mutable struct CHELSA_monthly{A <: AxisArray} <: AbstractClimate
    array::A
    function CHELSA_monthly(array::A) where {A <: AxisArray}
        size(array, 3) == 12 ||
            error("There should be 12 months of data for CHELSA")
        return new{A}(array)
    end
end
