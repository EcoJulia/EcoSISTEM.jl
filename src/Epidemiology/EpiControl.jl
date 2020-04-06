using Unitful

"""
    AbstractControl

Abstract type for all control strategies.
"""
abstract type AbstractControl end

"""
    NoControl <: AbstractControl

Default strategy, which implements no control measures.
"""
mutable struct NoControl <: AbstractControl end
