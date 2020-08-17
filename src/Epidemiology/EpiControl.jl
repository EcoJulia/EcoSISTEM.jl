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

"""
    Lockdown <: AbstractControl


"""
mutable struct Lockdown <: AbstractControl
    lockdown_date::Unitful.Time
    current_date::Unitful.Time
end

function Lockdown(lockdown_date::Unitful.Time)
    return Lockdown(lockdown_date, 0.0days)
end
