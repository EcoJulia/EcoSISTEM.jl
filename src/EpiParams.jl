using Unitful

mutable struct EpiGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      beta::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function EpiGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, beta::TimeUnitType{U}, sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, beta, sigma)
    end
end
