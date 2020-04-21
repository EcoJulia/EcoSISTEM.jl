using Unitful

mutable struct SIRGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      beta::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SIRGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, beta::TimeUnitType{U}, sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, beta, sigma)
    end
end

mutable struct SEI2HRDGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      beta::TimeUnitType{U}
      sigma_1::TimeUnitType{U}
      sigma_2::TimeUnitType{U}
      sigma_hospital::TimeUnitType{U}
      mu_1::TimeUnitType{U}
      mu_2::TimeUnitType{U}
      hospitalisation::TimeUnitType{U}
      death_home::TimeUnitType{U}
      death_hospital::TimeUnitType{U}
    function SEI2HRDGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, beta::TimeUnitType{U}, sigma_1::TimeUnitType{U}, sigma_2::TimeUnitType{U}, sigma_hospital::TimeUnitType{U}, mu_1::TimeUnitType{U}, mu_2::TimeUnitType{U}, hospitalisation::TimeUnitType{U}, death_home::TimeUnitType{U}, death_hospital::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, beta, sigma_1, sigma_2, sigma_hospital, mu_1, mu_2, hospitalisation, death_home, death_hospital)
    end
end

"""
SEI2HRDGrowth(birth::Vector{TimeUnitType{U}},
    death::Vector{TimeUnitType{U}}, beta::TimeUnitType{U}, prob_sym::Float64, prob_hosp::Float64, case_fatality_ratio::Float64, T_lat::Unitful.Time, T_asym::Unitful.Time, T_sym::Unitful.Time, T_hosp::Unitful.Time, T_rec::Unitful.Time) where {U <: Unitful.Units}

Function to calculate SEI2HRD parameters from initial probabilities of developing symptoms and needing hospitalisation, case fatality ratio, and time for transition between different categories.
"""
function SEI2HRDGrowth(birth::Vector{TimeUnitType{U}},
    death::Vector{TimeUnitType{U}}, beta::TimeUnitType{U}, prob_sym::Float64, prob_hosp::Float64, case_fatality_ratio::Float64, T_lat::Unitful.Time, T_asym::Unitful.Time, T_sym::Unitful.Time, T_hosp::Unitful.Time, T_rec::Unitful.Time) where {U <: Unitful.Units}
    # Prob of death at hospital
    prob_hosp_death = case_fatality_ratio/prob_hosp
    # Prob of death at home
    prob_death = case_fatality_ratio/(1 - prob_hosp)
    # Exposed -> asymptomatic
    mu_1 = 1/T_lat
    # Asymptomatic -> symptomatic
    mu_2 = prob_sym * 1/T_asym
    # Symptomatic -> hospital
    hospitalisation = prob_hosp * 1/T_sym
    # Asymptomatic -> recovered
    sigma_1 = (1 - prob_sym) * 1/T_asym
    # Symptomatic -> recovered
    sigma_2 = (1 - prob_hosp) * (1 - prob_death) * 1/T_rec
    # Hospital -> recovered
    sigma_hospital = (1 - prob_hosp_death) * 1/T_hosp
    # Symptomatic -> death
    death_home = (1 - prob_hosp) * prob_death * 2/T_hosp
    # Hospital -> death
    death_hospital = prob_hosp_death * 1/T_hosp
    return SEI2HRDGrowth{U}(birth, death, beta, sigma_1, sigma_2, sigma_hospital, mu_1, mu_2, hospitalisation, death_home, death_hospital)
end
