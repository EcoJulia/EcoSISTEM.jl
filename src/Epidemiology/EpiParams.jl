using Unitful

mutable struct SIRGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SIRGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U},
        virus_decay::TimeUnitType{U}, beta::TimeUnitType{U},
        sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, virus_growth, virus_decay, beta, sigma)
    end
end
mutable struct SISGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SISGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U},
        virus_decay::TimeUnitType{U}, beta::TimeUnitType{U},
        sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, virus_growth, virus_decay, beta, sigma)
    end
end
mutable struct SEIRGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta::TimeUnitType{U}
      mu::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SEIRGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U},
        virus_decay::TimeUnitType{U}, beta::TimeUnitType{U}, mu::TimeUnitType{U}, sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, virus_growth, virus_decay, beta, mu, sigma)
    end
end
mutable struct SEIRSGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta::TimeUnitType{U}
      mu::TimeUnitType{U}
      epsilon::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SEIRSGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U},
        virus_decay::TimeUnitType{U}, beta::TimeUnitType{U},
        mu::TimeUnitType{U}, sigma::TimeUnitType{U}, epsilon::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, virus_growth, virus_decay, beta, mu, sigma, epsilon)
    end
end

mutable struct SEI2HRDGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::Vector{TimeUnitType{U}}
      virus_decay::Vector{TimeUnitType{U}}
      beta::Vector{TimeUnitType{U}}
      sigma_1::TimeUnitType{U}
      sigma_2::TimeUnitType{U}
      sigma_hospital::TimeUnitType{U}
      mu_1::TimeUnitType{U}
      mu_2::TimeUnitType{U}
      hospitalisation::TimeUnitType{U}
      death_home::TimeUnitType{U}
      death_hospital::TimeUnitType{U}
    function SEI2HRDGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::Vector{TimeUnitType{U}},
        virus_decay::Vector{TimeUnitType{U}}, beta::Vector{TimeUnitType{U}},
        sigma_1::TimeUnitType{U}, sigma_2::TimeUnitType{U}, sigma_hospital::TimeUnitType{U}, mu_1::TimeUnitType{U},
        mu_2::TimeUnitType{U}, hospitalisation::TimeUnitType{U}, death_home::TimeUnitType{U},
        death_hospital::TimeUnitType{U}) where {U <: Unitful.Units}
        new{U}(birth, death, virus_growth, virus_decay, beta, sigma_1,
        sigma_2, sigma_hospital, mu_1, mu_2, hospitalisation, death_home, death_hospital)
    end
end

"""
SEI2HRDGrowth(birth::Vector{TimeUnitType{U}},
    death::Vector{TimeUnitType{U}}, beta::TimeUnitType{U}, prob_sym::Float64, prob_hosp::Float64, case_fatality_ratio::Float64, T_lat::Unitful.Time, T_asym::Unitful.Time, T_sym::Unitful.Time, T_hosp::Unitful.Time, T_rec::Unitful.Time) where {U <: Unitful.Units}

Function to calculate SEI2HRD parameters from initial probabilities of developing symptoms and needing hospitalisation, case fatality ratio, and time for transition between different categories.
"""
function SEI2HRDGrowth(birth::Vector{TimeUnitType{U}},
    death::Vector{TimeUnitType{U}}, virus_growth::Vector{TimeUnitType{U}},
    virus_decay::Vector{TimeUnitType{U}}, beta::Vector{TimeUnitType{U}},
    prob_sym::Float64, prob_hosp::Float64, case_fatality_ratio::Float64,
    T_lat::Unitful.Time, T_asym::Unitful.Time, T_sym::Unitful.Time,
    T_hosp::Unitful.Time, T_rec::Unitful.Time) where {U <: Unitful.Units}
    # Prob of death at hospital
    prob_hosp_death = (case_fatality_ratio * prob_hosp)/prob_hosp
    # Prob of death at home
    prob_death = (case_fatality_ratio * (1 - prob_hosp))/(1 - prob_hosp)
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
    return SEI2HRDGrowth{U}(birth, death, virus_growth, virus_decay, beta,
    sigma_1, sigma_2, sigma_hospital, mu_1, mu_2, hospitalisation, death_home, death_hospital)
end

mutable struct EpiParams{U <: Unitful.Units} <: AbstractParams
    births::Vector{TimeUnitType{U}}
    virus_growth::Vector{TimeUnitType{U}}
    virus_decay::Vector{TimeUnitType{U}}
    transition::Matrix{TimeUnitType{U}}
    transition_virus::Matrix{TimeUnitType{U}}
    function EpiParams{U}(births::Vector{TimeUnitType{U}}, virus_growth::Vector{TimeUnitType{U}},
        virus_decay::Vector{TimeUnitType{U}}, transition::Matrix{TimeUnitType{U}}, transition_virus::Matrix{TimeUnitType{U}}) where {U <: Unitful.Units}
        new{U}(births, virus_growth, virus_decay, transition, transition_virus)
    end
end

function transition(params::SISGrowth)
    tmat = zeros(typeof(params.beta), 4, 4)
    tmat[end, :] .+= params.death
    tmat[2, 3] = params.sigma
    vmat = zeros(typeof(params.beta), 4, 4)
    vmat[3, 2] = params.beta
    v_growth = v_decay = fill(0.0 * unit(params.virus_growth), 5)
    v_growth[3] = params.virus_growth
    v_decay[3] = params.virus_decay
  return EpiParams{typeof(unit(params.beta))}(params.birth, v_growth, v_decay, tmat, vmat)
end

function transition(params::SIRGrowth)
    tmat = zeros(typeof(params.beta), 5, 5)
    tmat[4, 3] = params.sigma
    tmat[end, :] .+= params.death
    vmat = zeros(typeof(params.beta), 5, 5)
    vmat[3, 2] = params.beta
    v_growth = v_decay = fill(0.0 * unit(params.virus_growth), 5)
    v_growth[3] = params.virus_growth
    v_decay[3] = params.virus_decay
  return EpiParams{typeof(unit(params.beta))}(params.birth, v_growth, v_decay, tmat, vmat)
end

function transition(params::SEIRGrowth)
    tmat = zeros(typeof(params.beta), 6, 6)
    ordered_transitions = [params.mu, params.sigma]
    from = [3, 4]
    to = [4, 5]
    for i in eachindex(to)
            tmat[to[i], from[i]] = ordered_transitions[i]
    end
    tmat[end, :] .+= params.death
    vmat = zeros(typeof(params.beta), 5, 5)
    vmat[3, 2] = params.beta
    v_growth = v_decay = fill(0.0 * unit(params.virus_growth), 5)
    v_growth[3] = params.virus_growth
    v_decay[3] = params.virus_decay
  return EpiParams{typeof(unit(params.beta))}(params.birth, v_growth, v_decay, tmat, vmat)
end

function transition(params::SEIRSGrowth)
    tmat = zeros(typeof(params.beta), 6, 6)
    ordered_transitions = [params.mu, params.sigma, params.epsilon]
    from = [3, 4, 5]
    to = [4, 5, 2]
    for i in eachindex(to)
            tmat[to[i], from[i]] = ordered_transitions[i]
    end
    tmat[end, :] .+= params.death
    vmat = zeros(typeof(params.beta), 5, 5)
    vmat[3, 2] = params.beta
    v_growth = v_decay = fill(0.0 * unit(params.virus_growth), 5)
    v_growth[3] = params.virus_growth
    v_decay[3] = params.virus_decay
  return EpiParams{typeof(unit(params.beta))}(params.birth, v_growth, v_decay, tmat, vmat)
end

function transition(params::SEI2HRDGrowth)
    ordered_transitions = [params.mu_1, params.mu_2, params.hospitalisation,
    params.sigma_1, params.sigma_2, params.sigma_hospital, params.death_home,
    params.death_hospital]
    from = [3, 4, 5, 4, 5, 6, 5, 6]
    to = [4, 5, 6, 7, 7, 7, 8, 8]
    tmat = zeros(typeof(params.beta[1]), 8, 8)
    for i in eachindex(to)
            tmat[to[i], from[i]] = ordered_transitions[i]
    end
    vmat = zeros(typeof(params.beta[1]), 8, 8)
    vmat[:, 2] .= params.beta
    tmat[end, :] .+= params.death
  return EpiParams{typeof(unit(params.beta[1]))}(params.birth, params.virus_growth, params.virus_decay, tmat, vmat)
end
