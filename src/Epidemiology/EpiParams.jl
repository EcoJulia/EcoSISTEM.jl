using Unitful

"""
    SIRGrowth{U <: Unitful.Units} <: AbstractParams

Parameter set that houses information on birth and death rates of different classes in an SIR model, as well as growth and decay of virus, and infection/recovery parameters.
"""
mutable struct SIRGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta_force::TimeUnitType{U}
      beta_env::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SIRGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U}, virus_decay::TimeUnitType{U}, beta_force::TimeUnitType{U}, beta_env::TimeUnitType{U}, sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        length(birth) == length(death) || ("Birth and death vector lengths differ")
        new{U}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
    end
end

"""
    SISGrowth{U <: Unitful.Units} <: AbstractParams

Parameter set that houses information on birth and death rates of different classes in an SIS model, as well as growth and decay of virus, and infection/recovery parameters.
"""
mutable struct SISGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta_force::TimeUnitType{U}
      beta_env::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SISGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U},
        virus_decay::TimeUnitType{U}, beta_force::TimeUnitType{U}, beta_env::TimeUnitType{U},
        sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        length(birth) == length(death) || ("Birth and death vector lengths differ")
        new{U}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma)
    end
end

"""
    SEIRGrowth{U <: Unitful.Units} <: AbstractParams

Parameter set that houses information on birth and death rates of different classes in an SEIR model, as well as growth and decay of virus, and infection/incubation/recovery parameters.
"""
mutable struct SEIRGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta_force::TimeUnitType{U}
      beta_env::TimeUnitType{U}
      mu::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SEIRGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U},
        virus_decay::TimeUnitType{U}, beta_force::TimeUnitType{U},
        beta_env::TimeUnitType{U}, mu::TimeUnitType{U}, sigma::TimeUnitType{U}) where {U <: Unitful.Units}
        length(birth) == length(death) || ("Birth and death vector lengths differ")
        new{U}(birth, death, virus_growth, virus_decay, beta_force, beta_env, mu, sigma)
    end
end

"""
    SEIRSGrowth{U <: Unitful.Units} <: AbstractParams

Parameter set that houses information on birth and death rates of different classes in an SEIRS model, as well as growth and decay of virus, and infection/incubation/recovery parameters.
"""
mutable struct SEIRSGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::TimeUnitType{U}
      virus_decay::TimeUnitType{U}
      beta_force::TimeUnitType{U}
      beta_env::TimeUnitType{U}
      mu::TimeUnitType{U}
      epsilon::TimeUnitType{U}
      sigma::TimeUnitType{U}
    function SEIRSGrowth{U}(birth::Vector{TimeUnitType{U}},
        death::Vector{TimeUnitType{U}}, virus_growth::TimeUnitType{U},
        virus_decay::TimeUnitType{U}, beta_force::TimeUnitType{U}, beta_env::TimeUnitType{U},
        mu::TimeUnitType{U}, sigma::TimeUnitType{U}, epsilon::TimeUnitType{U}) where {U <: Unitful.Units}
        length(birth) == length(death) || ("Birth and death vector lengths differ")
        new{U}(birth, death, virus_growth, virus_decay, beta_force, beta_env, mu, sigma, epsilon)
    end
end

"""
    SEI2HRDGrowth{U <: Unitful.Units} <: AbstractParams

Parameter set that houses information on birth and death rates of different classes in an SEI2HRD model, as well as growth and decay of virus, and infection/incubation/hospitalisation/recovery parameters.
"""
mutable struct SEI2HRDGrowth{U <: Unitful.Units} <: AbstractParams
      birth::Vector{TimeUnitType{U}}
      death::Vector{TimeUnitType{U}}
      virus_growth::Vector{TimeUnitType{U}}
      virus_decay::Vector{TimeUnitType{U}}
      beta_force::TimeUnitType{U}
      beta_env::TimeUnitType{U}
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
        virus_decay::Vector{TimeUnitType{U}}, beta_force::TimeUnitType{U}, beta_env::TimeUnitType{U},
        sigma_1::TimeUnitType{U}, sigma_2::TimeUnitType{U}, sigma_hospital::TimeUnitType{U}, mu_1::TimeUnitType{U},
        mu_2::TimeUnitType{U}, hospitalisation::TimeUnitType{U}, death_home::TimeUnitType{U},
        death_hospital::TimeUnitType{U}) where {U <: Unitful.Units}
        length(birth) == length(death) || ("Birth and death vector lengths differ")
        length(virus_growth) == length(virus_decay) || ("Virus growth and decay vector lengths differ")
        (beta_force != 0/day) & (beta_env != 0/day) || warning("Transmission rate is zero.")
        new{U}(birth, death, virus_growth, virus_decay, beta_force, beta_env, sigma_1, sigma_2, sigma_hospital, mu_1, mu_2, hospitalisation, death_home, death_hospital)
    end
end

"""
SEI2HRDGrowth(birth::Vector{TimeUnitType{U}},
    death::Vector{TimeUnitType{U}}, beta::TimeUnitType{U}, prob_sym::Float64, prob_hosp::Float64, case_fatality_ratio::Float64, T_lat::Unitful.Time, T_asym::Unitful.Time, T_sym::Unitful.Time, T_hosp::Unitful.Time, T_rec::Unitful.Time) where {U <: Unitful.Units}

Function to calculate SEI2HRD parameters from initial probabilities of developing symptoms and needing hospitalisation, case fatality ratio, and time for transition between different categories.
"""
function SEI2HRDGrowth(birth::Vector{TimeUnitType{U}},
    death::Vector{TimeUnitType{U}}, virus_growth::Vector{TimeUnitType{U}},
    virus_decay::Vector{TimeUnitType{U}}, beta_force::TimeUnitType{U}, beta_env::TimeUnitType{U},
    prob_sym::Float64, prob_hosp::Float64, cfr_home::Float64, cfr_hosp::Float64,
    T_lat::Unitful.Time, T_asym::Unitful.Time, T_sym::Unitful.Time,
    T_hosp::Unitful.Time, T_rec::Unitful.Time) where {U <: Unitful.Units}
    # Exposed -> asymptomatic
    mu_1 = 1/T_lat
    # Asymptomatic -> symptomatic
    mu_2 = prob_sym * 1/T_asym
    # Symptomatic -> hospital
    hospitalisation = prob_hosp * 1/T_sym
    # Asymptomatic -> recovered
    sigma_1 = (1 - prob_sym) * 1/T_asym
    # Symptomatic -> recovered
    sigma_2 = (1 - prob_hosp) * (1 - cfr_home) * 1/T_rec
    # Hospital -> recovered
    sigma_hospital = (1 - cfr_hosp) * 1/T_hosp
    # Symptomatic -> death
    death_home = cfr_home * 2/T_hosp
    # Hospital -> death
    death_hospital = cfr_hosp * 1/T_hosp
    return SEI2HRDGrowth{U}(birth, death, virus_growth, virus_decay, beta_force, beta_env,
    sigma_1, sigma_2, sigma_hospital, mu_1, mu_2, hospitalisation, death_home, death_hospital)
end

"""
    EpiParams{U <: Unitful.Units} <: AbstractParams

Parameter set for any epi model type, which stores information on birth, virus generation and decay probabilities, as well as matrices for transitions between different states. `transition` houses straightforward transition probabilities between classes, whereas `transition_virus` houses probabilities that should be multiplied by the amount of virus in the system, such as infection transitions.
"""
mutable struct EpiParams{U <: Unitful.Units} <: AbstractParams
    births::Vector{TimeUnitType{U}}
    virus_growth::Vector{TimeUnitType{U}}
    virus_decay::Vector{TimeUnitType{U}}
    transition::Matrix{TimeUnitType{U}}
    transition_force::Matrix{TimeUnitType{U}}
    transition_virus::Matrix{TimeUnitType{U}}
    function EpiParams{U}(births::Vector{TimeUnitType{U}}, virus_growth::Vector{TimeUnitType{U}},
        virus_decay::Vector{TimeUnitType{U}}, transition::Matrix{TimeUnitType{U}}, transition_force::Matrix{TimeUnitType{U}}, transition_virus::Matrix{TimeUnitType{U}}) where {U <: Unitful.Units}
        size(transition, 1) == size(transition, 2) || error("Transition matrix should be square.")
        size(transition, 1) == size(transition_virus, 1) || error("Transition matrices should match dimensions.")
        size(transition_force, 1) == size(transition_virus, 1) || error("Transition matrices should match dimensions.")
        new{U}(births, virus_growth, virus_decay, transition, transition_force, transition_virus)
    end
end

"""
    transition(params::SISGrowth)

Function to create transition matrix from SIS parameters and return an `EpiParams` type that can be used by the model update.
"""
function transition(params::SISGrowth)
    tmat = zeros(typeof(params.beta_force), 4, 4)
    tmat[end, :] .+= params.death
    tmat[2, 3] = params.sigma
    vmat = zeros(typeof(params.beta_force), 4, 4)
    vfmat = zeros(typeof(params.beta_force), 4, 4)
    vfmat[3, 2] = params.beta_force
    vmat[3, 2] = params.beta_env
    v_growth = fill(0.0 * unit(params.virus_growth), 5)
    v_decay = fill(0.0 * unit(params.virus_growth), 5)
    v_growth[3] = params.virus_growth
    v_decay[3] = params.virus_decay
  return EpiParams{typeof(unit(params.beta_force))}(params.birth, v_growth, v_decay, tmat, vfmat, vmat)
end

"""
    transition(params::SIRGrowth)

Function to create transition matrix from SIR parameters and return an `EpiParams` type that can be used by the model update.
"""
function transition(params::SIRGrowth)
    tmat = zeros(typeof(params.beta_force), 5, 5)
    tmat[4, 3] = params.sigma
    tmat[end, :] .+= params.death
    vmat = zeros(typeof(params.beta_force), 5, 5)
    vfmat = zeros(typeof(params.beta_force), 5, 5)
    vfmat[3, 2] = params.beta_force
    vmat[3, 2] = params.beta_env
    v_growth = fill(0.0 * unit(params.virus_growth), 5)
    v_decay = fill(0.0 * unit(params.virus_growth), 5)
    v_growth[3] = params.virus_growth
    v_decay[3] = params.virus_decay
  return EpiParams{typeof(unit(params.beta_force))}(params.birth, v_growth, v_decay, tmat, vfmat, vmat)
end

"""
    transition(params::SEIRGrowth)

Function to create transition matrix from SEIR parameters and return an `EpiParams` type that can be used by the model update.
"""
function transition(params::SEIRGrowth)
    tmat = zeros(typeof(params.beta_force), 6, 6)
    ordered_transitions = [params.mu, params.sigma]
    from = [3, 4]
    to = [4, 5]
    for i in eachindex(to)
            tmat[to[i], from[i]] = ordered_transitions[i]
    end
    tmat[end, :] .+= params.death
    vmat = zeros(typeof(params.beta_force), 6, 6)
    vfmat = zeros(typeof(params.beta_force), 6, 6)
    vfmat[3, 2] = params.beta_force
    vmat[3, 2] = params.beta_env
    v_growth = fill(0.0 * unit(params.virus_growth), 6)
    v_decay = fill(0.0 * unit(params.virus_growth), 6)
    v_growth[4] = params.virus_growth
    v_decay[4] = params.virus_decay
  return EpiParams{typeof(unit(params.beta_force))}(params.birth, v_growth, v_decay, tmat, vfmat, vmat)
end

"""
    transition(params::SEIRSGrowth)

Function to create transition matrix from SEIRS parameters and return an `EpiParams` type that can be used by the model update.
"""
function transition(params::SEIRSGrowth)
    tmat = zeros(typeof(params.beta_force), 6, 6)
    ordered_transitions = [params.mu, params.sigma, params.epsilon]
    from = [3, 4, 5]
    to = [4, 5, 2]
    for i in eachindex(to)
            tmat[to[i], from[i]] = ordered_transitions[i]
    end
    tmat[end, :] .+= params.death
    vmat = zeros(typeof(params.beta_force), 6, 6)
    vfmat = zeros(typeof(params.beta_force), 6, 6)
    vfmat[3, 2] = params.beta_force
    vmat[3, 2] = params.beta_env
    v_growth = fill(0.0 * unit(params.virus_growth), 6)
    v_decay = fill(0.0 * unit(params.virus_growth), 6)
    v_growth[4] = params.virus_growth
    v_decay[4] = params.virus_decay
  return EpiParams{typeof(unit(params.beta_force))}(params.birth, v_growth, v_decay, tmat, vfmat, vmat)
end

"""
    transition(params::SEI2HRDGrowth)

Function to create transition matrix from SEI2HRD parameters and return an `EpiParams` type that can be used by the model update.
"""
function transition(params::SEI2HRDGrowth)
    ordered_transitions = (incubation_period = params.mu_1, symptoms_develop = params.mu_2, symptoms_worsen = params.hospitalisation, recovery_asymptomatic = params.sigma_1, recovery_symptomatic = params.sigma_2, recovery_hospital = params.sigma_hospital, death_symptomatic = params.death_home,
    death_hospitalised = params.death_hospital)
    from = [3, 4, 5, 4, 5, 6, 5, 6]
    to = [4, 5, 6, 7, 7, 7, 8, 8]
    tmat = zeros(typeof(params.beta_force), 8, 8)
    for i in eachindex(to)
            tmat[to[i], from[i]] = ordered_transitions[i]
    end
    vmat = zeros(typeof(params.beta_force), 8, 8)
    vfmat = zeros(typeof(params.beta_force), 8, 8)
    vfmat[3, 2] = params.beta_force
    vmat[3, 2] = params.beta_env
    tmat[end, :] .+= params.death
  return EpiParams{typeof(unit(params.beta_force))}(params.birth, params.virus_growth, params.virus_decay, tmat, vfmat, vmat)
end
