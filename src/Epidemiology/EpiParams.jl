using Unitful
using LinearAlgebra
using DataFrames


"""
    EpiParams{U <: Unitful.Units} <: AbstractParams

Parameter set for any epi model type, which stores information on birth, virus generation and decay probabilities, as well as matrices for transitions between different states. `transition` houses straightforward transition probabilities between classes, whereas `transition_virus` houses probabilities that should be multiplied by the amount of virus in the system, such as infection transitions.
"""
mutable struct EpiParams{U <: Unitful.Units} <: AbstractParams
    births::Vector{TimeUnitType{U}}
    virus_growth::Vector{TimeUnitType{U}}
    virus_decay::TimeUnitType{U}
    transition::Matrix{TimeUnitType{U}}
    transition_force::Matrix{TimeUnitType{U}}
    transition_virus::Matrix{TimeUnitType{U}}
    freq_vs_density_force::Float64
    freq_vs_density_env::Float64
    age_mixing::Matrix{Float64}
    env_virus_scale::Float64
    function EpiParams{U}(births::Vector{TimeUnitType{U}}, virus_growth::Vector{TimeUnitType{U}},
        virus_decay::TimeUnitType{U}, transition::Matrix{TimeUnitType{U}}, transition_force::Matrix{TimeUnitType{U}}, transition_virus::Matrix{TimeUnitType{U}}, freq_vs_density_force::Float64, freq_vs_density_env::Float64, age_mixing::Matrix{Float64}, env_virus_scale::Float64 = 1.0) where {U <: Unitful.Units}
        size(transition, 1) == size(transition, 2) || error("Transition matrix should be square.")
        size(transition, 1) == size(transition_virus, 1) || error("Transition matrices should match dimensions.")
        size(transition_force, 1) == size(transition_virus, 1) || error("Transition matrices should match dimensions.")
        new{U}(births, virus_growth, virus_decay, transition, transition_force, transition_virus, freq_vs_density_force, freq_vs_density_env, age_mixing, env_virus_scale)
    end
end


function create_transition_matrix(params::NamedTuple, paramDat::DataFrame, age_categories::Int64, nclasses::Int64)
    # Set up number of classes etc
    tm_size = age_categories * nclasses
    cat_idx = reshape(1:tm_size, age_categories, nclasses)

    # Set up transition matrix
    tmat = zeros(eltype(params.beta_env), tm_size, tm_size)

    # Death
    for i in 1:nclasses
        dmat = @view tmat[cat_idx[:, end], cat_idx[:, i]]
        dmat[diagind(dmat)].= params.death[cat_idx[:, i]]
    end

    # Other transitions
    ordered_transitions = paramDat[!, :prob]
    from = paramDat[!, :from_ind]
    to = paramDat[!, :to_ind]
    for i in eachindex(to)
            view_mat = @view tmat[cat_idx[:, to[i]], cat_idx[:, from[i]]]
            view_mat[diagind(view_mat)] .= ordered_transitions[i]
    end
    return tmat
end

function create_virus_matrix(beta::Vector{TimeUnitType{U}}, age_categories::Int64, nclasses::Int64) where U <: Unitful.Units
    vm_size = age_categories * nclasses
    cat_idx = reshape(1:vm_size, age_categories, nclasses)
    vmat = zeros(eltype(beta), vm_size, vm_size)
    bmat = @view vmat[cat_idx[:, 2], cat_idx[:, 1]]
    bmat[diagind(bmat)] .= beta
    return vmat
end
function create_virus_matrix(beta::TimeUnitType{U}, age_categories::Int64, nclasses::Int64) where U <: Unitful.Units
    return create_virus_matrix([beta], age_categories, nclasses)
end

function create_virus_vector(virus_growth::Vector{TimeUnitType{U}}, age_categories::Int64, nclasses::Int64, inf_cat::Vector{Int64}) where U <: Unitful.Units
    vm_size = age_categories * nclasses
    cat_idx = reshape(1:vm_size, age_categories, nclasses)
    v_growth = zeros(eltype(virus_growth), vm_size)
    v_growth[cat_idx[:, inf_cat]] .= virus_growth
    return v_growth
end

function create_virus_vector(virus_growth::Matrix{TimeUnitType{U}}, age_categories::Int64, nclasses::Int64, inf_cat::Vector{Int64}) where U <: Unitful.Units
    vm_size = age_categories * nclasses
    cat_idx = reshape(1:vm_size, age_categories, nclasses)
    v_growth = zeros(eltype(virus_growth), vm_size)
    v_growth[cat_idx[:, inf_cat]] .= virus_growth
    return v_growth
end

function create_virus_vector(virus_growth::TimeUnitType{U}, age_categories::Int64, nclasses::Int64, inf_cat::Vector{Int64}) where U <: Unitful.Units
    return create_virus_vector([virus_growth], age_categories, nclasses, inf_cat)
end

"""
    transition(params::NamedTuple, paramDat::DataFrame, nclasses::Int64, inf_cat = [2], age_categories = 1)

Function to create transition matrix from SIS parameters and return an `EpiParams` type that can be used by the model update.
"""
function transition(params::NamedTuple, paramDat::DataFrame, nclasses::Int64, inf_cat = [2], age_categories = 1)
    # Check for missing params
    if !haskey(params, :freq_vs_density_force)
        params = (; params..., freq_vs_density_force = 1.0)
    end
    if !haskey(params, :freq_vs_density_env)
        params = (; params..., freq_vs_density_env = 1.0)
    end
    if !haskey(params, :env_scale)
        params = (; params..., env_scale = 1.0)
    end
    if !haskey(params, :age_mixing)
        params = (; params..., age_mixing = fill(1.0, 1, 1))
    end

    tmat = create_transition_matrix(params, paramDat, age_categories, nclasses)

    # Env virus matrix
    vmat = create_virus_matrix(params.beta_env, age_categories, nclasses)

    # Force matrix
    vfmat = create_virus_matrix(params.beta_force, age_categories, nclasses)

    # Virus growth and decay
    v_growth = create_virus_vector(params.virus_growth, age_categories, nclasses, inf_cat)

  return EpiParams{typeof(unit(params.beta_force[1]))}(params.birth[1:end], v_growth, params.virus_decay, tmat, vfmat, vmat, params.freq_vs_density_force, params.freq_vs_density_env, params.age_mixing, params.env_scale)
end
