using Phylo
using Calculus

import Base.show
function varcovar(tree::AbstractTree)
    tips = collect(nodenamefilter(isleaf, tree))
    root = collect(nodenamefilter(isroot, tree))[1]
    V = zeros(Float64, length(tips), length(tips))
    for i in 1:(length(tips) - 1)
        for j in (i+1):length(tips)
            V[i, i] =  distance(tree, root, tips[i])
            V[j, j]= V[i,i]
            inter = getancestors(tree, tips[i]) ∩ getancestors(tree, tips[j])
            common = indmax(map(x-> distance(tree, root, x), inter))
            V[i, j] = distance(tree, root, inter[common])
            V[j, i] = V[i, j]
        end
    end
    return V
end

mutable struct Brownian
    optimum::AbstractArray
    se::AbstractArray
    H::AbstractMatrix
    LL::Float64
end

function show(io::IO, m::Brownian)
    roundedopts = round.(m.optimum, 2)
    roundedses = round.(m.se, 2)
    roundedLL = round(m.LL, 2)
    return print(io, "σ² = $(roundedopts[1]) ($(roundedopts[1] - 2*roundedses[1]) - $(roundedopts[1] + 2*roundedses[1]))", "\n",
    "z̄₀ = $(roundedopts[2]) ($(roundedopts[2] - 2*roundedses[2]) - $(roundedopts[2] + 2*roundedses[2]))","\n",
    "log-likelihood = $roundedLL")
end

mutable struct Lambda
    optimum::AbstractArray
    se::AbstractArray
    H::AbstractMatrix
    LL::Float64
end

function show(io::IO, m::Lambda)
    roundedopts = round.(m.optimum, 2)
    roundedses = round.(m.se, 2)
    roundedLL = round(m.LL, 2)
    return print(io, "σ² = $(roundedopts[1]) ($(roundedopts[1] - 2*roundedses[1]) - $(roundedopts[1] + 2*roundedses[1]))", "\n",
    "z̄₀ = $(roundedopts[2]) ($(roundedopts[2] - 2*roundedses[2]) - $(roundedopts[2] + 2*roundedses[2]))", "\n",
    "λ = $(roundedopts[3]) ($(roundedopts[3] - 2*roundedses[3]) - $(roundedopts[3] + 2*roundedses[3]))","\n",
    "log-likelihood = $roundedLL")
end

@init @require Optim="429524aa-4258-5aef-a3af-852621145aeb" @eval begin
    using .Optim

    function fitBrownian(tree::AbstractTree, traits::Vector{F} where F <: AbstractFloat)
        tips= collect(nodenamefilter(isleaf, tree))
        n = length(tips)
        V = varcovar(tree)
        O = ones(n)
        LL(x) = 1/2 * (n * log(2π) + log(abs(det(x[1] * V))) +
        transpose(traits - x[2] * O) * inv(x[1] * V) * (traits - x[2] * O))
        result = optimize(LL, [0.1, 0.1])
        opts = Optim.minimizer(result)
        H = Calculus.hessian(LL, opts)
        se = sqrt.(diag(abs.(inv(H))))
        logL = -LL(opts)
        return Brownian(opts, se, H, logL)
    end

    function fitLambda(tree::AbstractTree, traits::Vector{F} where F <: AbstractFloat)
        tips= collect(nodenamefilter(isleaf, tree))
        n = length(tips)
        V = varcovar(tree)
        function LL(x, n, V, traits)
            O = ones(n)
            dV = diagm(diag(V))
            V = (x[3] * (V - dV) + dV)
            return 1/2 * (n * log(2π) + log(abs(det(x[1] * V))) +
            transpose(traits - x[2] * O) * inv(x[1] * V) * (traits - x[2] * O))
        end
        result = optimize(x -> LL(x, n, V, traits), [0.1, 0.1, 0.2],
        [exp(-100), -Inf, 0.0],
        [Inf, Inf, 1.0])
        opts = Optim.minimizer(result)
        H = Calculus.hessian(x -> LL(x, n, V, traits), opts)
        se = sqrt.(diag(abs.(inv(H))))
        logL = -LL(opts, n, V, traits)
        return Lambda(opts, se, H, logL)
    end

    export fitBrownian, fitLambda
end
