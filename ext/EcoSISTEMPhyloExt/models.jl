# SPDX-License-Identifier: LGPL-3.0-or-later

# --- Phylogenetic trait models ------------------------------------------------------------------------
#
# The whole of the phylogenetics: the variance–covariance matrix of a tree, the Brownian/Lambda
# fits over it, and the two fitted-model types themselves. `src/extensions.jl` keeps the function
# stubs, the two abstract supertypes and all the docstrings — `api.md`'s `@autodocs` cannot see into
# an extension, so documentation has to live in the parent. Moving code into a fresh module is what
# compiles it: `varcovar` sat dead since Julia 0.7, calling a long-removed `indmax`, until the move
# forced it through the compiler.

# **The concrete fitted-model types.** Their supertypes are declared in `src/extensions.jl`,
# carrying the documentation, because an extension cannot add a binding to the parent module and a
# struct cannot be stubbed the way a function can. The names are reused deliberately: this module
# refers to the parent qualified (`EcoSISTEM.Brownian`), so the bare name is free here -- exactly
# what `EcoSISTEMMPIExt` does for `MPIEcosystem`.

struct Brownian <: EcoSISTEM.Brownian
    optimum::AbstractArray
    se::AbstractArray
    H::AbstractMatrix
    LL::Float64
end

function Base.show(io::IO, m::Brownian)
    # The digits are a KEYWORD: `round(x, 2)` was Julia 0.6 syntax and has been a `MethodError`
    # since 0.7, so this method threw on every call until it was first tested.
    roundedopts = round.(m.optimum, digits = 2)
    roundedses = round.(m.se, digits = 2)
    roundedLL = round(m.LL, digits = 2)
    return print(io, "Brownian\n",
                 "σ² = $(roundedopts[1]) ($(roundedopts[1] - 2*roundedses[1]) - $(roundedopts[1] + 2*roundedses[1]))",
                 "\n",
                 "z̄₀ = $(roundedopts[2]) ($(roundedopts[2] - 2*roundedses[2]) - $(roundedopts[2] + 2*roundedses[2]))",
                 "\n",
                 "log-likelihood = $roundedLL")
end

function EcoSISTEM.varcovar(tree::AbstractTree)
    tips = collect(nodenamefilter(isleaf, tree))
    root = collect(nodenamefilter(isroot, tree))[1]
    V = zeros(Float64, length(tips), length(tips))
    for i in 1:(length(tips) - 1)
        for j in (i + 1):length(tips)
            V[i, i] = distance(tree, root, tips[i])
            V[j, j] = V[i, i]
            inter = getancestors(tree, tips[i]) ∩ getancestors(tree, tips[j])
            # **Was `indmax`, which has not existed since Julia 0.7** — so `varcovar` threw an
            # `UndefVarError` on every call, and nothing noticed because nothing tests it. Found by
            # moving the function, not caused by it. `argmax` is the documented replacement and
            # returns the index, which is what `inter[common]` below wants.
            common = argmax(map(x -> distance(tree, root, x), inter))
            V[i, j] = distance(tree, root, inter[common])
            V[j, i] = V[i, j]
        end
    end
    return V
end

function EcoSISTEM.fitbrownian(tree::AbstractTree,
                               traits::Vector{F} where {F <: AbstractFloat})
    tips = collect(nodenamefilter(isleaf, tree))
    n = length(tips)
    V = varcovar(tree)
    O = ones(n)
    LL(x) = 1 / 2 * (n * log(2π) + log(abs(det(x[1] * V))) +
             transpose(traits - x[2] * O) * inv(x[1] * V) * (traits - x[2] * O))
    result = optimize(LL, [0.1, 0.1])
    opts = Optim.minimizer(result)
    H = Calculus.hessian(LL, opts)
    se = sqrt.(diag(abs.(inv(H))))
    logL = -LL(opts)
    return Brownian(opts, se, H, logL)
end
