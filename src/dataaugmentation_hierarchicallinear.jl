#
# Date created: 2021-07-27
# Author: aradclif
#
#
############################################################################################
#### p. 59-65, 2021-07-27
## Functions for augmenting groupings into a generic framework for hierarchical linear
## models; also works for hierarchical generalized linear models.
############################################################################################
"""
Construct a `jj` from `K_beta` groups of hierarchically modeled parameters,
thereby permitting their treatment as a single augmented parameter vector.
The indexing achieved associates each element of the augmented parameter vector
to its respective hierarchical variance, such that the covariance matrix
can be constructed as: diag(σ_beta[jj_beta]).
"""
function to_jj_beta(Js::Vector{<:Integer})
    K_beta = length(Js)
    pieces = Vector{Vector{Int}}(undef, K_beta)
    for k = 1:K_beta
        pieces[k] = k * ones(Int, Js[k])
    end
    jj_beta = vcat(pieces...)
    return jj_beta
end
to_jj_beta(BGs::Vector{<:BasicGroup}) = to_jj_beta(map(x -> getproperty(x, :J), BGs))

"""
Construct a `W_beta` from `K_beta` groups of hierarchically modeled parameters,
thereby permitting their treatment as a single augmented parameter vector.
The indexing achieved associates each element of the augmented parameter vector
to its respective hierarchical variance, such that the covariance matrix
can be constructed as: diag(W_beta * σ_beta).
"""
function to_W_beta(Js::Vector{<:Integer})
    K_beta = length(Js)
    J = sum(Js)
    W_beta = zeros(J, K_beta)
    pos = 1
    for k = 1:K_beta
        W_beta[pos:(pos + Js[k] - 1), k] .= 1.0
        pos += Js[k]
    end
    return W_beta
end
to_W_beta(BGs::Vector{<:BasicGroup}) = to_W_beta(map(x -> getproperty(x, :J), BGs))

"""
Construct the augmented matrix which is used to represent the hierarchical
groupings for linear algebra purposes.
"""
function to_x_beta(BGs::Vector{<:BasicGroup})
    K_beta = length(BGs)
    mats = Vector{Matrix{Int}}(undef, K_beta)
    for k = 1:K_beta
        mats[k] = onehot_dense(BGs[k])
    end
    x = hcat(convert.(Matrix{Float64}, mats)...)
    return x
end
