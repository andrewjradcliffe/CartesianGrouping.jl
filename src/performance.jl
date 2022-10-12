#
# Date created: 2022-10-12
# Author: aradclif
#
#
############################################################################################
# Performance tests
function rand_onehot(::Type{T}, M::Integer, N::Integer) where {T<:Real}
    r = 1:N
    A = zeros(T, M, N)
    for i ∈ axes(A, 1)
        A[i, rand(r)] = one(T)
    end
    A
end
rand_onehot(M::Integer, N::Integer) = rand_onehot(Float64, M, N)
# Obviously, this could be done in a more efficient way
sprand_onehot(::Type{T}, M::Integer, N::Integer) where {T<:Real} = sparse(rand_onehot(T, M, N))
sprand_onehot(M::Integer, N::Integer) = sprand_onehot(Float64, M, N)

x = rand_onehot(100, 100)
@timev onehot_toΓ(x);
@timev onehot_to(x);
