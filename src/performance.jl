#
# Date created: 2022-10-12
# Author: aradclif
#
#
############################################################################################
# Performance tests
function rand_onehot(::Type{T}, m::Integer, n::Integer) where {T<:Real}
    r = 1:n
    A = zeros(T, m, n)
    for i ∈ axes(A, 1)
        A[i, rand(r)] = one(T)
    end
    A
end
rand_onehot(m::Integer, n::Integer) = rand_onehot(Float64, m, n)

function sprand_onehot(::Type{T}, m::Integer, n::Integer) where {T<:Real}
    I = 1:m
    J = rand(1:n, m)
    V = ones(T, m)
    sparse(I, J, V, m, n)
end
sprand_onehot(m::Integer, n::Integer) = sprand_onehot(Float64, m, n)

x = rand_onehot(100, 100)
@timev onehot_toΓ(x);
@timev onehot_to(x);

y1 = rand_onehot(Int, 10000, 10);
Γ1 = onehot_toΓ(y1);
y2 = rand_onehot(Int, 10000, 10);
Γ2 = onehot_toΓ(y2);
y3 = rand_onehot(Int, 10000, 10);
Γ3 = onehot_toΓ(y3);

jj1 = onehot_tojj(y1);
jj2 = onehot_tojj(y2);
jj3 = onehot_tojj(y3);

function getall_cis(Ŝ::AbstractVector, ξ)
    Γ = [Vector{Int}() for _ = 1:length(Ŝ)]
    for i ∈ eachindex(ξ)
        for j ∈ eachindex(Ŝ)
            isequal(Ŝ[j], ξ[i]) && (push!(Γ[j], i); break)
        end
    end
    Γ
end

ξ = rand(1:100, 100000);
Ŝ = collect(1:100);
@timev getall_cartesianindexsets(Ŝ, ξ);
@timev getall_cis(Ŝ, ξ);
getall_cartesianindexsets(Ŝ, ξ) == getall_cis(Ŝ, ξ)

bg1 = BasicGroup(jj1);
bg2 = BasicGroup(jj2);
bg3 = BasicGroup(jj3);

bg1 × bg2 × bg3
@timev ncross(bg1, bg2, bg3, bg1, bg2, bg3, bg3, bg2, bg3);
