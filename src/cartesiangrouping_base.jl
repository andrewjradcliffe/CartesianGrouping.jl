#
# Date created: 2021-07-07
# Author: aradclif
#
#
############################################################################################
#### Order of sections
# - construct_
# - Cartesian index sets
# - Conversions to/from Γ, jj, onehot
# - Types
# - misc
################################################################
# using SparseArrays
################################################################
#### Construction of Gⱼₘ, Γₘ, allΓ
# construct_Gⱼₘ(Gⱼ::Vector{Int}, P::Int, m::Int) = [p + (m - 1)P for p in Gⱼ]
function construct_Gⱼₘ(Gⱼ::Vector{Int}, P::Int, m::Int)
    Gⱼₘ = [p + (m - 1)P for p in Gⱼ]
    return Gⱼₘ
end

function construct_Γₘ(Γ::Vector{Vector{Int}}, P::Int, m::Int)
    Γₘ = [construct_Gⱼₘ(Gⱼ, P, m) for Gⱼ in Γ]
    return Γₘ
end

function construct_allΓ(Γ::Vector{Vector{Int}}, P::Int, M::Int)
    allΓ = [construct_Γₘ(Γ, P, m) for m = 1:M]
    return allΓ
end

# Construction of G̃ⱼ, Γ̃ from Gⱼₘ's and Γₘ's
function construct_Γⱼ(allΓ::Vector{Vector{Vector{Int}}}, j::Int)
    Γⱼ = [Γₘ[j] for Γₘ in allΓ]
    return Γⱼ
end

function construct_G̃ⱼ(allΓ::Vector{Vector{Vector{Int}}}, j::Int)
    G̃ⱼ = union(construct_Γⱼ(allΓ, j)...)
    return G̃ⱼ
end

function construct_Γ̃(allΓ::Vector{Vector{Vector{Int}}}, J::Int)
    Γ̃ = [construct_G̃ⱼ(allΓ, j) for j = 1:J]
    return Γ̃
end
#### Special functions for augmenting non-simple cases: p. 52-54, 2021-07-13
# ρ = [P₁, P₂, …, Pₘ], Pₘ = sum(length.(Γ)) for the mᵗʰ Γ being augmented
## Discussion
# The other, much simpler, way to faithfully construct a BasicGroup is to augment
# together the respective jj's, then construct BasicGroup(jj).
# Though the result is the same, the 3 `constructspecial` functions may be useful
# if one desires allΓ. Incidentally, one reaches allΓ, one uses the same machinery as
# the non-special case for constructing G̃ⱼ, Γ̃.
function constructspecial_Gⱼₘ(Gⱼ::Vector{<:Integer}, ρ::Vector{<:Integer}, m::Integer)
    Gⱼₘ = [p + sum(ρ[1:(m - 1)]) for p in Gⱼ]
    return Gⱼₘ
end
function constructspecial_Γₘ(Γ::Vector{<:Vector{<:Integer}}, ρ::Vector{<:Integer}, m::Integer)
    Γₘ = [constructspecial_Gⱼₘ(Gⱼ, ρ, m) for Gⱼ in Γ]
    return Γₘ
end
function constructspecial_allΓ(catΓ::Vector{<:Vector{<:Vector{<:Integer}}}, ρ)
    allΓ = [constructspecial_Γₘ(Γ, ρ, m) for (m, Γ) in enumerate(catΓ)]
    return allΓ
end
ρ(catΓ) = sum.(map(x -> length.(x), catΓ))
# Testing
# cells₁ = ["a_v1", "a_w1", "a_v2", "a_v3", "b_v1", "b_v2", "b=v2", "c_v1"]
# Γ₁ = getall_cartesianindexsets(regs_family, cells₁)
# jj₁ = indexsets_tojj(Γ₁)
# cells₂ = ["b_v1", "b_v2", "b_v3", "c_v1", "c=v1", "c_v2"]
# Γ₂ = getall_cartesianindexsets(regs_family, cells₂)
# jj₂ = indexsets_tojj(Γ₂)
# cells₃ = ["a=v1", "a=w1", "b=v1", "b=v2", "b=v3"]
# Γ₃ = getall_cartesianindexsets(regs_family, cells₃)
# jj₃ = indexsets_tojj(Γ₃)
# catΓ = [Γ₁, Γ₂, Γ₃]
# allΓ = constructspecial_allΓ(catΓ, sum.(map(x -> length.(x), catΓ)))
# Ĵ = 3  # number of groups which appear in Ŝ -- known from length(regs_family)
# jj = indexsets_tojj(construct_Γ̃(allΓ, Ĵ))
# catjj = [jj₁; jj₂; jj₃]
# jj == catjj
################################################################
#### Some default grouping methods: p.22
function indexsets_groupbym(P::Int, M::Int)
    # Note! Γ for grouping by m only is already a Γ̃
    Γ = [[p + (m - 1)P for p = 1:P] for m = 1:M]
    return Γ
end
function indexsets_groupbyp(P::Int)
    Γ = [[p] for p = 1:P]
    return Γ
end
################################################################
#### Essential functions: p. 25-26
# p. 33, 2021-07-08: ξ ≡ cartpairs, Ŝ ≡ cartset
# function construct_cartpairs(jj₁::AbstractVector{<:Integer}, jj₂::AbstractVector{<:Integer})
#     cartpairs = [(jj₁[p], jj₂[p]) for p = 1:length(jj₁)] # p in eachindex(jj₁)
#     # cartpairs = CartesianIndex.(jj₁, jj₂)
#     return cartpairs
# end
function construct_cartpairs(jj₁::AbstractVector{<:Integer}, jj₂::AbstractVector{<:Integer})
    Tₒ = promote_type(S, T)
    ξ = Vector{Tuple{Tₒ, Tₒ}}(undef, length(jj₁))
    @inbounds for i ∈ eachindex(jj₁, jj₂, ξ)
        ξ[i] = (jj₁[i], jj₂[i])
    end
    ξ
end

# 2021-07-26: varargs version -- albeit, it is actually faster in most cases to just
# compute the crosses via repeated binary operations.
function construct_cartpairs(jj₁::AbstractVector{<:Integer}, jj₂::AbstractVector{<:Integer},
                             jj₃...)
    Tₒ = promote_type(S, T, W)
    ξ = Vector{NTuple{N+2, Tₒ}}(undef, length(jj₁))
    @inbounds for i ∈ eachindex(jj₁, jj₂, ξ, jjs...)
        ξ[i] = (jj₁[i], jj₂[i], ntuple(k -> jjs[k][i], Val(N))...)
    end
    ξ
end
# each Gⱼ : {p | ξₚ = Ŝⱼ}, Ŝⱼ ≡ the jᵗʰ element in Ŝ
# function getall_cartesianindexsets(cartset, cartpairs)
#     Γ = [findall(x -> isequal(id, x), cartpairs) for id in cartset]
# end
# # in fact, this is the conversion: jj_to_Γ -- see note p. 40 for elaboration
# getall_cartesianindexsets(cartpairs) =
#     getall_cartesianindexsets(unique!(sort(cartpairs)), cartpairs)
#### p. 55-57, 2021-07-13
function getall_cartesianindexsets(F̂::Vector{T} where {T<:Function}, ξ)
    Γ = [findall(fⱼ, ξ) for fⱼ in F̂]
end
function getall_cartesianindexsets(Ŝ::AbstractVector, ξ)
    Γ = [findall(ξₚ -> isequal(ξₚ, Ŝⱼ), ξ) for Ŝⱼ in Ŝ] # or, simply: findall(isequal(Ŝⱼ), ξ)
end
function getall_cartesianindexsets(Ŝ::Vector{Regex}, ξ)
    Γ = [findall(ξₚ -> occursin(Ŝⱼ, ξₚ), ξ) for Ŝⱼ in Ŝ]
end
# The case of Sⱼ ⊂ Ŝⱼ
# getall_cartesianindexsets(ξ) = getall_cartesianindexsets(sort!(unique(ξ)), ξ)
getall_cartesianindexsets(ξ) = getall_cartesianindexsets(unique(ξ), ξ)
#### Checks of necessary conditions
# If it meets both, then it can be usd at will. If it does not meet both
# necessary conditions, then either re-specifty F̂ or Ŝ; or construct the
# subset Ξ = ξ[union(Γ...)] and repeat the procedure on Ξ
function necessarycond1(ξ, Γ::Vector{<:Vector{<:Integer}})
    issetequal(LinearIndices(ξ), union(Γ...))
end
function necessarycond1(ξ, Γ::Vector{<:Vector{<:CartesianIndex}})
    issetequal(CartesianIndices(ξ), union(Γ...))
end
function necessarycond2(Γ)
    Ĵ = length(Γ)
    for jₛ = 1:Ĵ
        Gⱼ = Γ[jₛ]
        for l in (j for j = 1:Ĵ if j != jₛ)
            if !isdisjoint(Gⱼ, Γ[l])
                return false
            end
        end
    end
    return true
end
# Testing
# ξ = [13:24;];
# ξnew = reshape(ξ, 3, 4)
# f1(x) = 13 ≤ x ≤ 15
# f2(x) = 16 ≤ x ≤ 18
# f3(x) = 19 ≤ x ≤ 21
# f4(x) = 22 ≤ x ≤ 24
# fvec = [f1, f2, f3, f4]
# # linear indexed
# Γ = getall_cartesianindexsets(fvec, ξ)
# necessarycond1(ξ, Γ)
# necessarycond2(Γ)
# # cartesian indexed
# Γnew = getall_cartesianindexsets(fvec, ξnew)
# necessarycond1(ξnew, Γnew)
# necessarycond2(Γnew)
################################################################
#### Essential functions: p. 18-20
function indexsets_tojj(Γ::Vector{Vector{Int}})
    J = length(Γ)
    jj = Vector{Int}(undef, sum(length, Γ))
    for j = 1:J
        for n in Γ[j]
            @inbounds jj[n] = j
        end
    end
    return jj
end
# indexsets_tojj(Γ) = indexsets_tojj(Γ, sum(length, Γ))
#### jj -- routes to/from one-hot matrices
function onehot_sparse(jj::AbstractVector{<:Integer}, P::Integer)
    return sparse([p for p = 1:P], jj, ones(Int, P))
end
onehot_sparse(jj::AbstractVector{<:Integer}) = onehot_sparse(jj, length(jj))
onehot_sparse(Γ::AbstractVector{<:AbstractVector{<:Integer}}) =
    onehot_sparse(indexsets_tojj(Γ))

function onehot_dense(jj::AbstractVector{<:Integer}, P::Integer, J::Integer)
    X = zeros(Int, P, J)
    for p = 1:P
        @inbounds X[p, jj[p]] = 1
    end
    return X
end
onehot_dense(jj::AbstractVector{<:Integer}) = onehot_dense(jj, length(jj), length(unique(jj)))
onehot_dense(Γ::AbstractVector{<:AbstractVector{<:Integer}}) =
    onehot_dense(indexsets_tojj(Γ))
#### Conversion of one-hot matrix to Γ
# dense matrices
function denseonehotcol_toGⱼ(xⱼ)
    Gⱼ = findall(!iszero, xⱼ)
    return Gⱼ
end
function onehot_toΓ(X::AbstractMatrix{<:Real})
    Γ = [denseonehotcol_toGⱼ(xⱼ) for xⱼ in eachcol(X)]
    return Γ
end
onehot_to(X::AbstractMatrix{<:Real}) = map(Base.Fix1(findall, !iszero), eachcol(X))

# sparse matrices
function sparseonehotcol_toGⱼ(rowvals, nzrange)
    Gⱼ = [rowvals[i] for i in nzrange]
    return Gⱼ
end
function onehot_toΓ(X::AbstractSparseMatrix{<:Real, <:Integer})
    Γ = [sparseonehotcol_toGⱼ(rowvals(X), nzrange(X, j)) for j = 1:X.n]
    # Alternative
    # Γ = [[rowvals[i] for i in nzrange(X, j)] for j = 1:X.n]
    return Γ
end
onehot_tojj(X) = indexsets_tojj(onehot_toΓ(X))
#### Conversion of Γ, jj to ELM_y, n_y: p. 44
function ELMy_ny(Γ::AbstractVector{<:AbstractVector{<:Integer}})
    J = length(Γ)
    ny = length.(Γ)
    N = maximum(ny)
    ELMy = zeros(Int, N, J)
    for j in 1:J
        @inbounds ELMy[1:ny[j], j] .= Γ[j]
    end
    return ELMy, ny
end
################################################################
#### BasicGroup, AugmentedGroup Types: p. 30-32
struct BasicGroup{T<:AbstractVector{<:AbstractVector{<:Integer}},
                  S<:AbstractVector{<:Integer},
                  R<:Integer, Q<:Integer}
    Γ::T
    jj::S
    J::R
    P::Q
end
# Constructor which will see major use -- the single interface
BasicGroup(Γ::AbstractVector{<:AbstractVector{<:Integer}}) =
    BasicGroup(Γ, indexsets_tojj(Γ), length(Γ), sum(length.(Γ)))
# Additional outer constructors
BasicGroup(X::AbstractMatrix{<:Real}) = BasicGroup(onehot_toΓ(X))
BasicGroup(X::AbstractSparseMatrix{<:Real, <:Integer}) = BasicGroup(onehot_toΓ(X))
# Additional constructor -- USE WITH CAUTION! (it can be abused)
BasicGroup(jj::AbstractVector{<:Integer}) = BasicGroup(getall_cartesianindexsets(jj))

# Type-Functions -- essential operation for composing new Cartesian Groups
function cross(a::BasicGroup, b::BasicGroup)
    # cartpairs = construct_cartpairs(a.jj, b.jj)
    # cartset = unique!(sort(cartpairs))
    # return BasicGroup(getall_cartesianindexsets(cartset, cartpairs))
    # 2021-07-13: implementation below is identical, just uses symbol convention
    ξ = construct_cartpairs(a.jj, b.jj)
    return BasicGroup(getall_cartesianindexsets(ξ))
end
×(a::BasicGroup, b::BasicGroup) = cross(a, b)
# 2021-07-26: ncross, for some special use cases
function ncross(groups::Vector{<:BasicGroup})
    K = length(groups)
    group = groups[1]
    for k = 2:K
        temp = group × groups[k]
        group = temp
    end
    return group
    # alas, slower, but uses less memory.
    # ξ = construct_cartpairs2(getproperty.(groups, :jj)...)
    # return BasicGroup(getall_cartesianindexsets(unique(ξ), ξ))
end
function ncross(a::BasicGroup, others...)
    group = a
    for b in others
        temp = group × b
        group = temp
    end
    return group
end
# Alternatively, foldl can be used to express the same idea: foldl(×, grps)
# Other operations
function Base.isequal(a::BasicGroup, b::BasicGroup)
    # all(getproperty.(Ref(a), propertynames(a)) == getproperty.(Ref(b), propertynames(b)))
    # syms = (:Γ, :jj, :J, :P)
    # state = true
    # for s in syms
    #     if getproperty(a, s) != getproperty(b, s)
    #         state = false
    #         break
    #     end
    # end
    # return state
    # all(s -> getproperty(a, s) == getproperty(b, s), syms)
    # all(i -> getfield(a, i) == getfield(b, i), 1:nfields(a))
    # for i = 1:nfields(a)
    #     getfield(a, i) == getfield(b, i) || return false
    # end
    # return true
    a.Γ == b.Γ && a.jj == b.jj && a.J == b.J && a.P == b.P
end
Base.:(==)(a::BasicGroup, b::BasicGroup) = Base.isequal(a, b)

# Conversion from of BasicGroup to one-hot and enumerated list matrix
onehot_dense(a::BasicGroup) = onehot_dense(a.jj, a.P, a.J)
onehot_sparse(a::BasicGroup) = onehot_sparse(a.jj, a.P)
ELMy_ny(a::BasicGroup) = ELMy_ny(a.Γ)

####
struct AugmentedGroup{T<:AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}},
                      S<:Integer, R<:Integer, Q<:Integer, O<:Integer}
    allΓ::T
    J::S
    P::R
    M::Q
    N::O
end
# Primary constructor
AugmentedGroup(Γ::AbstractVector{<:AbstractVector{<:Integer}}, P::Int, M::Int) =
    AugmentedGroup(construct_allΓ(Γ, P, M), length(Γ), P, M, M * P)
# Convenience constructor which simplifies use
AugmentedGroup(a::BasicGroup, M::Int) =
    AugmentedGroup(construct_allΓ(a.Γ, a.P, M), a.J, a.P, M, M * a.P)
# Other operations
function Base.isequal(a::AugmentedGroup, b::AugmentedGroup)
    a.allΓ == b.allΓ && a.J == b.J && a.P == b.P && a.M == b.M && a.N == b.N
end
Base.:(==)(a::AugmentedGroup, b::AugmentedGroup) = Base.isequal(a, b)
# Convenience constructor for conversion back from AugmentedGroup
BasicGroup(a::AugmentedGroup) = BasicGroup(construct_Γ̃(a.allΓ, a.J))
#### Deficient Group: p. 41, 42 -- Xᵍʳ
struct DeficientGroup{T<:AbstractVector{<:AbstractVector{<:Integer}}, S<:Integer, R<:Integer}
    Γ::T
    J::S
    P::R
end
function onehot_dense(a::DeficientGroup)
    X = zeros(Int, a.P, a.J)
    for j = 1:J
        for p in Γ[j]
            X[p, j] = 1
        end
    end
    return X
end
################
# function getall_indexsets(regexes::Vector{Regex}, v::AbstractVector{<:AbstractString})
#     Γ = [findall(x -> occursin(r, x), v) for r in regexes]
#     return Γ
# end
# # Test example: p.23-24
# cells = ["a_v1", "a_v2", "a_v3", "b_v1", "b_v2", "c_v1"]
# Γ = getall_indexsets(regs_family, cells)
# jj = indexsets_tojj(Γ)
# J = length(Γ)
# M = 3
# P = length(cells)
# N = M * P
# allΓ = construct_allΓ(Γ, P, M)
# Γ̃ = construct_Γ̃(allΓ, J)
# jj̃ = indexsets_tojj(Γ̃)
# # Test examples: p. 27-29, 31
# cells = ["a_v1", "a_w1", "a_v2", "a_v3", "b_v1", "b_v2", "b=v2", "c_v1"]
# regs_family = [r"^a"i, r"^b"i, r"^c"i]
# regs_drive = [r"1$"i, r"2$"i, r"3$"i]
# regs_metal = [r"v"i, r"w"i]
# # family × drive for 8-cell
# Γ_family = getall_indexsets(regs_family, cells)
# jj_family = indexsets_tojj(Γ_family)
# Γ_drive = getall_indexsets(regs_drive, cells)
# jj_drive = indexsets_tojj(Γ_drive)
# cartpairs_familycrossdrive = construct_cartpairs(jj_family, jj_drive)
# cartset_familycrossdrive = unique(cartpairs_familycrossdrive)
# Γ_familycrossdrive = getall_cartesianindexsets(cartset_familycrossdrive, cartpairs_familycrossdrive)
# jj_familycrossdrive = indexsets_tojj(Γ_familycrossdrive)
# # (family × drive) × metal
# Γ_metal = getall_indexsets(regs_metal, cells)
# jj_metal = indexsets_tojj(Γ_metal)
# cartpairs_3cross = construct_cartpairs(jj_familycrossdrive, jj_metal)
# cartset_3cross = unique(cartpairs_3cross)
# Γ_3cross = getall_cartesianindexsets(cartset_3cross, cartpairs_3cross)
# jj_3cross = indexsets_tojj(Γ_3cross)
# # family × family
# cartpairs_familycrossfamily = construct_cartpairs(jj_family, jj_family)
# cartset_familycrossfamily = unique(cartpairs_familycrossfamily)
# Γ_familycrossfamily = getall_cartesianindexsets(cartset_familycrossfamily,
#                                                 cartpairs_familycrossfamily)
# jj_familycrossfamily = indexsets_tojj(Γ_familycrossfamily)
# jj_familycrossfamily == jj_family
# #### with Types: p.31
# byfamily = BasicGroup(getall_indexsets(regs_family, cells))
# bydrive = BasicGroup(getall_indexsets(regs_drive, cells))
# bymetal = BasicGroup(getall_indexsets(regs_metal, cells))
# #
# familybydrive = byfamily × bydrive
# familybydrive.Γ == Γ_familycrossdrive
# familybydrive.jj == jj_familycrossdrive
# #
# familybydrivebymetal = byfamily × bydrive × bymetal
# familybydrivebymetal.Γ == Γ_3cross
# familybydrivebymetal.jj == jj_3cross
# familybydrivebymetal.jj == BasicGroup(familybydrivebymetal.jj).jj
# #### AugmentedGroup
# aug_byfamily = AugmentedGroup(byfamily, M)
# basic_aug_byfamily = BasicGroup(aug_byfamily)
# # need to re-eval allΓ for compares to be valid
# aug_byfamily.allΓ == allΓ
# basic_aug_byfamily.jj == jj̃
# #### Test from p. 40, 2021-07-12
# xd = onehot_dense(familybydrivebymetal)
# xs = onehot_sparse(familybydrivebymetal)
# onehot_toΓ(xd) == onehot_toΓ(xs)
# familybydrivebymetal_rev = BasicGroup(xd)
# familybydrivebymetal_rev == familybydrivebymetal # as expected since isequal is undefined
# familybydrivebymetal == familybydrive
# ELMy_ny(familybydrivebymetal)
# # BasicGroup constructor testing
# x1 = BasicGroup(jj_3cross)
# x1 == familybydrivebymetal
# x2 = BasicGroup(xd)
# x2 == familybydrivebymetal
# x3 = BasicGroup(xs)
# x3 == familybydrivebymetal
# x4 = BasicGroup(Γ_3cross)
# x4 == familybydrivebymetal
