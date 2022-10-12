using CartesianGrouping
using Test

@testset "CartesianGrouping.jl" begin
    # Test examples: p. 27-29, 31
    cells = ["a_v1", "a_w1", "a_v2", "a_v3", "b_v1", "b_v2", "b=v2", "c_v1"]
    regs_family = [r"^a"i, r"^b"i, r"^c"i]
    regs_drive = [r"1$"i, r"2$"i, r"3$"i]
    regs_metal = [r"v"i, r"w"i]
    @testset "construct" begin
        # family × drive for 8-cell
        Γ_family = getall_cartesianindexsets(regs_family, cells)
        jj_family = indexsets_tojj(Γ_family)
        @test jj_family == [1, 1, 1, 1, 2, 2, 2, 3]
        Γ_drive = getall_cartesianindexsets(regs_drive, cells)
        jj_drive = indexsets_tojj(Γ_drive)
        @test jj_drive == [1, 1, 2, 3, 1, 2, 2, 1]
        cartpairs_familycrossdrive = construct_cartpairs(jj_family, jj_drive)
        cartset_familycrossdrive = unique(cartpairs_familycrossdrive)
        Γ_familycrossdrive = getall_cartesianindexsets(cartset_familycrossdrive,
                                                       cartpairs_familycrossdrive)
        jj_familycrossdrive = indexsets_tojj(Γ_familycrossdrive)
        @test jj_familycrossdrive == [1, 1, 2, 3, 4, 5, 5, 6]
        # (family × drive) × metal
        Γ_metal = getall_cartesianindexsets(regs_metal, cells)
        jj_metal = indexsets_tojj(Γ_metal)
        @test jj_metal == [1, 2, 1, 1, 1, 1, 1, 1]
        cartpairs_3cross = construct_cartpairs(jj_familycrossdrive, jj_metal)
        cartset_3cross = unique(cartpairs_3cross)
        Γ_3cross = getall_cartesianindexsets(cartset_3cross, cartpairs_3cross)
        jj_3cross = indexsets_tojj(Γ_3cross)
        # family × family
        cartpairs_familycrossfamily = construct_cartpairs(jj_family, jj_family)
        cartset_familycrossfamily = unique(cartpairs_familycrossfamily)
        Γ_familycrossfamily = getall_cartesianindexsets(cartset_familycrossfamily,
                                                        cartpairs_familycrossfamily)
        jj_familycrossfamily = indexsets_tojj(Γ_familycrossfamily)
        @test jj_familycrossfamily == jj_family
    end
    @testset "constructspecial" begin
        # Testing
        cells₁ = ["a_v1", "a_w1", "a_v2", "a_v3", "b_v1", "b_v2", "b=v2", "c_v1"]
        Γ₁ = getall_cartesianindexsets(regs_family, cells₁)
        jj₁ = indexsets_tojj(Γ₁)
        cells₂ = ["b_v1", "b_v2", "b_v3", "c_v1", "c=v1", "c_v2"]
        Γ₂ = getall_cartesianindexsets(regs_family, cells₂)
        jj₂ = indexsets_tojj(Γ₂)
        cells₃ = ["a=v1", "a=w1", "b=v1", "b=v2", "b=v3"]
        Γ₃ = getall_cartesianindexsets(regs_family, cells₃)
        jj₃ = indexsets_tojj(Γ₃)
        catΓ = [Γ₁, Γ₂, Γ₃]
        allΓ = constructspecial_allΓ(catΓ, sum.(map(x -> length.(x), catΓ)))
        Ĵ = 3  # number of groups which appear in Ŝ -- known from length(regs_family)
        jj = indexsets_tojj(construct_Γ̃(allΓ, Ĵ))
        catjj = [jj₁; jj₂; jj₃]
        @test jj == catjj
    end
    @testset "cartesian_indexsets" begin
        ξ = [13:24;];
        ξnew = reshape(ξ, 3, 4);
        f1(x) = 13 ≤ x ≤ 15
        f2(x) = 16 ≤ x ≤ 18
        f3(x) = 19 ≤ x ≤ 21
        f4(x) = 22 ≤ x ≤ 24
        fvec = [f1, f2, f3, f4]
        # linear indexed
        Γ = getall_cartesianindexsets(fvec, ξ)
        @test necessarycond1(ξ, Γ)
        @test necessarycond2(Γ)
        # cartesian indexed
        Γnew = getall_cartesianindexsets(fvec, ξnew)
        @test necessarycond1(ξnew, Γnew)
        @test necessarycond2(Γnew)
    end
    @testset "Types" begin
        # family × drive for 8-cell
        Γ_family = getall_cartesianindexsets(regs_family, cells)
        jj_family = indexsets_tojj(Γ_family)
        Γ_drive = getall_cartesianindexsets(regs_drive, cells)
        jj_drive = indexsets_tojj(Γ_drive)
        cartpairs_familycrossdrive = construct_cartpairs(jj_family, jj_drive)
        cartset_familycrossdrive = unique(cartpairs_familycrossdrive)
        Γ_familycrossdrive = getall_cartesianindexsets(cartset_familycrossdrive,
                                                       cartpairs_familycrossdrive)
        jj_familycrossdrive = indexsets_tojj(Γ_familycrossdrive)
        # (family × drive) × metal
        Γ_metal = getall_cartesianindexsets(regs_metal, cells)
        jj_metal = indexsets_tojj(Γ_metal)
        cartpairs_3cross = construct_cartpairs(jj_familycrossdrive, jj_metal)
        cartset_3cross = unique(cartpairs_3cross)
        Γ_3cross = getall_cartesianindexsets(cartset_3cross, cartpairs_3cross)
        jj_3cross = indexsets_tojj(Γ_3cross)
        #### with Types: p.31
        byfamily = BasicGroup(getall_cartesianindexsets(regs_family, cells))
        bydrive = BasicGroup(getall_cartesianindexsets(regs_drive, cells))
        bymetal = BasicGroup(getall_cartesianindexsets(regs_metal, cells))
        #
        familybydrive = byfamily × bydrive
        @test familybydrive.Γ == Γ_familycrossdrive
        @test familybydrive.jj == jj_familycrossdrive
        #
        familybydrivebymetal = byfamily × bydrive × bymetal
        @test familybydrivebymetal.Γ == Γ_3cross
        @test familybydrivebymetal.jj == jj_3cross
        @test familybydrivebymetal.jj == BasicGroup(familybydrivebymetal.jj).jj
        #### AugmentedGroup
        Γ = getall_cartesianindexsets(regs_family, cells)
        jj = indexsets_tojj(Γ)
        J = length(Γ)
        M = 3
        P = length(cells)
        N = M * P
        allΓ = construct_allΓ(Γ, P, M)
        Γ̃ = construct_Γ̃(allΓ, J)
        jj̃ = indexsets_tojj(Γ̃)
        # construct
        aug_byfamily = AugmentedGroup(byfamily, M)
        basic_aug_byfamily = BasicGroup(aug_byfamily)
        @test aug_byfamily.allΓ == allΓ
        @test basic_aug_byfamily.jj == jj̃
        #### Test from p. 40, 2021-07-12
        xd = onehot_dense(familybydrivebymetal)
        xs = onehot_sparse(familybydrivebymetal)
        @test onehot_toΓ(xd) == onehot_toΓ(xs)
        familybydrivebymetal_rev = BasicGroup(xd)
        @test familybydrivebymetal_rev == familybydrivebymetal
        @test familybydrivebymetal_rev !== familybydrivebymetal
        @test familybydrivebymetal != familybydrive
        # ELMy_ny(familybydrivebymetal)
        # BasicGroup constructor testing
        @test BasicGroup(jj_3cross) == familybydrivebymetal
        @test BasicGroup(xd) == familybydrivebymetal
        @test BasicGroup(xs) == familybydrivebymetal
        @test BasicGroup(Γ_3cross) == familybydrivebymetal
    end

end
