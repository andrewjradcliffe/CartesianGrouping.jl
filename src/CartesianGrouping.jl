module CartesianGrouping

using SparseArrays
import LinearAlgebra: ×

# cartesiangrouping_base.jl
export construct_Gⱼₘ, construct_Γₘ, construct_allΓ, construct_Γⱼ, construct_G̃ⱼ, construct_Γ̃,
    constructspecial_Gⱼₘ, constructspecial_Γₘ, constructspecial_allΓ

export construct_cartpairs, getall_cartesianindexsets, necessarycond1, necessarycond2

export indexsets_tojj, onehot_sparse, onehot_dense, onehot_toΓ, ELMy_ny

export BasicGroup, cross, ×, ncross, AugmentedGroup, DeficientGroup

include("cartesiangrouping_base.jl")

# dataaugmentation_hierarchicallinear.jl
export to_jj_beta, to_W_beta, to_x_beta

include("dataaugmentation_hierarchicallinear.jl")

end
