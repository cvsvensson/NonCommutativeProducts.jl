module NonCommutativeProducts

using Base.ScopedValues
using LinearAlgebra: LinearAlgebra, UniformScaling
using BangBang: push!!, pushfirst!!, setindex!!, append!!
using BangBang.Extras: modify!!
using VectorInterface: VectorInterface, One
using PrecompileTools: @setup_workload, @compile_workload

public NCMul, NCAdd, @nc, @commutative, mul_effect, enable_autosort!, disable_autosort!, autosort, with_autosort,
    ncmap, ncmapreduce, scalar, add!!, terms, factors, coefficient

include("mul.jl")
include("add.jl")
include("muladd.jl")
include("accessors.jl")
include("promotion.jl")
include("sorting.jl")
include("grading.jl")
include("traversal.jl")
include("vectorinterface.jl")
include("precompile.jl")

end
