module NonCommutativeProducts

using Base.ScopedValues
using LinearAlgebra: LinearAlgebra, UniformScaling
using TestItems: @testitem
using BangBang: push!!, pushfirst!!, setindex!!, append!!
using BangBang.Extras: modify!!
using VectorInterface: VectorInterface, One
using PrecompileTools: @setup_workload, @compile_workload

include("mul.jl")
include("add.jl")
include("muladd.jl")
include("promotion.jl")
include("sorting.jl")
include("traversal.jl")
include("vectorinterface.jl")
include("precompile.jl")

end
