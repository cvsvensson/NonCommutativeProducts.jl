module NonCommutativeProducts

using Base.ScopedValues
using LinearAlgebra: LinearAlgebra, UniformScaling
using TestItems: @testitem
using BangBang: push!!, pushfirst!!, setindex!!, append!!
using BangBang.Extras: modify!!
using VectorInterface: VectorInterface, One

include("mul.jl")
include("add.jl")
include("muladd.jl")
include("sorting.jl")
include("traversal.jl")
include("vectorinterface.jl")

end
