module NonCommutativeProducts

using Base.ScopedValues
using LinearAlgebra
using TestItems
using BangBang: push!!, pushfirst!!, setindex!!, append!!, mergewith!!
using BangBang.Extras: modify!!
import VectorInterface
using VectorInterface: One

include("mul.jl")
include("add.jl")
include("muladd.jl")
include("sorting.jl")
include("traversal.jl")
include("vectorinterface.jl")

end
