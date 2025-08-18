module NonCommutativeProducts

using LinearAlgebra
using TestItems
import BangBang: push!!, pushfirst!!, setindex!!, append!!, mergewith!!

include("mul.jl")
include("add.jl")
include("muladd.jl")
include("sorting.jl")

end
