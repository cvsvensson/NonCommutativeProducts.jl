module NonCommutativeProducts

using LinearAlgebra
using TermInterface
using TestItems
import BangBang: push!!, pushfirst!!, setindex!!, append!!
export @fermions, @majoranas

include("mul.jl")
include("add.jl")
include("muladd.jl")
include("symbolic_fermions.jl")
include("symbolic_majoranas.jl")
include("sorting.jl")

end
