module NonCommutativeProducts

using LinearAlgebra
using TermInterface
using TestItems
# using BangBang
export @fermions, @majoranas

include("muladd.jl")
include("symbolic_fermions.jl")
include("symbolic_majoranas.jl")
include("sorting.jl")

end
