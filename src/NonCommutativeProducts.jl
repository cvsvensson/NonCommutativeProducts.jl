"""
    NonCommutativeProducts

Symbolic arithmetic with atoms that do not commute. Atoms are registered with `@nc`. Products are `NCMul` values and
sums are `NCAdd` values. The order of factors is set by defining `mul_effect` for the atom types, and autosorting
during multiplication is controlled by `enable_autosort!`, `disable_autosort!` and `with_autosort`.
"""
module NonCommutativeProducts

using Base.ScopedValues
using LinearAlgebra: LinearAlgebra, UniformScaling
using BangBang: push!!, pushfirst!!, setindex!!, append!!
using BangBang.Extras: modify!!
using VectorInterface: VectorInterface, One
using PrecompileTools: @setup_workload, @compile_workload

public NCMul, NCAdd, @nc, @commutative, mul_effect, enable_autosort!, disable_autosort!, autosort, with_autosort,
    ncmap, ncmapreduce, scalar, add!!, terms, nterms, termpairs, factors, coefficient,
    map_coefficients

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
