"""
    NonCommutativeProducts

Symbolic arithmetic with atoms that do not commute. Atoms are registered with `@nc`. Products are `NCMul` values and
sums are `NCAdd` values. The order of factors is set by defining `mul_effect` for the atom types. Products are sorted
as they are formed, unless inside `with_autosort(false)`.
"""
module NonCommutativeProducts

using Base.ScopedValues
using LinearAlgebra: LinearAlgebra, UniformScaling
using VectorInterface: VectorInterface, One
using PrecompileTools: @setup_workload, @compile_workload

public NCMul, NCAdd, @nc, @commutative, mul_effect, Swap, AddTerms, autosort, with_autosort,
    ncmap, ncmapreduce, scalar, add!!, terms, nterms, factors, coefficient

include("mul.jl")
include("add.jl")
include("muladd.jl")
include("accessors.jl")
include("promotion.jl")
include("sorting.jl")
include("grading.jl")
include("traversal.jl")
include("segments.jl")
include("vectorinterface.jl")
include("truncation.jl")
include("precompile.jl")

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        exc.f === mul_effect || return
        print(io, "\n\nProducts are sorted as they are formed and by `inner`, so define `mul_effect` for these types, or, if you don't need `inner`, form the product inside `with_autosort(false)`.")
    end
end

end
