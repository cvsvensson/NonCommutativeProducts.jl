using NonCommutativeProducts
import NonCommutativeProducts: @nc

"""
    SpinOp(label, op)

General spin operator using raising/lowering basis.
- `label`: identifies the site/particle
- `op`: :z for Sz, :+ for S+, :- for S-
"""
struct SpinOp{L}
    label::L
    op::Symbol
end

function Base.adjoint(x::SpinOp)
    new_op = x.op == :+ ? :- : x.op == :- ? :+ : :z
    return SpinOp(x.label, new_op)
end

SpinOp(label) = SpinOp(label, :z)

function Base.show(io::IO, x::SpinOp)
    print(io, "S", x.op, "[$(x.label)]")
end

@nc SpinOp

const OP_ORDER = Dict(:z => 1, :+ => 2, :- => 3)

function should_swap(a::SpinOp, b::SpinOp)
    a.label != b.label ? a.label > b.label : OP_ORDER[a.op] > OP_ORDER[b.op]
end

function NonCommutativeProducts.mul_effect(a::SpinOp, b::SpinOp)
    if a.label != b.label
        # Different sites commute
        return should_swap(a, b) ? b * a : nothing
    end

    # Commutation relations: [Sz, S+] = S+, [Sz, S-] = -S-, [S+, S-] = 2Sz
    # Canonical order: :z < :+ < :-
    if a.op == :+ && b.op == :z
        # S+ Sz = Sz S+ - S+
        return b * a - SpinOp(a.label, :+)
    elseif a.op == :- && b.op == :z
        # S- Sz = Sz S- + S-
        return b * a + SpinOp(a.label, :-)
    elseif a.op == :- && b.op == :+
        # S- S+ = S+ S- - 2Sz
        return b * a - 2 * SpinOp(a.label, :z)
    end

    return nothing
end