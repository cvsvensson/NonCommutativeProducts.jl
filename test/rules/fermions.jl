using NonCommutativeProducts
import NonCommutativeProducts: @nc
struct Fermion{L}
    label::L
    creation::Bool
end
Base.adjoint(x::Fermion) = Fermion(x.label, !x.creation)
Fermion(k) = Fermion(k, false)
Base.show(io::IO, x::Fermion) = print(io, "c", x.creation ? "†" : "", "[", x.label, "]")
@nc Fermion

function should_swap(a::Fermion, b::Fermion)
    if a.creation == b.creation
        return a.label > b.label
    else
        a.creation < b.creation
    end
end

function NonCommutativeProducts.mul_effect(a::Fermion, b::Fermion)
    a == b && return 0 # a*b => 0
    if should_swap(a, b)
        a.label == b.label && xor(a.creation, b.creation) && return -b * a + 1
        return -b * a # a*b => -b*a
    else
        return nothing
    end
end

struct FState{S}
    occ::Bool
    label::S
    adj::Bool
end
State(occ, f::Fermion) = FState(occ, f.label, false)
Base.adjoint(s::FState) = FState(s.occ, s.label, !s.adj)
function NonCommutativeProducts.mul_effect(f::Fermion, s::FState)
    f.label == s.label || throw(ArgumentError("Cannot multiply Fermion by FState with different space"))
    s.adj && throw(ArgumentError("Cannot multiply Fermion by adjoint FState"))
    s.occ == f.creation && return 0
    return FState(!s.occ, s.label, s.adj)
end
function NonCommutativeProducts.mul_effect(s::FState, f::Fermion)
    f.label == s.label || throw(ArgumentError("Cannot multiply FState by Fermion with different space"))
    !s.adj && throw(ArgumentError("Cannot multiply Ket by fermion from the right"))
    s.occ == !f.creation && return 0
    return FState(!s.occ, s.label, s.adj)
end
function NonCommutativeProducts.mul_effect(s1::FState, s2::FState)
    if s1.label == s2.label
        s1.adj && !s2.adj && return (s1.occ == s2.occ)
        throw(ArgumentError("Cannot multiply two FStates with the same space and same adj"))
    end
    s1.adj < s2.adj && return NonCommutativeProducts.Swap(1)
    s2.adj < s1.adj && return nothing
    s1.label < s2.label && return nothing
    return NonCommutativeProducts.Swap(1)
end
NonCommutativeProducts.@nc FState Fermion

