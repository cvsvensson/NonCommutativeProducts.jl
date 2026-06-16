using NonCommutativeProducts
import NonCommutativeProducts: @nc
struct Boson
    exp::Int
end
Base.adjoint(x::Boson) = Boson(-x.exp)
Boson() = Boson(-1)
Base.show(io::IO, x::Boson) = print(io, "b", x.exp > 0 ? "†" : "", abs(x.exp) > 1 ? "^($(x.exp))" : "")
@nc Boson

function NonCommutativeProducts.mul_effect(a::Boson, b::Boson)
    sign(a.exp) == sign(b.exp) && return Boson(a.exp + b.exp)
    if a.exp < 0 && b.exp > 0
        return b * a + 1
    else
        return nothing
    end
end

struct BState{S}
    occ::Int
    adj::Bool
end
State(occ) = BState(occ, false)
Base.adjoint(s::BState) = BState(s.occ, s.label, !s.adj)
function NonCommutativeProducts.mul_effect(b::Boson, s::BState)
    s.adj && throw(ArgumentError("Cannot multiply Boson by adjoint BState"))
    newocc = s.occ + b.exp
    newocc < 0 && return 0
    return BState(newocc, s.space, s.adj)
end
function NonCommutativeProducts.mul_effect(s::BState, b::Boson)
    !s.adj && throw(ArgumentError("Cannot multiply Ket by boson from the right"))
    newocc = s.occ + !b.exp
    newocc < 0 && return 0
    return BState(newocc, s.space, s.adj)
end
function NonCommutativeProducts.mul_effect(s1::BState, s2::BState)
    s1.adj && !s2.adj && return (s1.occ == s2.occ)
    throw(ArgumentError("Cannot multiply two BStates with the same space and same adj"))
end
NonCommutativeProducts.@nc BState Boson

