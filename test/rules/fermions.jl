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
