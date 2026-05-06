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
