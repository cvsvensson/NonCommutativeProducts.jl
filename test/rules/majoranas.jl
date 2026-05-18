using NonCommutativeProducts
import NonCommutativeProducts: mul_effect, @nc
struct Majorana{L}
    label::L
end
Base.adjoint(x::Majorana) = Majorana(x.label)
Base.show(io::IO, x::Majorana) = print(io, "γ[", x.label, "]")

@nc Majorana
function mul_effect(a::Majorana, b::Majorana)
    if a.label == b.label
        return 1 # a*b => 1
    elseif a.label < b.label
        return nothing # a*b => a*b
    elseif a.label > b.label
        return -b * a # a*b => -b*a
    else
        throw(ArgumentError("Don't know how to multiply $a * $b"))
    end
end