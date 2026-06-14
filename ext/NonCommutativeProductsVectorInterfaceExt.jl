module NonCommutativeProductsVectorInterfaceExt

using NonCommutativeProducts
import VectorInterface

const NC = NonCommutativeProducts

_scalartype(x) = VectorInterface.scalartype(typeof(x))

function _set_ncadd!(y::NC.NCAdd, x::NC.NCAdd)
    y.coeff = x.coeff
    empty!(y.dict)
    merge!(y.dict, x.dict)
    return y
end

# Required by VectorInterface for custom vector-like types:
# - scalartype(::Type{T})
# - zerovector(x, S)
# - zerovector!!(x)
# - scale(x, α), scale!(...), scale!!(...)
# - add(y, x, α, β), add!(...), add!!(...)
# - inner(x, y)

VectorInterface.scalartype(::Type{<:NC.NCMul{C}}) where {C<:Number} = C
VectorInterface.scalartype(::Type{<:NC.NCAdd{C}}) where {C<:Number} = C

function VectorInterface.zerovector(x::NC.NCMul, ::Type{S}) where {S<:Number}
    K = typeof(NC.NCMul(1, x.factors))
    return NC.NCAdd(zero(S), Dict{K,S}())
end
function VectorInterface.zerovector(x::NC.NCAdd, ::Type{S}) where {S<:Number}
    K = keytype(typeof(x.dict))
    return NC.NCAdd(zero(S), Dict{K,S}())
end

function VectorInterface.zerovector!(x::NC.NCAdd)
    x.coeff = zero(_scalartype(x))
    empty!(x.dict)
    return x
end

function VectorInterface.zerovector!!(x::NC.NCMul)
    return VectorInterface.zerovector(x, _scalartype(x))
end
function VectorInterface.zerovector!!(x::NC.NCAdd)
    return VectorInterface.zerovector!(x)
end

function VectorInterface.scale(x::NC.NCMul, α::Number)
    return α * x
end
function VectorInterface.scale(x::NC.NCAdd, α::Number)
    return α * x
end

function VectorInterface.scale!(x::NC.NCMul, α::Number)
    throw(ArgumentError("NCMul is immutable; use scale or scale!!"))
end
function VectorInterface.scale!(x::NC.NCAdd, α::Number)
    x.coeff *= α
    for (k, v) in x.dict
        x.dict[k] = v * α
    end
    return x
end
function VectorInterface.scale!(y::NC.NCMul, x::NC.NCMul, α::Number)
    throw(ArgumentError("NCMul is immutable; use scale or scale!!"))
end
function VectorInterface.scale!(y::NC.NCAdd, x::NC.NCAdd, α::Number)
    return _set_ncadd!(y, VectorInterface.scale(x, α))
end

function VectorInterface.scale!!(x::NC.NCMul, α::Number)
    return VectorInterface.scale(x, α)
end
function VectorInterface.scale!!(x::NC.NCAdd, α::Number)
    try
        return VectorInterface.scale!(x, α)
    catch
        return VectorInterface.scale(x, α)
    end
end
function VectorInterface.scale!!(y::NC.NCMul, x::NC.NCMul, α::Number)
    return VectorInterface.scale(x, α)
end
function VectorInterface.scale!!(y::NC.NCAdd, x::NC.NCAdd, α::Number)
    try
        return VectorInterface.scale!(y, x, α)
    catch
        return VectorInterface.scale(x, α)
    end
end

function VectorInterface.add(y::NC.NCMul, x::NC.NCMul, α::Number, β::Number)
    return β * y + α * x
end
function VectorInterface.add(y::NC.NCAdd, x::NC.NCAdd, α::Number, β::Number)
    return β * y + α * x
end

function VectorInterface.add!(y::NC.NCMul, x::NC.NCMul, α::Number, β::Number)
    throw(ArgumentError("NCMul is immutable; use add or add!!"))
end
function VectorInterface.add!(y::NC.NCAdd, x::NC.NCAdd, α::Number, β::Number)
    return _set_ncadd!(y, VectorInterface.add(y, x, α, β))
end

function VectorInterface.add!!(y::NC.NCMul, x::NC.NCMul, α::Number, β::Number)
    return VectorInterface.add(y, x, α, β)
end
function VectorInterface.add!!(y::NC.NCAdd, x::NC.NCAdd, α::Number, β::Number)
    try
        return VectorInterface.add!(y, x, α, β)
    catch
        return VectorInterface.add(y, x, α, β)
    end
end

# function VectorInterface.inner(x::NC.NCMul, y::NC.NCMul)
#     return VectorInterface.inner(x + 0, y + 0)
# end
# function VectorInterface.inner(x::NC.NCAdd, y::NC.NCAdd)
#     return _inner_terms(x, y)
# end
# function _inner_terms(x::NC.NCAdd, y::NC.NCAdd)
#     s = conj(x.coeff) * y.coeff
#     for (k, vx) in x.dict
#         vy = get(y.dict, k, zero(vx))
#         s += conj(vx) * vy
#     end
#     return s
# end


end
