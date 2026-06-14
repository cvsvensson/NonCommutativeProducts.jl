_scalartype(x) = VectorInterface.scalartype(typeof(x))

function _set_ncadd!(y::NCAdd, x::NCAdd)
    y.coeff = x.coeff
    empty!(y.dict)
    merge!(y.dict, x.dict)
    return y
end

VectorInterface.scalartype(::Type{<:NCMul{C}}) where {C<:Number} = C
VectorInterface.scalartype(::Type{<:NCAdd{C}}) where {C<:Number} = C

function VectorInterface.zerovector(x::NCMul, ::Type{S}) where {S<:Number}
    K = typeof(NCMul(1, x.factors))
    return NCAdd(zero(S), Dict{K,S}())
end
function VectorInterface.zerovector(x::NCAdd, ::Type{S}) where {S<:Number}
    K = keytype(typeof(x.dict))
    return NCAdd(zero(S), Dict{K,S}())
end

function VectorInterface.zerovector!(x::NCAdd)
    x.coeff = zero(_scalartype(x))
    empty!(x.dict)
    return x
end

function VectorInterface.zerovector!!(x::NCMul)
    return VectorInterface.zerovector(x, _scalartype(x))
end
function VectorInterface.zerovector!!(x::NCAdd)
    return VectorInterface.zerovector!(x)
end

VectorInterface.scale(x::MulAdd, α::Number) = α * x

function VectorInterface.scale!(x::NCMul, α::Number)
    throw(ArgumentError("NCMul is immutable; use scale or scale!!"))
end
function VectorInterface.scale!(x::NCAdd, α::Number)
    x.coeff *= α
    for (k, v) in x.dict
        x.dict[k] = v * α
    end
    return x
end
function VectorInterface.scale!(y::NCMul, x::NCMul, α::Number)
    throw(ArgumentError("NCMul is immutable; use scale or scale!!"))
end
function VectorInterface.scale!(y::NCAdd, x::MulAdd, α::Number)
    return _set_ncadd!(y, VectorInterface.scale(x, α))
end

function VectorInterface.scale!!(x::NCMul, α::Number)
    return VectorInterface.scale(x, α)
end
VectorInterface.scale!!(x::NCAdd, α::Number) = scale!!(x, α)
function VectorInterface.scale!!(y::NCMul, x::MulAdd, α::Number)
    return VectorInterface.scale(x, α)
end
function VectorInterface.scale!!(y::NCAdd, x::MulAdd, α::Number)
    try
        return VectorInterface.scale!(y, x, α)
    catch
        return VectorInterface.scale(x, α)
    end
end

function VectorInterface.add(y::MulAdd, x::MulAdd, α::Number, β::Number)
    return β * y + α * x
end

function VectorInterface.add!(y::NCMul, x::NCMul, α::Number, β::Number)
    throw(ArgumentError("NCMul is immutable; use add or add!!"))
end
function VectorInterface.add!(y::NCAdd, x::NCAdd, α::Number, β::Number)
    return _set_ncadd!(y, VectorInterface.add(y, x, α, β))
end

function VectorInterface.add!!(y::MulAdd, x::MulAdd, α::Number, β::Number)
    return add!!(y, x, α, β)
end

VectorInterface.inner(x::MulAdd, y::MulAdd) = _inner(x, y)
_inner(x, y) = scalar(x' * y)

LinearAlgebra.norm(x::MulAdd) = sqrt(VectorInterface.inner(x, x))
