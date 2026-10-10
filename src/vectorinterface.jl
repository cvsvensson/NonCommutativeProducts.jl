_scalartype(x) = VectorInterface.scalartype(typeof(x))
# The in-place methods below compute the result out of place and then store it with _set_ncadd! (see add.jl), which
# converts everything before it mutates `y`, so that `y` is left unchanged if the result doesn't fit in it.

VectorInterface.scalartype(::Type{<:NCMul{C}}) where {C<:Number} = C
VectorInterface.scalartype(::Type{<:NCAdd{C}}) where {C<:Number} = C

VectorInterface.zerovector(x::MulAdd, ::Type{S}) where {S<:Number} = zero(NCAdd{S,to_add_dict_type(typeof(x))})

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

# A product isn't closed under addition, so as a vector it is represented by a sum. Solvers like KrylovKit store
# vectors in containers typed by scale's result, so scale must return the same type as zerovector and add.
VectorInterface.scale(x::NCAdd, α::Number) = α * x
VectorInterface.scale(x::NCMul, α::Number) = α * NCAdd(x)

function VectorInterface.scale!(x::NCMul, α::Number)
    throw(ArgumentError("NCMul is immutable; use scale or scale!!"))
end
VectorInterface.scale!(x::NCAdd, α::Number) = scale!(x, α)
function VectorInterface.scale!(y::NCMul, x::MulAdd, α::Number)
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
    scale!!(y, x, α)
end

function VectorInterface.add(y::MulAdd, x::MulAdd, α::Number, β::Number)
    return β * y + α * x
end

function VectorInterface.add!(y::NCMul, x::MulAdd, α::Number, β::Number)
    throw(ArgumentError("NCMul is immutable; use add or add!!"))
end
function VectorInterface.add!(y::NCAdd, x::MulAdd, α::Number, β::Number)
    return _set_ncadd!(y, VectorInterface.add(y, x, α, β))
end

VectorInterface.add!!(y::MulAdd, x, α::Number, β::Number) = add!!(y, x, α, β)
VectorInterface.add!!(y, x::MulAdd, α::Number, β::Number) = add!!(y, x, α, β)
VectorInterface.add!!(y::MulAdd, x::MulAdd, α::Number, β::Number) = add!!(y, x, α, β)

VectorInterface.inner(x::MulAdd, y::MulAdd) = _inner(x, y)
function _inner(x, y)
    autosort() && return _inner_sorted(_as_muladd(x), _as_muladd(y))
    z = x' * y
    _reduces_to_scalar(z) || throw(_inner_error())
    return scalar(z)
end
_inner_error() = ArgumentError("inner(x, y) requires x' * y to reduce to a scalar; check that adjoint and mul_effect are defined and autosort is enabled")
_as_muladd(x::MulAdd) = x
_as_muladd(x) = NCMul(x)

# With autosort, x' * y is the sum of the sorted products of each term of x' with each term of y. Sorting the products
# as one list, instead of collecting them in a sum first, avoids hashing every unsorted product, and the products that
# reduce to scalars are added up directly. Only the rest is collected in a sum, to check that it cancels.
function _inner_sorted(x::MulAdd, y::MulAdd)
    init = zero(_scalartype(x))' * zero(_scalartype(y))
    xterms, yterms = _inner_terms(x), _inner_terms(y)
    (isempty(xterms) || isempty(yterms)) && return init
    sorted = with(_autosort => false) do
        xterms_adjoint = map(adjoint, xterms)
        products = vec([catenate(a, b) for a in xterms_adjoint, b in yterms])
        # inference may leave the element type abstract, e.g. when x and y have different types
        __bubble_sort!(eltype(products) <: NCMul ? products : identity.(products))
    end
    return _sum_inner_terms(sorted, init)
end
# The terms of x as products, with the additive coefficient as a product without factors
_inner_terms(x::NCMul) = [x]
function _inner_terms(x::NCAdd{C,K}) where {C,K}
    terms = [NCMul(v, k.factors) for (k, v) in x.dict]
    iszero(additive_coeff(x)) || push!(terms, NCMul(additive_coeff(x), one(K).factors))
    return terms
end
function _sum_inner_terms(terms::Vector{T}, init) where {T<:NCMul}
    s = init
    rest = T[]
    for term in terms
        isscalar(term) ? (s += prefactor(term)) : push!(rest, term)
    end
    isempty(rest) && return s
    remainder = _sum_sorted_terms(rest, typeof(s))
    isscalar(remainder) || throw(_inner_error())
    return s + scalar(remainder)
end
_reduces_to_scalar(::Number) = true
_reduces_to_scalar(z::MulAdd) = isscalar(z)
_reduces_to_scalar(z) = false

LinearAlgebra.norm(x::MulAdd) = sqrt(real(VectorInterface.inner(x, x)))
# p-norms depend on a choice of basis, which a sum of products doesn't have
LinearAlgebra.norm(x::MulAdd, p::Real) = p == 2 ? LinearAlgebra.norm(x) : throw(ArgumentError("norm(x, p) is only supported for p = 2, the norm induced by inner"))
LinearAlgebra.dot(x::MulAdd, y::MulAdd) = VectorInterface.inner(x, y)

# Whether scalartype(T) is defined by a method other than VectorInterface's generic (warning) fallback
function _has_scalartype(T::Type)
    m = which(VectorInterface.scalartype, Tuple{Type{T}})
    return m.sig !== Tuple{typeof(VectorInterface.scalartype),Type}
end

# The VectorInterface methods for an atom type T, which act on it as the product with it as the only factor.
# Called by @nc_common.
macro nc_vectorinterface(T)
    T = esc(T)
    quote
        VectorInterface.inner(x::MulAdd, y::$T) = _inner(x, y)
        VectorInterface.inner(x::$T, y::MulAdd) = _inner(x, y)
        VectorInterface.inner(x::$T, y::$T) = _inner(x, y)
        LinearAlgebra.norm(x::$T) = sqrt(real(VectorInterface.inner(x, x)))
        LinearAlgebra.dot(x::MulAdd, y::$T) = VectorInterface.inner(x, y)
        LinearAlgebra.dot(x::$T, y::MulAdd) = VectorInterface.inner(x, y)
        LinearAlgebra.dot(x::$T, y::$T) = VectorInterface.inner(x, y)

        # an atom acts as 1 * atom, but an existing scalartype method for T (e.g. if T is an AbstractArray) is kept.
        # The methods whose VectorInterface defaults call scalartype are defined directly, so they work either way.
        if !NonCommutativeProducts._has_scalartype($T)
            VectorInterface.scalartype(::Type{<:$T}) = Int
        end
        VectorInterface.scale(x::$T, α::Number) = VectorInterface.scale(NCMul(x), α)
        VectorInterface.scale!!(x::$T, α::Number) = VectorInterface.scale(x, α)
        VectorInterface.scale!!(a::NCAdd, x::$T, α::Number) = add!!(a, NCMul(x), α, VectorInterface.Zero())
        VectorInterface.scale!!(y::$T, x::Union{$T,MulAdd}, α::Number) = VectorInterface.scale(x, α)
        VectorInterface.scale!!(y::NCMul, x::$T, α::Number) = VectorInterface.scale(x, α)
        # an atom is immutable, so like an NCMul it can't be the destination of the in-place methods
        VectorInterface.scale!(x::$T, α::Number) = throw(ArgumentError(string(typeof(x), " is immutable; use scale or scale!!")))
        VectorInterface.scale!(y::$T, x::Union{$T,MulAdd}, α::Number) = throw(ArgumentError(string(typeof(y), " is immutable; use scale or scale!!")))
        VectorInterface.scale!(y::MulAdd, x::$T, α::Number) = VectorInterface.scale!(y, NCMul(x), α)
        VectorInterface.add(y::$T, x::Union{$T,MulAdd}, α::Number, β::Number) = VectorInterface.add(NCMul(y), x, α, β)
        VectorInterface.add(y::MulAdd, x::$T, α::Number, β::Number) = VectorInterface.add(y, NCMul(x), α, β)
        VectorInterface.add!!(y::$T, x::$T, α::Number, β::Number) = add!!(y, x, α, β)
        VectorInterface.add!(y::$T, x::Union{$T,MulAdd}, α::Number, β::Number) = throw(ArgumentError(string(typeof(y), " is immutable; use add or add!!")))
        VectorInterface.add!(y::NCMul, x::$T, α::Number, β::Number) = VectorInterface.add!(y, NCMul(x), α, β)
        VectorInterface.add!(y::NCAdd, x::$T, α::Number, β::Number) = VectorInterface.add!(y, NCAdd(x), α, β)
        VectorInterface.zerovector(x::$T) = VectorInterface.zerovector(NCMul(x))
        VectorInterface.zerovector(x::$T, ::Type{S}) where {S<:Number} = VectorInterface.zerovector(NCMul(x), S)
        VectorInterface.zerovector!!(x::$T) = VectorInterface.zerovector(x)
        VectorInterface.zerovector!!(x::$T, ::Type{S}) where {S<:Number} = VectorInterface.zerovector(x, S)
    end
end
