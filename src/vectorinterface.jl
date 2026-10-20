_scalartype(x) = VectorInterface.scalartype(typeof(x))
# The in-place methods store an out-of-place result with _set_ncadd!, which leaves `y` unchanged if it doesn't fit.

VectorInterface.scalartype(::Type{<:NCMul{C}}) where {C<:Number} = C
VectorInterface.scalartype(::Type{<:NCAdd{C}}) where {C<:Number} = C

# The zero vector has the grades of x (see grading.jl)
function VectorInterface.zerovector(x::MulAdd, ::Type{S}) where {S<:Number}
    return _set_zero_grades!(zero(NCAdd{S,to_add_dict_type(typeof(x))}), _grades(x))
end

function VectorInterface.zerovector!(x::NCAdd)
    gs = _grades(x)
    x.coeff = zero(_scalartype(x))
    empty!(x.dict)
    return _set_zero_grades!(x, gs)
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
    result = β * y + α * x
    # as in add!!, the grades of y are dropped when β isa Zero
    β isa VectorInterface.Zero || return result
    return _with_zero_grades(result, () -> _grades(x))
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
# Autosort is enabled so the result doesn't depend on the setting. The result is promoted to the product of the scalar
# types, so its type doesn't depend on which terms vanish.
function _inner(x, y)
    z = with_autosort(() -> x' * y)
    _reduces_to_scalar(z) || throw(_inner_error())
    return zero(_scalartype(x))' * zero(_scalartype(y)) + scalar(z)
end
_inner_error() = ArgumentError("inner(x, y) requires x' * y to reduce to a scalar; check that adjoint and mul_effect are defined")
_reduces_to_scalar(::Number) = true
_reduces_to_scalar(z::MulAdd) = isscalar(z)
_reduces_to_scalar(z) = false

LinearAlgebra.norm(x::MulAdd) = _norm_from_inner(VectorInterface.inner(x, x))
# p-norms depend on a choice of basis, which a sum of products doesn't have
LinearAlgebra.norm(x::MulAdd, p::Real) = _norm_p(x, p)
_norm_p(x, p::Real) = p == 2 ? LinearAlgebra.norm(x) : throw(ArgumentError("norm(x, p) is only supported for p = 2, the norm induced by inner"))
function _norm_from_inner(z::Number)
    r = real(z)
    return sqrt(max(r, zero(r)))
end
LinearAlgebra.dot(x::MulAdd, y::MulAdd) = VectorInterface.inner(x, y)

# Whether scalartype(T) is defined by a method other than VectorInterface's generic (warning) fallback
function _has_scalartype(T::Type)
    m = which(VectorInterface.scalartype, Tuple{Type{T}})
    return m.sig !== Tuple{typeof(VectorInterface.scalartype),Type}
end

# An atom x as a single-factor product with coefficient one(scalartype(x)), so that the methods below return the same
# type as zerovector(x)
_atom_as_muladd(x::MulAdd) = x
_atom_as_muladd(x) = NCMul(one(_scalartype(x)), [x])

# VectorInterface methods for an atom type T, treating an atom as a single-factor product. Methods with an atom as
# the first argument accept any second argument, so they cover atoms of other types (e.g. inner(A(), B()) after
# @nc A B) without ambiguities.
macro nc_vectorinterface(T)
    T = esc(T)
    quote
        VectorInterface.inner(x::MulAdd, y::$T) = _inner(x, _atom_as_muladd(y))
        VectorInterface.inner(x::$T, y) = _inner(_atom_as_muladd(x), _atom_as_muladd(y))
        LinearAlgebra.norm(x::$T) = _norm_from_inner(VectorInterface.inner(x, x))
        LinearAlgebra.norm(x::$T, p::Real) = _norm_p(x, p)
        LinearAlgebra.dot(x::MulAdd, y::$T) = VectorInterface.inner(x, y)
        LinearAlgebra.dot(x::$T, y) = VectorInterface.inner(x, y)

        # an atom acts as 1 * atom, but an existing scalartype method for T (e.g. if T is an AbstractArray) is kept.
        # The methods whose VectorInterface defaults call scalartype are defined directly, so they work either way.
        if !NonCommutativeProducts._has_scalartype($T)
            VectorInterface.scalartype(::Type{<:$T}) = Int
        end
        VectorInterface.scale(x::$T, α::Number) = VectorInterface.scale(_atom_as_muladd(x), α)
        VectorInterface.scale!!(x::$T, α::Number) = VectorInterface.scale(x, α)
        VectorInterface.scale!!(a::NCAdd, x::$T, α::Number) = add!!(a, _atom_as_muladd(x), α, VectorInterface.Zero())
        VectorInterface.scale!!(y::$T, x, α::Number) = VectorInterface.scale(_atom_as_muladd(x), α)
        VectorInterface.scale!!(y::NCMul, x::$T, α::Number) = VectorInterface.scale(x, α)
        # an atom is immutable, so like an NCMul it can't be the destination of the in-place methods
        VectorInterface.scale!(x::$T, α::Number) = throw(ArgumentError(string(typeof(x), " is immutable; use scale or scale!!")))
        VectorInterface.scale!(y::$T, x, α::Number) = throw(ArgumentError(string(typeof(y), " is immutable; use scale or scale!!")))
        VectorInterface.scale!(y::MulAdd, x::$T, α::Number) = VectorInterface.scale!(y, _atom_as_muladd(x), α)
        VectorInterface.add(y::$T, x, α::Number, β::Number) = VectorInterface.add(_atom_as_muladd(y), _atom_as_muladd(x), α, β)
        VectorInterface.add(y::MulAdd, x::$T, α::Number, β::Number) = VectorInterface.add(y, _atom_as_muladd(x), α, β)
        # (y::T, x::MulAdd) avoids an ambiguity with add!!(y, x::MulAdd, α, β)
        VectorInterface.add!!(y::$T, x, α::Number, β::Number) = add!!(_atom_as_muladd(y), _atom_as_muladd(x), α, β)
        VectorInterface.add!!(y::$T, x::MulAdd, α::Number, β::Number) = add!!(_atom_as_muladd(y), x, α, β)
        VectorInterface.add!(y::$T, x, α::Number, β::Number) = throw(ArgumentError(string(typeof(y), " is immutable; use add or add!!")))
        VectorInterface.add!(y::NCMul, x::$T, α::Number, β::Number) = VectorInterface.add!(y, _atom_as_muladd(x), α, β)
        VectorInterface.add!(y::NCAdd, x::$T, α::Number, β::Number) = VectorInterface.add!(y, _atom_as_muladd(x), α, β)
        VectorInterface.zerovector(x::$T) = VectorInterface.zerovector(_atom_as_muladd(x), _scalartype(x))
        VectorInterface.zerovector(x::$T, ::Type{S}) where {S<:Number} = VectorInterface.zerovector(_atom_as_muladd(x), S)
        VectorInterface.zerovector!!(x::$T) = VectorInterface.zerovector(x)
        VectorInterface.zerovector!!(x::$T, ::Type{S}) where {S<:Number} = VectorInterface.zerovector(x, S)
    end
end
