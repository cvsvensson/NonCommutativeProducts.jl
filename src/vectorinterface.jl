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
# The sorted sum of a' * b over the pairs of terms a of x and b of y that can contract. The result is promoted to the
# product of the scalar types, so its type doesn't depend on which terms vanish or on the autosort setting.
function _inner(x, y)
    xterms, yterms = _collect_terms(x), _collect_terms(y)
    bras = [NCMul(adjoint(prefactor(t)), adjoint_factors(t.factors)) for t in xterms]
    S = _factor_eltype(eltype(bras), eltype(yterms))
    C = promote_type(_scalartype(x), _scalartype(y))
    acc = zero(NCAdd{Base.promote_op(+, C, C),NCMul{Int,S,Vector{S}}})
    buffer = S[]
    z = with_autosort(false) do
        fold_retyping!!(acc, _contracting_pairs(xterms, yterms)) do acc, (i, j, p)
            _inner_pair!!(acc, bras[i], yterms[j], p, buffer)
        end
    end
    _reduces_to_scalar(z) || throw(_inner_error())
    return zero(_scalartype(x))' * zero(_scalartype(y)) + scalar(z)
end
# The terms of x in a Vector with a concrete element type
_collect_terms(x) = fold_terms!!(push!!, Union{}[], x)
"""
    contraction_key(atom)

A key that lets `inner(x, y)` skip products that vanish, or `nothing` (the default). The *keyed run* of a term is the
keys of its leading factors, up to the first factor without a key. `inner` forms `a' * b` for a term `a` of `x` and a
term `b` of `y` only if one of their keyed runs is a prefix of the other, and contracts the factors of the common
prefix directly. This gives the same result as sorting `x' * y` when, for all atoms `k` and `l` with keys:
- `mul_effect(k', l)` is zero whenever their keys differ (by `isequal`),
- `mul_effect(k', l)` is a `Number` whenever their keys are equal, and
- `mul_effect(u, k')` is `nothing` for every atom `u`.

Then once `a₁' * b₁` reduces to a scalar, `a₂'` meets `b₂`, and so on along the common prefix.

Kets of an orthonormal basis can have keys when nothing acts on a bra from the left. Bras can't, since
`⟨1|' * ⟨2| = |1⟩⟨2|` is nonzero:
```julia
NonCommutativeProducts.contraction_key(k::Ket) = k.bra ? nothing : k.n
```
"""
contraction_key(atom) = nothing
# A trie over the keyed runs of the kets. Node 1 is the root, and the kets whose runs end at node n are
# order[own[n]], and those in its subtree are order[sub[n]].
struct _KeyTrie
    children::Vector{Dict{Any,Int}}
    order::Vector{Int}
    own::Vector{UnitRange{Int}}
    sub::Vector{UnitRange{Int}}
end
function _KeyTrie(kets)
    children, ends = [Dict{Any,Int}()], [Int[]]
    for (j, t) in pairs(kets)
        node = 1
        for f in t.factors
            key = contraction_key(f)
            isnothing(key) && break
            node = get!(children[node], key) do
                push!(children, Dict{Any,Int}()), push!(ends, Int[])
                length(children)
            end
        end
        push!(ends[node], j)
    end
    n = length(children)
    return _number_subtrees!(_KeyTrie(children, Int[], fill(1:0, n), fill(1:0, n)), ends)
end
# Lists the kets of each subtree contiguously in trie.order, those whose runs end at its root first. Depth first with
# an explicit stack, since a path is as long as a product.
function _number_subtrees!(trie, ends)
    stack = [(1, false)]
    while !isempty(stack)
        node, finished = pop!(stack)
        if finished
            trie.sub[node] = first(trie.own[node]):length(trie.order)
        else
            start = length(trie.order) + 1
            append!(trie.order, ends[node])
            trie.own[node] = start:length(trie.order)
            push!(stack, (node, true))
            foreach(child -> push!(stack, (child, false)), values(trie.children[node]))
        end
    end
    return trie
end
# The triples (i, j, p) such that xterms[i] and kets[j] can contract, with p the length of the common prefix of their
# keyed runs. Each pair appears once.
function _contracting_pairs(xterms, kets)
    trie = _KeyTrie(kets)
    isempty(trie.children[1]) && return ((i, j, 0) for (i, j) in Iterators.product(eachindex(xterms), eachindex(kets)))
    segments = Tuple{Int,UnitRange{Int},Int}[]
    for (i, t) in pairs(xterms)
        _push_segments!(segments, trie, i, t.factors)
    end
    return ((i, trie.order[k], p) for (i, range, p) in segments for k in range)
end
# Walks the trie along the keyed run of `factors`, pairing term i with the kets whose runs end on the way, and with
# the whole subtree where its own run ends
function _push_segments!(segments, trie, i, factors)
    node = 1
    for d in 0:length(factors)
        children = trie.children[node]
        key = d == length(factors) || isempty(children) ? nothing : contraction_key(factors[d+1])
        isnothing(key) && return push!(segments, (i, trie.sub[node], d))
        isempty(trie.own[node]) || push!(segments, (i, trie.own[node], d))
        node = get(children, key, 0)
        node == 0 && return segments
    end
    return segments
end
# acc + bra * ket sorted, where the first p factors of ket contract with the last p factors of bra to a scalar. The
# rest of the product is formed in `buffer`.
function _inner_pair!!(acc, bra, ket, p, buffer)
    p == 0 && return _mul_term_pair!!(acc, bra, ket, buffer)
    coeff = prefactor(bra) * prefactor(ket)
    iszero(coeff) && return acc
    bf, kf = bra.factors, ket.factors
    nb = length(bf)
    for d in 1:p
        effect = _mul_effect(buffer, bf[nb-d+1], kf[d])
        effect isa Number || throw(_contraction_key_error(bf[nb-d+1], kf[d]))
        coeff *= effect
        iszero(coeff) && return acc
    end
    append!(empty!(buffer), Iterators.take(bf, nb - p), Iterators.drop(kf, p))
    return sort_into!!(acc, coeff, buffer)
end
function _contraction_key_error(a, b)
    return ArgumentError(string("mul_effect(", a, ", ", b, ") must be a Number, since contraction_key(", a,
        "') and contraction_key(", b, ") are equal; see the contraction_key docstring"))
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
