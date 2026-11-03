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
function VectorInterface.add!(y::NCAdd{C}, x::MulAdd, α::Number, β::Number) where {C}
    isone(β) && _keeps_type(C, typeof(β)) && _add_terms!(y, x, α) && return y
    return _set_ncadd!(y, VectorInterface.add(y, x, α, β))
end
# Sets y to y + α * x by updating only the terms of x, and returns whether it did. It doesn't handle a scalar product,
# or a zero result, whose grades include those of y before the update. Everything is computed before y is mutated.
function _add_terms!(y::NCAdd{C,K}, x::MulAdd, α) where {C,K}
    x isa NCMul && isempty(x.factors) && return false
    newcoeff = convert(C, x isa NCAdd ? additive_coeff(y) + α * additive_coeff(x) : additive_coeff(y))
    updates = Pair{K,C}[]
    newlength = length(y.dict)
    for t in NCterms(x)
        c = α * prefactor(t)
        iszero(c) && continue
        k = term_key(K, t)
        old = get(y.dict, k, nothing)
        v = convert(C, isnothing(old) ? c : old + c)
        newlength += isnothing(old) ? !iszero(v) : -iszero(v)
        push!(updates, k => v)
    end
    newlength == 0 && iszero(newcoeff) && return false
    y.coeff = newcoeff
    for (k, v) in updates
        iszero(v) ? delete!(y.dict, k) : (y.dict[k] = v)
    end
    y.zero_grades = nothing
    return true
end

VectorInterface.add!!(y::MulAdd, x, α::Number, β::Number) = add!!(y, x, α, β)
VectorInterface.add!!(y, x::MulAdd, α::Number, β::Number) = add!!(y, x, α, β)
VectorInterface.add!!(y::MulAdd, x::MulAdd, α::Number, β::Number) = add!!(y, x, α, β)

VectorInterface.inner(x::MulAdd, y::MulAdd) = _inner(x, y)
# The sorted sum of a' * b over the pairs of terms a of x and b of y that can contract, promoted to at least the
# product of the scalar types.
function _inner(x, y)
    C = promote_type(_scalartype(x), _scalartype(y))
    T = Base.promote_op(+, C, C)
    (iszero(x) || iszero(y)) && return zero(_scalartype(x))' * zero(_scalartype(y)) + zero(T)
    xterms, yterms = _collect_terms(x), _collect_terms(y)
    candidates, needed = _contracting_pairs(xterms, yterms)
    bras = _bras(xterms, needed)
    S = _factor_eltype(eltype(bras), eltype(yterms))
    acc = zero(NCAdd{T,NCMul{Int,S,Vector{S}}})
    buffer = S[]
    z = with_autosort(false) do
        fold_retyping!!(acc, candidates) do acc, (i, j, d)
            _inner_pair!!(acc, xterms[i].factors, bras[i], yterms[j], d, buffer)
        end
    end
    _reduces_to_scalar(z) || throw(_inner_error())
    return zero(_scalartype(x))' * zero(_scalartype(y)) + scalar(z)
end
# The terms of x in a Vector with a concrete element type
_collect_terms(x) = fold_terms!!(push!!, Union{}[], x)
"""
    contraction_key(atom)

Return a key that lets `inner(x, y)` skip orthogonal pairs of terms, or `nothing` (the default).

`inner` treats two terms as orthogonal if their keys differ at any position before the first
factor without a key, without sorting them first. Terms must list subsystems in the same order:
`|a⟩|b⟩` and `|b⟩|a⟩` count as orthogonal. Keys must differ between subsystems.

For atoms `k`, `l` with keys and any atom `u`:
- `mul_effect(k', l)` must be a `Number` when the keys are equal.
- `mul_effect(u, k')` must be `nothing` or a `Swap`.

Kets of an orthonormal basis can have keys; bras can't. See [Fast inner products](@ref).
"""
contraction_key(atom) = nothing
# The bras of the terms of x, formed for the terms in the mask `needed` (all of them if it is `nothing`), with the
# others left undefined
function _bras(xterms, needed)
    isnothing(needed) && return [NCMul(adjoint(prefactor(t)), adjoint_factors(t.factors)) for t in xterms]
    is = findall(needed)
    formed = [NCMul(adjoint(prefactor(xterms[i])), adjoint_factors(xterms[i].factors)) for i in is]
    length(is) == length(xterms) && return formed
    bras = similar(formed, length(xterms))
    bras[is] = formed
    return bras
end
_key_at(t, d) = d > length(t.factors) ? nothing : contraction_key(t.factors[d])
# An index over the keyed runs of the kets (the keys of their leading factors, up to the first factor without a key):
# a tree with root node 1, where only the root and the nodes with more than one ket below have children
struct _KeyIndex{K}
    children::Dict{Tuple{Int,UInt},Int} # (node, hash of key) => first child with that hash
    next_sibling::Vector{Int} # next child whose key has the same hash, or 0
    keys::Vector{K} # key that leads to each node
    order::Vector{Int} # ket indices, grouped by node
    ending::Vector{UnitRange{Int}} # order[ending[n]]: kets whose runs end at split node n
    below::Vector{UnitRange{Int}} # order[below[n]]: kets whose runs pass through node n
end
_issplit(idx::_KeyIndex, node) = node == 1 || length(idx.below[node]) > 1
# Splits the nodes by the keys at their depth, with an explicit stack since a path is as long as a product
function _KeyIndex(kets)
    K = Base.promote_op(contraction_key, _factor_eltype(eltype(kets)))
    n = length(kets)
    idx = _KeyIndex{K}(Dict{Tuple{Int,UInt},Int}(), [0], Vector{K}(undef, 1), collect(1:n), [1:0], [1:n])
    # scratch space: the child that each ket goes to (0 if its run ends), the numbers of kets of the new children, and
    # the kets in their new order
    child_of, counts, sorted = zeros(Int, n), Int[], similar(idx.order)
    stack = [(1, 0)]
    while !isempty(stack)
        node, d = pop!(stack)
        for c in _split_node!(idx, kets, node, d, child_of, counts, sorted)
            length(idx.below[c]) > 1 && push!(stack, (c, d + 1))
        end
    end
    return idx
end
# Adds the children of `node`, at depth d, by the keys at depth d + 1, and sorts the kets below it stably: first those
# whose runs end at node, then those of each child in turn. Returns the range of the children.
function _split_node!(idx, kets, node, d, child_of, counts, sorted)
    r = idx.below[node]
    first_child = length(idx.keys) + 1
    empty!(counts)
    nending = 0
    for k in r
        key = _key_at(kets[idx.order[k]], d + 1)
        if isnothing(key)
            child_of[k] = 0
            nending += 1
        else
            c = _child!(idx, node, key)
            child_of[k] = c
            if c - first_child + 1 > length(counts)
                push!(counts, 0)
            end
            counts[c-first_child+1] += 1
        end
    end
    idx.ending[node] = first(r):first(r)+nending-1
    children = first_child:length(idx.keys)
    start = first(r) + nending
    for (m, c) in enumerate(children)
        idx.below[c] = start:start+counts[m]-1
        start += counts[m]
    end
    # the next free slots of the kets whose runs end here and of each child (counted back from its end)
    ending_slot = first(r)
    for k in r
        c = child_of[k]
        if c == 0
            sorted[ending_slot] = idx.order[k]
            ending_slot += 1
        else
            m = c - first_child + 1
            sorted[last(idx.below[c])-counts[m]+1] = idx.order[k]
            counts[m] -= 1
        end
    end
    copyto!(idx.order, first(r), sorted, first(r), length(r))
    return children
end
# The child of `node` with the key `key`, which is added if there is none
function _child!(idx, node, key)
    h = hash(key)
    c = get(idx.children, (node, h), 0)
    if c == 0
        c = _add_node!(idx, key)
        idx.children[(node, h)] = c
        return c
    end
    while !isequal(idx.keys[c], key)
        if idx.next_sibling[c] == 0
            idx.next_sibling[c] = _add_node!(idx, key)
            return idx.next_sibling[c]
        end
        c = idx.next_sibling[c]
    end
    return c
end
function _add_node!(idx, key)
    push!(idx.keys, key)
    push!(idx.next_sibling, 0)
    push!(idx.ending, 1:0)
    push!(idx.below, 1:0)
    return length(idx.keys)
end
# The child of `node` with the key `key`, or 0 if there is none
function _find_child(idx, node, key)
    c = get(idx.children, (node, hash(key)), 0)
    while c != 0 && !isequal(idx.keys[c], key)
        c = idx.next_sibling[c]
    end
    return c
end
# The triples (i, j, d) such that xterms[i] and yterms[j] may contract, where the first d keys of their runs are equal,
# and a mask of the terms of x that appear in them, or `nothing` if all do. The index is over the side with fewer
# terms.
function _contracting_pairs(xterms, yterms)
    index_x = length(xterms) < length(yterms)
    indexed, walking = index_x ? (xterms, yterms) : (yterms, xterms)
    if all(t -> isnothing(_key_at(t, 1)), indexed)
        return Iterators.map(((i, j),) -> (i, j, 0), Iterators.product(eachindex(xterms), eachindex(yterms))), nothing
    end
    idx = _KeyIndex(indexed)
    # (term of the walking side, range of idx.order, d)
    segments = Tuple{Int,UnitRange{Int},Int}[]
    for (w, t) in enumerate(walking)
        _push_segments!(segments, idx, w, t.factors)
    end
    needed = falses(length(xterms))
    if index_x
        for (_, r, _) in segments, k in r
            needed[idx.order[k]] = true
        end
        return ((idx.order[k], j, d) for (j, r, d) in segments for k in r), needed
    end
    for (i, _, _) in segments
        needed[i] = true
    end
    return ((i, idx.order[k], d) for (i, r, d) in segments for k in r), needed
end
# Walks the index along the keyed run of `factors`, pairing term w with the kets whose runs end on the way, and with
# all kets below the node where its run ends or the index stops splitting.
function _push_segments!(segments, idx, w, factors)
    node, d = 1, 0
    while true
        key = _issplit(idx, node) && d < length(factors) ? contraction_key(factors[d+1]) : nothing
        if isnothing(key)
            push!(segments, (w, idx.below[node], d))
            return segments
        end
        isempty(idx.ending[node]) || push!(segments, (w, idx.ending[node], d))
        node = _find_child(idx, node, key)
        node == 0 && return segments
        d += 1
    end
end
# acc + bra * ket sorted, contracting the common prefix of the keyed runs of the term of x (with factors xfactors)
# and ket, of which the first d keys are known to be equal. The rest of the product is formed in `buffer`.
function _inner_pair!!(acc, xfactors, bra, ket, d, buffer)
    bf, kf = bra.factors, ket.factors
    # Known at compile time: factors whose type has no keys contract nothing
    Base.promote_op(contraction_key, eltype(kf)) === Nothing && return _mul_term_pair!!(acc, bra, ket, buffer)
    nb = length(bf)
    effect = true
    p = 0
    while p < min(nb, length(kf))
        if p >= d
            kk = contraction_key(kf[p+1])
            isnothing(kk) && break
            kx = contraction_key(xfactors[p+1])
            isnothing(kx) && break
            isequal(kx, kk) || return acc
        end
        e = _mul_effect(buffer, bf[nb-p], kf[p+1])
        e isa Number || throw(_contraction_key_error(bf[nb-p], kf[p+1]))
        effect *= e
        iszero(effect) && return acc
        p += 1
    end
    p == 0 && return _mul_term_pair!!(acc, bra, ket, buffer)
    append!(empty!(buffer), Iterators.take(bf, nb - p), Iterators.drop(kf, p))
    return sort_into!!(acc, prefactor(bra) * prefactor(ket) * effect, buffer)
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
