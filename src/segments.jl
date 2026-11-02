"""
    atom_group(atom)

The group that `atom` belongs to, used to split products into segments with [`foreach_segment`](@ref),
[`split_by_group`](@ref) and [`segment_ranges`](@ref). Two factors are in the same group if their groups are `isequal`.
Defaults to `typeof(atom)`, so parametric variants of one type are different groups, e.g. `Fermion{Int}` and
`Fermion{Int32}`.

Define methods for your types to group atoms by value, e.g. by the Hilbert space they act on, or to put several types
in one group. The group can also encode a kind, e.g. `(space_id, :ket)` and `(space_id, :op)`, to keep kets and
operators of the same subsystem in separate segments. Groups do not change how products are sorted: the order of
factors is set by `mul_effect` alone.

# Example
```julia
NonCommutativeProducts.atom_group(f::Fermion) = f.site # fermions on the same site form one group
```
"""
atom_group(atom) = typeof(atom)

# Call h(g, range) for each maximal run of factors with isequal groups, from left to right. The group of each factor
# is computed once: the group that ends a run is carried over as the group of the next run.
function _foreach_run(h, factors)
    i, last = firstindex(factors), lastindex(factors)
    i > last && return nothing
    g = atom_group(factors[i])
    while i <= last
        j = i
        gnext = g
        while j < last
            gnext = atom_group(factors[j+1])
            isequal(gnext, g) || break
            j += 1
        end
        h(g, i:j)
        i = j + 1
        g = gnext
    end
    return nothing
end

# An atom is a product with one factor, as for `factors` and `terms`
function _as_product(x, name)
    _is_atom(typeof(x)) || throw(ArgumentError("$name takes a product or an atom, not a $(typeof(x))"))
    return NCMul(x)
end

"""
    segment_ranges(x)

The segments of `x` (see [`foreach_segment`](@ref)) as a `Vector` of pairs `g => range`, from left to right, where `g`
is the [`atom_group`](@ref) of the segment and `range` the indices of its factors in `factors(x)`.

# Example
```julia
segment_ranges(NCMul(2, [a1, a2, b, a3])) # [A => 1:2, B => 3:3, A => 4:4], with A the group of `a` and B that of `b`
```
"""
function segment_ranges(x::NCMul)
    out = Pair[]
    _foreach_run((g, r) -> push!(out, g => r), x.factors)
    return isempty(out) ? out : [p for p in out]
end
segment_ranges(x::NCAdd) = segment_ranges(_only_term(x, "segment_ranges"))
segment_ranges(x) = segment_ranges(_as_product(x, "segment_ranges"))

"""
    foreach_segment(f, x)

Call `f(g, seg)` for each segment of the product `x`, from left to right. A segment is a maximal run of consecutive
factors whose [`atom_group`](@ref)s are `isequal`, `g` is that group and `seg` is an `AbstractVector` of the factors.
`x` is a product or an atom; the coefficient of a product is ignored, and a product without factors has no segments.
A sum with one term counts as that term; more terms throw an `ArgumentError`, so apply it to each of its `terms`.

Interleaved groups give several segments.

`seg` has a concrete element type when its factors share one. It may alias the factors of `x`, so do not mutate it.

# Example
```julia
foreach_segment((g, seg) -> println(g, " => ", seg), NCMul(2, [a1, a2, b, a3]))
# prints the group of `a` with [a1, a2], the group of `b` with [b], and the group of `a` with [a3]
```
"""
foreach_segment(f, x::NCMul) = _foreach_segment(f, factors_vector(x.factors))
foreach_segment(f, x::NCAdd) = foreach_segment(f, _only_term(x, "foreach_segment"))
foreach_segment(f, x) = foreach_segment(f, _as_product(x, "foreach_segment"))
function _foreach_segment(f, factors::AbstractVector)
    _foreach_run((g, r) -> f(g, _segment(factors, r)), factors)
    return nothing
end
# The factors in `range`, with a concrete element type if they share a concrete type
_segment(factors::AbstractVector{T}, range) where {T} = isconcretetype(T) ? view(factors, range) : _narrowed_copy(factors, range)
function _narrowed_copy(factors, range)
    U = typeof(factors[first(range)])
    all(k -> factors[k] isa U, range) || return factors[range]
    return _copy_as(U, factors, range)
end
# A function barrier, so that the copy is specialized on U
_copy_as(::Type{U}, factors, range) where {U} = U[factors[k] for k in range]

"""
    split_by_group(x)

The segments of `x` as a `Vector` of products, from left to right. The first has the coefficient of `x`, the others
one, so their unsorted product equals `x`. A product without factors gives `[x]`. Inferable for factors of a concrete
element type; otherwise prefer [`foreach_segment`](@ref) in hot paths.

# Example
```julia
split_by_group(2 * a1 * a2 * b) # [2 * a1 * a2, 1 * b]
```
"""
split_by_group(x::NCAdd) = split_by_group(_only_term(x, "split_by_group"))
split_by_group(x) = split_by_group(_as_product(x, "split_by_group"))
function split_by_group(x::NCMul{C}) where {C}
    factors = factors_vector(x.factors)
    return _split_by_group(x, C, prefactor(x), factors, Val(isconcretetype(eltype(factors))))
end
# The parts are Vector-backed whatever the factors of `x` are, so that the result type depends only on their eltype
function _split_by_group(_, ::Type{C}, c, factors::AbstractVector{T}, ::Val{true}) where {C,T}
    parts = Vector{NCMul{C,T,Vector{T}}}()
    _foreach_run(factors) do _, r
        push!(parts, NCMul{C}(isempty(parts) ? c : one(c), _owned_copy(view(factors, r))))
    end
    isempty(parts) && push!(parts, NCMul{C}(c, T[]))
    return parts
end
function _split_by_group(x, ::Type{C}, c, factors, ::Val{false}) where {C}
    isempty(factors) && return [x]
    parts = NCMul[]
    _foreach_segment(factors) do _, seg
        push!(parts, NCMul{C}(isempty(parts) ? c : one(c), _owned_copy(seg)))
    end
    return [p for p in parts]
end
# `seg` as a Vector that can be stored in a product: the copy made by foreach_segment, or a copy of anything else
_owned_copy(seg::Vector) = seg
_owned_copy(seg::AbstractVector) = collect(seg)
