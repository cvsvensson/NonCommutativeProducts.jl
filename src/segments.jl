"""
    atom_group(atom)

The group that `atom` belongs to, used to split products into segments with [`foreach_segment`](@ref),
[`split_by_group`](@ref) and [`segment_ranges`](@ref). Two factors are in the same group if their groups are `isequal`.
Defaults to `typeof(atom)`, so parametric variants of one type are different groups, e.g. `Fermion{Int}` and
`Fermion{Int32}`.

Define methods for your types to group atoms by value, e.g. by the Hilbert space they act on, or to put several types
in one group. The group can also encode a kind, e.g. `(space_id, :ket)` and `(space_id, :op)`, to keep kets and
operators of the same subsystem in separate segments. Groups do not change how products are sorted: the order of
factors is set by `mul_effect` alone. Declaring how groups commute, with [`group_commutation`](@ref), makes sorting
faster.

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

_segments_of_sum_error(name) = ArgumentError("$name is not defined for a sum; apply it to each term, e.g. via `terms(x)`")
# An atom is a product with one factor, as for `factors` and `terms`
function _as_product(x, name)
    _is_atom(typeof(x)) || throw(ArgumentError("$name takes a product, an atom, or a vector or tuple of factors, not a $(typeof(x))"))
    return NCMul(x)
end

"""
    segment_ranges(x)

The segments of the product `x` (see [`foreach_segment`](@ref)) as a `Vector` of pairs `g => range`, from left to
right, where `g` is the [`atom_group`](@ref) of the segment and `range` the indices of its factors in `factors(x)`; the
vector has the narrowest pair type that holds all pairs. `x` is an `NCMul`, an atom, or a vector or tuple of factors.
A group appears once for each run of its factors. Throws an `ArgumentError` for a sum; apply it to each of its
`terms`.

# Example
```julia
segment_ranges(2 * a1 * a2 * b * a3) # [A => 1:2, B => 3:3, A => 4:4], with A the group of `a` and B that of `b`
```
"""
segment_ranges(x::NCMul) = segment_ranges(x.factors)
segment_ranges(::NCAdd) = throw(_segments_of_sum_error("segment_ranges"))
segment_ranges(x) = segment_ranges(_as_product(x, "segment_ranges"))
function segment_ranges(factors::Union{AbstractVector,Tuple})
    out = Pair[]
    _foreach_run((g, r) -> push!(out, g => r), factors)
    return [p for p in out]
end

"""
    foreach_segment(f, x)

Call `f(g, seg)` for each segment of the product `x`, from left to right. A segment is a maximal run of consecutive
factors whose [`atom_group`](@ref)s are `isequal`, `g` is that group and `seg` is an `AbstractVector` of the factors.
`x` is an `NCMul`, an atom, or a vector or tuple of factors; the coefficient of a product is ignored, and a product
without factors has no segments. Returns `nothing`. Throws an `ArgumentError` for a sum; apply it to each of its
`terms`. Use [`segment_ranges`](@ref) for the positions of the segments.

Segments are runs, so interleaved groups give several runs, and a group appears once for each run of its factors:
the factors of `a1*b*a2`, with `a1` and `a2` in one group and `b` in another, give the three segments `[a1]`, `[b]`
and `[a2]`. Groups are contiguous only if `mul_effect` sorts them so.

`seg` has a concrete element type whenever all its factors have the same concrete type:
- if the element type of the factors of `x` is concrete, `seg` is a view of them;
- otherwise `seg` is a new `Vector`, with the type of the factors as element type if they share one, and the element
  type of the factors of `x` if not.

`f` is called through a function barrier, so it is compiled for the type of each segment, also when the factors of `x`
have an abstract element type. A view aliases `x.factors`, which by the convention of this package are never mutated,
so `seg` stays valid after `f` returns. `seg` must not be mutated; copy it if you intend to mutate it.

# Example
```julia
foreach_segment((g, seg) -> println(g, " => ", seg), 2 * a1 * a2 * b * a3)
# prints the group of `a` with [a1, a2], the group of `b` with [b], and the group of `a` with [a3]
```
"""
foreach_segment(f, x::NCMul) = foreach_segment(f, x.factors)
foreach_segment(f, ::NCAdd) = throw(_segments_of_sum_error("foreach_segment"))
foreach_segment(f, x) = foreach_segment(f, _as_product(x, "foreach_segment"))
foreach_segment(f, factors::Tuple) = foreach_segment(f, factors_vector(factors))
function foreach_segment(f, factors::AbstractVector)
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

The segments of `x` (see [`foreach_segment`](@ref)) as a `Vector` of products, from left to right. `x` is an `NCMul`
or an atom; a sum throws an `ArgumentError`. The first product has the coefficient of `x` and the others have
coefficient `one(coefficient(x))`, with the coefficient type of `x`, so the products of the parts, in order and
unsorted, equals `x`. A zero product with factors is split as any other, with the zero coefficient on the first part.
A product without factors gives the vector `[x]`. The factors of each part are a `Vector` copy, narrowed to a concrete
element type when they share a concrete type.

If the element type of the factors of `x` is concrete, the result is a concretely typed vector and `split_by_group`
is inferable. Otherwise the vector has the narrowest element type that holds all parts, which is in general abstract
and not known at compile time; in hot paths, use [`foreach_segment`](@ref) instead.

The group of a part `p` is `atom_group(first(factors(p)))`. A group appears once for each run of its factors, as in
`foreach_segment`.

# Example
```julia
split_by_group(2 * a1 * a2 * b) # [2 * a1 * a2, 1 * b]
```
"""
split_by_group(::NCAdd) = throw(_segments_of_sum_error("split_by_group"))
split_by_group(x) = split_by_group(_as_product(x, "split_by_group"))
function split_by_group(x::NCMul{C}) where {C}
    isempty(x.factors) && return [x]
    factors = factors_vector(x.factors)
    return _split_by_group(C, prefactor(x), factors, Val(isconcretetype(eltype(factors))))
end
function _split_by_group(::Type{C}, c, factors::AbstractVector{T}, ::Val{true}) where {C,T}
    parts = Vector{NCMul{C,T,Vector{T}}}()
    _foreach_run(factors) do _, r
        push!(parts, NCMul{C}(isempty(parts) ? c : one(c), _owned_copy(view(factors, r))))
    end
    return parts
end
function _split_by_group(::Type{C}, c, factors, ::Val{false}) where {C}
    parts = NCMul[]
    foreach_segment(factors) do _, seg
        push!(parts, NCMul{C}(isempty(parts) ? c : one(c), _owned_copy(seg)))
    end
    return [p for p in parts]
end
# `seg` as a Vector that can be stored in a product: the copy made by foreach_segment, or a copy of anything else
_owned_copy(seg::Vector) = seg
_owned_copy(seg::AbstractVector) = collect(seg)
