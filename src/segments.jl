"""
    group(atom)

The group that `atom` belongs to, used to split products into segments with [`foreach_segment`](@ref) and
[`split_by_group`](@ref). Two factors are in the same group if their groups are `isequal`. Defaults to `typeof(atom)`.

Define methods for your types to group atoms by value, e.g. by the Hilbert space they act on, or to put several types
in one group. Groups do not change how products are sorted: the order of factors is set by `mul_effect` alone.

# Example
```julia
NonCommutativeProducts.group(f::Fermion) = f.site # fermions on the same site form one group
```
"""
group(atom) = typeof(atom)

"""
    foreach_segment(f, x)

Call `f(g, seg)` for each segment of the product `x`, from left to right. A segment is a maximal run of consecutive
factors whose [`group`](@ref)s are `isequal`, `g` is that group and `seg` is an `AbstractVector` of the factors. `x`
is an `NCMul` or a vector or tuple of factors; the coefficient of a product is ignored, and a product without factors
has no segments. Returns `nothing`.

Segments are runs, so a group appears once for each run of its factors: the factors of `a1*b*a2`, with `a1` and
`a2` in one group and `b` in another, give the three segments `[a1]`, `[b]` and `[a2]`. Groups are contiguous only if
`mul_effect` sorts them so.

`seg` has a concrete element type whenever all its factors have the same concrete type:
- if the element type of the factors of `x` is concrete, `seg` is a view of them;
- otherwise `seg` is a new `Vector`, with the type of the factors as element type if they share one, and the element
  type of the factors of `x` if not.

`f` is called through a function barrier, so it is compiled for the type of each segment, also when the factors of `x`
have an abstract element type. `seg` must not be mutated, and should be copied to be kept after `f` returns.

# Example
```julia
foreach_segment((g, seg) -> println(g, " => ", seg), 2 * a1 * a2 * b * a3)
# prints the group of `a` with [a1, a2], the group of `b` with [b], and the group of `a` with [a3]
```
"""
foreach_segment(f, x::NCMul) = foreach_segment(f, x.factors)
foreach_segment(f, factors::Tuple) = foreach_segment(f, factors_vector(factors))
function foreach_segment(f, factors::AbstractVector)
    i, last = firstindex(factors), lastindex(factors)
    while i <= last
        g = group(factors[i])
        j = i
        while j < last && isequal(group(factors[j+1]), g)
            j += 1
        end
        f(g, _segment(factors, i:j))
        i = j + 1
    end
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
    split_by_group(x::NCMul)

The segments of `x` (see [`foreach_segment`](@ref)) as a `Vector` of products, from left to right. The first product
has the coefficient of `x` and the others have coefficient `one(coefficient(x))`, so the products of the parts, in
order and unsorted, equals `x`. A product without factors gives the vector `[x]`. The factors of each part are a copy,
narrowed to a concrete element type when they share a concrete type, and the vector has the narrowest element type
that holds all parts.

The group of a part `p` is `group(first(factors(p)))`. A group appears once for each run of its factors, as in
`foreach_segment`.

# Example
```julia
split_by_group(2 * a1 * a2 * b) # [2 * a1 * a2, 1 * b]
```
"""
function split_by_group(x::NCMul)
    isempty(x.factors) && return [x]
    c = prefactor(x)
    parts = NCMul[]
    foreach_segment(x) do _, seg
        push!(parts, NCMul(isempty(parts) ? c : one(c), _owned_copy(seg)))
    end
    return [p for p in parts]
end
# `seg` as a Vector that can be stored in a product: a copy of a view, or the copy made by foreach_segment
_owned_copy(seg::SubArray) = collect(seg)
_owned_copy(seg::Vector) = seg
