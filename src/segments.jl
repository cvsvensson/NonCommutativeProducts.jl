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
function atom_group(atom)
    _is_atom(typeof(atom)) || throw(ArgumentError("$(typeof(atom)) has no group: it is not declared with `@nc`, and `atom_group` is not defined for it"))
    return typeof(atom)
end

# An atom is a product with one factor, as for `factors` and `terms`
function _as_product(x, name)
    _is_atom(typeof(x)) || throw(ArgumentError("$name takes a product or an atom, not a $(typeof(x))"))
    return NCMul(x)
end

# Iterates `g => range` over the maximal runs of factors with isequal groups. The group of each factor is computed once:
# the group that ends a run is in the state as the group of the next run.
struct SegmentRanges{F}
    factors::F
end
Base.IteratorSize(::Type{<:SegmentRanges}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:SegmentRanges}) = Base.EltypeUnknown()
@inline function Base.iterate(s::SegmentRanges)
    factors = s.factors
    isempty(factors) && return nothing
    i = firstindex(factors)
    return iterate(s, (i, atom_group(factors[i])))
end
@inline function Base.iterate(s::SegmentRanges, (i, g))
    factors = s.factors
    last = lastindex(factors)
    i > last && return nothing
    j = i
    gnext = g
    while j < last
        gnext = atom_group(factors[j+1])
        isequal(gnext, g) || break
        j += 1
    end
    return g => i:j, (j + 1, gnext)
end

"""
    segment_ranges(x)

A lazy iterator over the segments of `x` (see [`foreach_segment`](@ref)), from left to right. It gives pairs
`g => range`, where `g` is the [`atom_group`](@ref) of the segment and `range` the indices of its factors in
`factors(x)`. Use `collect` to get a `Vector`.

# Example
```julia
collect(segment_ranges(NCMul(2, [a1, a2, b, a3]))) # [A => 1:2, B => 3:3, A => 4:4], with A the group of `a` and B that of `b`

# Whether a product with factors has them all in one group, and which
g, r = first(segment_ranges(x))
length(r) == length(factors(x))
```
"""
segment_ranges(x::NCMul) = SegmentRanges(x.factors)
segment_ranges(x::NCAdd) = segment_ranges(_only_term(x, "segment_ranges"))
function segment_ranges(x)
    _is_atom(typeof(x)) || throw(ArgumentError("segment_ranges takes a product or an atom, not a $(typeof(x))"))
    return SegmentRanges((x,))
end

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
    for (g, r) in SegmentRanges(factors)
        f(g, _segment(factors, r))
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
    split_by_group(x)

The segments of `x` (see [`foreach_segment`](@ref)) as a `Vector` of pairs `g => part`, from left to right, where `g` is
the [`atom_group`](@ref) of the segment and `part` a product of its factors. The first part has the coefficient of `x`,
the others one, so the unsorted product of the parts equals `x`. A product without factors gives an empty `Vector`.

Segments are runs, not groups: interleaved factors give several parts with the same `g`. A part may share its factors
with `x`. For factors of an abstract element type, prefer [`segment_ranges`](@ref) or [`foreach_segment`](@ref) in hot
paths.

# Example
```julia
split_by_group(2 * a1 * a2 * b) # [A => 2 * a1 * a2, B => 1 * b], with A the group of `a` and B that of `b`
```
"""
split_by_group(x::NCAdd) = split_by_group(_only_term(x, "split_by_group"))
split_by_group(x) = split_by_group(_as_product(x, "split_by_group"))
function split_by_group(x::NCMul{C}) where {C}
    factors = factors_vector(x.factors)
    T = eltype(factors)
    # Fixed pair types, as a comprehension or `g => part` is 2-4x slower for abstract T, where pair types are made at runtime
    G = Base.promote_op(atom_group, T)
    P = isconcretetype(T) ? NCMul{C,T,Vector{T}} : NCMul{C}
    parts = Pair{G,P}[]
    for (g, r) in SegmentRanges(factors)
        push!(parts, Pair{G,P}(g, _part(x, factors, r)))
    end
    return parts
end
# Without a copy if `x` is already what the only part would be
function _part(x::NCMul{C,T,Vector{T}}, factors, r) where {C,T}
    isconcretetype(T) && length(r) == length(factors) && return x
    return _copied_part(x, factors, r)
end
_part(x, factors, r) = _copied_part(x, factors, r)
function _copied_part(x::NCMul{C}, factors, r) where {C}
    c = prefactor(x)
    return NCMul{C}(first(r) == firstindex(factors) ? c : one(c), _owned_copy(_segment(factors, r)))
end
# `seg` as a Vector that can be stored in a product: the copy made by _segment, or a copy of anything else
_owned_copy(seg::Vector) = seg
_owned_copy(seg::AbstractVector) = collect(seg)
