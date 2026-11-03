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
# Fixed G, so that the state and the pairs are not built with runtime types if the group type is unknown.
struct SegmentState{G}
    i::Int
    g::G
end
@inline function Base.iterate(s::SegmentRanges)
    factors = s.factors
    isempty(factors) && return nothing
    i = firstindex(factors)
    return iterate(s, SegmentState{Base.promote_op(atom_group, eltype(factors))}(i, atom_group(factors[i])))
end
@inline function Base.iterate(s::SegmentRanges, st::SegmentState{G}) where {G}
    factors = s.factors
    last = lastindex(factors)
    i, g = st.i, st.g
    i > last && return nothing
    j, gnext = _run_end(factors, i, g, last)
    return Pair{G,UnitRange{Int}}(g, i:j), SegmentState{G}(j + 1, gnext)
end
# The last index of the run that starts at i with group g, and the group of the next factor
@inline function _run_end(factors, i, g, last)
    j = i
    isconcretetype(eltype(factors)) || return _dynamic_run_end(factors, j, g, last)
    gnext = g
    while j < last
        gnext = atom_group(factors[j+1])
        isequal(gnext, g) || break
        j += 1
    end
    return j, gnext
end
@inline function _dynamic_run_end(factors, j, g, last)
    while j < last
        r = _scan(factors, j, g, last, factors[j+1])
        r isa Int || return r[1]::Int, r[2]
        j = r
    end
    return j, g
end
# A function barrier: one dispatch for each run of factors of one type, not two for each factor. Gives the new j if the
# run ended by the end of the factors or by another factor type, or (j, group of the next factor) if it ended by its group.
function _scan(factors, j, g, last, f::F) where {F}
    while true
        h = atom_group(f)
        isequal(h, g) || return (j, h)
        j += 1
        j == last && return j
        next = factors[j+1]
        next isa F || return j
        f = next
    end
end

"""
    segment_ranges(x)

A lazy iterator over the segments of `x` (see [`foreach_segment`](@ref)), from left to right. It gives pairs
`g => range`, where `g` is the [`atom_group`](@ref) of the segment and `range` the indices of its factors in
`factors(x)`.

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
segment_ranges(x) = segment_ranges(_as_product(x, "segment_ranges"))

"""
    foreach_segment(f, x)

Call `f(g, seg)` for each segment of the product `x`, from left to right. A segment is a maximal run of consecutive
factors whose [`atom_group`](@ref)s are `isequal`, `g` is that group and `seg` is an `AbstractVector` of the factors.
`x` is a product or an atom; the coefficient of a product is ignored, and a product without factors has no segments.
A sum with one term counts as that term; any other sum, including zero, throws an `ArgumentError`, so apply it to each
of its `terms`.

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
_segment(factors::AbstractVector{T}, range) where {T} = isconcretetype(T) ? view(factors, range) : _narrowed_copy(factors, range, factors[first(range)])
# A function barrier on the type of the first factor
function _narrowed_copy(factors, range, ::U) where {U}
    out = Vector{U}(undef, length(range))
    for (n, k) in enumerate(range)
        f = factors[k]
        f isa U || return factors[range]
        out[n] = f
    end
    return out
end

"""
    split_by_group(x)

The segments of `x` (see [`foreach_segment`](@ref)) as a `Vector` of pairs `g => part`, from left to right, where `g` is
the [`atom_group`](@ref) of the segment and `part` a product of its factors. The first part has the coefficient of `x`,
the others one, so the unsorted product of the parts equals `x`, except that a product without factors or a scalar-only
sum gives an empty `Vector`. A sum with several terms, or a zero sum, throws an `ArgumentError`.

Segments are runs, not groups: interleaved factors give several parts with the same `g`. A part may share its factors
with `x`, so do not mutate it. For factors of an abstract element type, prefer [`segment_ranges`](@ref) or [`foreach_segment`](@ref) in hot
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
    PP = Pair{G,P}
    segments = SegmentRanges(factors)
    y = iterate(segments)
    y === nothing && return PP[]
    (g, r), state = y
    first_part = PP(g, _part(x, factors, r))
    # One segment is the common case: `push!` would grow a Vector with more capacity, by a slower path
    if last(r) == lastindex(factors)
        parts = Vector{PP}(undef, 1)
        parts[1] = first_part
        return parts
    end
    parts = PP[first_part]
    y = iterate(segments, state)
    while y !== nothing
        (g, r), state = y
        push!(parts, PP(g, _part(x, factors, r)))
        y = iterate(segments, state)
    end
    return parts
end
# Without a copy if `x` is already what the only part would be
function _part(x::NCMul{C}, factors::AbstractVector{T}, r) where {C,T}
    c = prefactor(x)
    c = first(r) == firstindex(factors) ? c : one(c)
    isconcretetype(T) || return _narrowed_part(C, c, factors, r, factors[first(r)])
    x.factors isa Vector{T} && length(r) == length(factors) && return x
    return NCMul{C}(c, collect(view(factors, r)))
end
# A function barrier, so that the constructor of the product is static
function _narrowed_part(::Type{C}, c, factors, r, f1::U) where {C,U}
    out = _narrowed_copy(factors, r, f1)
    return out isa Vector{U} ? NCMul{C,U,Vector{U}}(c, out) : NCMul{C}(c, out)
end
