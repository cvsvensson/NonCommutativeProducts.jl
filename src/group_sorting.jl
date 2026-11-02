"""
    group_commutation(g1, g2)

How the factors of the distinct [`atom_group`](@ref)s `g1` and `g2` commute, which lets `sort` order them by group
without calling [`mul_effect`](@ref) on each pair of their factors:
- `:commute`: `a*b == b*a` for all atoms `a` in `g1` and `b` in `g2`;
- `:graded`: `a*b == (-1)^(p(a)*p(b)) * b*a`, with `p` the [`atom_parity`](@ref) of an atom;
- `nothing` (the default): unknown, and products with factors of both groups are sorted by `mul_effect` alone.

It must be symmetric, and is called with the groups in either order. A declaration is a promise about `mul_effect`,
which must still be defined for these pairs: for `a` in `g1` and `b` in `g2`, either `mul_effect(a, b)` is `nothing`
for all such pairs, or it is the swap (`Swap(1)`, or `Swap(-1)` for two odd atoms of graded groups) for all of them.
So the order of the two groups is the same for all their atoms, as with [`@commutative`](@ref). Sorting then gives
the same result as without the declaration, including the order of the factors, and only gets faster; see
[`group_isless`](@ref) for the order of groups.

The declarations are used to sort products whose factors have an abstract element type, have at least 16 factors, and
whose groups all have a declared commutation with each other; other products are sorted by `mul_effect` alone. The
results for groups that are types are cached until methods are next defined, so they must depend only on the groups.

# Example
```julia
# Fermions of two species anticommute, and both commute with bosons
NonCommutativeProducts.group_commutation(::Type{<:Union{FermionA,FermionB}}, ::Type{<:Union{FermionA,FermionB}}) = :graded
NonCommutativeProducts.group_commutation(::Type{Boson}, ::Type{<:Union{FermionA,FermionB}}) = :commute
NonCommutativeProducts.group_commutation(::Type{<:Union{FermionA,FermionB}}, ::Type{Boson}) = :commute
NonCommutativeProducts.atom_parity(::Union{FermionA,FermionB}) = true
```
"""
group_commutation(g1, g2) = nothing

"""
    group_isless(g1, g2)

Whether the factors of the group `g1` come before those of the group `g2` in a sorted product, for groups whose
[`group_commutation`](@ref) is not `nothing`. Returns `true`, `false`, or `nothing` (the default) to take the order
from `mul_effect`: `g1` comes first if `mul_effect(a, b)` is `nothing` for an atom `a` of `g1` and an atom `b` of `g2`.

Define it to avoid these calls. It must give the same order as `mul_effect`, which by the promise of
`group_commutation` depends only on the groups, so the default never reorders products differently from `mul_effect`.
It must be a strict total order on the groups that occur in one product.

# Example
```julia
NonCommutativeProducts.group_isless(::Type{FermionA}, ::Type{FermionB}) = true
NonCommutativeProducts.group_isless(::Type{FermionB}, ::Type{FermionA}) = false
```
"""
group_isless(g1, g2) = nothing

"""
    atom_parity(atom)

The parity of `atom`, `true` for odd and `false` for even, used for groups whose [`group_commutation`](@ref) is
`:graded`: swapping two odd atoms of such groups gives a factor `-1`. It is called only for atoms of these groups, and
has no default. Parity is separate from the [`grade`](@ref) of an atom, which can be any grading, e.g. a charge; for
the charge of fermions it is `isodd(grade(atom))`.
"""
atom_parity(atom) = throw(ArgumentError("atom_parity is not defined for $(typeof(atom)), whose group is declared \
    `:graded` by group_commutation; define `NonCommutativeProducts.atom_parity(::$(typeof(atom)))`"))

# For factors of type Any, the hooks are called behind an inference barrier, as mul_effect is (see _mul_effect): the
# default methods would otherwise be inlined into the precompiled sorting code, and downstream methods would
# invalidate it
_atom_group(::Vector{Any}, x) = Base.inferencebarrier(atom_group)(x)
_atom_group(factors, x) = atom_group(x)
_group_commutation(g1, g2) = _check_commutation(Base.inferencebarrier(group_commutation)(g1, g2), g1, g2)
_group_isless(g1, g2) = Base.inferencebarrier(group_isless)(g1, g2)
_atom_parity(::Vector{Any}, x) = Base.inferencebarrier(atom_parity)(x)::Bool
_atom_parity(factors, x) = atom_parity(x)::Bool
function _check_commutation(c, g1, g2)
    (isnothing(c) || c === :commute || c === :graded) && return c
    throw(ArgumentError("group_commutation($(repr(g1)), $(repr(g2))) returned $(repr(c)), but it must return \
        :commute, :graded or nothing"))
end

# Sorts `coeff * factors` into `acc` as bubble!! does, using `factors` as scratch space, with the same result. When the
# group hooks allow it, the factors are first moved to the order of their groups. Each group is then sorted in place,
# on a copy with a concrete element type when its factors share one, as long as its effects are swaps. From the first
# other effect on, the product is sorted by bubble!!, which then finds the factors before that effect sorted.
#
# Only products of at least `_group_sort_min_length` factors with an abstract element type take this path: for a
# concrete one, mul_effect is already resolved at compile time, and for fewer factors, finding the groups and their
# order takes about as long as it saves (with three groups of fermions, sorting 24 factors is about 6 times faster,
# while products of 12 factors whose groups give AddTerms are slower).
const _group_sort_min_length = 16
_sort_by_groups!!(acc, coeff, factors) = bubble!!(acc, coeff, factors, nothing, 1)
function _sort_by_groups!!(acc, coeff, factors::Vector{T}) where {T}
    (isconcretetype(T) || length(factors) < _group_sort_min_length) && return bubble!!(acc, coeff, factors, nothing, 1)
    _first_pair_commutes(factors) || return bubble!!(acc, coeff, factors, nothing, 1)
    plan = _group_plan(factors)
    isnothing(plan) && return bubble!!(acc, coeff, factors, nothing, 1)
    groups, ids, rank, order, graded = plan
    # A single group: sort a copy with a concrete element type, behind a function barrier
    length(groups) == 1 && return bubble!!(acc, coeff, _group_segment(factors, eachindex(factors), factors[1]), nothing, 1)
    _crossing_sign(factors, ids, rank, graded) && (coeff = -coeff)
    sorted, ranges = _permute_by_group(factors, ids, order)
    start = 0 # where bubble!! continues, or 0 if all groups are sorted
    for r in ranges
        length(r) == 1 && continue
        λ, i = _swap_sort_group!(sorted, r, sorted[first(r)])
        coeff *= λ
        iszero(coeff) && return acc
        (i == 0 || start != 0) || (start = first(r) + i - 1)
    end
    start == 0 && return add_sorted!!(acc, coeff, sorted)
    return bubble!!(acc, coeff, sorted, nothing, start)
end

# Whether the first two groups of `factors` have a declared commutation, or there is one group. A check without
# allocations, so that products of undeclared groups fall back to bubble!! quickly.
function _first_pair_commutes(factors)
    g = _atom_group(factors, factors[1])
    for i in 2:length(factors)
        h = _atom_group(factors, factors[i])
        _same_group(g, h) || return !isnothing(first(_group_hooks(g, h)))
    end
    return true
end

# The factors in the order of `order`, the group indices from first to last, keeping the order within each group, and
# the range of each group in the result
function _permute_by_group(factors, ids, order)
    k = length(order)
    counts = zeros(Int, k)
    foreach(id -> counts[id] += 1, ids)
    starts = Vector{Int}(undef, k) # the next position of each group
    pos = 1
    for id in order
        starts[id] = pos
        pos += counts[id]
    end
    ranges = [starts[id]:starts[id]+counts[id]-1 for id in order]
    sorted = similar(factors)
    for (x, id) in zip(factors, ids)
        sorted[starts[id]] = x
        starts[id] += 1
    end
    return sorted, ranges
end

# The distinct groups of `factors`, the group index of each factor, the rank of each group in the sorted order, the
# groups in that order, and the pairs of groups that are graded. `nothing` if some pair of groups has unknown commutation.
function _group_plan(factors)
    groups = Any[]
    ids = Vector{Int}(undef, length(factors))
    reps = Int[] # the index of the first factor of each group
    graded = Tuple{Int,Int}[]
    current = 0
    for (i, x) in enumerate(factors)
        g = _atom_group(factors, x)
        if current == 0 || !_same_group(g, groups[current])
            k = _find_group(groups, g)
            if k == 0
                for (j, h) in enumerate(groups)
                    c = first(_group_hooks(h, g))
                    isnothing(c) && return nothing
                    c === :graded && push!(graded, (j, length(groups) + 1))
                end
                push!(groups, g)
                push!(reps, i)
                k = length(groups)
            end
            current = k
        end
        ids[i] = current
    end
    # Insertion sort of the few groups by the order of group_isless or mul_effect
    order = collect(1:length(groups))
    for m in 2:length(order)
        for n in m:-1:2
            a, b = order[n-1], order[n]
            _group_before(groups[b], groups[a], factors[reps[b]], factors[reps[a]], factors) || break
            order[n-1], order[n] = b, a
        end
    end
    return groups, ids, invperm(order), order, graded
end
function _find_group(groups, g)
    for (k, h) in enumerate(groups)
        _same_group(h, g) && return k
    end
    return 0
end
# isequal(g, h), which for two types is a slow type comparison. DataTypes are unique, so they are equal only if they
# are identical.
@inline _same_group(g, h) = g === h || (!(g isa DataType && h isa DataType) && isequal(g, h))
# Whether the group g1 of the atom `a` comes before the group g2 of the atom `b`
function _group_before(g1, g2, a, b, factors)
    before = last(_group_hooks(g1, g2))
    isnothing(before) || return before::Bool
    return isnothing(_mul_effect(factors, a, b))
end

# group_commutation(g1, g2) and group_isless(g1, g2). Hooks that dispatch on types are slow to call at runtime, since
# their methods are found by subtyping, so their results for pairs of types are cached, in a short list that is
# searched by identity: hashing types is slow too. The list is never mutated, but replaced by a longer one, so that it
# is read without a lock. It is cleared whenever methods are defined, since the hooks may then have changed. Other
# groups are not cached, since they can be many.
struct _HookCache
    world::UInt
    entries::Vector{Tuple{DataType,DataType,Union{Nothing,Symbol},Union{Nothing,Bool}}}
end
mutable struct _HookCacheRef
    @atomic cache::_HookCache
end
const _hook_cache = _HookCacheRef(_HookCache(typemax(UInt), []))
const _hook_cache_lock = ReentrantLock()
const _hook_cache_size = 256
_call_group_hooks(g1, g2) = (_group_commutation(g1, g2), _group_isless(g1, g2)::Union{Nothing,Bool})
# A single method, so that calls with groups of unknown type are not dispatched at runtime: dispatch on the value of a
# type is slow
function _group_hooks(g1, g2)
    (g1 isa DataType && g2 isa DataType) || return _call_group_hooks(g1, g2)
    cache = @atomic _hook_cache.cache
    world = Base.get_world_counter()
    if cache.world == world
        for (h1, h2, c, before) in cache.entries
            (h1 === g1 && h2 === g2) && return c, before
        end
    end
    hooks = _call_group_hooks(g1, g2)
    @lock _hook_cache_lock begin
        cache = @atomic _hook_cache.cache
        entries = cache.world == world ? cache.entries : eltype(cache.entries)[]
        if length(entries) < _hook_cache_size
            @atomic _hook_cache.cache = _HookCache(world, [entries; (g1, g2, hooks...)])
        end
    end
    return hooks
end

# Whether moving the factors to the order of their groups crosses an odd number of pairs of odd atoms in graded groups.
# The relative order of factors in one group is kept.
function _crossing_sign(factors, ids, rank, graded)
    isempty(graded) && return false
    k = length(rank)
    isgraded = falses(k, k)
    for (i, j) in graded
        isgraded[i, j] = isgraded[j, i] = true
    end
    hasparity = [any(@view isgraded[:, j]) for j in 1:k]
    odd = [hasparity[ids[i]] && _atom_parity(factors, factors[i]) for i in eachindex(factors)]
    oddcount = zeros(Int, k) # the number of odd factors so far in each group
    negate = false
    for i in eachindex(factors)
        odd[i] || continue
        g = ids[i]
        for h in 1:k
            rank[h] > rank[g] && isgraded[h, g] && isodd(oddcount[h]) && (negate = !negate)
        end
        oddcount[g] += 1
    end
    return negate
end

# A copy of the factors at `idxs`, with the type U of `x` as element type if they all have that type
function _group_segment(factors, idxs, ::U) where {U}
    all(i -> factors[i] isa U, idxs) || return factors[idxs]
    return U[factors[i] for i in idxs]
end

# Sorts sorted[r] in place by _swap_sort!, on a copy with the type of its first factor `x` as element type if all its
# factors have that type, and returns the result of _swap_sort!. A function barrier, so that sorting is specialized on
# the type of the factors.
function _swap_sort_group!(sorted, r, x)
    seg = _group_segment(sorted, r, x)
    λ, i = _swap_sort!(seg)
    copyto!(sorted, first(r), seg, 1, length(r))
    return λ, i
end
# The steps of bubble!! on `seg` while the effects are `nothing` or a Swap with an Int factor, which are applied in
# place. Returns the product of the factors, and 0 if `seg` is then sorted, or else the position of the first other
# effect. The factors before that position are sorted.
function _swap_sort!(seg)
    λ, i = 1, 1
    while i < length(seg)
        effect = _mul_effect(seg, seg[i], seg[i+1])
        isnothing(effect) && (i += 1; continue)
        effect isa Swap{Int} || return λ, i
        seg[i], seg[i+1] = seg[i+1], seg[i]
        λ *= effect.λ
        iszero(λ) && return λ, 0
        i = max(i - 1, 1)
    end
    return λ, 0
end
