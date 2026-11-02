"""
    group_commutation(g1, g2)

How the factors of the distinct [`atom_group`](@ref)s `g1` and `g2` commute, which lets `sort` order them by group
without calling [`mul_effect`](@ref) on each pair of their factors:
- `:commute`: `a*b == b*a` for all atoms `a` in `g1` and `b` in `g2`;
- `:graded`: `a*b == (-1)^(p(a)*p(b)) * b*a`, with `p` the [`atom_parity`](@ref) of an atom;
- `nothing` (the default): unknown, and products with factors of both groups are sorted by `mul_effect` alone.

A declaration is a promise that `mul_effect` orders the two groups by group alone: for all `a` in `g1` and `b` in
`g2`, either `mul_effect(a, b)` is `nothing` and `mul_effect(b, a)` is the swap, or the other way around, where the
swap is `Swap(-1)` for two odd atoms of `:graded` groups and `Swap(1)` otherwise, as with [`@commutative`](@ref). The
order of the groups is taken from `mul_effect` on one atom of each, and the groups of a product must then be in a
strict total order. Sorting gives the same result as without the declaration, including the order of the factors, if
`mul_effect` gives unique sorted forms, and only gets faster. See [Commuting groups](@ref) for what is checked.
"""
group_commutation(g1, g2) = nothing

"""
    atom_parity(atom)

The parity of `atom`, `true` for odd and `false` for even, used for groups whose [`group_commutation`](@ref) is
`:graded`: swapping two odd atoms of such groups gives a factor `-1`. It is called only for atoms of these groups, and
has no default. Parity is separate from the [`grade`](@ref) of an atom, which can be any grading, e.g. a charge; for
the charge of fermions it is `isodd(grade(atom))`.
"""
atom_parity(atom) = throw(ArgumentError("atom_parity is not defined for $(typeof(atom)), whose group is declared \
    `:graded` by group_commutation; define `NonCommutativeProducts.atom_parity(::$(typeof(atom)))`"))

# The hooks are called behind an inference barrier, as mul_effect is (see _mul_effect), so that downstream methods
# don't invalidate the precompiled sorting code
_atom_group(::Vector{Any}, x) = Base.inferencebarrier(atom_group)(x)
_atom_group(factors, x) = atom_group(x)
_atom_parity(::Vector{Any}, x) = Base.inferencebarrier(atom_parity)(x)::Bool
_atom_parity(factors, x) = atom_parity(x)::Bool
function _group_commutation(g1, g2)
    c = Base.inferencebarrier(group_commutation)(g1, g2)
    (isnothing(c) || c === :commute || c === :graded) || throw(ArgumentError("group_commutation($(repr(g1)), \
        $(repr(g2))) returned $(repr(c)), but it must return :commute, :graded or nothing"))
    c2 = Base.inferencebarrier(group_commutation)(g2, g1)
    c2 === c || throw(ArgumentError("group_commutation is not symmetric: it returned $(repr(c)) for \
        ($(repr(g1)), $(repr(g2))) but $(repr(c2)) for ($(repr(g2)), $(repr(g1))); define it for both orders"))
    return c
end

# Sorts `coeff * factors` into `acc` as bubble!! does, using `factors` as scratch space, with the same result. When all
# groups of the factors have a declared commutation, the factors are first moved to the order of their groups, keeping
# the order within each group. Each group is then sorted in place, on a copy with a concrete element type when its
# factors share one, as long as its effects are swaps. From the first other effect on, bubble!! sorts the product.
#
# Only factors of type Any or of an abstract type that is not a Union take this path: for others, mul_effect is resolved
# at compile time. Nor do products shorter than `_group_sort_min_length`: finding their groups and their order takes
# about as long as it saves (for three groups of fermions and spins, with eltype Any and abstract: 0.85x and 0.77x as
# fast for 8 factors, 1.4x and 1.0x for 12, 1.7x and 1.1x for 16, 2.7x and 1.9x for 24).
const _group_sort_min_length = 12
_sort_by_groups!!(acc, coeff, factors) = bubble!!(acc, coeff, factors, nothing, 1)
function _sort_by_groups!!(acc, coeff, factors::Vector{T}) where {T}
    (T isa Union || isconcretetype(T) || length(factors) < _group_sort_min_length) &&
        return bubble!!(acc, coeff, factors, nothing, 1)
    return _group_sort!!(acc, coeff, factors)
end
function _group_sort!!(acc, coeff, factors)
    plan = _group_plan(factors)
    isnothing(plan) && return bubble!!(acc, coeff, factors, nothing, 1)
    _crossing_sign(factors, plan) && (coeff = -coeff)
    sorted, ranges = _permute_by_group(factors, plan)
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

# The groups of `factors`, numbered by first appearance: the group of each factor (`ids`), the groups in sorted order
# (`order`), and the pairs (g, h) of graded groups where g comes after h (`flips`). `nothing` if there is one group, or
# a pair of groups with unknown commutation. The order of each pair of groups is found by mul_effect on their first
# atoms, which is checked against the declared commutation.
function _group_plan(factors)
    groups, reps = Any[], Int[] # the groups, and the index of the first factor of each
    nbefore = Int[] # the number of groups that come before each group
    flips = Tuple{Int,Int}[]
    ids = Vector{Int}(undef, length(factors))
    current = 0
    for (i, x) in enumerate(factors)
        g = _atom_group(factors, x)
        if current == 0 || !_same_group(g, groups[current])
            current = 0
            for (j, h) in enumerate(groups)
                _same_group(h, g) && (current = j; break)
            end
            if current == 0
                push!(nbefore, 0)
                for (j, h) in enumerate(groups)
                    c = _group_commutation(h, g)
                    isnothing(c) && return nothing
                    swapped = _swaps_groups(factors, factors[reps[j]], x, c, h, g)
                    nbefore[swapped ? j : end] += 1
                    c === :graded && push!(flips, swapped ? (j, length(groups) + 1) : (length(groups) + 1, j))
                end
                push!(groups, g)
                push!(reps, i)
                current = length(groups)
            end
        end
        ids[i] = current
    end
    length(groups) == 1 && return nothing
    # The order is transitive if and only if the groups have distinct numbers of groups before them
    allunique(nbefore) || throw(ArgumentError("mul_effect does not order the groups $(join(repr.(groups), ", ")) \
        consistently: define it so that the order of these groups is transitive"))
    return (; ids, order = invperm(nbefore .+ 1), flips)
end
# isequal(g, h), which for two types is a slow type comparison. DataTypes are unique, so they are equal only if they
# are identical.
@inline _same_group(g, h) = g === h || (!(g isa DataType && h isa DataType) && isequal(g, h))
# Whether the group g of `a` comes after the group h of `b`, from mul_effect(a, b), which must be `nothing` or the swap
# of commutation `c`, and if it is `nothing`, from mul_effect(b, a), which must be that swap. The order of the groups
# then doesn't depend on which comes first in the product.
function _swaps_groups(factors, a, b, c, g, h)
    λ = c === :graded && _atom_parity(factors, a) && _atom_parity(factors, b) ? -1 : 1
    effect = _mul_effect(factors, a, b)
    if !isnothing(effect)
        effect isa Swap && effect.λ == λ && return true
        throw(_swap_error(c, g, h, a, b, λ, effect))
    end
    effect = _mul_effect(factors, b, a)
    effect isa Swap && effect.λ == λ && return false
    isnothing(effect) && throw(ArgumentError("groups $(repr(g)) and $(repr(h)) are declared to commute but \
        mul_effect does not order them; define mul_effect to return Swap($λ) in one order"))
    throw(_swap_error(c, h, g, b, a, λ, effect))
end
_swap_error(c, g, h, a, b, λ, effect) = ArgumentError("group_commutation($(repr(g)), $(repr(h))) is $(repr(c)), so \
    mul_effect($(repr(a)), $(repr(b))) must be nothing or Swap($λ), but it is $(repr(effect)); define mul_effect to \
    order these groups, or don't declare their commutation")

# Whether moving the factors to the order of their groups crosses an odd number of pairs of odd atoms in graded groups
function _crossing_sign(factors, (; ids, flips))
    isempty(flips) && return false
    odd = [any(p -> id in p, flips) && _atom_parity(factors, x) for (x, id) in zip(factors, ids)]
    crossings = 0
    for (g, h) in flips # count the odd atoms of g before odd atoms of h
        seen = 0
        for (o, id) in zip(odd, ids)
            o || continue
            id == g && (seen += 1)
            id == h && (crossings += seen)
        end
    end
    return isodd(crossings)
end

# The factors in the order of their groups, keeping the order within each group, and the range of each group in it
function _permute_by_group(factors, (; ids, order))
    counts = zeros(Int, length(order))
    foreach(id -> counts[id] += 1, ids)
    starts = similar(counts) # the next position of each group
    pos = 1
    for g in order
        starts[g] = pos
        pos += counts[g]
    end
    ranges = [starts[g]:starts[g]+counts[g]-1 for g in order]
    sorted = similar(factors)
    for (x, id) in zip(factors, ids)
        sorted[starts[id]] = x
        starts[id] += 1
    end
    return sorted, ranges
end

# Sorts sorted[r] in place by _swap_sort!, on a copy with the type U of its first factor as element type if all its
# factors have that type, and returns the result of _swap_sort!. A function barrier, so that sorting is specialized on
# the type of the factors.
function _swap_sort_group!(sorted, r, ::U) where {U}
    seg = all(i -> sorted[i] isa U, r) ? U[sorted[i] for i in r] : sorted[r]
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
