"""
    AddTerms(terms)

Return value for [`mul_effect`](@ref) that replaces `a*b` by the sum of `terms`. `terms` is an iterable such as a
vector or a tuple, and its elements can be any other allowed return value of `mul_effect` except `nothing`, including
sums and nested `AddTerms`.
"""
struct AddTerms{T}
    terms::T
end
"""
    Swap(λ)

Return value for [`mul_effect`](@ref) that replaces `a*b` by `λ*b*a`.
"""
struct Swap{T}
    λ::T
end

function splice!!(ncmul::NCMul, i::Integer, swap::Swap)
    factors = ncmul.factors
    a = factors[i]
    factors[i] = factors[i+1]
    factors[i+1] = a
    coeff = prefactor(ncmul) * swap.λ
    return NCMul(coeff, factors)
end
function splice!!_and_add(ncmul::NCMul, i, terms::AddTerms)
    # resolved at compile time for tuples, which keeps the common case type stable
    any(term -> term isa Union{NCAdd,AddTerms}, terms.terms) && return _splice!!_and_add_flatten(ncmul, i, terms.terms)
    firstterm, state = iterate(terms.terms)
    newterms = map(term -> splice!!(copy(ncmul), i, term), Base.rest(terms.terms, state))
    ncmul2 = splice!!(ncmul, i, firstterm)
    return ncmul2, newterms
end
# Elements that are sums expand to several terms. The first element is spliced last, since that mutates ncmul.
function _splice!!_and_add_flatten(ncmul::NCMul, i, terms)
    firstterm, state = iterate(terms)
    newterms = Any[]
    for term in Base.rest(terms, state)
        ncmul2, newterms2 = splice!!_and_add(copy(ncmul), i, term)
        push!(newterms, ncmul2)
        append!(newterms, newterms2)
    end
    ncmul2, newterms2 = splice!!_and_add(ncmul, i, firstterm)
    append!(newterms, newterms2)
    # narrow the element type, since the sorting loop needs a Vector of NCMuls
    return ncmul2, identity.(newterms)
end

"""
    sort(x::Union{NCMul,NCAdd})

Sort `x` by repeatedly applying [`mul_effect`](@ref) to neighbouring factors. The result is always an `NCAdd`, also
when `x` is a product that sorts to a single term.
"""
Base.sort(a::NCMul) = bubble_sort(a) |> _as_ncadd
Base.sort(a::NCAdd) = bubble_sort(a)
_as_ncadd(a::NCAdd) = a
_as_ncadd(a::NCMul) = NCAdd(a)
function bubble_sort(a::NCMul)
    return _with_zero_grades(bubble_sort!(copy(a)), () -> _grades(a))
end
# Sorts `a` in place, so its factors must not be shared.
function bubble_sort!(a::NCMul{C}) where C
    if length(a.factors) <= 1
        return a
    end
    return _bubble_sort!([NCMul(prefactor(a), factors_vector(a.factors))], C)
end
function bubble_sort(ncadd::NCAdd{C}) where C
    length(ncadd.dict) == 0 && return ncadd
    terms = collect(NCterms(ncadd, copy_factors))
    sorted = _add_scaled!!(_bubble_sort!(terms, C), additive_coeff(ncadd), One())
    return _with_zero_grades(sorted, () -> _grades(ncadd))
end
function _bubble_sort!(terms::Vector{T}, ::Type{C}=Int) where {T<:NCMul,C}
    sorted_terms = with(_autosort => false) do
        __bubble_sort!(terms)
    end
    # function barrier: the element type of sorted_terms is only known at runtime
    return _sum_sorted_terms(sorted_terms, C)
end
function _sum_sorted_terms(sorted_terms::Vector{T}, ::Type{C}) where {T<:NCMul,C}
    if length(sorted_terms) == 0
        return NCAdd(zero(C), Dict{to_add_dict_type(T),C}())
    end
    newadd = NCAdd(zero(C), to_add_dict(first(sorted_terms)))
    for term in Iterators.drop(sorted_terms, 1)
        newadd = _add_scaled!!(newadd, term, One())
    end
    sizehint!(newadd.dict, length(newadd.dict))
    return newadd
end

# Sorts terms[n:end]. The optional arguments let the recursive call resume: `start` is the factor position in
# terms[n] and `done` whether terms[n] is sorted.
function __bubble_sort!(terms::Vector{T}, n::Int=1, start::Int=1, done::Bool=false) where {T<:NCMul}
    while n <= length(terms)
        while !done && n <= length(terms)
            newterms, done, start = __bubble_sort_step!(terms, n, start)
            # new terms can widen the vector's element type; continue in a call specialized on the new type
            newterms isa Vector{T} || return __bubble_sort!(newterms, n, start, done)
            terms = newterms
        end
        n += 1
        done = false
        start = 1
    end
    return terms
end

function __bubble_sort_step!(terms, index, start)
    ncmul = terms[index]
    factors = ncmul.factors
    N::Int = length(factors)
    if N <= 1
        done = true
        return terms, done, start
    end
    i::Int = max(0, start - 1)
    while i < N - 1
        i += 1
        a, b = factors[i], factors[i+1]
        effect = _mul_effect(factors, a, b)
        isnothing(effect) && continue
        return _apply_effect!!(terms, index, ncmul, i, effect)
    end
    done = true
    newstart = i - 1
    return terms, done, newstart
end
# Replaces the factors i and i+1 of `ncmul == terms[index]` according to `effect`, adding any new terms to `terms`
function _apply_effect!!(terms, index, ncmul, i, effect)
    # When effect is only known at runtime, these results infer as Any, and code compiled for Any is invalidated by
    # downstream convert methods. The type assertions here and in _add_newterms!! prevent that.
    newncmul, newterms = splice!!_and_add(ncmul, i, effect)
    newncmul::NCMul

    terms_with_newterms = _add_newterms!!(terms, newterms)
    if iszero(newncmul)
        deleteat!(terms_with_newterms, index)
        return terms_with_newterms, false, 1
    end
    return setindex!!(terms_with_newterms, newncmul, index), false, i - 1
end
# Mixed-type products dispatch mul_effect at runtime anyway. The inference barrier keeps the compiled sorting loop
# independent of the mul_effect methods, so downstream methods don't invalidate the precompiled code.
_mul_effect(::Vector{Any}, a, b) = Base.inferencebarrier(mul_effect)(a, b)
_mul_effect(factors, a, b) = mul_effect(a, b)

function _add_newterms!!(terms::Vector{T}, newterms) where {T}
    isempty(newterms) && return terms
    newterms isa Union{Vector{T},Tuple{Vararg{T}}} && return append!(terms, newterms)
    # the new terms widen the element type. The inference barrier keeps append!! from being compiled for Any, which
    # downstream convert methods would invalidate.
    return Base.inferencebarrier(append!!)(terms, newterms)::Vector
end


function mysplice!!(v::V, i::Integer, replacement::W) where {V,W}
    T = promote_type(eltype(W), eltype(V))
    if T <: eltype(V)
        _mysplice!(v, i, replacement)
        return v
    end
    # a new vector of the n factors before v[i:i+1], the replacement, and the factors after v[i:i+1]
    n = i - firstindex(v)
    m = length(replacement)
    nafter = lastindex(v) - i - 1
    newv = Vector{T}(undef, n + m + nafter)
    copyto!(newv, 1, v, firstindex(v), n)
    for (k, x) in enumerate(replacement)
        newv[n+k] = x
    end
    copyto!(newv, n + m + 1, v, i + 2, nafter)
    return newv
end

function _mysplice!(a::Vector, i::Integer, ins)
    m = length(ins)
    if m == 0
        deleteat!(a, i:(i+1))
        return nothing
    end
    left = 2i < length(a)
    if m == 1
        deleteat!(a, left ? i : i + 1)
    elseif m > 2
        Base._growat!(a, left ? i : i + 2, m - 2)
    end
    # m == 2 needs no resize
    for (k, x) in enumerate(ins)
        a[i+k-1] = x
    end
    return nothing
end


function splice!!(ncmul::NCMul, i, term::NCMul)
    coeff = prefactor(term) * prefactor(ncmul)
    newfactors = mysplice!!(ncmul.factors, i, term.factors)
    return NCMul(coeff, newfactors)
end
# an atom of the factor type. For mixed products (S == Any) this matches any value, so non-factors are rejected here.
function splice!!(ncmul::NCMul{C,S}, i, term::S) where {C,S}
    _check_atom_effect(term)
    splice!!(ncmul, i, NCMul(one(C), (term,)))
end
# an atom of another type than the factors. NCMul(term) is defined by @nc, so only registered atoms are accepted.
function splice!!(ncmul::NCMul, i, term)
    _check_atom_effect(term)
    splice!!(ncmul, i, NCMul(term)::NCMul)
end
_check_atom_effect(term) = nothing
# `nothing` inside AddTerms would keep a*b as a term, which gives the same effect again and never terminates
_check_atom_effect(::Nothing) = throw(ArgumentError("unsupported mul_effect value: `nothing` can only be returned as the whole result of mul_effect, not inside AddTerms."))
function splice!!(ncmul::NCMul, i::Integer, coeff::Number)
    deleteat!(ncmul.factors, i:(i+1))
    return NCMul(coeff * prefactor(ncmul), ncmul.factors)
end
function splice!!_and_add(ncmul::NCMul, i, add::NCAdd)
    if iszero(additive_coeff(add))
        terms = NCterms(add)
        firstterm, state = iterate(terms)
        newterms = map(term -> splice!!(copy(ncmul), i, term), Base.rest(terms, state))
        ncmul2 = splice!!(ncmul, i, firstterm)
    else
        newterms = map(term -> splice!!(copy(ncmul), i, term), NCterms(add))
        ncmul2 = splice!!(ncmul, i, additive_coeff(add))
    end
    return ncmul2, newterms
end
function splice!!_and_add(ncmul::T, i, effect) where T<:NCMul
    ncmul2 = splice!!(ncmul, i, effect)
    return ncmul2, T[]
end

"""
    mul_effect(a, b)

When sorting a non-commutative expression, `mul_effect(a, b)` gives the replacement for `a*b`. 
It must have a branch that returns `nothing` in order for the sorting to terminate. 

This package defines no methods for `mul_effect`, so define them for your types.

The return value can be `nothing`, a number, an object, a sum or product of these, or one of the helper types `Swap(λ)`
(replace `a*b` by `λ*b*a`) and `AddTerms(terms)` (replace `a*b` by a sum of `terms`). `Swap` and `AddTerms` are
available as `NonCommutativeProducts.Swap` and `NonCommutativeProducts.AddTerms`.

# Example
```julia
# Anticommuting Grassmann numbers
struct θ id::Int end
NonCommutativeProducts.@nc θ
NonCommutativeProducts.mul_effect(a::θ, b::θ) = a.id < b.id ? nothing : -b * a * (a.id==b.id ? 0 : 1)
NonCommutativeProducts.enable_autosort!()
θ(2)*θ(1) # = -θ(1)*θ(2)
θ(1)*θ(1) # = 0
```
"""
function mul_effect end
