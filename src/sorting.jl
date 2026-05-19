struct AddTerms{T}
    terms::T
end
struct Swap{T}
    λ::T
end

function splice!!(ncmul::NCMul, i::UnitRange, swap::Swap)
    i1, i2 = first(i), last(i)
    i2 == i1 + 1 && length(i) == 2 || throw(ArgumentError("Invalid index for swap"))
    factors = ncmul.factors
    a = factors[i1]
    b = factors[i2]
    factors[i1] = b
    factors[i2] = a
    coeff = prefactor(ncmul) * swap.λ
    return NCMul(coeff, factors)
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
    firstterm, state = iterate(terms.terms)
    newterms = map(term -> splice!!(copy(ncmul), i, term), Base.rest(terms.terms, state))
    ncmul2 = splice!!(ncmul, i, firstterm)
    return ncmul2, newterms
end

Base.sort(a::NCMul) = bubble_sort(a)
Base.sort(a::NCAdd) = bubble_sort(a)
Base.sort!(a::NCMul) = bubble_sort!(a)
Base.sort!(a::NCAdd) = bubble_sort(a)
function bubble_sort(a::NCMul)
    return bubble_sort!(copy(a))
end
function bubble_sort!(a::NCMul)
    if length(a.factors) <= 1
        return a
    end
    res = _bubble_sort!([a])
    return res
end
function bubble_sort(ncadd::NCAdd)
    terms = collect(NCMul(v, copy(k.factors)) for (k, v) in pairs(ncadd.dict))
    res = add!!(_bubble_sort!(terms), additive_coeff(ncadd))
end
function _bubble_sort!(terms::Vector{T}) where {T<:NCMul}
    sorted_terms = with(_autosort => false) do
        __bubble_sort!(terms)
    end
    if length(sorted_terms) == 0
        return NCAdd(0, Dict{T,Int}())
    end
    newadd = NCAdd(0, to_add_dict(first(sorted_terms)))
    for term in Iterators.drop(sorted_terms, 1)
        newadd = add!!(newadd, term)
    end
    return filter_scalars!(filter_zeros!(newadd))
end

function __bubble_sort!(terms::Vector{T}) where {T<:NCMul}
    n = 1
    while n <= length(terms)
        done = false
        start = 1
        while !done && n <= length(terms)
            terms, done, start = __bubble_sort!(terms, n, start)
        end
        n += 1
    end
    return terms
end

function __bubble_sort!(terms, index, start)
    no_effect = true
    ncmul = terms[index]
    factors = ncmul.factors
    N::Int = length(factors)
    if N <= 1
        done = true
        return terms, done, start
    end
    i::Int = max(0, start - 1)
    while no_effect && i < N - 1
        i += 1
        a, b = factors[i], factors[i+1]
        effect = mul_effect(a, b)
        isnothing(effect) && continue

        no_effect = false
        newncmul, newterms = splice!!_and_add(ncmul, i, effect)

        terms_with_newterms = _add_newterms!!(terms, newterms)
        if iszero(prefactor(newncmul))
            deleteat!(terms_with_newterms, index)
            return terms_with_newterms, false, 1
        end
        return setindex!!(terms_with_newterms, newncmul, index), false, i - 1
    end
    done = true #no_effect
    newstart = i - 1
    return terms, done, newstart
end
function _add_newterms!!(terms, newterms)
    Nnew = length(newterms)::Int
    if Nnew > 0
        return append!!(terms, newterms)
    else
        return terms
    end
end


function mysplice!!(v::V, i::UnitRange, replacement::W) where {V,W}
    T = promote_type(eltype(W), eltype(V))
    if T <: eltype(V)
        _mysplice!(v, i, replacement)
        return v
    else
        return Vector{T}(vcat(v[1:first(i)-1], replacement, v[last(i)+1:end]))
    end
end
function mysplice!!(v::V, i::Integer, replacement::W) where {V,W}
    T = promote_type(eltype(W), eltype(V))
    if T <: eltype(V)
        _mysplice!(v, i, replacement)
        return v
    else
        return Vector{T}(vcat(v[1:i-1], replacement, v[i+2:end]))
    end
end

function _mysplice!(a::Vector, r::AbstractUnitRange{<:Integer}, ins)
    # this is a copy of the splice! in base, except that it does not return the deleted elements,
    # in order to avoid allocations
    m = length(ins)
    if m == 0
        deleteat!(a, r)
        return nothing
    end

    n = length(a)
    f = first(r)
    l = last(r)
    d = length(r)

    if m < d
        delta = d - m
        Base._deleteat!(a, (f - 1 < n - l) ? f : (l - delta + 1), delta)
    elseif m > d
        Base._growat!(a, (f - 1 < n - l) ? f : (l + 1), m - d)
    end

    k = 1
    for x in ins
        a[f+k-1] = x
        k += 1
    end
    return nothing
end
function _mysplice!(a::Vector, i::Integer, ins)
    m = length(ins)
    if m == 0
        Base._deleteat!(a, i, 2)
        return nothing
    end
    left = 2i < length(a)
    if m == 1
        Base._deleteat!(a, left ? i : i + 1, 1)
    elseif m > 2
        Base._growat!(a, left ? i : i + 2, m - 2)
    end
    # m == 2: no resize, just overwrite
    for (k, x) in enumerate(ins)
        a[i+k-1] = x
    end
    return nothing
end


function splice!!(ncmul::NCMul, i::UnitRange)
    splice!(ncmul.factors, i)
    return ncmul
end
function splice!!(ncmul::NCMul, i, term::NCMul)
    coeff = prefactor(term) * prefactor(ncmul)
    newfactors = mysplice!!(ncmul.factors, i, term.factors)
    return NCMul(coeff, newfactors)
end
function splice!!(ncmul::NCMul{C,S}, i, term::S) where {C,S}
    splice!!(ncmul, i, NCMul(one(C), (term,)))
end
function splice!!(ncmul::NCMul, i::UnitRange, coeff::Number)
    deleteat!(ncmul.factors, i)
    return NCMul(coeff * prefactor(ncmul), ncmul.factors)
end
function splice!!(ncmul::NCMul, i::Integer, coeff::Number)
    deleteat!(ncmul.factors, i:i+1)
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

When sorting a non-commutative expression, this function defines what a*b should be replaced by. 
It must have a branch that returns `nothing` in order for the sorting to terminate. 

No methods for `mul_effect` are defined in this package, so you need to define them for your types.

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

@testitem "Collecting powers, signed swap" begin
    import NonCommutativeProducts: bubble_sort, @nc, NCMul, mul_effect
    NonCommutativeProducts.disable_autosort!()
    struct NCInt
        n::Int
        exp::Int
    end
    NCInt(n::Int) = NCInt(n, 1)
    Base.show(io::IO, x::NCInt) = print(io, "[$(x.n)^$(x.exp)]")
    Base.:(==)(a::NCInt, b::NCInt) = a.n == b.n
    @nc NCInt
    function mul_effect(a::NCInt, b::NCInt)
        (a.n == b.n && return NCInt(a.n, a.exp + b.exp)) # Collect powers
        a.n > b.n && return -b * a # swap and multiply by -1
        return nothing # Do nothing
    end

    a, b, c, d = NCInt.(1:4)
    ab = a * b
    @test ab == NCMul(1, [a, b])
    ab2 = bubble_sort(ab)
    ab3 = bubble_sort(b * a)
    @test ab == ab2 == -ab3
    @test hash(ab) == hash(ab2)
    @test bubble_sort(ab * a) == -1 * bubble_sort((a * ab))

    op = bubble_sort(prod([a, b, c, d, a, b, c, d, a, b]))
    @test length(op.dict) == 1
    @test map(x -> x.n, only(keys(op.dict)).factors) == 1:4
    @test sum(x -> x.exp, only(keys(op.dict)).factors) == 10
end