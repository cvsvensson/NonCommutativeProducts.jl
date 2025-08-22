
struct AddTerms{T}
    terms::T
end
struct Swap{T}
    λ::T
end


function bubble_sort(a::NCMul, ordering)
    return bubble_sort!(copy(a), ordering)
end
function bubble_sort!(a::NCMul, ordering)
    if length(a.factors) <= 1
        return a
    end
    res = _bubble_sort!([a], ordering)
    return res
end
function bubble_sort(ncadd::NCAdd, ordering)
    terms = collect(NCMul(v, copy(k.factors)) for (k, v) in pairs(ncadd.dict))
    res = add!!(_bubble_sort!(terms, ordering), ncadd.coeff)
end
function _bubble_sort!(terms::Vector{T}, ordering) where {T<:NCMul}
    sorted_terms = __bubble_sort!(terms, ordering)
    if length(sorted_terms) == 0
        return NCAdd(0, Dict{T,Int}())
    end
    newadd = NCAdd(0, to_add_dict(first(sorted_terms)))
    for term in Iterators.drop(sorted_terms, 1)
        newadd = add!!(newadd, term)
    end
    return filter_scalars!(filter_zeros!(newadd))
end

function __bubble_sort!(terms::Vector{T}, ordering) where {T<:NCMul}
    n = 1
    while n <= length(terms)
        done = false
        start = 1
        while !done && n <= length(terms)
            terms, done, start = __bubble_sort!(terms, n, ordering, start)
            if iszero(terms[n].coeff)
                done = false
                start = 1
                deleteat!(terms, n)
            end
        end
        n += 1
    end
    return terms
end

function __bubble_sort!(terms::Vector{NCMul{C,S,F}}, index, ordering, start) where {C,S,F}
    no_effect = true
    ncmul = terms[index]
    if length(ncmul.factors) <= 1
        done = true
        return terms, done, start
    end
    i::Int = max(0, start - 1)
    N::Int = length(ncmul.factors)
    while no_effect && i < N - 1
        i += 1
        a, b = ncmul.factors[i]::S, ncmul.factors[i+1]::S
        effect = mul_effect(a, b, ordering)
        isnothing(effect) && continue

        no_effect = false
        ncmul, newterms = splice!!_and_add(ncmul, i:i+1, effect)
        if length(newterms) > 0
            terms = append!!(terms, newterms)
        end
        terms = setindex!!(terms, ncmul, index)
    end
    done = no_effect
    newstart = i - 1
    return terms::Vector{NCMul{C,S,F}}, done, newstart
end


function mysplice!!(v::V, i, replacement::W) where {V,W}
    T = promote_type(eltype(W), eltype(V))
    if T <: eltype(V)
        splice!(v, i, replacement)
        return v
    else
        return Vector{T}(vcat(v[1:first(i)-1], replacement, v[last(i)+1:end]))
    end
end
function splice!!(ncmul::NCMul, i)
    splice!(ncmul.factors, i)
    return ncmul
end
function splice!!(ncmul::NCMul, i, term::NCMul)
    ncmul.coeff *= term.coeff
    newfactors = mysplice!!(ncmul.factors, i, term.factors)
    newfactors === ncmul.factors && return ncmul
    return NCMul(ncmul.coeff, newfactors)
end
function splice!!(ncmul::NCMul{C,S}, i, term::S) where {C,S}
    splice!!(ncmul, i, NCMul(1, (term,)))
end
function splice!!(ncmul::NCMul, i, coeff::Number)
    deleteat!(ncmul.factors, i)
    ncmul.coeff *= coeff
    return ncmul
end
function splice!!(ncmul::NCMul, i, swap::Swap)
    length(i) == 2 && i[2] == i[1] + 1 || throw(ArgumentError("Invalid index for swap"))
    a, b = ncmul.factors[i]
    splice!!(ncmul, i, NCMul(swap.λ, [b, a]))
end
function splice!!_and_add(ncmul::NCMul, i, terms::AddTerms)
    newterms = [splice!!(copy(ncmul), i, term) for term in Iterators.drop(terms.terms, 1)]
    ncmul = splice!!(ncmul, i, first(terms.terms))
    return ncmul, newterms
end
function splice!!_and_add(ncmul::T, i, effect) where T<:NCMul
    ncmul = splice!!(ncmul, i, effect)
    return ncmul, T[]
end

bubble_sort(a::Number, ordering; kwargs...) = a
function mul_effect end

@testitem "Collecting powers, signed swap" begin
    import NonCommutativeProducts: bubble_sort, Swap, @nc, NCMul, mul_effect, AddTerms

    struct NCInt
        n::Int
        exp::Int
    end
    NCInt(n::Int) = NCInt(n, 1)
    Base.show(io::IO, x::NCInt) = print(io, "[$(x.n)^$(x.exp)]")
    Base.:(==)(a::NCInt, b::NCInt) = a.n == b.n
    @nc NCInt
    struct IntOrder end
    function mul_effect(a::NCInt, b::NCInt, ::IntOrder)
        (a.n == b.n && return NCInt(a.n, a.exp + b.exp)) # Collect powers
        a.n > b.n && return Swap(-1) # swap and multiply by -1
        return nothing # Do nothing
    end

    a, b, c, d = NCInt.(1:4)
    ab = a * b
    @test ab == NCMul(1, [a, b])
    ab2 = bubble_sort(ab, IntOrder())
    ab3 = bubble_sort(b * a, IntOrder())
    @test ab == ab2 == -ab3
    @test hash(ab) == hash(ab2)
    @test bubble_sort(ab * a, IntOrder()) == -1 * bubble_sort((a * ab), IntOrder())

    op = bubble_sort(prod([a, b, c, d, a, b, c, d, a, b]), IntOrder())
    @test length(op.dict) == 1
    @test map(x -> x.n, only(keys(op.dict)).factors) == 1:4
    @test sum(x -> x.exp, only(keys(op.dict)).factors) == 10
end