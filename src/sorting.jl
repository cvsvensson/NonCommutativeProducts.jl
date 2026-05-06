Base.sort(a::NCMul) = bubble_sort(a)
Base.sort(a::NCAdd) = bubble_sort(a)
Base.sort!(a::NCMul) = bubble_sort!(a)
Base.sort!(a::NCAdd) = bubble_sort!(a)
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
            if iszero(prefactor(terms[n]))
                done = false
                start = 1
                deleteat!(terms, n)
            end
        end
        n += 1
    end
    return terms
end

function __bubble_sort!(terms::Vector, index, start)
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
        a, b = ncmul.factors[i], ncmul.factors[i+1]
        effect = mul_effect(a, b)
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
    return terms, done, newstart
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
    coeff = prefactor(term) * prefactor(ncmul)
    newfactors = mysplice!!(ncmul.factors, i, term.factors)
    return NCMul(coeff, newfactors)
end
function splice!!(ncmul::NCMul{C,S}, i, term::S) where {C,S}
    splice!!(ncmul, i, NCMul(1, (term,)))
end
function splice!!(ncmul::NCMul, i, coeff::Number)
    deleteat!(ncmul.factors, i)
    return coeff * ncmul
end
function splice!!_and_add(ncmul::NCMul, i, add::NCAdd)
    if iszero(additive_coeff(add))
        terms = NCterms(add)
        newterms = [splice!!(copy(ncmul), i, term) for term in Iterators.drop(terms, 1)]
        ncmul2 = splice!!(ncmul, i, first(terms))
    else
        newterms = [splice!!(copy(ncmul), i, term) for term in NCterms(add)]
        ncmul2 = splice!!(ncmul, i, additive_coeff(add))
    end
    return ncmul2, newterms
end
function splice!!_and_add(ncmul::T, i, effect) where T<:NCMul
    ncmul2 = splice!!(ncmul, i, effect)
    return ncmul2, T[]
end

bubble_sort(a::Number; kwargs...) = a
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