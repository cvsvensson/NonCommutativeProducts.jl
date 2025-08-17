
struct ScalarMul{T}
    λ::T
end
scalar(s::ScalarMul) = s.λ

struct AddTerms{T}
    terms::T
end
struct NCExp{S}
    term::S
    exp::Int
end
Base.show(io::IO, e::NCExp) = print(io, "$(e.term)^$(e.exp)")


struct Nilpotent end
struct MaybeSwap end

should_swap(a, b, ::NaiveOrdering) = false


## Normal ordering
struct NCExpEffect{E}
    effect::E
end
function mul_effect(a::NCExp, b::NCExp)
    term_effect = mul_effect(a.term, b.term)
    if term_effect isa NCExp
        return NCExp(a.term, a.exp + b.exp)
    elseif term_effect isa MaybeSwap
        return MaybeSwap()
    else
        return NCExpEffect(term_effect)
    end
end
function mul_effect(a::NCExp, b)
    term_effect = mul_effect(a.term, b)
    if term_effect isa NCExp
        return NCExp(a.term, a.exp + 1)
    elseif term_effect isa MaybeSwap
        return MaybeSwap()
    else
        return NCExpEffect(term_effect)
    end
end
function mul_effect(a, b::NCExp)
    term_effect = mul_effect(a, b.term)
    if term_effect isa NCExp
        return NCExp(b.term, b.exp + 1)
    elseif term_effect isa MaybeSwap
        return MaybeSwap()
    else
        return NCExpEffect(term_effect)
    end
end
should_swap(a::NCExp, b, ordering) = should_swap(a.term, b, ordering)
should_swap(a, b::NCExp, ordering) = should_swap(a, b.term, ordering)
should_swap(a::NCExp, b::NCExp, ordering) = should_swap(a.term, b.term, ordering)
function swap_effect(a::NCExp, b::NCExp)
    term_effect = swap_effect(a.term, b.term)
    if term_effect isa ScalarMul
        return ScalarMul(term_effect.λ^(a.exp * b.exp))
    else
        throw(ArgumentError("Don't know how to swap $a and $b"))
    end
end
swap_effect(a::NCExp, b) = swap_effect(a, NCExp(b, 1))
swap_effect(a, b::NCExp) = swap_effect(NCExp(a, 1), b)


function bubble_sort(a::NCMul, ordering; kwargs...)
    bubble_sort(NCAdd(0, to_add(a)), ordering; kwargs...)
end
function bubble_sort(ncadd::NCAdd, ordering)
    terms = collect(prod(p) for p in pairs(ncadd.dict))
    sorted_terms = bubble_sort!(terms, ordering)
    newadd = zero(ncadd)
    for term in sorted_terms
        newadd = tryadd!(newadd, term)
    end
    return ncadd.coeff + newadd
end

function bubble_sort!(terms::AbstractVector{<:NCMul}, ordering)
    # ncmul = copy(terms)
    # done = false
    # count = 0 # Just to avoid infinite loops for now. Remove later
    # start = 1
    # terms_to_sort = [(ncmul, false, 1)]
    n = 1
    while n <= length(terms)
        done = false
        count = 0
        start = 1
        while !done && count < 10
            count += 1
            terms, done, start = bubble_sort!(terms, n, ordering; start)
        end
        n += 1
    end
    return terms
    # while !done && count < 10
    #     count += 1
    #     # newterms = []
    #     ncmul, done, start = bubble_sort!(ncmul, ordering, terms_to_sort; start)
    # end
end

function bubble_sort!(terms::AbstractVector{<:NCMul}, index, ordering; start=1)
    no_effect = true
    ncmul = terms[index]
    if length(ncmul.factors) <= 1
        done = true
        return terms, done, start
    end
    # muloraddvec::Union{Number,MulAdd} = a
    # newterms = []
    i = max(0, start - 1)
    # triple_prod(a, b, c) = ordered_product(ordered_product(a, b, NaiveOrdering()), c, NaiveOrdering())
    while no_effect && i < length(eachindex(ncmul.factors)) - 1
        i += 1
        a, b = ncmul.factors[i], ncmul.factors[i+1]
        effect = mul_effect(a, b)
        if effect isa NCMul
            no_effect = false
            ncmul = splice!!(ncmul, i:i+1, effect.factors)
        elseif effect isa NCExp
            no_effect = false
            ncmul = splice!!(ncmul, i:i+1, [effect])
        elseif effect isa ScalarMul
            no_effect = false
            ncmul.coeff *= effect.λ
            ncmul = splice!!(ncmul, i:i+1, ncmul.factors[i+2:end])
        elseif effect isa MaybeSwap
            if should_swap(a, b, ordering)
                no_effect = false
                effect = swap_effect(a, b)
                if effect isa ScalarMul
                    ncmul.coeff *= effect.λ
                    ncmul = splice!!(ncmul, i:i+1, NCMul(1, [b, a]))
                elseif effect isa AddTerms
                    no_effect = false
                    ncmul, newterms = splice!!(ncmul, i:i+1, effect)
                    terms = append!!(terms, newterms)
                end
            end
        elseif !isnothing(effect)
            throw(ArgumentError("Don't know how to multiply $a * $b"))
        end
        if !no_effect
            terms = setindex!!(terms, ncmul, index)
        end
    end
    done = no_effect
    newstart = i - 1
    return terms, done, newstart
end

# function BangBang.splice!!(ncmul::NCMul{C,S}, i, replacement::V) where {C,S,V}
#     possible = promote_type(S, eltype(V)) <: S
#     if possible
#         splice!(ncmul.factors, i, replacement)
#     else
#         # println(ncmul.factors |> typeof)
#         println(S)
#         println(replacement |> typeof)
#         T = Union{S,eltype(replacement)}
#         ncmul = NCMul{C,T}(ncmul.coeff, Vector{T}(vcat(ncmul.factors[1:first(i)-1], replacement, ncmul.factors[last(i)+1:end])))
#         # println(ncmul.factors |> typeof)
#     end
#     return ncmul#NCMul(ncmul.coeff, newfactors)
# end

function mysplice!!(v::V, i, replacement::W) where {V,W}
    T = promote_type(eltype(W), eltype(V))
    if T <: eltype(V)
        splice!(v, i, replacement)
        return v
    else
        return Vector{T}(vcat(v[1:first(i)-1], replacement, v[last(i)+1:end]))
    end
end
function BangBang.splice!!(ncmul::NCMul, i, term::NCMul)
    ncmul.coeff *= term.coeff
    newfactors = mysplice!!(ncmul.factors, i, term.factors)
    newfactors === ncmul.factors && return ncmul
    return NCMul(ncmul.coeff, newfactors)
end
function BangBang.splice!!(ncmul::NCMul, i, coeff::Number)
    deleteat!(ncmul.factors, i)
    ncmul.coeff *= coeff
    return ncmul
end
function BangBang.splice!!(ncmul::NCMul, i, terms::AddTerms)
    newterms = [splice!!(copy(ncmul), i, term) for term in Iterators.drop(terms.terms, 1)]
    ncmul = splice!!(ncmul, i, first(terms.terms))
    return ncmul, newterms
    # terms = append!!(terms, newterms)
end
# function append_terms!(terms_with_sort_info, effect::AddTerms)
#     for term in effect.newterms
#         push!(terms_with_sort_info, (term, false, 1))
#     end
# end

# normal_order(a::MulAdd) = bubble_sort(a)
# normal_order(a::Number) = a
bubble_sort(a::Number, ordering; kwargs...) = a
# function mul_effect end

@testitem "Signed permutation" begin
    using NonCommutativeProducts
    import NonCommutativeProducts: ordered_product, NaiveOrdering, should_swap, NCMul, NCAdd, bubble_sort, swap_effect, ScalarMul, NCExp, MaybeSwap, mul_effect, AddTerms

    struct NCInt
        n::Int
    end
    Base.show(io::IO, x::NCInt) = print(io, "̲$(x.n)")
    Base.:(==)(a::NCInt, b::NCInt) = a.n == b.n
    struct IntOrder end
    Base.:*(a::NCInt, b::NCInt) = NCMul(1, [a, b])
    mul_effect(a::NCInt, b::NCInt) = (a == b && return NCExp(a, 2)) || return MaybeSwap()
    should_swap(a::NCInt, b::NCInt, ::IntOrder) = a.n > b.n
    swap_effect(::NCInt, ::NCInt) = ScalarMul(-1)


    a, b, c, d = NCInt.(1:4)
    ab = a * b
    # ab = ordered_product(a, b, NaiveOrdering())
    @test ab == NCMul(1, [a, b])
    ab2 = bubble_sort(ab, IntOrder())
    ab3 = -bubble_sort(b * a, IntOrder())
    @test ab == ab2 == ab3
    @test hash(ab) == hash(ab2)
    @test bubble_sort(ab * a, IntOrder()) == bubble_sort(-1 * (a * ab), IntOrder())

    bubble_sort(prod(rand((a, b, d, c)) for k in 1:10), IntOrder())
end