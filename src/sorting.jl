
struct ScalarMul{T}
    λ::T
end
scalar(s::ScalarMul) = s.λ

struct AddTerm end
struct NCExp{S}
    term::S
    exp::Int
end
Base.show(io::IO, e::NCExp) = print(io, "$(e.term)^$(e.exp)")


struct Nilpotent end
struct MaybeSwap end
self_product_effect_type(::Fermion) = Nilpotent()

should_swap(a, b, ::NaiveOrdering) = false
should_swap(a::Fermion, b::Fermion, ::NormalOrdering) = !a.creation && b.creation

function swap_effect(a::Fermion, b::Fermion)
    if a == b
        return self_product_effect(a)
    end
    a.label == b.label && xor(a.creation, b.creation) && return AddTerm()
    return ScalarMul(-1)
end

## Normal ordering
function bubble_sort(a::NCAdd, ordering; start=1)
    c = zero(a)
    c.coeff = a.coeff
    for term in NCterms(a)
        add!(c, bubble_sort(term, ordering; start))
    end
    return c
end
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

function bubble_sort(ncmul::NCMul, ordering; start=1)
    if length(ncmul.factors) == 1
        return ncmul
    end
    # muloraddvec::Union{Number,MulAdd} = a
    # newterms = []
    swapped = false
    i = max(0, start - 1)
    # triple_prod(a, b, c) = ordered_product(ordered_product(a, b, NaiveOrdering()), c, NaiveOrdering())
    while !swapped && i < length(eachindex(ncmul.factors)) - 1
        i += 1
        a, b = ncmul.factors[i], ncmul.factors[i+1]
        effect = mul_effect(a, b)
        if effect isa NCMul
            swapped = true
            ncmul = NCMul(ncmul.coeff, vcat(ncmul.factors[1:i-1], effect.factors, ncmul.factors[i+2:end]))
        elseif effect isa NCExp
            swapped = true
            ncmul = NCMul(ncmul.coeff, vcat(ncmul.factors[1:i-1], [effect], ncmul.factors[i+2:end]))
        elseif effect isa MaybeSwap
            if should_swap(a, b, ordering)
                swapped = true
                effect = swap_effect(a, b)
                if effect isa ScalarMul
                    ncmul = NCMul(scalar(effect) * ncmul.coeff, vcat(ncmul.factors[1:i-1], [b, a], ncmul.factors[i+2:end]))
                end
                # product = ordered_product(ncmul.factors[i], ncmul.factors[i+1], ordering)
                # muloraddvec = triple_prod(left_factors, product, right_factors)
            end
        elseif !isnothing(effect)
            throw(ArgumentError("Don't know how to multiply $a * $b"))
        end
    end
    if !swapped
        return ncmul
    end
    bubble_sort(ncmul, ordering; start=i - 1)
end

# normal_order(a::MulAdd) = bubble_sort(a)
# normal_order(a::Number) = a
bubble_sort(a::Number, ordering; kwargs...) = a
function mul_effect end

@testitem "Signed permutation" begin
    using NonCommutativeProducts
    import NonCommutativeProducts: ordered_product, NaiveOrdering, should_swap, NCMul, NCAdd, bubble_sort, swap_effect, ScalarMul, NCExp, MaybeSwap, mul_effect

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
    ab = ordered_product(a, b, NaiveOrdering())
    @test ab == NCMul(1, [a, b])
    @test bubble_sort(ab, IntOrder()) == ab == -bubble_sort(b * a, IntOrder())
    @test bubble_sort(ab * a, IntOrder()) == bubble_sort(-1 * (a * ab), IntOrder())

    bubble_sort(prod(rand((a, b, d, c)) for k in 1:10), IntOrder())
end