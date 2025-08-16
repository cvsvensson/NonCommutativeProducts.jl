
struct ScalarMul{T}
    λ::T
end
scalar(s::ScalarMul) = s.λ

struct AddTerm end
# struct Power{S}
#     sym::S
#     exponent::Int
# end
# struct Nilpotent end
# self_product_effect_type(::Fermion) = Nilpotent()

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

function bubble_sort(ncmul::NCMul, ordering; start=1)
    if length(ncmul.factors) == 1
        return ncmul
    end
    # muloraddvec::Union{Number,MulAdd} = a
    newterms = []
    swapped = false
    i = max(0, start - 1)
    triple_prod(a, b, c) = ordered_product(ordered_product(a, b, NaiveOrdering()), c, NaiveOrdering())
    while !swapped && i < length(eachindex(ncmul.factors)) - 1
        i += 1
        a, b = ncmul.factors[i], ncmul.factors[i+1]
        if should_swap(a, b, ordering)
            swapped = true
            # left_factors = NCMul(ncmul.coeff, ncmul.factors[1:i-1])
            # right_factors = NCMul(1, ncmul.factors[i+2:end])
            effect = swap_effect(a, b)
            if effect isa ScalarMul
                ncmul = NCMul(scalar(effect) * ncmul.coeff, vcat(ncmul.factors[1:i-1], [b, a], ncmul.factors[i+2:end]))
            end
            # product = ordered_product(ncmul.factors[i], ncmul.factors[i+1], ordering)
            # muloraddvec = triple_prod(left_factors, product, right_factors)
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


@testitem "Signed permutation" begin
    using NonCommutativeProducts
    import NonCommutativeProducts: ordered_product, NaiveOrdering, should_swap, NCMul, NCAdd, bubble_sort, swap_effect, ScalarMul

    struct NCInt
        n::Int
    end
    Base.:==(a::NCInt, b::NCInt) = a.n == b.n
    struct IntOrder end
    Base.:*(a::NCInt, b::NCInt) = NCMul(1, [a, b])
    should_swap(a::NCInt, b::NCInt, ::IntOrder) = a.n > b.n
    swap_effect(::NCInt, ::NCInt) = ScalarMul(-1)


    a, b, c, d = NCInt.(1:4)
    ab = ordered_product(a, b, NaiveOrdering())
    @test ab == NCMul(1, [a, b])
    @test bubble_sort(ab, IntOrder()) == ab
    @test bubble_sort(ab * a, IntOrder()) == -1 * (a * ab)

end