
struct ScalarMul end
struct AddTerm end

# struct Power{S}
#     sym::S
#     exponent::Int
# end
# struct Nilpotent end
# self_product_effect_type(::Fermion) = Nilpotent()

should_swap(a, b, ::NaiveOrdering) = false
should_swap(a::Fermion, b::Fermion, ::NormalOrdering) = !a.creation && b.creation

function swap_effect_type(a::Fermion, b::Fermion)
    if a == b
        return self_product_effect(a)
    end
    a.label == b.label && xor(a.creation, b.creation) && return AddTerm()
    return ScalarMul()
end

## Normal ordering
function bubble_sort(a::NCAdd; start=1)
    c = zero(a)
    c.coeff = a.coeff
    for term in NCterms(a)
        add!(c, bubble_sort(term; start))
    end
    return c
end

function bubble_sort(a::NCMul; start=1)
    if length(a.factors) == 1
        return a
    end
    muloraddvec::Union{Number,MulAdd} = a

    swapped = false
    i = max(0, start - 1)
    triple_prod(a, b, c) = ordered_product(ordered_product(a, b, NaiveOrdering()), c, NaiveOrdering())
    while !swapped && i < length(eachindex(a.factors)) - 1
        i += 1
        if a.factors[i] > a.factors[i+1] || isequal(a.factors[i], a.factors[i+1])
            swapped = true
            product = ordered_product(a.factors[i], a.factors[i+1], NormalOrdering())
            left_factors = NCMul(a.coeff, a.factors[1:i-1])
            right_factors = NCMul(1, a.factors[i+2:end])
            muloraddvec = triple_prod(left_factors, product, right_factors)
        end
    end
    if !swapped
        return a
    end
    bubble_sort(muloraddvec; start=i - 1)
end

# normal_order(a::MulAdd) = bubble_sort(a)
# normal_order(a::Number) = a
bubble_sort(a::Number; kwargs...) = a
