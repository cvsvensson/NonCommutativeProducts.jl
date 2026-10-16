"""
    terms(x)

An iterator over the terms of `x`, each an `NCMul` that includes its coefficient. A nonzero scalar part of a sum
appears as a product without factors, so `sum(terms(x)) == x`. Zero has no terms: `terms(0)`, the terms of a zero
sum and the terms of a zero product are empty, and `sum` of an empty iterator throws.

Works for sums, products, numbers and atoms registered with `@nc`.

# Example
```julia
collect(terms(2 * a * b + 3 * c + 1)) # the products 2*a*b, 3*c and 1, in no particular order
terms(2 * a) # the single product 2*a
```
"""
function terms end

terms(x::NCAdd) = Iterators.flatten((_scalar_terms(x), Iterators.filter(!iszero, NCterms(x))))
terms(x::NCMul) = iszero(x) ? () : NCterms(x)
terms(x::Number) = iszero(x) ? () : (NCMul(x, Any[]),)

function _scalar_terms(x::NCAdd{C,<:NCMul{<:Any,S}}) where {C,S}
    return Iterators.filter(!iszero, (NCMul(additive_coeff(x), S[]),))
end

"""
    factors(m::NCMul)

The factors of the product `m` as an `AbstractVector`, in order and without the coefficient. The result is a copy,
so modifying it does not change `m`, and it is a `Vector` whether or not `m` stores its factors in a `Tuple`.

# Example
```julia
factors(2 * a * b) == [a, b]
```
"""
factors(m::NCMul) = collect(eltype(m.factors), m.factors)

"""
    coefficient(m::NCMul)

The scalar coefficient of the product `m`.

# Example
```julia
coefficient(2 * a * b) == 2
```
"""
coefficient(m::NCMul) = prefactor(m)
