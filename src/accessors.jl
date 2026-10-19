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
    factors(x)

The factors of `x` as an `AbstractVector`, in order and without the coefficient. Works for products, numbers
and atoms registered with `@nc`: a number has no factors, and an atom is its own single factor. A sum with exactly one
term has the factors of that term. A sum with more terms throws an `ArgumentError`; use `terms` to get its products.
The result is a copy, so modifying it does not change `x`.

# Example
```julia
factors(2 * a * b) == [a, b]
factors(2) == []
factors(a) == [a]
```
"""
factors(m::NCMul) = collect(eltype(m.factors), m.factors)
factors(::Number) = Any[]
factors(x::NCAdd) = factors(_only_term(x, "factors"))

"""
    coefficient(x)

The scalar coefficient of `x`. Works for products, numbers and atoms registered with `@nc`: a number is its own
coefficient, and an atom has coefficient `1`. A sum with exactly one term has the coefficient of that term. A sum with
more terms throws an `ArgumentError`; use `terms` to get its products.

# Example
```julia
coefficient(2 * a * b) == 2
coefficient(2) == 2
coefficient(a) == 1
```
"""
coefficient(m::NCMul) = prefactor(m)
coefficient(x::Number) = x
coefficient(x::NCAdd) = coefficient(_only_term(x, "coefficient"))

# The only term of a sum, as a product that includes its coefficient
function _only_term(x::NCAdd, name::AbstractString)
    ts = collect(terms(x))
    length(ts) == 1 || throw(ArgumentError("`$name` needs a sum with exactly one term, but this sum has $(length(ts)) terms: $x. Use `terms(x)` to get its products."))
    return only(ts)
end
