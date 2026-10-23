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

_scalar_terms(x::NCAdd) = Iterators.filter(!iszero, (NCMul(additive_coeff(x), _scalar_factors(x)),))
# The factors of the scalar part of a sum
_scalar_factors(::NCAdd{C,<:NCMul{<:Any,S}}) where {C,S} = S[]

"""
    nterms(x)

The number of terms of `x`, `length(collect(terms(x)))`, computed in constant time.

# Example
```julia
nterms(2 * a * b + 3 * c + 1) == 3
nterms(2 * a) == 1
nterms(0) == 0
```
"""
function nterms end

# the terms of a sum never have zero coefficients
nterms(x::NCAdd) = length(x.dict) + !iszero(additive_coeff(x))
nterms(x::NCMul) = Int(!iszero(x))
nterms(x::Number) = Int(!iszero(x))

"""
    termpairs(x)

An iterator over the terms of `x`, as in `terms(x)`, as pairs `factors => coefficient` without creating a product for
each term. The factors may be the internal storage of `x`, so they must not be mutated.

# Example
```julia
Dict(termpairs(2 * a * b + 3 * c + 1)) == Dict([a, b] => 2, [c] => 3, [] => 1)
collect(termpairs(2 * a)) == [[a] => 2]
```
"""
function termpairs end

termpairs(x::NCAdd) = Iterators.flatten((_scalar_termpairs(x), (k.factors => v for (k, v) in pairs(x.dict))))
termpairs(x::NCMul) = iszero(x) ? () : (factors_vector(x.factors) => prefactor(x),)
termpairs(x::Number) = iszero(x) ? () : (Any[] => x,)

_scalar_termpairs(x::NCAdd) = Iterators.filter(!iszero ∘ last, (_scalar_factors(x) => additive_coeff(x),))

"""
    map_coefficients(f, x)

`x` with the coefficient `c` of each term in `termpairs(x)` replaced by `f(factors, c)`. Terms that become zero are
dropped, and a zero result keeps the grades of `x`. The factors must not be mutated.

# Example
```julia
map_coefficients((fs, c) -> length(fs) * c, 2 * a * b + 3 * c + 1) == 4 * a * b + 3 * c
map_coefficients((fs, c) -> 2c, 3 * a) == 6 * a
```
"""
function map_coefficients end

function map_coefficients(f, x::NCAdd{C,K}) where {C,K}
    coeff = additive_coeff(x)
    newcoeff = iszero(coeff) ? nothing : f(_scalar_factors(x), coeff)
    newterms = [k => f(k.factors, v) for (k, v) in pairs(x.dict)]
    T = mapreduce(typeof ∘ last, promote_type, newterms; init=isnothing(newcoeff) ? Union{} : typeof(newcoeff))
    T = T === Union{} ? C : T
    result = NCAdd(isnothing(newcoeff) ? zero(T) : newcoeff, Dict{K,T}(newterms))
    return _with_zero_grades(result, () -> _grades(x))
end
map_coefficients(f, x::NCMul) = iszero(x) ? x : NCMul(f(factors_vector(x.factors), prefactor(x)), x.factors)
map_coefficients(f, x::Number) = iszero(x) ? x : f(Any[], x)

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
    n = count(Returns(true), terms(x))
    n == 1 || throw(ArgumentError("`$name` needs a sum with exactly one term, but this sum has $n terms: $x. Use `terms(x)` to get its products."))
    return first(terms(x))
end
