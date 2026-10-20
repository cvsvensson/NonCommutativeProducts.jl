"""
    ncmapreduce(f, ops, expr; scalarmap=identity)

Map `f` over the noncommutative factors of `expr`, then reduce them with the binary operations `ops = (add, mul)`.
`mul` combines the factors of each term, starting from the term's prefactor. `add` combines the terms, starting from
the additive constant. Each reduction is a left fold. `scalarmap` is applied to each prefactor, to the additive
constant, and to a bare number `expr`. Returns the single value produced by the reduction.

# Example 
(using infix notation for add and mul)
```julia
ncmapreduce(f, (add, mul), 2*a*b + 1; scalarmap=g) == g(1) add ((g(2) mul f(a)) mul f(b))
```
"""
function ncmapreduce(f, (add, mul), x::NCMul; scalarmap=identity)
    return mapfoldl(f, mul, x.factors; init=scalarmap(prefactor(x)))
end

function ncmapreduce(f, (add, mul), x::NCAdd; scalarmap=identity)
    terms = (ncmapreduce(f, (add, mul), term; scalarmap=scalarmap) for term in NCterms(x))
    return foldl(add, terms; init=scalarmap(additive_coeff(x)))
end

ncmapreduce(f, (add, mul), x::Number; scalarmap=identity) = scalarmap(x)

"""
    ncmap(f, expr; scalarmap=identity)

Map `f` over the noncommutative factors in `expr`, keeping the sums and products of `expr`. `scalarmap` is applied to
the scalar coefficients: the prefactor of each term and the additive constant. It defaults to `identity`.
Equivalent to `ncmapreduce(f, (+, *), expr; scalarmap)`.

# Example
If `a` and `b` are of types registered with @nc then 
```julia
ncmap(f, 2*a*b + 1) == 2*f(a)*f(b) + 1
```
See `ncmapreduce` for more general mapping and reduction.
"""
ncmap(f, x; scalarmap=identity) = ncmapreduce(f, (+, *), x; scalarmap)
