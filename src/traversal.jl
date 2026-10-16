"""
    ncmapreduce(f, ops, expr; scalarmap=identity)

Map `f` over the noncommutative factors in `expr` and combine terms with
custom binary operations in `ops = (add, mul)`. Scalars are transformed with `scalarmap`. Each term is reduced with `mul` and then the terms are reduced with `add`. Uses foldl for the ordering of operations.

# Example 
(using infix notation for op1 and op2)
```julia
ncmapreduce(f, (op1, op2), 2*a*b + 1; scalarmap=g) == ((g(2) op1 f(a)) op1 f(b)) op2 g(1)
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
    ncmap(f, expr)

Map `f` over the noncommutative factors in `expr`.

# Example
If `a` and `b` are of types registered with @nc then 
```julia
ncmap(f, 2*a*b + 1) == 2*f(a)*f(b) + 1
```
See `ncmapreduce` for more general mapping and reduction.
"""
ncmap(f, x; scalarmap=identity) = ncmapreduce(f, (+, *), x; scalarmap)
