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

@testitem "ncmap" setup = [Fermions, Majoranas, Bosons] begin
    import NonCommutativeProducts: ncmap, ncmapreduce
    NonCommutativeProducts.disable_autosort!()
    f = Fermion(:a)
    change_label(x::Fermion) = Fermion(:b, x.creation)
    @test ncmap(change_label, f) == Fermion(:b)
    @test ncmap(change_label, 1) == 1
    expr = 3 - 4 * f + 2 * f' * f
    fb = Fermion(:b)
    expected = 3 - 4 * fb + 2 * fb' * fb
    @test ncmap(change_label, expr) == expected
    @test ncmap(identity, expr) == expr

    scalar_scale(x::Number) = 2 * x
    @test ncmap(identity, expr; scalarmap=scalar_scale) == 6 - 8 * f + 4 * f' * f
    @test ncmap(change_label, 1; scalarmap=scalar_scale) == 2

    # Test ncmapreduce with custom operations, here by interchangin + and *.
    fermion_weight(x::Fermion) = x.creation ? 2 : 1
    @test ncmapreduce(fermion_weight, (*, +), 2 * f' * f + 6; scalarmap=x -> x + 1) == (3 + 2 + 1) * 7

    scale_fermion_and_boson(x::Fermion) = 2 * x
    scale_fermion_and_boson(x::Boson) = 3 * x
    b = Boson()
    expr = 3 - 4 * f + 2 * f' * b * f + b
    expected = 3 - 8 * f + 24 * f' * b * f + 3 * b
    @test ncmap(scale_fermion_and_boson, expr) == expected

    γ = Majorana.(1:2)
    ferm_to_maj(x::Fermion) = x.creation ? γ[1] - im * γ[2] : γ[1] + im * γ[2]
    ferm_to_maj(x::Majorana) = x
    expr = 1 + f + 2 * f' * f + γ[1]
    expected = 1 + (γ[1] + im * γ[2]) + 2 * (γ[1] - im * γ[2]) * (γ[1] + im * γ[2]) + γ[1]
    @test ncmap(ferm_to_maj, expr) == expected

    @test ncmap(adjoint, f'f) == f * f'
    NonCommutativeProducts.enable_autosort!()
    @test ncmap(adjoint, f'f) == 1 - f' * f
end
