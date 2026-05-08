"""
    ncmap(f, expr)

Map `f` over the noncommutative factors in `expr`.
# Example
If `a` and `b` are of types registered with @nc then 
```julia
ncmap(f, 2*a*b + 1) == 2*f(a)*f(b) + 1
"""
ncmap(f, x::NCMul) = prod(f(factor) for factor in x.factors; init=prefactor(x))
ncmap(f, x::NCAdd) = sum(ncmap(f, term) for term in NCterms(x); init=additive_coeff(x))
ncmap(f, x::Number) = x

@testitem "ncmap" setup = [Fermions, Majoranas, Bosons] begin
    import NonCommutativeProducts: ncmap
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
