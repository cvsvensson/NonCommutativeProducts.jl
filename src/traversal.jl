"""
    ncmap(f, expr)

Map the user-defined leaf function `f` over the noncommutative symbolic
factors in `expr`.

Scalar coefficients and constant terms are left unchanged.
"""
ncmap(f, x::NCMul) = prod(f(factor) for factor in x.factors; init=prefactor(x))
ncmap(f, x::NCAdd) = sum(ncmap(f, term) for term in NCterms(x); init=additive_coeff(x))
ncmap(f, x::Number) = x

@testitem "ncmap" setup = [Fermions, Majoranas] begin
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

    scale_fermion(x::Fermion) = 2 * x
    expected = 3 - 8 * f + 8 * f' * f
    @test ncmap(scale_fermion, expr) == expected

    γ = Majorana.(1:2)
    ferm_to_maj(x::Fermion) = x.creation ? γ[1] - im * γ[2] : γ[1] + im * γ[2]
    ferm_to_maj(x::Majorana) = x
    expr2 = 1 + f + 2 * f' * f + γ[1]
    expected2 = 1 + (γ[1] + im * γ[2]) + 2 * (γ[1] - im * γ[2]) * (γ[1] + im * γ[2]) + γ[1]
    @test ncmap(ferm_to_maj, expr2) == expected2

    @test ncmap(adjoint, f'f) == f * f'
    NonCommutativeProducts.enable_autosort!()
    @test ncmap(adjoint, f'f) == 1 - f' * f
end
