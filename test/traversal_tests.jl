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
