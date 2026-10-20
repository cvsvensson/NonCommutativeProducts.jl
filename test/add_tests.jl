@testitem "Consistency between + and add!!" setup = [Fermions] begin
    import NonCommutativeProducts: add!!
    using LinearAlgebra: I
    NonCommutativeProducts.disable_autosort!()
    f = Fermion.(1:2)
    a = 1.0 * f[2] * f[1] + 1 + f[1]
    for b in [1.0, 1, f[1], 1.0 * f[1], f[2] * f[1], a]
        a2 = copy(a)
        a3 = add!!(a2, b) # Should mutate
        @test a + b == a3
        @test a2 == a3
        anew = add!!(a, 1im * b) #Should not mutate
        @test a2 !== anew
    end
    @test a == 1.0 * f[2] * f[1] + 1 + f[1]

    NonCommutativeProducts.disable_autosort!()
    @test add!!(f[1] + 1, 1) == f[1] + 2
    @test add!!(f[1] + 1, 2I) == f[1] + 3
    @test add!!(f[1] + 1, f[1]) == 2f[1] + 1
    @test add!!(f[1] + 1, 5 * f[1]) == 6f[1] + 1
    @test add!!(f[1] + 1, 1, 2, 5) == 5f[1] + 7
    @test add!!(f[1] + 1, f[1], 2, 5) == 7f[1] + 5
    @test add!!(f[1] + 1, 5 * f[1], 2, 5) == 15f[1] + 5
end

@testitem "oneunit of NCAdd" setup = [Fermions] begin
    f1 = Fermion(1)
    x = 2 * f1 + 3
    @test oneunit(x) == 1
    @test oneunit(x) isa typeof(x)
    @test oneunit(typeof(x)) == oneunit(x)
    @test oneunit(x) * x == x
end
