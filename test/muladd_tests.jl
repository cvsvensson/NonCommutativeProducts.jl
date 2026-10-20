@testitem "Addition and multiplication with different number types" setup = [Fermions] begin
    # This tests that addition and multiplication with different number types does not throw errors due to promotion issues and in place operations
    f1 = Fermion(1)
    f2 = Fermion(2)
    nums = (1, 1/2, 1//3, 0.1+1im, big(4))
    for x in nums
        @test x + f1 == f1 + x
        for y in nums
            @test x + y*f1 == y*f1 + x
            for z in nums
                @test z*(x+y*f1) isa NonCommutativeProducts.NCAdd
            end
        end
    end
end

@testitem "Multiplication with UniformScaling" setup = [Fermions] begin
    using LinearAlgebra: I
    a, b = Fermion(1), Fermion(2)
    @test a * I == a
    @test I * a == a
    @test (2a) * I == 2a
    @test I * (2a) == 2a
    @test a * (2I) == 2a
    @test (2I) * a == 2a
    @test (a * b) * I == a * b
    @test I * (a + b) == a + b
    @test (a + b) * I == a + b
    @test (a + 1) * (3I) == 3a + 3
    @test iszero((0I) * a)
    @test (2I) + a * b == 2 + a * b
end

@testitem "Multiplication by a zero UniformScaling keeps the grade" setup = [GradedKets] begin
    using LinearAlgebra: I
    import NonCommutativeProducts as NC
    z = (0I) * ket(1)
    @test iszero(z)
    @test NC.grades(z) == Set([KET])
    @test NC.grades(ket(1) * (0I)) == Set([KET])
    @test NC.grades((0I) * (ket(1) + ket(2))) == Set([KET])
end

@testitem "oneunit and one coefficient types" setup = [Fermions] begin
    f = Fermion.(1:2)
    m = 1.0 * f[1] * f[2]
    s = 1.0 * f[1] + f[2]
    @test oneunit(typeof(m)).coeff isa Float64
    @test oneunit(typeof(s)).coeff isa Float64
    @test oneunit(typeof(m)) == 1
    @test oneunit(typeof(s)) == 1
    @test one(typeof(m)).coeff isa Float64
    @test one(typeof(s)) isa typeof(s)
end

@testitem "Powers with any Integer exponent" setup = [Bosons] begin
    b = Boson()
    for x in (b, 2 * b' * b, b + b')
        for p in (Int8(3), Int32(3), UInt8(3), big(3))
            @test x^p == x^3
        end
        @test x^2 == x * x
        @test x^0 == 1
        n = -1
        @test_throws DomainError x^n
        @test_throws DomainError x^big(-2)
        @test_throws DomainError x^-1
        @test_throws DomainError x^-2
    end
end
