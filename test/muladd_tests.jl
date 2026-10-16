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
                @test z*(x+y*f1) isa NonCommutativeProducts.NCAdd #throws error on v0.4.4
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
