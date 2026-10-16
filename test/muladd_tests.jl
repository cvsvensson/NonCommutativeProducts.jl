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
