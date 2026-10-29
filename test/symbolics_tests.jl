@testitem "Symbolics: add!! with the default weights" setup = [Fermions] begin
    using Symbolics
    import VectorInterface
    @variables x y
    c = Fermion.(1:2)
    a = x * c[1] + y * c[2]
    @test isequal(NonCommutativeProducts.add!!(copy(a), a), 2x * c[1] + 2y * c[2])
    @test isequal(VectorInterface.add!!(copy(a), a, x), (x + x^2) * c[1] + (y + x * y) * c[2])
end

@testitem "Symbolics: isequal" setup = [Fermions] begin
    using Symbolics
    @variables x y
    c = Fermion.(1:2)
    @test isequal(x * c[1], x * c[1])
    @test !isequal(x * c[1], y * c[1])
    @test isequal(x * c[1] + y, y + x * c[1])
    @test isequal(x * c[1] + 0, x * c[1])
    @test isequal(1 * c[1] + 0, c[1])
    @test length(Set([x * c[1] + y, y + x * c[1], x * c[2]])) == 2
end

@testitem "Symbolics: show parenthesizes a symbolic scalar term" setup = [Fermions] begin
    using Symbolics
    @variables x y
    c = Fermion(1)
    @test sprint(show, x + y + c) == "(x + y)I + c[1]"
end
