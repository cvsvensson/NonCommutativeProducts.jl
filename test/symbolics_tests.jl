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

@testitem "Symbolics: substitute acts on coefficients" setup = [Fermions] begin
    using Symbolics
    @variables x y
    c = Fermion.(1:2)
    a = x * c[1] * c[2] + y * c[2] + x
    @test isequal(substitute(a, Dict(x => 2)), 2 * c[1] * c[2] + y * c[2] + 2)
    @test isequal(substitute(x * c[1], Dict(x => y)), y * c[1])
    @test iszero(substitute(a, Dict(x => 0, y => 0)))
end

@testitem "Symbolics: truncation keeps terms with symbolic coefficients" begin
    using Symbolics
    import NonCommutativeProducts as NC
    struct BoundedAtom
        n::Int
    end
    NC.@nc BoundedAtom
    NC.mul_effect(::BoundedAtom, ::BoundedAtom) = nothing
    NC.opnorm_bound(::BoundedAtom) = 1
    @variables x
    c = BoundedAtom.(1:3)
    # constants have the bound of their value
    y = x * c[1] + 0.001 * c[2] + 2 * c[3]
    result = NC.truncate_terms(y; atol=1)
    @test isequal(result.expr, x * c[1] + 2 * c[3])
    @test result.error == 0.001
    @test NC.opnorm_bound(y) === nothing
end
