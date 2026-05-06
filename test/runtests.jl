using TestItemRunner

@testmodule Fermions begin
    export Fermion
    include("rules/fermions.jl")
end

@testmodule Bosons begin
    export Boson
    include("rules/bosons.jl")
end

@testitem "Fermions" setup = [Fermions] begin
    import NonCommutativeProducts: bubble_sort
    using Symbolics, LinearAlgebra
    using Random: seed!
    seed!(1)
    @variables a::Real z::Complex
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    f3 = Fermion((1, :↑))
    NonCommutativeProducts.disable_autosort!()
    ord(op) = bubble_sort(op)
    ord_equals(a, b) = iszero(ord(a - b))
    ord_equals(a, b, c, xs...) = (ord_equals(a, b) && ord_equals(b, c, xs...))

    @test 1 * f1 == f1
    @test 1 * f1 + 0 == f1
    @test 1 * f1 + 0 == 1 * f1
    @test 0 * f1 == 0
    @test 1 * f1 !== 1
    @test 0 * f1 + 1 == 1
    @test f1 / 2 + 1 / 2 == (1 + f1) / 2 == (f1 * 1) / 2 + 1 / 2 == (f1 * 1 + 1) / 2
    @test hash(f1) == hash(1 * f1) == hash(1 * f1 + 0)

    # Test canonical commutation relations
    @test ord_equals(f1' * f1 + f1 * f1', 1)
    @test ord_equals(f1 * f2 + f2 * f1, 0)
    @test ord_equals(f1' * f2 + f2 * f1', 0)

    @test_nowarn display(f1)
    @test_nowarn display(f3)
    @test_nowarn display(1 * f1)
    @test_nowarn display(2 * f3)
    @test_nowarn display(1 + f1)
    @test_nowarn display(1 + f3)
    @test_nowarn display(1 + a * f2 - 5 * f1 + 2 * z * f1 * f2)

    @test iszero(f1 - f1)
    @test ord_equals(f1 * f1, 0)
    @test f1 * f2 isa NonCommutativeProducts.NCMul
    @test iszero(2 * f1 - 2 * f1)
    @test iszero(0 * f1)
    @test 2 * f1 isa NonCommutativeProducts.NCMul
    @test iszero(f1 * 0)
    @test ord_equals(f1^2, 0)
    @test ord_equals(0 * (f1 + f2), 0)
    @test iszero((f1 + f2) * 0)
    @test ord_equals(f1 * f2 * f1, 0)
    f12 = f1 * f2
    @test iszero(f12'' - f12)
    @test ord_equals(f12 * f12, 0)
    @test ord_equals(f12' * f12', 0)
    nf1 = f1' * f1
    @test ord_equals(nf1^2, nf1)
    @test f1' * f1 isa NonCommutativeProducts.NCMul
    @test ord(f1 * f1') isa NonCommutativeProducts.NCAdd

    @test ord_equals(1 + (f1 + f2), 1 + f1 + f2, f1 + f2 + 1, f1 + 1 + f2, 1 * f1 + f2 + 1, f1 + 0.5 * f2 + 1 + (0 * f1 + 0.5 * f2), (0.5 + 0.5 * f1 + 0.2 * f2) + 0.5 + (0.5 * f1 + 0.8 * f2), (1 + f1' + (1 * f2)')')
    @test 1 + (f1 + f2) == 1 + f1 + f2 == f1 + f2 + 1 == f1 + 1 + f2 == 1 * f1 + f2 + 1 == f1 + 0.5 * f2 + 1 + (0 * f1 + 0.5 * f2) == (0.5 + 0.5 * f1 + 0.2 * f2) + 0.5 + (0.5 * f1 + 0.8 * f2) == (1 + f1' + (1 * f2)')'
    @test ord_equals((2 * f1) * (2 * f1), 0)
    @test ord_equals((2 * f1)^2, 0)
    @test ord_equals((2 * f2) * (2 * f1), -4 * f1 * f2)
    @test ord_equals(f1, (f1 * (f1 + 1)), (f1 + 1) * f1)
    @test ord_equals(f1 * (f1 + f2) * f1, 0)
    @test ord_equals((f1 * (f1 + f2)), f1 * f2)
    @test ord_equals((2nf1 - 1) * (2nf1 - 1), 1)

    @test ord_equals((1 * f1) * f2, f1 * f2)
    @test ord_equals((1 * f1) * (1 * f2), f1 * f2)
    @test ord_equals(f1 * f2, f1 * (1 * f2), f1 * f2)
    @test ord_equals(f1 - 1, (1 * f1) - 1, (0.5 + f1) - 1.5)

    for op in [f1, f1 + f2, f1' * f2, f1 * f2 + 1]
        @test ord_equals(op + 1.0I, op + 1.0)
        @test ord_equals(op - 1.0I, op - 1.0, -(1.0I - op), -(1.0 - op))
    end

    bigprod = prod(Fermion(rand(1:15), rand(Bool)) + Fermion(rand(1:15), rand(Bool)) for k in 1:10)
    @test ord(bigprod) == ord(ord(bigprod))
    bigprodhc = bigprod + bigprod'
    @test ord(bigprodhc) == ord(bigprodhc')
end


@testmodule Majoranas begin
    export Majorana
    include("rules/majoranas.jl")
end

@testitem "Majoranas" setup = [Majoranas] begin
    using Random: seed!
    seed!(1)
    NonCommutativeProducts.enable_autosort!()
    γ = Majorana.(1:4)
    @test 1 * γ[1] == γ[1]
    @test 1 * γ[1] + 0 == γ[1]
    @test 1 * γ[1] + 0 == 1 * γ[1]
    @test hash(γ[1]) == hash(1 * γ[1]) == hash(1 * γ[1] + 0)

    #test canonical anticommutation relations
    @test γ[1] * γ[1] == 1
    @test γ[1] * γ[2] == -γ[2] * γ[1]

    @test γ[1] * γ[2] * γ[1] == -γ[2]
    @test γ[1] * γ[2] * γ[3] == -γ[3] * γ[2] * γ[1]

    f1 = (γ[1] + 1im * γ[2]) / 2
    f2 = (γ[3] + 1im * γ[4]) / 2

    @test 1 * f1 == f1
    @test 1 * f1 + 0 == f1
    @test hash(f1) == hash(1 * f1) == hash(1 * f1 + 0)


    @test iszero(f1 - f1)
    @test iszero(f1 * f1)
    @test iszero(2 * f1 - 2 * f1)
    @test iszero(0 * f1)
    @test iszero(f1 * 0)
    @test iszero(f1^2)
    @test iszero(0 * (f1 + f2))
    @test iszero((f1 + f2) * 0)
    @test iszero(f1 * f2 * f1)
    f12 = f1 * f2
    @test iszero(f12'' - f12)
    @test iszero(f12 * f12)
    @test iszero(f12' * f12')
    nf1 = f1' * f1
    @test nf1^2 == nf1
    @test 1 + (f1 + f2) == 1 + f1 + f2 == f1 + f2 + 1 == f1 + 1 + f2 == 1 * f1 + f2 + 1 == f1 + 0.5 * f2 + 1 + (0 * f1 + 0.5 * f2) == (0.5 + 0.5 * f1 + 0.2 * f2) + 0.5 + (0.5 * f1 + 0.8 * f2) == (1 + f1' + (1 * f2)')'
    @test iszero((2 * f1) * (2 * f1))
    @test iszero((2 * f1)^2)
    @test (2 * f2) * (2 * f1) == -4 * f1 * f2
    @test f1 == (f1 * (f1 + 1)) == (f1 + 1) * f1
    @test iszero(f1 * (f1 + f2) * f1)
    @test (f1 * (f1 + f2)) == f1 * f2
    @test (2nf1 - 1) * (2nf1 - 1) == 1

    @test (1 * f1) * f2 == f1 * f2
    @test (1 * f1) * (1 * f2) == f1 * f2
    @test f1 * f2 == f1 * (1 * f2) == f1 * f2
    @test f1 - 1 == (1 * f1) - 1 == (0.5 + f1) - 1.5

    bigprod = prod(Majorana(rand(1:15)) for k in 1:20)
    @test NonCommutativeProducts.bubble_sort(bigprod) == bigprod
    bigprod2 = prod(Majorana(rand(1:15)) + Majorana(rand(1:15)) for k in 1:10)
    @test NonCommutativeProducts.bubble_sort(bigprod2) == bigprod2
    @test NonCommutativeProducts.bubble_sort(bigprod2') == bigprod2'
end


@testitem "Fermions+Bosons" setup = [Fermions, Bosons] begin
    import NonCommutativeProducts: bubble_sort, @nc, add!!
    using LinearAlgebra
    using Random: seed!
    NonCommutativeProducts.enable_autosort!()
    NonCommutativeProducts.@commutative Fermion Boson
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    b = Boson()

    @test b * b == b^2
    @test b * b' == b' * b + 1
    @test (b + b' + f1)^2 == b^2 + b'^2 + f1^2 + b * b' + b' * b + 2 * b * f1 + 2 * b' * f1

    @test 1 * f1 * b == b * f1
    @test 1 * f1 + b == f1 + b
    @test b * 1 * f1 + 0 == b * 1 * f1
    @test 0 * (f1 * b) == 0
    @test hash(f1 * b) == hash(1 * f1 * b) == hash(1 * b * f1 + 0)

    base = 1.0 + f1 + b
    for term in (f2, b', f2 * b, f1 * b')
        x = copy(base)
        y = add!!(x, 1im * term) # should not mutate
        @test y !== x
        @test y == base + 1im * term
        y = add!!(x, term) # should not mutate
        @test y === x == base + term
    end
end

@testitem "Zeros, promotion, conversion" setup = [Fermions] begin
    f1 = Fermion(:a)
    f1mul = 1 * f1
    f1add = 1 * f1 + 0
    @test zero(f1) == zero(f1mul) == zero(f1add) == 0
    @test zero(typeof(f1)) == zero(typeof(f1mul)) == zero(typeof(f1add)) == 0

    @test promote_rule(typeof(f1), typeof(f1)) == typeof(f1)
    @test promote_rule(typeof(f1), typeof(1 * f1)) == typeof(1 * f1)
    @test promote_rule(typeof(f1), typeof(f1 + 0)) == typeof(f1 + 0)

    @test promote_rule(typeof(f1mul), typeof(f1mul)) == typeof(f1mul)
    @test promote_rule(typeof(f1mul), typeof(1im * f1mul)) == typeof(1im * f1mul)

    @test promote_rule(typeof(f1mul), typeof(f1add)) == typeof(f1add)
    @test promote_rule(typeof(f1add), typeof(f1mul)) == typeof(f1add)
    @test promote_rule(typeof(f1add), typeof(f1add)) == typeof(f1add)

    @test isconcretetype(eltype([f1, f1 * 1]))
    @test isconcretetype(eltype([f1, f1 + 1]))
    @test isconcretetype(eltype([f1 * 2, f1 + 1]))
end

@testmodule WrappedRules begin
    using NonCommutativeProducts
    export Wrapped, Sym
    import NonCommutativeProducts: @nc, mul_effect, ncmap

    struct Wrapped{S}
        sym::S
    end
    Base.show(io::IO, x::Wrapped) = print(io, "W(", x.sym, ")")
    @nc Wrapped

    function NonCommutativeProducts.mul_effect(a::Wrapped, b::Wrapped)
        effect = mul_effect(a.sym, b.sym)
        effect isa Union{Number,Nothing} && return effect
        return ncmap(Wrapped, effect)
    end
    Base.adjoint(x::Wrapped) = Wrapped(adjoint(x.sym))
end

@testitem "Wrapped mul_effect delegation" setup = [WrappedRules, Fermions] begin
    import NonCommutativeProducts: ncmap

    NonCommutativeProducts.enable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    w1, w2 = Wrapped.((f1, f2))

    @test w1 * w2 == ncmap(Wrapped, f1 * f2)
    @test w2 * w2 == ncmap(Wrapped, f2 * f2)
    @test w2 * w1 == ncmap(Wrapped, f2 * f1)
    @test iszero(w1 * w1)
    @test w1' * w1 + w1 * w1' == 1
end

@run_package_tests verbose = true
