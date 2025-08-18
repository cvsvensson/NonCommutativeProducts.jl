using TestItemRunner

@testmodule Fermions begin
    export Fermion, Ordering
    using NonCommutativeProducts
    import NonCommutativeProducts: AddTerms, Swap, @nc
    struct Fermion{L}
        label::L
        creation::Bool
    end
    Base.adjoint(x::Fermion) = Fermion(x.label, !x.creation)
    Fermion(k) = Fermion(k, false)
    # Base.iszero(x::Fermion) = false
    Base.show(io::IO, x::Fermion) = print(io, "c", x.creation ? "†" : "", "[", x.label, "]")
    Base.:(==)(a::Fermion, b::Fermion) = a.creation == b.creation && a.label == b.label
    @nc Fermion

    struct Ordering end
    function should_swap(a::Fermion, b::Fermion, ::Ordering)
        if a.creation == b.creation
            return a.label > b.label
        else
            a.creation < b.creation
        end
    end

    function NonCommutativeProducts.mul_effect(a::Fermion, b::Fermion, ordering::Ordering)
        a == b && return 0
        if should_swap(a, b, ordering)
            a.label == b.label && xor(a.creation, b.creation) && return AddTerms((Swap(-1), 1))
            return Swap(-1)
        else
            return nothing
        end
    end
end

@testitem "Fermions" setup = [Fermions] begin
    import NonCommutativeProducts: filter_zeros!, bubble_sort
    using Symbolics, LinearAlgebra

    @variables a::Real z::Complex
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    f3 = Fermion((1, :↑))
    ord(op) = bubble_sort(op, Fermions.Ordering())
    # ord(op) = bubble_sort(op, Ordering())
    ord_equals(a, b) = (iszero(filter_zeros!(ord(a - b))))
    ord_equals(a, b, c, xs...) = (ord_equals(a, b) && ord_equals(b, c, xs...))

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
end


@testmodule Majoranas begin
    using NonCommutativeProducts
    export Majorana
    import NonCommutativeProducts: AddTerms, Swap, mul_effect, @nc_eager
    struct Majorana{L}
        label::L
    end
    Base.:(==)(a::Majorana, b::Majorana) = a.label == b.label
    Base.adjoint(x::Majorana) = Majorana(x.label)
    Base.show(io::IO, x::Majorana) = print(io, "γ[", x.label, "]")
    struct Ordering end
    @nc_eager Majorana Ordering()
    function NonCommutativeProducts.mul_effect(a::Majorana, b::Majorana, ::Ordering)
        if a == b
            return 1
        elseif a.label < b.label
            return nothing
        elseif a.label > b.label
            a.label == b.label && return AddTerms(Swap(-1), 1)
            return Swap(-1)
        else
            throw(ArgumentError("Don't know how to multiply $a * $b"))
        end
    end
end

@testitem "Majoranas" setup = [Majoranas] begin
    using Symbolics
    @variables a::Real z::Complex
    γ = Majorana.(1:4)

    #test canonical anticommutation relations
    @test γ[1] * γ[1] == 1
    @test γ[1] * γ[2] == -γ[2] * γ[1]

    @test γ[1] * γ[2] * γ[1] == -γ[2]
    @test γ[1] * γ[2] * γ[3] == -γ[3] * γ[2] * γ[1]

    f1 = (γ[1] + 1im * γ[2]) / 2
    f2 = (γ[3] + 1im * γ[4]) / 2
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
    # @test ord_equals(1 + (f1 + f2), 1 + f1 + f2, f1 + f2 + 1, f1 + 1 + f2, 1 * f1 + f2 + 1, f1 + 0.5 * f2 + 1 + (0 * f1 + 0.5 * f2), (0.5 + 0.5 * f1 + 0.2 * f2) + 0.5 + (0.5 * f1 + 0.8 * f2), (1 + f1' + (1 * f2)')')
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
end
@run_package_tests
