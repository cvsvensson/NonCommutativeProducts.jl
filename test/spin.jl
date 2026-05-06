@testmodule Paulis begin
    export Pauli
    include("rules/paulis.jl")
end

@testitem "Paulis" setup = [Paulis] begin
    import NonCommutativeProducts: bubble_sort
    using Symbolics

    @variables a::Real
    σx = Pauli(1)
    σy = Pauli(2)
    σz = Pauli(3)

    ord(op) = bubble_sort(op)
    ord_equals(x, y) = iszero(ord(x - y))

    # Basic algebraic properties
    @test ord_equals(σx * σx, 1)
    @test ord_equals(σy * σy, 1)
    @test ord_equals(σz * σz, 1)

    # Anti-commutation: σᵢ σⱼ = -σⱼ σᵢ
    @test ord_equals(σx * σy + σy * σx, 0)
    @test ord_equals(σy * σz + σz * σy, 0)
    @test ord_equals(σz * σx + σx * σz, 0)

    # Multi-site operators commute
    σx1 = Pauli(:A, 1)
    σy1 = Pauli(:A, 2)
    σx2 = Pauli(:B, 1)
    σy2 = Pauli(:B, 2)

    @test ord_equals(σx1 * σx2, σx2 * σx1)
    @test ord_equals(σx1 * σy2, σy2 * σx1)
    @test ord_equals(σx1 * σy1 + σy1 * σx1, 0)

    # Display tests
    @test_nowarn display(σx)
    @test_nowarn display(σx1)
    @test_nowarn display(2 * σx)
    @test_nowarn display(σx * σy)

    # Coefficient handling
    @test 1 * σx == σx
    @test 1 * σx + 0 == σx
    @test hash(σx) == hash(1 * σx)

    # Products
    @test ord_equals((2 * σx) * (2 * σx), 4)
    @test ord_equals(σx * σy * σx, -σy)
    @test ord_equals(σx * σy * σz, σz * σx * σy)

    # Addition
    @test ord_equals(σx + σx, 2 * σx)
    @test iszero(σx - σx)

    # Addition and subtraction between NCMul and Pauli
    @test (σx + σy * σz) isa NonCommutativeProducts.NCAdd
    @test (σx - σy * σz) isa NonCommutativeProducts.NCAdd
    @test (σy * σz + σx) isa NonCommutativeProducts.NCAdd
    @test (σy * σz - σx) isa NonCommutativeProducts.NCAdd
end


@testmodule SpinOperators begin
    export SpinOp
    include("rules/spin.jl")
end

@testitem "SpinOperators" setup = [SpinOperators] begin
    import NonCommutativeProducts: bubble_sort
    using Symbolics

    Sz = SpinOp(0, :z)
    Sp = SpinOp(0, :+)
    Sm = SpinOp(0, :-)

    ord(op) = bubble_sort(op)
    ord_equals(x, y) = iszero(ord(x - y))

    # Commutation relations
    # [Sz, S+] = S+
    @test ord_equals(Sz * Sp - Sp * Sz, Sp)
    # [Sz, S-] = -S-
    @test ord_equals(Sz * Sm - Sm * Sz, -Sm)
    # [S+, S-] = 2Sz
    @test ord_equals(Sp * Sm - Sm * Sp, 2 * Sz)

    # Multi-site operators commute
    Sz1 = SpinOp(:A, :z)
    Sp1 = SpinOp(:A, :+)
    Sz2 = SpinOp(:B, :z)
    Sp2 = SpinOp(:B, :+)

    @test ord_equals(Sz1 * Sz2, Sz2 * Sz1)
    @test ord_equals(Sz1 * Sp2, Sp2 * Sz1)

    # Basic properties
    @test 1 * Sz == Sz
    @test iszero(Sz - Sz)
    @test ord_equals(Sz + Sz, 2 * Sz)
end
