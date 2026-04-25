@testmodule Paulis begin
    export Pauli
    using NonCommutativeProducts
    import NonCommutativeProducts: Swap, @nc

    """
        Pauli(label, idx)

    Pauli matrix operator for spin-1/2 systems.
    - `label`: identifies the site/particle
    - `idx`: 1 for σx, 2 for σy, 3 for σz
    """
    struct Pauli{L}
        label::L
        idx::Int
        Pauli(label::L, idx::Int=1) where L = new{L}(label, idx)
    end

    Base.adjoint(x::Pauli) = x  # Hermitian
    Pauli(idx::Int) = Pauli(0, idx)

    function Base.show(io::IO, x::Pauli)
        syms = ("x", "y", "z")
        base = "σ" * syms[x.idx]
        print(io, x.label == 0 ? base : "$base[$(x.label)]")
    end

    @nc Pauli

    function should_swap(a::Pauli, b::Pauli)
        a.label != b.label ? a.label > b.label : a.idx > b.idx
    end

    function NonCommutativeProducts.mul_effect(a::Pauli, b::Pauli)
        if a.label != b.label
            # Different sites commute
            return should_swap(a, b) ? Swap(1) : nothing
        elseif a.idx == b.idx
            # σᵢ² = 1
            return 1
        elseif should_swap(a, b)
            # σᵢ σⱼ = -σⱼ σᵢ for i ≠ j (anti-commute)
            return Swap(-1)
        else
            return nothing
        end
    end
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
    using NonCommutativeProducts
    import NonCommutativeProducts: AddTerms, Swap, @nc

    """
        SpinOp(label, op)

    General spin operator using raising/lowering basis.
    - `label`: identifies the site/particle
    - `op`: :z for Sz, :+ for S+, :- for S-
    """
    struct SpinOp{L}
        label::L
        op::Symbol
    end

    function Base.adjoint(x::SpinOp)
        new_op = x.op == :+ ? :- : x.op == :- ? :+ : :z
        return SpinOp(x.label, new_op)
    end

    SpinOp(label) = SpinOp(label, :z)

    function Base.show(io::IO, x::SpinOp)
        print(io, "S", x.op, "[$(x.label)]")
    end

    @nc SpinOp

    const OP_ORDER = Dict(:z => 1, :+ => 2, :- => 3)

    function should_swap(a::SpinOp, b::SpinOp)
        a.label != b.label ? a.label > b.label : OP_ORDER[a.op] > OP_ORDER[b.op]
    end

    function NonCommutativeProducts.mul_effect(a::SpinOp, b::SpinOp)
        if a.label != b.label
            # Different sites commute
            return should_swap(a, b) ? Swap(1) : nothing
        end

        # Commutation relations: [Sz, S+] = S+, [Sz, S-] = -S-, [S+, S-] = 2Sz
        # Canonical order: :z < :+ < :-
        if a.op == :+ && b.op == :z
            # S+ Sz = Sz S+ - S+
            return AddTerms((Swap(1), -SpinOp(a.label, :+)))
        elseif a.op == :- && b.op == :z
            # S- Sz = Sz S- + S-
            return AddTerms((Swap(1), SpinOp(a.label, :-)))
        elseif a.op == :- && b.op == :+
            # S- S+ = S+ S- - 2Sz
            return AddTerms((Swap(1), -2 * SpinOp(a.label, :z)))
        end

        return nothing
    end
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
