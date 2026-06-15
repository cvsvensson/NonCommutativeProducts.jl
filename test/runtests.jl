using TestItemRunner
@run_package_tests verbose = true

@testmodule Fermions begin
    export Fermion
    include("rules/fermions.jl")
end

@testmodule Majoranas begin
    export Majorana
    include("rules/majoranas.jl")
end

@testmodule Bosons begin
    export Boson
    include("rules/bosons.jl")
end

@testmodule Paulis begin
    export Pauli
    include("rules/paulis.jl")
end

@testmodule SpinOperators begin
    export SpinOp
    include("rules/spin.jl")
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

@testitem "Autosort local override with ScopedValue" setup = [Majoranas] begin
    NonCommutativeProducts.enable_autosort!()
    γ1, γ2 = Majorana.(1:2)

    # Global autosort enabled: product is immediately canonicalized.
    canonical = γ2 * γ1
    @test canonical == -γ1 * γ2

    # Local override: multiplication is left unsorted only in this dynamic scope.
    local_unsorted = Base.ScopedValues.with(NonCommutativeProducts._autosort => false) do
        @test !NonCommutativeProducts.autosort()
        γ2 * γ1
    end
    @test local_unsorted != canonical
    @test NonCommutativeProducts.bubble_sort(local_unsorted) == canonical

    # Global setting remains unchanged outside the local override.
    @test NonCommutativeProducts.autosort()
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

    @test one(f1) == one(f1mul) == one(f1add) == 1
    @test oneunit(f1) == oneunit(f1mul) == oneunit(f1add) == 1
    @test one(typeof(f1)) == one(typeof(f1mul)) == one(typeof(f1add)) == 1

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

@testitem "VectorInterface extension" setup = [Fermions] begin
    using VectorInterface

    NonCommutativeProducts.disable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)

    x = 1 + 2 * f1 + 3 * f2
    y = 4 - f1

    @test VectorInterface.scalartype(typeof(x)) == Int

    z = VectorInterface.zerovector(x)
    @test z == 0

    z2 = copy(x)
    @test VectorInterface.zerovector!(z2) === z2
    @test z2 == 0

    z3 = copy(x)
    @test VectorInterface.zerovector!!(z3) === z3
    @test z3 == 0

    @test VectorInterface.scale(x, 2) == 2 * x
    s = copy(x)
    @test VectorInterface.scale!(s, 2) === s
    @test s == 2 * x
    @test VectorInterface.scale!!(copy(x), 2) == 2 * x

    @test VectorInterface.add(y, x, 2, 3) == 3 * y + 2 * x
    a = copy(y)
    @test VectorInterface.add!(a, x, 2, 3) === a
    @test a == 3 * y + 2 * x
    @test VectorInterface.add!!(copy(y), x, 2, 3) == 3 * y + 2 * x
end

@testitem "KrylovKit compatibility" setup = [Fermions] begin
    using LinearAlgebra, KrylovKit
    using VectorInterface

    NonCommutativeProducts.enable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    s1 = Fermions.State(false, f1)
    s2 = Fermions.State(false, f2)
    @test inner(s1, s1) == 1
    @test inner(s2, s2) == 1
    @test_throws ArgumentError inner(s1, s2)
    @test inner(f1' * s1, f1' * s1) == 1
    @test inner(s1, f1' * s1) == 0

    ham = f1' * f1 + f1 + f1'
    basis = [s1, f1' * s1]
    hammat = [inner(b, ham * a) for b in basis, a in basis]
    numvals, numvecs = eigen(hammat)
    diagbasis = numvecs' * basis
    diaghammat = [inner(b, ham * a) for b in diagbasis, a in diagbasis]
    @test diaghammat ≈ Diagonal(numvals)

    _vals, _vecs, _ = KrylovKit.eigsolve(v -> ham * v, s1 + f1' * s1, 2; ishermitian=true)
    perm = sortperm(_vals)
    vals = _vals[perm]
    vecs = _vecs[perm]
    @test vals ≈ numvals
    @test map(abs ∘ inner, diagbasis, vecs) ≈ [1, 1]

    numlinsol = hammat \ [1, 0]
    sol, _ = KrylovKit.linsolve(v -> ham * v, basis[1])
    @test map(b -> inner(b, sol), basis) ≈ numlinsol


    sol2, _ = KrylovKit.lssolve((v -> ham * v, v -> ham' * v), basis[1])
    @test inner(sol, sol2) ≈ inner(sol, sol)

    svals, svecsl, svecsr, _ = KrylovKit.svdsolve((v -> ham * v, v -> ham' * v), 0 + basis[1])
    svecsl2, svals2, svecsr2 = svd(hammat)
    @test svals ≈ svals2
    @test abs(dot(map(b -> inner(b, svecsl[1]), basis), svecsl2[:, 1])) ≈ 1
    @test abs(dot(map(b -> inner(b, svecsl[2]), basis), svecsl2[:, 2])) ≈ 1
    @test abs(dot(map(b -> inner(b, svecsr[1]), basis), svecsr2[:, 1])) ≈ 1
    @test abs(dot(map(b -> inner(b, svecsr[2]), basis), svecsr2[:, 2])) ≈ 1

    ## try exponentiate

end

@testitem "KrylovKit apply extension" setup = [Fermions] begin
    using KrylovKit
    import NonCommutativeProducts: add!!

    NonCommutativeProducts.enable_autosort!()
    f = Fermion(:a)
    s = Fermions.State(false, f)
    x = s + f' * s
    op = f' * f + f + f'

    @test KrylovKit.apply(op, x) == op * x
    @test KrylovKit.apply(op, x, 2, 3) == add!!(op * x, x, 2, 3)
end

@testitem "Regression: add!! and mul!! weighted accumulation" setup = [Fermions] begin
    import NonCommutativeProducts: add!!, scale!, mul!!

    NonCommutativeProducts.disable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)

    # scale! should scale both term coefficients and additive coefficient.
    a0 = f1 + 1
    a1 = copy(a0)
    @test scale!(a1, 5) == 5 * a0

    # add!!(a, b::NCMul, α, β) should scale every term in a by β, not just overlapping keys.
    a2 = f1 + f2 + 1
    got_mul = add!!(copy(a2), f1, 2, 5)
    expected_mul = 5 * a2 + 2 * f1
    @test got_mul == expected_mul

    # add!!(a, b::NCAdd, α, β) should also scale all terms in a by β.
    a3 = f1 + f2 + 1
    b3 = f1 + 2
    got_add = add!!(copy(a3), b3, 2, 5)
    expected_add = 5 * a3 + 2 * b3
    @test got_add == expected_add

    # mul!! relies on weighted add!! in accumulation and should match distributive expansion.
    a4 = f1 + 1
    b4 = f2 + 3
    got_prod = mul!!(zero(a4), a4, b4)
    expected_prod = a4 * b4
    @test got_prod == expected_prod
end

