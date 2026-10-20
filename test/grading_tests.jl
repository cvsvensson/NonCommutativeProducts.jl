@testmodule GradedKets begin
    # Kets and bras with the bicyclic grade k^m b^n, stored as (m, n). The product is
    # (m,n)(p,q) = (m-n+max(n,p), q-p+max(n,p)), and the identity (0,0) is normalized to `nothing`.
    import NonCommutativeProducts as NC
    export GKet, ket, bra, Bicyclic, KET, BRA, KETBRA

    struct Bicyclic
        m::Int
        n::Int
    end
    bicyclic(m, n) = (m, n) == (0, 0) ? nothing : Bicyclic(m, n)
    NC.grading(::Type{Bicyclic}) = Bicyclic
    function NC.mul_grade(g::Bicyclic, h::Bicyclic)
        t = max(g.n, h.m)
        bicyclic(g.m - g.n + t, h.n - h.m + t)
    end
    NC.adjoint_grade(g::Bicyclic) = bicyclic(g.n, g.m)
    const KET = Bicyclic(1, 0)
    const BRA = Bicyclic(0, 1)
    const KETBRA = Bicyclic(1, 1)

    struct GKet
        n::Int
        bra::Bool
    end
    ket(n) = GKet(n, false)
    bra(n) = GKet(n, true)
    Base.adjoint(k::GKet) = GKet(k.n, !k.bra)
    Base.show(io::IO, k::GKet) = print(io, k.bra ? "⟨$(k.n)|" : "|$(k.n)⟩")
    NC.@nc GKet
    NC.mul_effect(a::GKet, b::GKet) = a.bra && !b.bra ? Int(a.n == b.n) : nothing
    NC.grade(k::GKet) = k.bra ? BRA : KET

    # Graded like GKet, but |n⟩⟨m| swaps to ⟨m||n⟩, which changes the grade from KETBRA to the identity
    struct SwapKet
        n::Int
        bra::Bool
    end
    NC.@nc SwapKet
    NC.mul_effect(a::SwapKet, b::SwapKet) = !a.bra && b.bra ? NC.Swap(1) : nothing
    NC.grade(k::SwapKet) = k.bra ? BRA : KET
    export SwapKet
end

@testmodule ChargedFermions begin
    # Fermions with a U(1) charge grade: creation +1, annihilation -1, and charge 0 normalized to `nothing`.
    import NonCommutativeProducts as NC
    export GFermion, BadFermion, Charge, Spin

    struct Charge
        q::Int
    end
    charge(q) = iszero(q) ? nothing : Charge(q)
    NC.grading(::Type{Charge}) = Charge
    NC.mul_grade(a::Charge, b::Charge) = charge(a.q + b.q)
    NC.adjoint_grade(c::Charge) = charge(-c.q)

    # A second abelian grading, for products of gradings
    struct Spin
        s::Int
    end
    spin(s) = iszero(s) ? nothing : Spin(s)
    NC.grading(::Type{Spin}) = Spin
    NC.mul_grade(a::Spin, b::Spin) = spin(a.s + b.s)
    NC.adjoint_grade(a::Spin) = spin(-a.s)

    struct GFermion
        label::Int
        creation::Bool
    end
    GFermion(label) = GFermion(label, false)
    Base.adjoint(f::GFermion) = GFermion(f.label, !f.creation)
    Base.show(io::IO, f::GFermion) = print(io, "c", f.creation ? "†" : "", "[", f.label, "]")
    NC.@nc GFermion
    NC.grade(f::GFermion) = Charge(f.creation ? 1 : -1)

    should_swap(a, b) = a.creation == b.creation ? a.label > b.label : a.creation < b.creation
    function NC.mul_effect(a::GFermion, b::GFermion)
        a == b && return 0
        if should_swap(a, b)
            a.label == b.label && return -b * a + 1
            return -b * a
        end
        return nothing
    end

    # Same grade, but with rules that break it: a*a = 1 has charge 0 instead of ∓2, and c*c† = -c†*c + c† adds a
    # term of charge +1 instead of 0. Distinct labels swap correctly.
    struct BadFermion
        label::Int
        creation::Bool
    end
    BadFermion(label) = BadFermion(label, false)
    Base.adjoint(f::BadFermion) = BadFermion(f.label, !f.creation)
    NC.@nc BadFermion
    NC.grade(f::BadFermion) = Charge(f.creation ? 1 : -1)
    function NC.mul_effect(a::BadFermion, b::BadFermion)
        a == b && return 1
        if should_swap(a, b)
            a.label == b.label && return NC.AddTerms((NC.Swap(-1), b))
            return NC.Swap(-1)
        end
        return nothing
    end
end

@testmodule GradeFamilies begin
    import NonCommutativeProducts as NC
    export U1, Undeclared

    # A U(1) grading spread over several types, which share the key U1. Charges of different types promote.
    struct U1{T}
        q::T
    end
    NC.grading(::Type{<:U1}) = U1
    NC.mul_grade(a::U1, b::U1) = iszero(a.q + b.q) ? nothing : U1(a.q + b.q)
    NC.adjoint_grade(a::U1) = U1(-a.q)
    # equal charges are equal grades, whatever their type
    Base.:(==)(a::U1, b::U1) = a.q == b.q
    Base.isequal(a::U1, b::U1) = isequal(a.q, b.q)
    Base.hash(a::U1, h::UInt) = hash(a.q, hash(:U1, h))

    # A grade type without a declared grading, which can only be multiplied with itself
    struct Undeclared
        n::Int
    end
    NC.mul_grade(a::Undeclared, b::Undeclared) = Undeclared(a.n + b.n)
    NC.adjoint_grade(a::Undeclared) = Undeclared(-a.n)
end

@testitem "Grading: all zeros are equal" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    # -0.0 is left out, since isequal(-0.0, 0.0) is false in Base
    zs = [0, 0.0, 0 * ket(1), 0 * ket(1) * bra(2), NC.NCMul(0, GKet[]), ket(1) - ket(1), zero(ket(1))]
    for a in zs, b in zs
        @test a == b
        @test isequal(a, b)
        @test hash(a) == hash(b)
    end
end

@testitem "Grading: grade and grades" setup = [GradedKets, ChargedFermions] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    f = GFermion(1)
    rows = Any[
        (2, Set([nothing])),
        (0, Set()),
        (ket(1), Set([KET])),
        (bra(1), Set([BRA])),
        (f', Set([Charge(1)])),
        (ket(1) * bra(2), Set([KETBRA])),
        (ket(0) * bra(1) * ket(1), Set([KET])),
        (bra(1) * ket(1), Set([nothing])),
        (ket(1) + 2 * ket(2), Set([KET])),
        (ket(1) + ket(2) * bra(3), Set([KET, KETBRA])),
        # zeros remember their grades
        (0 * ket(1), Set([KET])),
        (ket(1) - ket(1), Set([KET])),
        (0 * (ket(1) + ket(2) * bra(3)), Set([KET, KETBRA])),
        # ⟨1|2⟩ = 0 is a zero scalar, which has the identity grade
        (bra(1) * ket(2), Set([nothing])),
    ]
    for (x, expected) in rows
        @test NC.grades(x) == expected
        if length(expected) == 1
            @test isequal(NC.grade(x), only(expected))
        else
            @test_throws ArgumentError NC.grade(x)
        end
    end
end

@testitem "Grading: zeros from products keep their grade" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    op = ket(0) * bra(1)
    z = op * ket(0)
    @test z == 0
    @test NC.grades(z) == Set([KET])
    @test NC.grades(bra(0) * op) == Set([BRA])
    @test NC.grades(z') == Set([BRA])

    # independent of bracketing
    @test NC.grades((z * bra(1)) * ket(1)) == Set([KET])
    @test NC.grades(z * (bra(1) * ket(1))) == Set([KET])
    @test NC.grades((0 * op) * ket(1)) == Set([KET])
    @test NC.grades(op * (0 * ket(1))) == Set([KET])

    # also when sorting explicitly
    NC.disable_autosort!()
    u = ket(0) * bra(1) * ket(0)
    @test !iszero(u)
    @test iszero(sort(u))
    @test NC.grades(sort(u)) == Set([KET])
    NC.enable_autosort!()
end

@testitem "Grading: linear operations keep the grade of zero" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using VectorInterface
    NC.enable_autosort!()
    x = 1.0 * ket(1) + 2.0 * ket(2)
    kb = 1.0 * ket(3) * bra(1)

    @test NC.grades(0 * x) == Set([KET])
    @test NC.grades(zero(x)) == Set([KET])
    @test NC.grades(zerovector(x)) == Set([KET])

    # adding zeros takes the union of grades
    z = x - x
    @test NC.grades(z + kb) == Set([KETBRA])
    @test NC.grades(z + (kb - kb)) == Set([KET, KETBRA])

    # bare zeros carry no grades and act as the unit
    @test NC.grades(0 + z) == Set([KET])
    @test NC.grades(0 * ket(1)) == Set([KET])

    # an overwritten destination doesn't contribute its grades, an added-to one does
    y = x' - x'
    @test NC.grades(scale!(y, x, 0)) == Set([KET])
    y = x' - x'
    @test NC.grades(add!(y, x, 0)) == Set([KET, BRA])

    # zeros don't keep grades from earlier values
    y = NC.add!!(x - x, kb)
    @test y == kb
    y = NC.add!!(y, kb, -1)
    @test iszero(y) && NC.grades(y) == Set([KETBRA])
end

@testitem "Grading: grade_part" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using LinearAlgebra: I
    NC.enable_autosort!()
    kb = ket(2) * bra(3)

    # a mixed sum splits into its graded parts
    x = ket(1) + 2 * ket(4) + kb + 3
    @test NC.grade_part(x, KET) == ket(1) + 2 * ket(4)
    @test NC.grade_part(x, KETBRA) == kb
    @test NC.grade_part(x, nothing) == 3
    @test sum(NC.grade_part(x, g) for g in NC.grades(x)) == x

    # a grade that isn't present gives a zero without grades
    z = NC.grade_part(x, BRA)
    @test iszero(z)
    @test isempty(NC.grades(z))

    # zeros keep the requested grade if they have it
    z0 = 0 * (ket(1) + kb)
    @test NC.grades(NC.grade_part(z0, KET)) == Set([KET])
    @test isempty(NC.grades(NC.grade_part(z0, BRA)))

    # a product has only its own grade, and scalars and uniform scalings have only the grade nothing
    @test NC.grade_part(kb, KETBRA) == kb
    @test iszero(NC.grade_part(kb, KET))
    @test NC.grade_part(3, nothing) == 3
    @test iszero(NC.grade_part(3, KET))
    @test iszero(NC.grade_part(2I, KET))
end

@testitem "Grading: fermion charge" setup = [ChargedFermions] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    f = GFermion.(1:3)

    @test NC.grade(f[1]' * f[1]) === nothing
    @test NC.grades(f[1]' * f[2]') == Set([Charge(2)])

    # a number-conserving Hamiltonian
    H = sum(f[i]' * f[i+1] + f[i+1]' * f[i] for i in 1:2) + 0.5 * f[1]' * f[1] * f[2]' * f[2]
    @test NC.grades(H) == Set([nothing])
    @test NC.grades(H * H) == Set([nothing])
    @test NC.grades(H + f[1]) == Set([nothing, Charge(-1)])

    # zeros remember their charge
    @test f[1] * f[1] == 0
    @test NC.grades(f[1] * f[1]) == Set([Charge(-2)])
    @test NC.grades((f[1] + f[2])^2) == Set([Charge(-2)])

    # (c† + c)^2 = 1: annihilation removes the charges ±2, so the result is charge-neutral
    x = (f[1]' + f[1])^2
    @test x == 1
    @test NC.grades(x) == Set([nothing])
end

@testitem "Grading: ungraded types" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    f1, f2 = Fermion(1), Fermion(2)

    @test NC.grade(f1) === nothing
    @test NC.grades(f1' * f1 + f2) == Set([nothing])
    # a zero of ungraded factors has the identity grade
    @test NC.grades(f1 * f1) == Set([nothing])

    # the grade of an ungraded product is inferred, so it costs nothing
    @test (@inferred NC.grade(NC.NCMul(1, [f1', f1]))) === nothing
end

@testitem "Grading: ProductGrade" setup = [GradedKets, ChargedFermions, GradeFamilies] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade

    # the order of the components doesn't matter
    @test ProductGrade(Charge(1), KET) == ProductGrade(KET, Charge(1))
    @test hash(ProductGrade(Charge(1), KET)) == hash(ProductGrade(KET, Charge(1)))
    @test ProductGrade(Charge(1), KET) != ProductGrade(Charge(1), BRA)
    @test ProductGrade(Charge(1), KET) != ProductGrade(Charge(1), KET, Spin(1))

    # invalid components
    @test_throws ArgumentError ProductGrade(Charge(1))
    @test_throws ArgumentError ProductGrade(Charge(1), nothing)
    @test_throws ArgumentError ProductGrade(ProductGrade(Charge(1), KET), Spin(1))
    @test_throws ArgumentError ProductGrade(Charge(1), Charge(2))
    @test_throws ArgumentError ProductGrade(U1(1), U1(Int8(2)))
    @test_throws ArgumentError ProductGrade(Charge(1), Undeclared(1))

    # grades of different gradings multiply to a ProductGrade, in either order
    @test NC._mul_grade(Charge(1), KET) == ProductGrade(Charge(1), KET)
    @test NC._mul_grade(KET, Charge(1)) == ProductGrade(Charge(1), KET)

    # components of the same grading are multiplied in order, since the bicyclic grade doesn't commute
    @test NC._mul_grade(ProductGrade(Charge(1), KET), ProductGrade(Charge(1), BRA)) == ProductGrade(Charge(2), KETBRA)
    @test NC._mul_grade(ProductGrade(Charge(1), BRA), ProductGrade(Charge(1), KET)) == Charge(2)

    # identity components are dropped, and a single remaining component is unwrapped
    @test NC._mul_grade(NC._mul_grade(Charge(1), KET), Charge(-1)) == KET
    @test NC._mul_grade(ProductGrade(Charge(1), BRA), ProductGrade(Charge(-1), KET)) === nothing

    # associativity, also with grades of one grading spread over several types
    gs = (Charge(1), Charge(-1), Spin(1), KET, BRA, KETBRA, ProductGrade(Charge(1), KET), ProductGrade(Spin(-1), BRA),
        U1(1), U1(Int8(-1)), U1(Int8(2)), ProductGrade(U1(1), Charge(1)), ProductGrade(U1(Int8(-1)), KET))
    for a in gs, b in gs, c in gs
        @test isequal(NC._mul_grade(NC._mul_grade(a, b), c), NC._mul_grade(a, NC._mul_grade(b, c)))
    end

    # adjoints act componentwise
    @test NC._adjoint_grade(ProductGrade(Charge(1), KET)) == ProductGrade(Charge(-1), BRA)

    # an undeclared grading can only be multiplied with itself
    @test NC._mul_grade(Undeclared(1), Undeclared(2)) == Undeclared(3)
    @test_throws ArgumentError NC._mul_grade(Charge(1), Undeclared(1))
end

@testitem "Grading: ProductGrade for mixed types" setup = [GradedKets, ChargedFermions] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade
    NC.@commutative GFermion GKet
    NC.enable_autosort!()
    f, k = GFermion(1), ket(1)

    # commuting factors of different types have the same grade in either order
    @test NC.grade(f * k) == NC.grade(k * f) == ProductGrade(Charge(-1), KET)
    @test NC.grade(f' * k * f) == KET
    @test NC.grades(f * k + k) == Set([KET, ProductGrade(Charge(-1), KET)])

    # zeros keep product grades, also through adjoints
    z = f * k * f
    @test z == 0
    @test NC.grades(z) == Set([ProductGrade(Charge(-2), KET)])
    @test NC.grades(z') == Set([ProductGrade(Charge(2), BRA)])

    for a in (f, f', GFermion(2)), b in (k, bra(1), ket(2))
        @test NC.check_grading(a, b)
        @test NC.check_grading(b, a)
    end
end

@testitem "Grading: check_grading" setup = [GradedKets, ChargedFermions] begin
    import NonCommutativeProducts as NC
    f1, f2 = GFermion(1), GFermion(2)
    for (a, b) in ((f1, f1), (f1, f1'), (f1', f1), (f2, f1), (f1, f2'))
        @test NC.check_grading(a, b)
    end
    for (a, b) in ((bra(1), ket(1)), (bra(1), ket(2)), (ket(1), bra(2)))
        @test NC.check_grading(a, b)
    end

    b1, b2 = BadFermion(1), BadFermion(2)
    @test NC.check_grading(b2, b1)
    @test !NC.check_grading(b1, b1)
    @test !NC.check_grading(b1, b1')

    # a swap has the grade of the swapped product, which differs from the original for a non-commutative grading
    @test !NC.check_grading(SwapKet(1, false), SwapKet(2, true))
    @test NC.check_grading(SwapKet(1, true), SwapKet(2, false))
end

@testitem "Grading: autosort disabled" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.disable_autosort!()
    x = 1.0 * ket(1) + 2.0 * ket(2)
    z = x - x
    # without sorting, products of nonzero terms stay nonzero, so zeros come from the coefficients
    @test NC.grades(z * bra(3)) == Set([KETBRA])
    @test NC.grades(z * (x' + 1)) == Set([KETBRA, KET])
    @test NC.grades(z') == Set([BRA])
    NC.enable_autosort!()
end

@testitem "Grading: unregistered objects" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    struct NotAnAtom end

    @test_throws ArgumentError NC.grade(NotAnAtom())
    @test_throws ArgumentError NC.grade("text")
    @test NC.grade(Fermion(1)) === nothing
end
