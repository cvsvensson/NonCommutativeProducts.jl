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
end

@testmodule ChargedFermions begin
    # Fermions with a U(1) charge grade: creation +1, annihilation -1, and charge 0 normalized to `nothing`.
    import NonCommutativeProducts as NC
    export GFermion, BadFermion, Charge

    struct Charge
        q::Int
    end
    charge(q) = iszero(q) ? nothing : Charge(q)
    NC.mul_grade(a::Charge, b::Charge) = charge(a.q + b.q)
    NC.adjoint_grade(c::Charge) = charge(-c.q)

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

@testitem "Grading: numbers and atoms" setup = [GradedKets, ChargedFermions] begin
    import NonCommutativeProducts as NC
    @test NC.grades(2) == Set([nothing])
    @test isempty(NC.grades(0))
    @test NC.grades(ket(1)) == Set([KET])
    @test NC.grades(bra(1)) == Set([BRA])
    @test NC.grades(GFermion(1)') == Set([Charge(1)])
    @test NC.grades(GFermion(1)) == Set([Charge(-1)])
    @test NC.grade(2) === nothing
    @test_throws ArgumentError NC.grade(0)
end

@testitem "Grading: all zeros are equal" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using VectorInterface
    NC.enable_autosort!()
    u = NC.NCAdd(0, Dict(1 * ket(1) => 0))
    # -0.0 is left out, since isequal(-0.0, 0.0) is false in Base
    zs = [0, 0.0, 0 * ket(1), 0 * bra(1), 0 * ket(1) * bra(2), NC.NCMul(0, GKet[]), ket(1) - ket(1), zero(ket(1)),
        zero(1.0 * ket(1) * bra(2)), u, zerovector(ket(1))]
    for a in zs, b in zs
        @test a == b
        @test isequal(a, b)
        @test hash(a) == hash(b)
    end
    @test length(Dict(z => i for (i, z) in enumerate(zs))) == 1
    @test length(Dict(0 * ket(1) => 1, 0 * bra(1) => 2, 0 => 3)) == 1
    for z in zs
        z isa NC.MulAdd || continue
        @test NC.isscalar(z)
        @test NC.scalar(z) == 0
    end

    # nonzero values still compare by their factors
    @test 2 * ket(1) != 2 * bra(1)
    @test NC.NCMul(2, GKet[]) == 2
end

@testitem "Grading: zeros don't keep grades from earlier values" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using VectorInterface
    NC.enable_autosort!()
    satisfies_Z(x) = !(x isa NC.NCAdd) || iszero(x) || x.zero_grades === nothing
    x = 1.0 * ket(1) + 2.0 * ket(2)
    kb = 1.0 * ket(3) * bra(1)

    # each row starts from a zero ket, makes it the nonzero kb and then zero again
    rows = (
        ("add!!", y -> NC.add!!(y, kb), y -> NC.add!!(y, kb, -1)),
        ("VectorInterface.add!", y -> add!(y, kb), y -> add!(y, kb, -1)),
        ("mul!!", y -> NC.mul!!(y, 1.0 * ket(3), 1 * bra(1)), y -> NC.mul!!(y, -1.0 * ket(3), 1 * bra(1))),
        ("scale!(y, 0)", y -> add!(y, kb), y -> scale!(y, 0)),
        ("scale!!(y, 0)", y -> NC.add!!(y, kb), y -> scale!!(y, 0)),
        ("VectorInterface.scale!(y, kb, 0)", y -> add!(y, kb), y -> VectorInterface.scale!(y, kb, 0)),
        ("scale!!(y, kb, 0)", y -> NC.add!!(y, kb), y -> scale!!(y, kb, 0)),
        ("zerovector!", y -> add!(y, kb), y -> zerovector!(y)),
    )
    for (name, make_nonzero, make_zero) in rows
        @testset "$name" begin
            y = x - x
            @test iszero(y) && NC.grades(y) == Set([KET])
            @test satisfies_Z(y)
            y = make_nonzero(y)
            @test y == kb
            @test satisfies_Z(y)
            y = make_zero(y)
            @test satisfies_Z(y)
            @test iszero(y) && NC.grades(y) == Set([KETBRA])
        end
    end

    # results of public operations satisfy invariant Z
    for r in (x + kb, x * bra(1), x', sort(x), copy(x), NC.anyadd(x), 2x, (x - x) + kb)
        @test satisfies_Z(r)
    end
end

@testitem "Grading: bare zeros are the unit" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using LinearAlgebra
    NC.enable_autosort!()
    x = 1.0 * ket(1) + 2.0 * ket(2)
    u = NC.NCAdd(0, Dict(1 * ket(1) => 0))
    mixed = ket(1) + ket(2) * bra(3)
    for b in Any[0, 0.0, NC.NCMul(0, GKet[]), u, zero(typeof(x)), zero(typeof(1 * ket(1)))]
        @test isempty(NC.grades(b))
        @test NC.grades(b + (x - x)) == Set([KET])
        @test NC.grades(b * ket(1)) == NC.grades(ket(1) * b) == Set([KET])
        @test NC.grades(b * mixed) == Set([KET, KETBRA])
        @test isempty(NC.grades(b'))
    end
    # UniformScaling * NCMul is not defined, so only grades and +
    @test isempty(NC.grades(0I))
    @test NC.grades(0I + (x - x)) == Set([KET])
end

@testitem "Grading: grade is the only element of grades" setup = [GradedKets, Fermions] begin
    import NonCommutativeProducts as NC
    using LinearAlgebra
    NC.enable_autosort!()
    u = NC.NCAdd(0, Dict(1 * ket(1) => 0))
    rows = Any[
        (2, Set([nothing])),
        (2I, Set([nothing])),
        (NC.NCMul(2, GKet[]), Set([nothing])),
        (0, Set()),
        (0I, Set()),
        (NC.NCMul(0, GKet[]), Set()),
        (u, Set()),
        (ket(1), Set([KET])),
        (0 * ket(1), Set([KET])),
        (zero(ket(1)), Set([KET])),
        (ket(1) - ket(1), Set([KET])),
        (zero(2 * ket(1) * bra(2)), Set([KETBRA])),
        (ket(1) + ket(2) * bra(3), Set([KET, KETBRA])),
        (0 * (ket(1) + ket(2) * bra(3)), Set([KET, KETBRA])),
    ]
    for (x, expected) in rows
        @test NC.grades(x) == expected
        if length(expected) == 1
            @test isequal(NC.grade(x), only(expected))
        else
            @test_throws ArgumentError NC.grade(x)
        end
    end

    # the grade of an ungraded product is inferred, also for a zero product and one without factors
    f1 = Fermion(1)
    @test (@inferred NC.grade(NC.NCMul(1, [f1', f1]))) === nothing
    @test (@inferred NC.grade(NC.NCMul(1, (f1', f1)))) === nothing
    @test (@inferred NC.grade(NC.NCMul(0, [f1', f1]))) === nothing
    @test (@inferred NC.grade(NC.NCMul(2, typeof(f1)[]))) === nothing
end

@testitem "Grading: zero of a product or atom" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    @test NC.grades(zero(ket(1))) == Set([KET])
    @test zero(ket(1)) isa typeof(zero(GKet))
    @test iszero(zero(ket(1)))
    @test zero(ket(1)) == 0

    p = NC.NCMul(2, [ket(1), bra(2)])
    @test p.factors isa Vector
    pt = NC.NCMul(2, (ket(1), bra(2)))
    @test pt.factors isa Tuple
    for q in (p, pt)
        @test typeof(zero(q)) == typeof(zero(typeof(q)))
        @test NC.grades(zero(q)) == Set([KETBRA])
        @test iszero(zero(q))
        @test zero(q) == 0
    end
end

@testitem "Grading: bicyclic kets and bras" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    op = ket(0) * bra(1)
    @test NC.grade(op) == KETBRA
    @test NC.grades(ket(0) * bra(1) * ket(1)) == Set([KET])
    @test NC.grades(bra(1) * ket(1)) == Set([nothing])
    # ⟨1|2⟩ = 0 is a zero scalar, which has the identity grade
    @test bra(1) * ket(2) == 0
    @test NC.grades(bra(1) * ket(2)) == Set([nothing])

    # annihilated kets and bras are zero kets and bras
    z = op * ket(0)
    @test NC.grades(z) == Set([KET])
    @test NC.grades(bra(0) * op) == Set([BRA])
    @test NC.grades(z') == Set([BRA])

    # the value is plain zero, with no phantom terms
    @test z == 0
    @test hash(z) == hash(0)
    @test iszero(z)
    @test NC.isscalar(z)
    @test isempty(z.dict)

    # mixed sums have several grades, and the single query throws
    x = ket(1) + ket(2) * bra(3)
    @test NC.grades(x) == Set([KET, KETBRA])
    @test_throws ArgumentError NC.grade(x)
    @test NC.grades(ket(1) + 2 * ket(2)) == Set([KET])
    @test NC.grades(op + 1) == Set([KETBRA, nothing])
end

@testitem "Grading: independent of bracketing" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    z = ket(0) * bra(1) * ket(0)
    @test NC.grades((z * bra(1)) * ket(1)) == Set([KET])
    @test NC.grades(z * (bra(1) * ket(1))) == Set([KET])
    @test NC.grades(z * bra(1)) == Set([KETBRA])

    op = ket(0) * bra(1)
    @test NC.grades((0 * op) * ket(1)) == Set([KET])
    @test NC.grades(op * (0 * ket(1))) == Set([KET])
    @test NC.grades((op + 2 * ket(3) * bra(1)) * (0 * ket(1))) == Set([KET])
end

@testitem "Grading: sorting keeps the grade of zero" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.disable_autosort!()
    u = ket(0) * bra(1) * ket(0)
    @test !iszero(u)
    @test iszero(sort(u))
    @test NC.grades(sort(u)) == Set([KET])
    @test NC.grades(sort(u + 2 * u)) == Set([KET])
    NC.enable_autosort!()
end

@testitem "Grading: linear operations keep the grade of zero" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using VectorInterface
    NC.enable_autosort!()
    x = 1.0 * ket(1) + 2.0 * ket(2)

    # scaling by 0 and cancellation
    @test NC.grades(0 * ket(1)) == Set([KET])
    @test NC.grades(0 * x) == Set([KET])
    @test NC.grades(ket(1) - ket(1)) == Set([KET])
    @test NC.grades(x - x) == Set([KET])

    # zerovector and zero
    @test NC.grades(zerovector(x)) == Set([KET])
    @test NC.grades(zerovector(x, ComplexF64)) == Set([KET])
    @test NC.grades(zero(x)) == Set([KET])

    # in place
    y = add!(copy(x), x, -1)
    @test iszero(y) && NC.grades(y) == Set([KET])
    y = scale!(copy(x), 0)
    @test iszero(y) && NC.grades(y) == Set([KET])
    y = zerovector!(copy(x))
    @test iszero(y) && NC.grades(y) == Set([KET])
    @test NC.grades(scale!!(copy(x), 0)) == Set([KET])
    @test NC.grades(add!!(copy(x), x, -1)) == Set([KET])

    # a zero added to something else keeps the union
    z = x - x
    @test NC.grades(z + ket(3) * bra(1)) == Set([KETBRA])
    @test NC.grades(z + (ket(3) * bra(1) - ket(3) * bra(1))) == Set([KET, KETBRA])

    # anyadd, convert and copy
    @test NC.grades(NC.anyadd(z)) == Set([KET])
    @test NC.grades(convert(NC.NCAdd{ComplexF64,keytype(z.dict)}, z)) == Set([KET])
    @test NC.grades(copy(z)) == Set([KET])
end

@testitem "Grading: in-place results don't keep stale grades" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using VectorInterface
    NC.enable_autosort!()
    x = 1.0 * ket(1) + 2.0 * ket(2)
    kb = 1.0 * ket(3) * bra(1)

    # zero ket → nonzero ketbra → zero: only the operands since the last nonzero value count
    y = x - x
    add!(y, kb)
    @test NC.grades(y) == Set([KETBRA])
    add!(y, kb, -1)
    @test iszero(y)
    @test NC.grades(y) == Set([KETBRA])

    # an overwritten destination doesn't contribute its grades
    y = x' - x'
    @test NC.grades(scale!(y, x, 0)) == Set([KET])
    y = x' - x'
    @test NC.grades(add!(y, x, 0, VectorInterface.Zero())) == Set([KET])
    y = x' - x'
    @test NC.grades(add!(y, x, 0)) == Set([KET, BRA])

    # in-place products, c + a * b
    c = x' - x'
    @test NC.grades(NC.mul!!(c, kb, 0 * ket(1))) == Set([KET, BRA])
    y = NC.mul!!(1.0 * kb, -1.0 * ket(3) * bra(2), ket(2) * bra(1))
    @test iszero(y) && NC.grades(y) == Set([KETBRA])

    # the grades of a zero can't be changed through the returned set
    z = x - x
    push!(NC.grades(z), BRA)
    @test NC.grades(z) == Set([KET])
end

@testitem "Grading: mixed and unknown zeros" setup = [GradedKets, ChargedFermions] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()

    # a zero computed from a mixed sum has all its grades
    z = 0 * (ket(1) + ket(2) * bra(3))
    @test NC.grades(z) == Set([KET, KETBRA])
    @test_throws ArgumentError NC.grade(z)
    f = GFermion(1)
    @test NC.grades(0 * (f' + f)) == Set([Charge(1), Charge(-1)])

    # zeros of unknown origin have no grades, and the single query throws
    @test isempty(NC.grades(NC.ncmap(identity, ket(1) - ket(1))))
    u = NC.NCAdd(0, Dict(1 * ket(1) => 0))
    @test iszero(u)
    @test isempty(NC.grades(u))
    @test_throws ArgumentError NC.grade(u)

    # equality and hashing ignore the grades
    zk, zb = ket(1) - ket(1), bra(1) - bra(1)
    @test NC.grades(zk) != NC.grades(zb)
    @test zk == zb
    @test isequal(zk, zb)
    @test hash(zk) == hash(zb)
    @test length(Dict(zk => 1, zb => 2)) == 1
end

@testitem "Grading: grade_part" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    kb = ket(2) * bra(3)

    # a mixed sum splits into its graded parts
    x = ket(1) + 2 * ket(4) + kb + 3
    @test NC.grade_part(x, KET) == ket(1) + 2 * ket(4)
    @test NC.grade_part(x, KETBRA) == kb
    @test NC.grade_part(x, nothing) == 3
    @test sum(NC.grade_part(x, g) for g in NC.grades(x)) == x
    @test x == ket(1) + 2 * ket(4) + kb + 3 # not mutated

    # atoms and products are kept or replaced by a zero
    @test NC.grade_part(ket(1), KET) === ket(1)
    @test iszero(NC.grade_part(ket(1), BRA))
    @test iszero(NC.grade_part(kb, KET))
    @test isempty(NC.grades(NC.grade_part(kb, KET)))

    # a grade that isn't present gives a zero without grades
    z = NC.grade_part(x, BRA)
    @test iszero(z) && z == 0
    @test isempty(NC.grades(z))

    # zeros keep the requested grade if they have it
    z0 = 0 * (ket(1) + kb)
    @test NC.grades(z0) == Set([KET, KETBRA])
    zk = NC.grade_part(z0, KET)
    @test iszero(zk) && zk == 0
    @test NC.grades(zk) == Set([KET])
    @test isempty(NC.grades(NC.grade_part(z0, BRA)))
    @test NC.grades(z0) == Set([KET, KETBRA])

    # numbers have the identity grade
    @test NC.grade_part(2, nothing) == 2
    @test iszero(NC.grade_part(2, KET))
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
    @test NC.grades(f[1]' * f[1] * f[1]) == Set([Charge(-1)])
    @test NC.grades((f[1] + f[2])^2) == Set([Charge(-2)])

    # (c† + c)^2 = 1: annihilation removes the charges ±2, so the result is charge-neutral
    x = (f[1]' + f[1])^2
    @test x == 1
    @test NC.grades(x) == Set([nothing])
end

@testitem "Grading: check_grading" setup = [GradedKets, ChargedFermions] begin
    import NonCommutativeProducts as NC
    f1, f2 = GFermion(1), GFermion(2)
    for (a, b) in ((f1, f1), (f1', f1'), (f1, f1'), (f1', f1), (f2, f1), (f2', f1), (f1, f2'))
        @test NC.check_grading(a, b)
    end
    for (a, b) in ((bra(1), ket(1)), (bra(1), ket(2)), (ket(1), bra(2)), (ket(1), ket(2)))
        @test NC.check_grading(a, b)
    end

    b1, b2 = BadFermion(1), BadFermion(2)
    @test NC.check_grading(b2, b1)
    @test NC.check_grading(b2', b1)
    @test !NC.check_grading(b1, b1)
    @test !NC.check_grading(b1, b1')
end

@testitem "Grading: ungraded types" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    f1, f2 = Fermion(1), Fermion(2)

    @test NC.grade(f1) === nothing
    @test NC.grades(f1' * f1) == Set([nothing])
    @test NC.grades(f1' * f1 + f2) == Set([nothing])
    @test NC.grades(f1' * f1 + 1) == Set([nothing])
    # the grades of an ungraded zero may be skipped, so they are either unknown or the identity
    @test NC.grades(f1 * f1) ⊆ Set([nothing])
    @test NC.grades(f1 - f1) ⊆ Set([nothing])
end
