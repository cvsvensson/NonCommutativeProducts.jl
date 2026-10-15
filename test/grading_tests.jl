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
    export GFermion, BadFermion, Charge

    struct Charge
        q::Int
    end
    charge(q) = iszero(q) ? nothing : Charge(q)
    NC.grading(::Type{Charge}) = Charge
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

    # Like GFermion, except that distinct annihilators swap with an extra term c†*c† of charge +2 instead of -2. The
    # term is zero once sorted, so only a check that doesn't sort the effect finds it.
    struct HiddenFermion
        label::Int
        creation::Bool
    end
    HiddenFermion(label) = HiddenFermion(label, false)
    Base.adjoint(f::HiddenFermion) = HiddenFermion(f.label, !f.creation)
    NC.@nc HiddenFermion
    NC.grade(f::HiddenFermion) = Charge(f.creation ? 1 : -1)
    function NC.mul_effect(a::HiddenFermion, b::HiddenFermion)
        a == b && return 0
        if should_swap(a, b)
            a.label == b.label && return -b * a + 1
            !a.creation && !b.creation && return -b * a + a' * a'
            return -b * a
        end
        return nothing
    end
    export HiddenFermion

    # Spinful fermions graded by charge and by spin (in units of 1/2), as a ProductGrade
    struct Spin
        s::Int
    end
    spin(s) = iszero(s) ? nothing : Spin(s)
    NC.grading(::Type{Spin}) = Spin
    NC.mul_grade(a::Spin, b::Spin) = spin(a.s + b.s)
    NC.adjoint_grade(a::Spin) = spin(-a.s)

    struct SpinFermion
        label::Tuple{Int,Bool} # (site, spin up)
        creation::Bool
    end
    SpinFermion(site, up) = SpinFermion((site, up), false)
    Base.adjoint(f::SpinFermion) = SpinFermion(f.label, !f.creation)
    NC.@nc SpinFermion
    function NC.grade(f::SpinFermion)
        sign = f.creation ? 1 : -1
        return NC.ProductGrade(Charge(sign), Spin(f.label[2] ? sign : -sign))
    end
    function NC.mul_effect(a::SpinFermion, b::SpinFermion)
        a == b && return 0
        if should_swap(a, b)
            a.label == b.label && return -b * a + 1
            return -b * a
        end
        return nothing
    end
    export Spin, SpinFermion
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
    # zero sums hold no zero terms
    for z in zs
        z isa NC.NCAdd && @test isempty(z.dict)
    end
    @test isempty((ket(1) + ket(2) - ket(1) - ket(2)).dict)
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
            y = make_nonzero(y)
            @test y == kb
            # a nonzero sum stores no grades
            @test y.zero_grades === nothing
            y = make_zero(y)
            @test iszero(y) && NC.grades(y) == Set([KETBRA])
        end
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

@testitem "Grading: grade is the only element of grades" setup = [GradedKets, ChargedFermions] begin
    import NonCommutativeProducts as NC
    using LinearAlgebra
    NC.enable_autosort!()
    u = NC.NCAdd(0, Dict(1 * ket(1) => 0))
    f = GFermion(1)
    rows = Any[
        (2, Set([nothing])),
        (2I, Set([nothing])),
        (NC.NCMul(2, GKet[]), Set([nothing])),
        (0, Set()),
        (0I, Set()),
        (NC.NCMul(0, GKet[]), Set()),
        (u, Set()),
        (NC.ncmap(identity, ket(1) - ket(1)), Set()),
        (ket(1), Set([KET])),
        (bra(1), Set([BRA])),
        (f', Set([Charge(1)])),
        (f, Set([Charge(-1)])),
        (0 * ket(1), Set([KET])),
        (zero(ket(1)), Set([KET])),
        (ket(1) - ket(1), Set([KET])),
        (zero(2 * ket(1) * bra(2)), Set([KETBRA])),
        (ket(1) + ket(2) * bra(3), Set([KET, KETBRA])),
        (0 * (ket(1) + ket(2) * bra(3)), Set([KET, KETBRA])),
        (0 * (f' + f), Set([Charge(1), Charge(-1)])),
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

@testitem "Grading: zero of a product" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    # products with vector and tuple factors
    for q in (NC.NCMul(2, [ket(1), bra(2)]), NC.NCMul(2, (ket(1), bra(2))))
        @test typeof(zero(q)) == typeof(zero(typeof(q)))
        @test NC.grades(zero(q)) == Set([KETBRA])
    end
end

@testitem "Grading: zero products without factors have no grades" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    # like numbers: grades(2) == {nothing} but grades(0 * 2) == ∅, since NCMul(0, []) can't remember a history
    k = 1 * ket(1)
    @test NC.grades(2 * one(k)) == Set([nothing])
    @test isempty(NC.grades(0 * one(k)))
    @test_throws ArgumentError NC.grade(0 * one(k))
    # zeros from sorting are sums and keep the identity grade, and so does a sorted power
    @test NC.grades(bra(1) * ket(2)) == Set([nothing])
    @test NC.grades(0 * k^0) == Set([nothing])
    # an unsorted power is a product without factors
    NC.disable_autosort!()
    @test NC.grades(k^0) == Set([nothing])
    @test isempty(NC.grades(0 * k^0))
    NC.enable_autosort!()
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

    @test NC.grades(0 * x) == Set([KET])
    @test NC.grades(zerovector(x)) == Set([KET])
    @test NC.grades(zerovector(x, ComplexF64)) == Set([KET])
    @test NC.grades(zero(x)) == Set([KET])

    # a zero added to something else keeps the union
    z = x - x
    @test NC.grades(z + ket(3) * bra(1)) == Set([KETBRA])
    @test NC.grades(z + (ket(3) * bra(1) - ket(3) * bra(1))) == Set([KET, KETBRA])

    # anyadd, convert and copy
    @test NC.grades(NC.anyadd(z)) == Set([KET])
    @test NC.grades(convert(NC.NCAdd{ComplexF64,keytype(z.dict)}, z)) == Set([KET])
    @test NC.grades(copy(z)) == Set([KET])
end

@testitem "Grading: in-place operations" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using VectorInterface
    NC.enable_autosort!()
    x = 1.0 * ket(1) + 2.0 * ket(2)
    kb = 1.0 * ket(3) * bra(1)

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

@testitem "Grading: ProductGrade" setup = [GradedKets, ChargedFermions, GradeFamilies] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade

    # the order of the components doesn't matter
    @test ProductGrade(Charge(1), KET) == ProductGrade(KET, Charge(1))
    @test hash(ProductGrade(Charge(1), KET)) == hash(ProductGrade(KET, Charge(1)))
    @test Set{Any}([ProductGrade(Charge(1), KET), ProductGrade(KET, Charge(1))]) == Set{Any}([ProductGrade(KET, Charge(1))])
    @test ProductGrade(Charge(1), KET) != ProductGrade(Charge(1), BRA)
    @test startswith(repr(ProductGrade(Charge(1), KET)), "ProductGrade(")

    # grades of different gradings multiply to a ProductGrade, in either order
    @test NC._mul_grade(Charge(1), KET) == ProductGrade(Charge(1), KET)
    @test NC._mul_grade(KET, Charge(1)) == ProductGrade(Charge(1), KET)
    @test NC._mul_grade(ProductGrade(Charge(1), KET), Spin(1)) == ProductGrade(Spin(1), KET, Charge(1))
    @test NC._mul_grade(Spin(1), ProductGrade(Charge(1), KET)) == ProductGrade(Spin(1), KET, Charge(1))

    # components of the same grading are multiplied in order, since the bicyclic grade doesn't commute
    @test NC._mul_grade(ProductGrade(Charge(1), KET), ProductGrade(Charge(1), BRA)) == ProductGrade(Charge(2), KETBRA)
    @test NC._mul_grade(ProductGrade(Charge(1), BRA), ProductGrade(Charge(1), KET)) == Charge(2)
    @test NC._mul_grade(ProductGrade(KET, Charge(1)), BRA) == ProductGrade(Charge(1), KETBRA)
    @test NC._mul_grade(BRA, ProductGrade(KET, Charge(1))) == Charge(1)

    # identity components are dropped, and a single remaining component is unwrapped
    @test NC._mul_grade(NC._mul_grade(Charge(1), KET), Charge(-1)) == KET
    @test NC._mul_grade(Charge(1), NC._mul_grade(KET, Charge(-1))) == KET
    @test NC._mul_grade(ProductGrade(Charge(1), BRA), ProductGrade(Charge(-1), KET)) === nothing
    @test NC._mul_grade(ProductGrade(Charge(1), KET), ProductGrade(Charge(1), KET)) == ProductGrade(Charge(2), Bicyclic(2, 0))

    # nothing is the unit
    @test NC._mul_grade(nothing, ProductGrade(Charge(1), KET)) == ProductGrade(Charge(1), KET)
    @test NC._mul_grade(ProductGrade(Charge(1), KET), nothing) == ProductGrade(Charge(1), KET)
    @test NC._mul_grade(nothing, nothing) === nothing

    # associativity, also with grades of one grading spread over several types
    gs = (Charge(1), Charge(-1), Spin(1), KET, BRA, KETBRA, ProductGrade(Charge(1), KET), ProductGrade(Spin(-1), BRA),
        U1(1), U1(Int8(-1)), U1(Int8(2)), ProductGrade(U1(1), Charge(1)), ProductGrade(U1(Int8(-1)), KET))
    for a in gs, b in gs, c in gs
        @test isequal(NC._mul_grade(NC._mul_grade(a, b), c), NC._mul_grade(a, NC._mul_grade(b, c)))
    end

    # adjoints act componentwise
    @test NC._adjoint_grade(ProductGrade(Charge(1), KET)) == ProductGrade(Charge(-1), BRA)
    @test NC._adjoint_grade(ProductGrade(Charge(2), KETBRA)) == ProductGrade(Charge(-2), KETBRA)
    @test NC._adjoint_grade(ProductGrade(Charge(1), KETBRA, Spin(-1))) == ProductGrade(Charge(-1), KETBRA, Spin(1))
    @test NC._adjoint_grade(nothing) === nothing

    # a grade type without its own mul_grade throws a MethodError
    struct Unmultipliable end
    NC.grading(::Type{Unmultipliable}) = Unmultipliable
    @test_throws MethodError NC._mul_grade(Unmultipliable(), Unmultipliable())
    @test_throws MethodError NC._mul_grade(ProductGrade(Unmultipliable(), KET), Unmultipliable())
    @test NC._mul_grade(Unmultipliable(), KET) == ProductGrade(KET, Unmultipliable())
end

@testitem "Grading: ProductGrade construction, equality and printing" setup = [GradedKets, ChargedFermions, GradeFamilies] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade

    # invalid components
    @test_throws ArgumentError ProductGrade()
    @test_throws ArgumentError ProductGrade(Charge(1))
    @test_throws ArgumentError ProductGrade(Charge(1), nothing)
    @test_throws ArgumentError ProductGrade(nothing, Charge(1))
    @test_throws ArgumentError ProductGrade(ProductGrade(Charge(1), KET), Spin(1))
    @test_throws ArgumentError ProductGrade(Spin(1), ProductGrade(Charge(1), KET))
    # undeclared grading
    @test_throws ArgumentError ProductGrade(Charge(1), Undeclared(1))
    @test_throws ArgumentError ProductGrade(Undeclared(1), Charge(1))
    # several components of one grading, also of different types in a family
    @test_throws ArgumentError ProductGrade(Charge(1), Charge(2))
    @test_throws ArgumentError ProductGrade(Charge(1), KET, Charge(2))
    @test_throws ArgumentError ProductGrade(U1(1), U1(Int8(2)))
    @test_throws ArgumentError ProductGrade(U1(1), Charge(1), U1(1))

    # ==, isequal, hash and show don't depend on the order of the components
    cs = (Charge(1), KET, Spin(-1))
    perms = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
    ps = [ProductGrade(map(i -> cs[i], p)...) for p in perms]
    for p in ps, q in ps
        @test p == q
        @test isequal(p, q)
        @test hash(p) == hash(q)
        @test repr(p) == repr(q)
        @test sprint(show, p) == sprint(show, q)
    end
    @test length(Set{Any}(ps)) == 1
    @test length(Dict(p => i for (i, p) in enumerate(ps))) == 1
    s = repr(first(ps))
    @test startswith(s, "ProductGrade(")
    for c in ("Charge(1)", "Bicyclic(1, 0)", "Spin(-1)")
        @test occursin(c, s)
    end

    # components are compared by grading
    @test ProductGrade(Charge(1), KET) != ProductGrade(Charge(1), KET, Spin(1))
    @test ProductGrade(Charge(1), KET, Spin(1)) != ProductGrade(Charge(1), KET)
    @test ProductGrade(Charge(1), KET) != ProductGrade(Charge(1), Spin(1))
    @test ProductGrade(Charge(1), KET) != ProductGrade(Charge(2), KET)
    @test !isequal(ProductGrade(Charge(1), KET), ProductGrade(Spin(1), KET))
    # with the components' own ==, isequal and hash
    @test ProductGrade(U1(1), Charge(1)) == ProductGrade(Charge(1), U1(Int8(1)))
    @test isequal(ProductGrade(U1(1), Charge(1)), ProductGrade(Charge(1), U1(Int8(1))))
    @test hash(ProductGrade(U1(1), Charge(1))) == hash(ProductGrade(Charge(1), U1(Int8(1))))
end

@testitem "Grading: mixing an undeclared grading throws" setup = [GradedKets, ChargedFermions, GradeFamilies] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade
    error_of(f) =
        try
            f()
            nothing
        catch e
            e
        end
    function test_guiding_error(f)
        e = error_of(f)
        @test e isa ArgumentError
        e isa ArgumentError && @test occursin("NonCommutativeProducts.grading(::Type{", e.msg)
    end

    # grades of one undeclared type need no declaration
    @test NC._mul_grade(Undeclared(1), Undeclared(2)) == Undeclared(3)
    @test NC._mul_grade(nothing, Undeclared(2)) == Undeclared(2)
    @test NC._mul_grade(Undeclared(2), nothing) == Undeclared(2)

    # but mixing it with another grading does
    test_guiding_error(() -> NC._mul_grade(Charge(1), Undeclared(1)))
    test_guiding_error(() -> NC._mul_grade(Undeclared(1), Charge(1)))
    test_guiding_error(() -> NC._mul_grade(KET, Undeclared(1)))
    test_guiding_error(() -> NC._mul_grade(ProductGrade(Charge(1), KET), Undeclared(1)))
    test_guiding_error(() -> NC._mul_grade(Undeclared(1), ProductGrade(Charge(1), KET)))

    # two undeclared types
    struct OtherUndeclared end
    test_guiding_error(() -> NC._mul_grade(Undeclared(1), OtherUndeclared()))

    # through expressions: atoms graded by an undeclared type, next to charged fermions
    struct UAtom
        n::Int
    end
    NC.@nc UAtom
    NC.mul_effect(::UAtom, ::UAtom) = nothing
    NC.grade(a::UAtom) = Undeclared(a.n)
    NC.@commutative GFermion UAtom
    NC.enable_autosort!()
    @test NC.grade(UAtom(1) * UAtom(2)) == Undeclared(3)
    test_guiding_error(() -> NC.grade(GFermion(1) * UAtom(1)))
    test_guiding_error(() -> NC.grade(UAtom(1) * GFermion(1)'))
end

@testitem "Grading: families of grade types" setup = [GradedKets, ChargedFermions, GradeFamilies] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade

    # members of a family are multiplied with the user's mul_grade
    @test NC._mul_grade(U1(1), U1(Int8(2))) === U1(3)
    @test NC._mul_grade(U1(Int8(2)), U1(1)) === U1(3)
    @test NC._mul_grade(U1(Int8(1)), U1(Int8(2))) === U1(Int8(3))
    @test NC._mul_grade(U1(1), U1(Int8(-1))) === nothing

    # also inside a product
    @test NC._mul_grade(ProductGrade(U1(1), Charge(1)), U1(Int8(-1))) == Charge(1)
    @test NC._mul_grade(U1(Int8(-1)), ProductGrade(U1(1), Charge(1))) == Charge(1)
    @test NC._mul_grade(ProductGrade(U1(1), Charge(1)), U1(Int8(1))) == ProductGrade(U1(2), Charge(1))
    @test NC._mul_grade(ProductGrade(U1(1), Charge(1)), ProductGrade(Charge(-1), U1(Int8(1)))) == U1(2)
    @test NC._mul_grade(ProductGrade(U1(Int8(1)), BRA), ProductGrade(KET, U1(-1))) === nothing
    @test NC._mul_grade(U1(1), KET) == ProductGrade(KET, U1(1))
    @test NC._adjoint_grade(ProductGrade(U1(Int8(2)), Charge(1))) == ProductGrade(U1(Int8(-2)), Charge(-1))

    # a family without a method for differing types gives a MethodError
    struct Parity{T}
        odd::T
    end
    NC.grading(::Type{<:Parity}) = Parity
    NC.mul_grade(a::Parity{T}, b::Parity{T}) where {T} = a.odd == b.odd ? nothing : Parity(one(T))
    @test NC._mul_grade(Parity(true), Parity(false)) == Parity(true)
    @test NC._mul_grade(Parity(true), Parity(true)) === nothing
    @test_throws MethodError NC._mul_grade(Parity(true), Parity(1))
    @test_throws MethodError NC._mul_grade(ProductGrade(Parity(true), Charge(1)), Parity(1))
    @test_throws MethodError NC._mul_grade(Parity(1), ProductGrade(Parity(true), Charge(1)))
    @test NC._mul_grade(Parity(true), Charge(1)) == ProductGrade(Charge(1), Parity(true))
end

@testitem "Grading: methods across gradings are never called" setup = [ChargedFermions] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade

    struct GA
        n::Int
    end
    struct GB
        n::Int
    end
    NC.grading(::Type{GA}) = GA
    NC.grading(::Type{GB}) = GB
    NC.mul_grade(a::GA, b::GA) = iszero(a.n + b.n) ? nothing : GA(a.n + b.n)
    NC.mul_grade(a::GB, b::GB) = iszero(a.n + b.n) ? nothing : GB(a.n + b.n)
    NC.mul_grade(::GA, ::GB) = error("should not be called")
    NC.mul_grade(::GB, ::GA) = error("should not be called")

    @test NC._mul_grade(GA(1), GB(1)) == ProductGrade(GA(1), GB(1))
    @test NC._mul_grade(GB(1), GA(1)) == ProductGrade(GA(1), GB(1))
    @test NC._mul_grade(ProductGrade(GA(1), Charge(1)), GB(1)) == ProductGrade(GA(1), GB(1), Charge(1))
    @test NC._mul_grade(GB(1), ProductGrade(GA(1), Charge(1))) == ProductGrade(GA(1), GB(1), Charge(1))
    @test NC._mul_grade(ProductGrade(GA(1), GB(1)), GA(-1)) == GB(1)

    # consistent at the top level and inside a product
    gs = (GA(1), GA(-1), GB(1), GB(2), Charge(1), ProductGrade(GA(1), Charge(-1)), ProductGrade(GB(-1), GA(2)))
    for a in gs, b in gs, c in gs
        @test isequal(NC._mul_grade(NC._mul_grade(a, b), c), NC._mul_grade(a, NC._mul_grade(b, c)))
    end
end

@testitem "Grading: closure of mul_grade and adjoint_grade is checked" setup = [ChargedFermions] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade
    error_of(f) =
        try
            f()
            nothing
        catch e
            e
        end

    # mul_grade and adjoint_grade return a grade of another grading
    struct Leaky
        n::Int
    end
    struct Target
        n::Int
    end
    NC.grading(::Type{Leaky}) = Leaky
    NC.grading(::Type{Target}) = Target
    NC.mul_grade(a::Leaky, b::Leaky) = Target(a.n + b.n)
    NC.adjoint_grade(a::Leaky) = Target(-a.n)

    e = error_of(() -> NC._mul_grade(ProductGrade(Leaky(1), Charge(1)), Leaky(1)))
    @test e isa ArgumentError
    e isa ArgumentError && @test occursin("mul_grade", e.msg)
    @test_throws ArgumentError NC._mul_grade(Leaky(1), ProductGrade(Leaky(1), Charge(1)))
    @test_throws ArgumentError NC._mul_grade(ProductGrade(Leaky(1), Charge(1)), ProductGrade(Leaky(1), Spin(1)))
    e = error_of(() -> NC._adjoint_grade(ProductGrade(Leaky(1), Charge(1))))
    @test e isa ArgumentError
    e isa ArgumentError && @test occursin("adjoint_grade", e.msg)

    # mul_grade and adjoint_grade return something without a grading
    struct Wrong
        n::Int
    end
    NC.grading(::Type{Wrong}) = Wrong
    NC.mul_grade(a::Wrong, b::Wrong) = a.n + b.n
    NC.adjoint_grade(a::Wrong) = -a.n
    e = error_of(() -> NC._mul_grade(ProductGrade(Wrong(1), Charge(1)), Wrong(1)))
    @test e isa ArgumentError
    e isa ArgumentError && @test occursin("mul_grade", e.msg)
    e = error_of(() -> NC._adjoint_grade(ProductGrade(Wrong(1), Charge(1))))
    @test e isa ArgumentError
    e isa ArgumentError && @test occursin("adjoint_grade", e.msg)
end

@testitem "Grading: ProductGrade for spinful fermions" setup = [ChargedFermions] begin
    import NonCommutativeProducts as NC
    using NonCommutativeProducts: ProductGrade
    NC.enable_autosort!()
    up, dn = SpinFermion(1, true), SpinFermion(1, false)
    @test NC.grade(up') == ProductGrade(Charge(1), Spin(1))
    @test NC.grade(dn) == ProductGrade(Charge(-1), Spin(1))
    @test NC.grade(up' * dn) == Spin(2)
    @test NC.grade(up' * dn') == Charge(2)
    @test NC.grade(up' * up) === nothing
    @test NC.grades(up' * up + up' * dn + dn) == Set([nothing, Spin(2), ProductGrade(Charge(-1), Spin(1))])

    # zeros keep their grades, also through adjoints
    z = up * up
    @test z == 0
    @test NC.grades(z) == Set([ProductGrade(Charge(-2), Spin(-2))])
    @test NC.grades(z') == Set([ProductGrade(Charge(2), Spin(2))])
    @test NC.grades(up' * dn * dn) == Set([ProductGrade(Charge(-1), Spin(3))])

    fs = (up, dn, SpinFermion(2, true))
    for a in (fs..., adjoint.(fs)...), b in (fs..., adjoint.(fs)...)
        @test NC.check_grading(a, b)
    end
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
    @test NC.grade(k * f' * bra(1)) == ProductGrade(Charge(1), KETBRA)
    @test NC.grades(f * k + k) == Set([KET, ProductGrade(Charge(-1), KET)])

    # zeros keep product grades, also through adjoints
    z = f * k * f
    @test z == 0
    @test NC.grades(z) == Set([ProductGrade(Charge(-2), KET)])
    @test NC.grades(z') == Set([ProductGrade(Charge(2), BRA)])
    @test NC.grades(bra(2) * f' * k) == Set([Charge(1)])

    for a in (f, f', GFermion(2)), b in (k, bra(1), ket(2))
        @test NC.check_grading(a, b)
        @test NC.check_grading(b, a)
    end
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

    # a swap has the grade of the swapped product, which differs from the original for a non-commutative grading
    @test !NC.check_grading(SwapKet(1, false), SwapKet(2, true))
    @test NC.check_grading(SwapKet(1, true), SwapKet(2, false))
end

@testitem "Grading: check_grading checks the effect as the sorter sees it" setup = [ChargedFermions] begin
    import NonCommutativeProducts as NC
    h1, h2 = HiddenFermion(1), HiddenFermion(2)
    for autosort in (false, true)
        autosort ? NC.enable_autosort!() : NC.disable_autosort!()
        @test !NC.check_grading(h2, h1)
        @test NC.check_grading(h2', h1')
        @test NC.check_grading(h1, h1')
    end
    NC.enable_autosort!()
end

@testitem "Grading: ungraded types" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    NC.enable_autosort!()
    f1, f2 = Fermion(1), Fermion(2)

    @test NC.grade(f1) === nothing
    @test NC.grades(f1' * f1) == Set([nothing])
    @test NC.grades(f1' * f1 + f2) == Set([nothing])
    @test NC.grades(f1' * f1 + 1) == Set([nothing])
    # a zero of ungraded factors has the identity grade
    @test NC.grades(f1 * f1) == Set([nothing])
    @test NC.grades(f1 - f1) == Set([nothing])

    # the grade of an ungraded product is inferred, so it costs nothing
    @test (@inferred NC.grade(NC.NCMul(1, [f1', f1]))) === nothing
    @test (@inferred NC.grade(NC.NCMul(1, (f1', f1)))) === nothing
    @test (@inferred NC.grade(NC.NCMul(0, [f1', f1]))) === nothing
    @test (@inferred NC.grade(NC.NCMul(2, typeof(f1)[]))) === nothing
end

@testitem "Grading: powers of zero products" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    for autosort in (true, false)
        autosort ? NC.enable_autosort!() : NC.disable_autosort!()
        # products with vector and tuple factors
        for z in (0 * ket(1), NC.NCMul(0, (ket(1),)))
            @test iszero(z^2)
            @test NC.grades(z^2) == Set([Bicyclic(2, 0)])
        end
        @test NC.grades((0 * ket(1) * bra(2))^3) == Set([KETBRA])
    end
    NC.enable_autosort!()
end

@testitem "Grading: operations on atoms" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    using VectorInterface
    NC.enable_autosort!()
    z = bra(1) - bra(1)
    @test NC.grades(NC.add!!(copy(z), ket(1), 0, 1)) == Set([BRA, KET])
    @test NC.grades(NC.add!!(copy(z), ket(1), 0, VectorInterface.Zero())) == Set([KET])
    @test NC.grades(scale!!(copy(z), ket(1), 0)) == Set([KET])
    @test NC.grades(add!(copy(z), ket(1), 0)) == Set([BRA, KET])
    @test NC.grades(VectorInterface.add(ket(1), bra(2), 0, VectorInterface.Zero())) == Set([BRA])
    @test NC.grades(VectorInterface.add(ket(1), bra(2), 0, 0)) == Set([KET, BRA])
    @test NC.grades(zerovector(ket(1))) == Set([KET])
    @test NC.grades(zerovector(bra(1), Float64)) == Set([BRA])
end

@testitem "Grading: autosort disabled" setup = [GradedKets] begin
    import NonCommutativeProducts as NC
    NC.disable_autosort!()
    x = 1.0 * ket(1) + 2.0 * ket(2)
    z = x - x
    # without sorting, products of nonzero terms stay nonzero, so zeros come from the coefficients
    @test NC.grades(z * bra(3)) == Set([KETBRA])
    @test NC.grades((0 * bra(3)) * x) == Set([nothing])
    @test NC.grades(z * (x' + 1)) == Set([KETBRA, KET])
    @test NC.grades(z') == Set([BRA])
    @test NC.grades((0 * ket(1) * bra(2))') == Set([KETBRA])
    c = x' - x'
    @test NC.grades(NC.mul!!(c, z, 1 * bra(3))) == Set([BRA, KETBRA])
    NC.enable_autosort!()
end
