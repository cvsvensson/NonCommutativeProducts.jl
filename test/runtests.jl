using TestItemRunner
@run_package_tests verbose = true

@testmodule Fermions begin
    export Fermion, State
    include("rules/fermions.jl")
end

@testmodule Majoranas begin
    export Majorana
    include("rules/majoranas.jl")
end

@testmodule Bosons begin
    export Boson, State
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
    using Symbolics, LinearAlgebra
    using Random: seed!
    seed!(1)
    @variables a::Real z::Complex
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    f3 = Fermion((1, :↑))
    NonCommutativeProducts.disable_autosort!()
    ord(op) = sort(op)
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
    @test sort(bigprod) == bigprod
    bigprod2 = prod(Majorana(rand(1:15)) + Majorana(rand(1:15)) for k in 1:10)
    @test sort(bigprod2) == bigprod2
    @test sort(bigprod2') == bigprod2'
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
    @test sort(local_unsorted) == canonical

    # Global setting remains unchanged outside the local override.
    @test NonCommutativeProducts.autosort()
end

@testitem "with_autosort overrides the global setting" setup = [Majoranas] begin
    γ1, γ2 = Majorana.(1:2)
    was_enabled = NonCommutativeProducts.autosort()
    try
        NonCommutativeProducts.disable_autosort!()
        @test !NonCommutativeProducts.autosort()

        # The default argument enables autosort inside the block.
        sorted_product = NonCommutativeProducts.with_autosort() do
            @test NonCommutativeProducts.autosort()
            γ2 * γ1
        end
        @test sorted_product == -γ1 * γ2
        @test !NonCommutativeProducts.autosort()

        NonCommutativeProducts.enable_autosort!()
        unsorted_product = NonCommutativeProducts.with_autosort(false) do
            @test !NonCommutativeProducts.autosort()
            γ2 * γ1
        end
        @test unsorted_product != sorted_product
        @test sort(unsorted_product) == sorted_product

        # The global setting is restored after the block.
        @test NonCommutativeProducts.autosort()
        NonCommutativeProducts.with_autosort(true) do
            @test NonCommutativeProducts.autosort()
        end
        @test NonCommutativeProducts.autosort()
    finally
        was_enabled ? NonCommutativeProducts.enable_autosort!() : NonCommutativeProducts.disable_autosort!()
    end
end


@testitem "Fermions+Bosons" setup = [Fermions, Bosons] begin
    import NonCommutativeProducts: @nc, add!!
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

    @test promote_type(typeof(f1), typeof(f1)) == typeof(f1)
    @test promote_type(typeof(f1), typeof(1 * f1)) == typeof(1 * f1)
    @test promote_type(typeof(f1), typeof(f1 + 0)) == typeof(f1 + 0)

    @test promote_type(typeof(f1mul), typeof(f1mul)) == typeof(f1mul)
    @test promote_type(typeof(f1mul), typeof(1im * f1mul)) == typeof(1im * f1mul)

    @test promote_type(typeof(f1mul), typeof(f1add)) == typeof(f1add)
    @test promote_type(typeof(f1add), typeof(f1mul)) == typeof(f1add)
    @test promote_type(typeof(f1add), typeof(f1add)) == typeof(f1add)

    @test isconcretetype(eltype([f1, f1 * 1]))
    @test isconcretetype(eltype([f1, f1 + 1]))
    @test isconcretetype(eltype([f1 * 2, f1 + 1]))
end

@testitem "Heterogeneous operator vectors" setup = [Fermions] begin
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    mul = f1 * f2
    complex_mul = (1 + 2im) * f1
    add = f1 + f2
    complex_add = 1im + f1 + f2

    # Array literals and vcat should accept every ordering of symbols, products,
    # and sums, while preserving the values and choosing one concrete eltype.
    combinations = (
        (f1, mul, add),
        (f1, add, mul),
        (mul, f1, add),
        (mul, add, f1),
        (add, f1, mul),
        (add, mul, f1),
        (f1, complex_mul, add),
        (complex_add, complex_mul, f1),
        (f1, 2 * f1, (1 + 2im) * f1, f1 + f2, 1im + f1 + f2)
    )
    for values in combinations
        vector = [values...]
        @test isconcretetype(eltype(vector))
        @test all(typeof(value) === eltype(vector) for value in vector)
    end

end

@testitem "Heterogeneous Fermion-Boson vectors" setup = [Fermions, Bosons] begin
    NonCommutativeProducts.enable_autosort!()
    NonCommutativeProducts.@commutative Fermion Boson

    f1 = Fermion(:a)
    f2 = Fermion(:b)
    b = Boson()
    fermion_mul = f1 * f2
    boson_mul = b * b'
    mixed_mul = f1 * b
    fermion_add = f1 + f2
    boson_add = b + b'
    mixed_add = f1 + b

    combinations = (
        (f1, b, fermion_mul, fermion_add),
        (b, f1, boson_mul, boson_add),
        (f1, b, mixed_mul, mixed_add),
        (mixed_add, mixed_mul, b, f1),
        (1im * f1, b, 1im * mixed_mul, 1im + mixed_add),
    )

    for values in combinations
        @test [values...] isa Vector
    end

    # Mixing species must give a concrete factors container (Vector{Any}, not the abstract Vector)
    factorstype(::Type{T}) where {T<:NonCommutativeProducts.NCMul} = fieldtype(T, :factors)
    factorstype(::Type{T}) where {T<:NonCommutativeProducts.NCAdd} = factorstype(T.parameters[2])
    for x in (fermion_mul + boson_mul, (1.0 * fermion_mul) + boson_mul, mixed_mul, mixed_add, fermion_mul + mixed_mul, f1 * f2 + b * b)
        @test isconcretetype(factorstype(typeof(x)))
    end
    fmul = NonCommutativeProducts.NCMul(1, [f1, f2])
    bmul = NonCommutativeProducts.NCMul(1.0, [b, b])
    @test factorstype(promote_type(typeof(fmul), typeof(bmul))) == Vector{Any}
    @test factorstype(typeof(fmul + bmul)) == Vector{Any}
end

@testitem "Promotion and conversion robustness" setup = [Fermions, Bosons] begin
    import NonCommutativeProducts: add!!, NCMul
    using VectorInterface
    NonCommutativeProducts.@commutative Fermion Boson
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    fi = Fermion(1)
    b = Boson()

    for sorted in (false, true)
        sorted ? NonCommutativeProducts.enable_autosort!() : NonCommutativeProducts.disable_autosort!()

        # Widening the additive coefficient widens the term coefficients too
        x = 1.0 + (f1 + f2)
        half = 0.5 + 0.5 * f1 + 0.5 * f2
        @test valtype(x.dict) == Float64
        @test 0.5 * x == half
        @test x / 2 == half
        # / divides each coefficient, as [c] / b does
        @test valtype(((2 * f1 + 2 * f2) / 2).dict) == Float64
        @test (0.1 * f1 + 0.1) / 7.0 == 0.1 / 7.0 * f1 + 0.1 / 7.0 != inv(7.0) * 0.1 * f1 + inv(7.0) * 0.1
        @test (2 * f1 + 2) // 2 == f1 + 1 && valtype(((2 * f1 + 2) // 2).dict) == Rational{Int}
        @test (2 * f1 * f2) // 2 == f1 * f2
        @test VectorInterface.scale!(copy(x), 0.5) == half
        @test VectorInterface.scale!!(copy(x), 0.5) == half
        @test VectorInterface.add(x, f1 + f2, 0.5, 0.5) == 0.5 + f1 + f2
        @test VectorInterface.add!!(copy(x), f1 + f2, 0.5, 0.5) == 0.5 + f1 + f2
        @test (f1 + f2) + 0.5 * f1 == 1.5 * f1 + f2 == 0.5 * f1 + (f1 + f2)
        @test add!!(f1 + f2, 0.5 * f1) == 1.5 * f1 + f2

        # Adding a number gives a new sum that does not share storage with the old one
        y = f1 + f2
        add!!(1 + y, f1)
        add!!(y + 1, f2)
        add!!(0 + y, f1)
        @test y == f1 + f2

        # convert returns exactly the requested type, or throws
        T = typeof(f1 + f2)
        @test convert(T, 3) isa T
        @test convert(T, 3) == 3
        @test_throws InexactError convert(T, 3.5)
        @test_throws InexactError push!([1.5 * f1 + 0], 1im)

        # Atoms, products and sums convert into containers of any compatible type
        for (container, value) in ((f1 * b, f1), (f1 + b, f1), (f1 + b, f1 * f2), (f1 * b, f1 * f2))
            v = push!([container], value)
            @test v[2] == value
        end

        # Heterogeneous vectors get one concrete element type
        fm, bm, mm = f1 * f2, b * b', f1 * b
        fa, ba, ma = f1 + f2, b + b', f1 + b
        combinations = ((fa, ba), (fa, ma), (fm, ma), (mm, fa), (f1, mm), (f1, ma), (fi, fa), (1 * fi, fa), (fm, f1, b, bm, fa, ma), (1im * f1, b, 1im * mm, 1im + ma))
        for values in combinations
            v = [values...]
            @test isconcretetype(eltype(v))
            @test all(x -> typeof(x) === eltype(v), v)
            @test all(splat(==), zip(v, values))
        end
        types = unique(typeof.((f1, fi, b, fm, bm, mm, fa, ba, ma, 1.0 * fm, 1im + fa)))
        for A in types, B in types
            @test promote_type(A, B) == promote_type(B, A)
        end

        # The result type of + does not depend on the order of the terms
        anyf1 = NCMul(1, Any[f1])
        @test typeof(anyf1 + 1.0 * f1) == typeof(1.0 * f1 + anyf1)
        @test typeof(fm + mm) == typeof(mm + fm)
    end
end

@testitem "Coefficient widening, Zero scaling, container types" setup = [Fermions] begin
    import NonCommutativeProducts: add!!, NCMul, NCAdd
    using VectorInterface
    NonCommutativeProducts.disable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)

    # Summing coefficients may widen beyond the promoted coefficient type (Bool + Bool isa Int)
    @test NCMul(true, [f1]) + NCMul(true, [f1]) == 2 * f1
    @test NCMul(true, [f1]) + (NCMul(true, [f1]) + NCMul(true, [f2])) == 2 * f1 + f2
    b = convert(NCAdd{Bool,NonCommutativeProducts.to_add_dict_type(typeof(NCMul(true, [f1])))}, NCMul(true, [f1]))
    @test b + b == 2 * f1
    @test b * b isa NCAdd{Int}

    # Tuple-backed and Vector-backed products promote to a concrete factors container
    tv = NCMul(1, (f1, f2)) + NCMul(1, [f1])
    K = keytype(tv.dict)
    @test isconcretetype(K)
    @test isconcretetype(fieldtype(K, :factors))
    @test eltype(fieldtype(K, :factors)) == K.parameters[2]
    @test tv == f1 * f2 + f1

    # Adding a multi-term sum into an atom
    @test add!!(f1, f2 + 3, 1, 1) == f1 + f2 + 3
    @test add!!(f1, f2 + 3, 2, 5) == 5 * f1 + 2 * f2 + 6
end

@testitem "Mixed factor containers, aliasing, Dict canonical" setup = [Fermions] begin
    import NonCommutativeProducts: add!!, NCMul, NCAdd
    using VectorInterface
    NonCommutativeProducts.disable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)

    # Tuple-backed and Vector-backed products are equal, hash equally, and share Dict keys
    @test NCMul(1, (f1, f2)) == NCMul(1, [f1, f2])
    @test hash(NCMul(2, (f1, f2))) == hash(NCMul(2, [f1, f2]))
    @test NCMul(1, (f1, f2)) != NCMul(1, [f1])

    # Scalar and zero products equal numbers, so they must hash like them
    for (x, n) in ((NCMul(2, typeof(f1)[]), 2), (NCMul(2, ()), 2), (NCMul(0, [f1]), 0),
        (NCMul(0, (f1, f2)), 0), (NCMul(0.0, [f1]), 0), (NCMul(0, typeof(f1)[]), 0))
        @test x == n
        @test hash(x) == hash(n)
        @test n in Set([x])
        @test x in Set([n])
    end
    s = f1 * f2 + 0
    @test add!!(copy(s), NCMul(1, (f1, f2))) == 2 * f1 * f2
    @test s + (NCMul(1, (f1, f2)) + 0) == 2 * f1 * f2
    @test length((s + (NCMul(1, (f1, f2)) + f1)).dict) == 2

    # Aliased destination and source
    x = f1 + 2 * f2 + 1
    @test VectorInterface.scale!!(x, x, 2) == 2 * f1 + 4 * f2 + 2
    x = f1 + 2 * f2 + 1
    @test add!!(x, x, 2, 3) == 5 * f1 + 10 * f2 + 5
end

@testitem "Sum keys own their factors, Tuple targets, scale!!" setup = [Fermions] begin
    import NonCommutativeProducts: add!!, NCMul, NCAdd
    using VectorInterface
    NonCommutativeProducts.disable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    f3 = Fermion(:c)
    f4 = Fermion(:d)

    # Sorting a product must not change the keys of sums built from it
    for build in (x -> x + f3 * f4, x -> x + (f3 * f4 + 1), x -> f3 * f4 + 1 + x, x -> add!!(f3 * f4 + 1, x), x -> NCAdd(x),
        x -> x + 1, x -> convert(typeof(1.0 + f1), x))
        x = f2 * f1
        s = build(x)
        s_ref = build(f2 * f1)
        sort(x)
        @test s == s_ref
        @test all(k -> haskey(s.dict, k), keys(s_ref.dict))
    end
    # nor products that share its factors, nor the product itself
    for share in (x -> 2 * x, x -> -x, x -> NCMul(x + 0), x -> convert(NCMul{Float64,eltype(x.factors),typeof(x.factors)}, x))
        x = f2 * f1
        y = share(x)
        y_ref = share(f2 * f1)
        @test sort(y) == sort(y_ref)
        @test y == y_ref
        @test x == f2 * f1
        s = x + 0
        sort(NCMul(s))
        @test haskey(s.dict, f2 * f1)
    end

    # Sums only hold Dicts, so that terms compare and are looked up by value
    K1 = typeof(NCMul(1, [f1]))
    a = NCAdd(1.0, IdDict{K1,Float64}(NCMul(1, [f1]) => 1.0))
    @test a.dict isa Dict{K1,Float64}
    @test a == 1.0 + f1
    @test add!!(a, f1) == 1.0 + 2 * f1
    @test length(a.dict) == 1

    # In-place scaling either scales everything or throws and leaves the sum unchanged
    x = 2 + f1 + 2 * f2
    @test_throws InexactError VectorInterface.scale!(x, 0.5)
    @test x == 2 + f1 + 2 * f2
    @test_throws InexactError NonCommutativeProducts.add!(x, 1, 1, 0.5)
    @test x == 2 + f1 + 2 * f2
    @test VectorInterface.scale!(x, 2.0) === x
    @test x == 4 + 2 * f1 + 4 * f2
    @test NonCommutativeProducts.add!(x, 1, 1, 0.5) == 3 + f1 + 2 * f2

    # Atoms and Vector-backed products convert to Tuple-backed product and sum types
    F = typeof(f1)
    T1 = NCMul{Int,F,Tuple{F}}
    @test convert(T1, f1) == NCMul(1, (f1,))
    @test convert(T1, f1).factors isa Tuple{F}
    @test convert(NCMul{Int,F,Tuple{F,F}}, f1 * f2) == f1 * f2
    S1 = NCAdd{Int,T1}
    @test convert(S1, f1) == f1 + 0
    @test keytype(convert(S1, f1).dict) == T1

    # scale!! scales in place when the coefficients can hold the result, and never scales twice
    x = f1 + 2 * f2 + 1
    @test VectorInterface.scale!!(x, 3) === x
    @test x == 3 * f1 + 6 * f2 + 3
    x = f1 + 2 * f2 + 1
    y = VectorInterface.scale!!(x, 1.5)
    @test y == 1.5 * f1 + 3 * f2 + 1.5
    @test x == f1 + 2 * f2 + 1
    x = f1 + 2 * f2 + 1
    @test VectorInterface.scale!!(x, VectorInterface.One()) == f1 + 2 * f2 + 1
    @test iszero(VectorInterface.scale!!(f1 + 2 * f2 + 1, VectorInterface.Zero()))
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
    using Symbolics

    @variables a::Real
    σx = Pauli(1)
    σy = Pauli(2)
    σz = Pauli(3)

    ord(op) = sort(op)
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
    using Symbolics

    Sz = SpinOp(0, :z)
    Sp = SpinOp(0, :+)
    Sm = SpinOp(0, :-)

    ord(op) = sort(op)
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
    # products are immutable, so scale! throws
    err = try
        VectorInterface.scale!(f1 * f2, 2)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("immutable", err.msg)

    @test VectorInterface.add(y, x, 2, 3) == 3 * y + 2 * x
    a = copy(y)
    @test VectorInterface.add!(a, x, 2, 3) === a
    @test a == 3 * y + 2 * x
    @test VectorInterface.add!!(copy(y), x, 2, 3) == 3 * y + 2 * x
end

@testitem "KrylovKit compatibility" setup = [Fermions] begin
    using LinearAlgebra, KrylovKit
    using VectorInterface
    using NonCommutativeProducts: anyadd

    NonCommutativeProducts.enable_autosort!()
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    s1 = State(false, f1)
    s2 = State(false, f2)
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

    _vals, _vecs, _ = KrylovKit.eigsolve(ham, s1 + f1' * s1, 2; ishermitian=true)
    perm = sortperm(_vals)
    vals = _vals[perm]
    vecs = _vecs[perm]
    @test vals ≈ numvals
    @test map(abs ∘ inner, diagbasis, vecs) ≈ [1, 1]

    numlinsol = hammat \ [1, 0]
    sol, _ = KrylovKit.linsolve(ham, basis[1])
    @test map(b -> inner(b, sol), basis) ≈ numlinsol

    sol2, _ = KrylovKit.lssolve(ham, basis[1])
    @test inner(sol, sol2) ≈ inner(sol, sol)

    svals, svecsl, svecsr, _ = KrylovKit.svdsolve(ham, anyadd(basis[1]))
    svecsl2, svals2, svecsr2 = svd(hammat)
    @test svals ≈ svals2
    @test abs(dot(map(b -> inner(b, svecsl[1]), basis), svecsl2[:, 1])) ≈ 1
    @test abs(dot(map(b -> inner(b, svecsl[2]), basis), svecsl2[:, 2])) ≈ 1
    @test abs(dot(map(b -> inner(b, svecsr[1]), basis), svecsr2[:, 1])) ≈ 1
    @test abs(dot(map(b -> inner(b, svecsr[2]), basis), svecsr2[:, 2])) ≈ 1

    ## try exponentiate
    expsol, info = exponentiate(ham, 2.0, anyadd(basis[1]))
    numsol = exp(2.0 * hammat) * [1, 0]
    @test map(b -> inner(b, expsol), basis) ≈ numsol
end

@testitem "VectorInterface: products and atoms are represented by sums" setup = [Fermions] begin
    using VectorInterface
    import NonCommutativeProducts: NCAdd, NCMul

    NonCommutativeProducts.disable_autosort!()
    f1 = Fermion(:a)
    p = 2 * f1

    # scale must return the same type as zerovector and add, so that solvers can store them in one container
    for x in (p, f1)
        for s in (VectorInterface.scale(x, 1.5), VectorInterface.scale!!(x, 1.5))
            @test s isa NCAdd
            @test typeof(s) == typeof(VectorInterface.zerovector(s))
            @test typeof(s) == typeof(VectorInterface.add(s, s, 1, 1))
        end
    end
    @test VectorInterface.scale(p, 1.5) == 3.0 * f1
    @test VectorInterface.scale(f1, 1.5) == 1.5 * f1
    @test VectorInterface.scalartype(typeof(f1)) == Int
    @test VectorInterface.scalartype(f1) == Int
    @test (@test_logs VectorInterface.zerovector(f1)) == 0
    @test VectorInterface.zerovector(f1) isa NCAdd{Int}
    @test VectorInterface.zerovector!!(f1, Float64) isa NCAdd{Float64}
    @test VectorInterface.zerovector!!(f1) isa NCAdd
    @test typeof(VectorInterface.zerovector!!(f1)) == typeof(VectorInterface.scale(f1, 1))

    # atoms work in every argument position, without hitting VectorInterface's (warning) fallbacks
    f2 = Fermion(:b)
    forms = (identity, x -> 1 * x, x -> x + 0) # atom, NCMul, NCAdd
    for fy in forms, fx in forms
        y, x = fy(f1), fx(f2)
        r = @test_logs VectorInterface.add(y, x, 2, 3)
        @test r isa NCAdd && r == 3 * f1 + 2 * f2
        r = @test_logs VectorInterface.add!!(deepcopy(y), x, 2, 3)
        @test r isa NCAdd && r == 3 * f1 + 2 * f2
        r = @test_logs VectorInterface.scale!!(deepcopy(y), x, 2)
        @test r isa NCAdd && r == 2 * f2
        if y isa NCAdd
            yc = deepcopy(y)
            r = @test_logs VectorInterface.add!(yc, x, 2, 3)
            @test r === yc && r == 3 * f1 + 2 * f2
            yc = deepcopy(y)
            r = @test_logs VectorInterface.scale!(yc, x, 2)
            @test r === yc && r == 2 * f2
        else # atoms and products are immutable
            @test_throws ArgumentError VectorInterface.add!(y, x, 2, 3)
            @test_throws ArgumentError VectorInterface.scale!(y, x, 2)
        end
    end
    @test_throws ArgumentError VectorInterface.scale!(f1, 2)
end

@testitem "KrylovKit: exponentiate from a product or an atom" begin
    using KrylovKit, LinearAlgebra
    using VectorInterface: inner
    import NonCommutativeProducts: @nc, mul_effect

    # A ket whose adjoint is a bra, with ⟨i|j⟩ = δᵢⱼ
    struct Ket
        n::Int
        bra::Bool
    end
    Ket(n) = Ket(n, false)
    Base.adjoint(k::Ket) = Ket(k.n, !k.bra)
    @nc Ket
    mul_effect(a::Ket, b::Ket) = a.bra && !b.bra ? Int(a.n == b.n) : nothing

    NonCommutativeProducts.enable_autosort!()
    H = Ket(0) * Ket(1)' + Ket(1) * Ket(0)' # σx

    for x0 in (Ket(0), 1.0 * Ket(0), Ket(0) + 0)
        x, info = exponentiate(H, 1.0, x0)
        @test info.converged > 0
        @test norm(x - (cosh(1.0) * Ket(0) + sinh(1.0) * Ket(1))) < 1e-10

        # real-time evolution: exp(-iπ/2 σx)|0⟩ = -i|1⟩
        x, info = exponentiate(H, -im * pi / 2, x0)
        @test info.converged > 0
        @test norm(x - (-im * Ket(1))) < 1e-10

        vals, _ = eigsolve(H, x0, 2, :SR; ishermitian=true)
        @test vals ≈ [-1, 1]
    end
    @test norm(im * Ket(0) + 0) isa Real

    # an operator that is a single product: exp(|1⟩⟨0|)|0⟩ = |0⟩ + |1⟩
    x, info = exponentiate(Ket(1) * Ket(0)', 1.0, Ket(0))
    @test info.converged > 0
    @test norm(x - (Ket(0) + Ket(1))) < 1e-10

    # (σx + 2)x = |0⟩ has the solution 2/3|0⟩ - 1/3|1⟩
    x, info = linsolve(H + 2, 1.0 * Ket(0))
    @test info.converged > 0
    @test norm(x - (2 / 3 * Ket(0) - 1 / 3 * Ket(1))) < 1e-10

    # inner agrees with computing x' * y as a sum
    vecs = (x, H * x, 1im * Ket(0) + 2 * Ket(1), Ket(1))
    @test all(inner(a, b) ≈ NonCommutativeProducts.scalar(a' * b) for a in vecs, b in vecs)

    # an atom whose overlap with itself is a complex number with zero imaginary part has a real norm
    struct CKet
        bra::Bool
    end
    Base.adjoint(k::CKet) = CKet(!k.bra)
    @nc CKet
    mul_effect(a::CKet, b::CKet) = a.bra && !b.bra ? 1 + 0im : nothing
    @test norm(CKet(false)) isa Real
    @test norm(CKet(false)) == 1
    NonCommutativeProducts.disable_autosort!()
end

@testitem "VectorInterface: dot and overlaps that don't reduce to scalars" begin
    using LinearAlgebra, VectorInterface
    import NonCommutativeProducts: @nc, mul_effect

    struct Ket
        n::Int
        bra::Bool
    end
    Ket(n) = Ket(n, false)
    Base.adjoint(k::Ket) = Ket(k.n, !k.bra)
    @nc Ket
    mul_effect(a::Ket, b::Ket) = a.bra && !b.bra ? Int(a.n == b.n) : nothing

    NonCommutativeProducts.enable_autosort!()
    k0, k1, k2 = Ket(0), Ket(1), Ket(2)
    x = k0 + 2im * k1
    y = 3 * k0 - k1

    # dot is inner, conjugate-linear in the first argument
    @test dot(x, y) == inner(x, y) == 3 + 2im
    @test dot(im * x, y) == -im * dot(x, y)
    @test norm(x, 2) == norm(x) ≈ sqrt(5)
    @test_throws ArgumentError norm(x, 1)


    # atoms of different types are converted to products, so they can be combined
    struct KetA
        bra::Bool
    end
    struct KetB
        bra::Bool
    end
    Base.adjoint(k::KetA) = KetA(!k.bra)
    Base.adjoint(k::KetB) = KetB(!k.bra)
    @nc KetA KetB
    mul_effect(a::Union{KetA,KetB}, b::Union{KetA,KetB}) = a.bra && !b.bra ? Int(typeof(a) == typeof(b)) : nothing
    a, b = KetA(false), KetB(false)
    @test inner(a + b, b) == 1
    @test add(a, b, 2, 1) == a + 2 * b

    # inner agrees with x' * y, and throws when that doesn't reduce to a scalar
    vecs = (x, k0, x + 1, zero(x), k0 * k1', k2 * k1' + 1im * k0 * k0')
    for a in vecs, b in vecs
        ref = a' * b
        if NonCommutativeProducts._reduces_to_scalar(ref)
            @test inner(a, b) ≈ NonCommutativeProducts.scalar(ref)
        else
            @test_throws ArgumentError inner(a, b)
        end
    end
    # products that don't reduce to scalars, but cancel: ⟨0|0⟩|2⟩ - ⟨1|1⟩|2⟩ = 0
    @test inner(k0 + k1, k0 * k2 - k1 * k2) == 0
    # inner sorts x' * y whether or not autosort is enabled
    NonCommutativeProducts.with_autosort(false) do
        @test norm(k0 + 2k1) ≈ sqrt(5)
        @test inner(k0 + k1, k1) == dot(k0 + k1, k1) == 1
        @test_throws ArgumentError inner(k0, k0 * k1')
    end
    NonCommutativeProducts.disable_autosort!()
end

@testitem "VectorInterface: the type of inner doesn't depend on autosort" setup = [Fermions] begin
    using VectorInterface: inner
    f = Fermion(1)
    x, y = zero(1.0im * f + 0), zero(f + 0)
    results = map((true, false)) do sorted
        try
            sorted ? NonCommutativeProducts.enable_autosort!() : NonCommutativeProducts.disable_autosort!()
            (inner(x, y), inner(y, x), inner(2 + 0 * f, 1.0 + 0 * f))
        finally
            NonCommutativeProducts.disable_autosort!()
        end
    end
    @test results[1] == results[2] == (0, 0, 2)
    @test map(typeof, results[1]) == map(typeof, results[2]) == (ComplexF64, ComplexF64, Float64)
end

@testitem "@nc keeps an existing scalartype method" begin
    using NonCommutativeProducts, VectorInterface
    import NonCommutativeProducts: @nc, NCAdd

    struct HasScalartype end
    VectorInterface.scalartype(::Type{HasScalartype}) = Float64
    @nc HasScalartype
    @test VectorInterface.scalartype(HasScalartype) == Float64

    struct NoScalartype end
    @nc NoScalartype
    @test VectorInterface.scalartype(NoScalartype) == Int
    @test VectorInterface.zerovector(HasScalartype()) isa NCAdd{Float64}
    @test VectorInterface.scalartype(VectorInterface.zerovector(HasScalartype())) == Float64
    # scale and add return the same type as zerovector, as KrylovKit requires
    x = HasScalartype()
    T = typeof(VectorInterface.zerovector(x))
    @test VectorInterface.scale(x, 1) isa T
    @test VectorInterface.scale!!(x, 1) isa T
    @test VectorInterface.add(x, x, 1, 1) isa T
    @test VectorInterface.add!!(x, x, 1, 1) isa T
end

@testitem "add!! and mul!! weighted accumulation" setup = [Fermions] begin
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


@testitem "Public API" begin
    using NonCommutativeProducts
    public_names = [:NCMul, :NCAdd, Symbol("@nc"), Symbol("@commutative"), :mul_effect, :enable_autosort!,
        :disable_autosort!, :autosort, :with_autosort, :ncmap, :ncmapreduce, :scalar, :add!!, :terms, :factors,
        :coefficient]
    for name in public_names
        @test Base.ispublic(NonCommutativeProducts, name)
        @test Base.Docs.hasdoc(NonCommutativeProducts, name)
    end
    for name in [:Swap, :AddTerms, :bubble_sort, :MulAdd, :_autosort, :grade, :mul!!, :anyadd, :isscalar]
        @test !Base.ispublic(NonCommutativeProducts, name)
    end
    @test !Base.isexported(NonCommutativeProducts, :NCMul)
end

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(NonCommutativeProducts)
end

@testitem "ExplicitImports" begin
    using ExplicitImports
    using KrylovKit: KrylovKit  # loads the extension so that it is checked too
    test_explicit_imports(NonCommutativeProducts;
        # Base internals without a public equivalent, KrylovKit operator hooks, and package internals used by the extension
        all_qualified_accesses_are_public=(; ignore=(:_growat!, :inferencebarrier, :literal_pow, :power_by_squaring, :promote_op,
            :apply, :apply_normal, :apply_adjoint, :MulAdd)),
        # qualified names in macro-generated code
        no_self_qualified_accesses=(; ignore=(:NCAdd, :NCMul, :_has_scalartype, :_is_atom, :add!!, :anyadd,
            :mul_effect, :ncmapreduce, :terms, :coefficient, :factors)),
    )
end
