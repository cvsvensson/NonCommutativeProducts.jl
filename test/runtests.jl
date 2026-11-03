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
    ord(op) = sort(op)
    ord_equals(a, b) = iszero(ord(a - b))
    ord_equals(a, b, c, xs...) = (ord_equals(a, b) && ord_equals(b, c, xs...))

    # explicit sorting of products left unsorted
    NonCommutativeProducts.with_autosort(false) do
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
end


@testitem "Majoranas" setup = [Majoranas] begin
    using Random: seed!
    seed!(1)
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

@testitem "Autosort is the default" setup = [Fermions] begin
    @test NonCommutativeProducts.autosort()
    f = Fermion(1)
    @test f * f == 0
    @test f * f' == 1 - f' * f
end

@testitem "A product without mul_effect asks for one" begin
    import NonCommutativeProducts: @nc, NCMul, with_autosort
    struct NoEffect end
    @nc NoEffect
    x = NoEffect()
    err = @test_throws MethodError x * x
    @test occursin("with_autosort(false)", sprint(showerror, err.value))
    @test with_autosort(() -> x * x, false) == NCMul(1, [x, x])
end

@testitem "Autosort local override with ScopedValue" setup = [Majoranas] begin
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
        NonCommutativeProducts.with_autosort(sorted) do
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
end

@testitem "Coefficient widening, Zero scaling, container types" setup = [Fermions] begin
    import NonCommutativeProducts: add!!, NCMul, NCAdd
    using VectorInterface
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
    f1 = Fermion(:a)
    f2 = Fermion(:b)
    f3 = Fermion(:c)
    f4 = Fermion(:d)

    NonCommutativeProducts.with_autosort(false) do
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
    @test convert(NCMul{Int,F,Tuple{F,F}}, NCMul(1, [f1, f2])) == f1 * f2
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
        VectorInterface.scale!(2 * f1, 2)
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
    # with Number effects of the same type, a zero sum on either side gives the type of a nonzero sum. Bool + Bool is an Int.
    for v in (NonCommutativeProducts.NCMul(true, [k0]), x)
        for z in (inner(zero(v), v), inner(v, zero(v)))
            @test iszero(z)
            @test typeof(z) == typeof(inner(v, v))
        end
    end
end

@testitem "VectorInterface: inner with contraction keys" begin
    using VectorInterface
    import NonCommutativeProducts: @nc, mul_effect, scalar, _reduces_to_scalar

    struct KKet
        n::Int
        bra::Bool
    end
    KKet(n) = KKet(n, false)
    Base.adjoint(k::KKet) = KKet(k.n, !k.bra)
    @nc KKet
    mul_effect(a::KKet, b::KKet) = a.bra && !b.bra ? Int(a.n == b.n) : nothing
    # |0⟩ and the bras have no key, so terms mix keyed and wildcard factors
    NonCommutativeProducts.contraction_key(k::KKet) = k.bra || k.n == 0 ? nothing : k.n

    k0, k1, k2, k3 = KKet.(0:3)
    x = k0 + 2im * k1
    y = 2k1 - k2 + 3k3
    vecs = (x, y, x + 1, k2, k1 * k2', zero(x))
    for a in vecs, b in vecs
        ref = a' * b
        if _reduces_to_scalar(ref)
            @test inner(a, b) == scalar(ref)
            @test typeof(inner(a, b)) == typeof(scalar(ref))
        else
            @test_throws ArgumentError inner(a, b)
        end
    end
end

@testitem "VectorInterface: inner with contraction keys on products" begin
    using LinearAlgebra, VectorInterface
    import NonCommutativeProducts: @nc, mul_effect, NCMul

    # a ket |s⟩ on a site, or its bra. State 0 has no key, so a keyed run can stop in the middle.
    struct SKet
        site::Int
        s::Int
        bra::Bool
    end
    SKet(site, s) = SKet(site, s, false)
    Base.adjoint(k::SKet) = SKet(k.site, k.s, !k.bra)
    @nc SKet
    mul_effect(a::SKet, b::SKet) = a.bra && !b.bra ? Int(a.site == b.site && a.s == b.s) : nothing
    NonCommutativeProducts.contraction_key(k::SKet) = k.bra || k.s == 0 ? nothing : (k.site, k.s)

    a1, a2, b1, b2, u2 = SKet(1, 1), SKet(2, 1), SKet(1, 2), SKet(2, 2), SKet(2, 0)
    @test inner(a1 * a2 + 2im * b1 * a2, a1 * a2 + b1 * a2) == 1 - 2im
    @test norm(a1 * a2 + 2b1 * a2 - a1 * b2) ≈ sqrt(6)
    @test inner(a1 * u2 + b1 * a2, a1 * u2) == 1
    # products backed by a Tuple
    @test inner(a1 * a2, NCMul(1, (a1, a2))) == 1
    # products of different lengths that don't reduce to a scalar
    @test_throws ArgumentError inner(a1 * a2, a1)

    # mul_effect between atoms with equal keys must return a Number
    struct BKet
        n::Int
        bra::Bool
    end
    Base.adjoint(k::BKet) = BKet(k.n, !k.bra)
    @nc BKet
    mul_effect(a::BKet, b::BKet) = a.bra && !b.bra ? (a.n == b.n ? 2 * BKet(0, false) : 0) : nothing
    NonCommutativeProducts.contraction_key(k::BKet) = k.bra ? nothing : k.n
    @test_throws "must be a Number" inner(BKet(1, false), BKet(1, false))
end

@testitem "VectorInterface: inner with keys over several groups" begin
    using VectorInterface
    import NonCommutativeProducts: @nc, mul_effect, with_autosort, scalar

    # kets |n⟩ of several groups, keyed by (group, n). A bra and a ket of different groups don't act on each other.
    struct GKet
        g::Int
        n::Int
        bra::Bool
    end
    GKet(g, n) = GKet(g, n, false)
    Base.adjoint(k::GKet) = GKet(k.g, k.n, !k.bra)
    @nc GKet
    mul_effect(a::GKet, b::GKet) = a.bra && !b.bra && a.g == b.g ? Int(a.n == b.n) : nothing
    NonCommutativeProducts.contraction_key(k::GKet) = k.bra ? nothing : (k.g, k.n)

    ψ = sum((m + 2im * n) * GKet(1, m) * GKet(2, n) for m in 1:5, n in 1:5)
    φ = GKet(1, 2) * GKet(2, 3) - 2 * GKet(1, 4) * GKet(2, 3) + GKet(1, 2) * GKet(2, 5)
    ket = GKet(1, 3) * GKet(2, 4)
    for (a, b) in ((ψ, ψ), (ψ, φ), (φ, ψ), (ket, ψ), (ψ, ket))
        @test inner(a, b) == scalar(with_autosort(() -> a' * b))
    end
    # keys define orthogonality: states of group 3 alone are orthogonal to those of groups 1 and 2, although x' * y
    # doesn't sort to a scalar
    @test inner(φ + 2 * GKet(3, 1), ψ + 3 * GKet(3, 1)) == inner(φ, ψ) + 6
end

@testitem "VectorInterface: inner with keys over signed swaps" begin
    using VectorInterface
    import NonCommutativeProducts: @nc, mul_effect, NCMul, Swap

    # kets of two sites whose factors anticommute: the signs cancel within each pair of terms
    struct SwapKet
        site::Int
        n::Int
        bra::Bool
    end
    SwapKet(site, n) = SwapKet(site, n, false)
    Base.adjoint(k::SwapKet) = SwapKet(k.site, k.n, !k.bra)
    @nc SwapKet
    function mul_effect(a::SwapKet, b::SwapKet)
        a.site == b.site || return a.site > b.site ? Swap(-1) : nothing
        return a.bra && !b.bra ? Int(a.n == b.n) : nothing
    end
    NonCommutativeProducts.contraction_key(k::SwapKet) = k.bra ? nothing : (k.site, k.n)

    state(a, b) = NCMul(1, [SwapKet(1, a), SwapKet(2, b)])
    x = 2 * state(0, 0) + 3 * state(0, 1) + 5 * state(1, 1)
    y = 7 * state(0, 1) + 11 * state(1, 1) + 13 * state(1, 0)
    @test inner(x, y) == 3 * 7 + 5 * 11
end

@testitem "VectorInterface: inner compares keys with isequal" begin
    using VectorInterface
    import NonCommutativeProducts: @nc, mul_effect, NCMul

    # keys whose hashes all collide
    struct HashKey
        n::Int
    end
    Base.hash(::HashKey, h::UInt) = h
    struct CKet
        n::Int
        bra::Bool
    end
    CKet(n) = CKet(n, false)
    Base.adjoint(k::CKet) = CKet(k.n, !k.bra)
    @nc CKet
    # nonzero for different keys, which define orthogonality, so inner never uses it
    mul_effect(a::CKet, b::CKet) = a.bra && !b.bra ? (a.n == b.n ? 2 : 5) : nothing
    NonCommutativeProducts.contraction_key(k::CKet) = k.bra ? nothing : HashKey(k.n)

    x = CKet(1) * CKet(2) + 3 * CKet(2) * CKet(1) + 5 * CKet(1) * CKet(1)
    y = CKet(1) * CKet(2) + CKet(2) * CKet(2) + 7 * CKet(1) * CKet(1)
    @test inner(x, y) == 4 + 35 * 4
    @test inner(NCMul(1, Any[CKet(1), CKet(2)]), x) == 4
end

@testitem "VectorInterface: inner with keys and an unkeyed ket that overlaps them" begin
    using LinearAlgebra, VectorInterface
    import NonCommutativeProducts: @nc, mul_effect

    # Keyed kets orthogonal to each other with ⟨n|n⟩ = 2, and an unkeyed u with ⟨u|u⟩ = 3, ⟨n|u⟩ = ⟨u|n⟩ = 10n.
    # The keys have different types.
    struct UKet
        n::Int
        bra::Bool
    end
    UKet(n) = UKet(n, false)
    Base.adjoint(k::UKet) = UKet(k.n, !k.bra)
    @nc UKet
    function mul_effect(x::UKet, y::UKet)
        x.bra && !y.bra || return nothing
        x.n == y.n && return x.n == 0 ? 3 : 2
        return 10 * max(x.n, y.n) * (x.n * y.n == 0)
    end
    NonCommutativeProducts.contraction_key(k::UKet) = k.bra || k.n == 0 ? nothing : isodd(k.n) ? k.n : string(k.n)

    u, a, b, c, d = UKet.(0:4)
    @test inner(a + b + c, u + a) == 62 # 10 + 20 + 30 + 2; y indexed
    @test inner(u + a, a + b + c) == 62 # x indexed
    @test inner(u + a + b, a + c) == 42 # 10 + 30 + 2
    @test inner(a + c, u + a + b) == 42
    @test inner(a * b + c * c, a * u + a * b) == 44 # 2 * 20 + 2 * 2; runs end at different depths
    @test inner(a * b + a * u, a * b + c * c + c * a) == 44
    @test inner(a * u + c * c, a * b + a * c) == 100 # 2 * 20 + 2 * 30; a run stops where the other side splits
    @test inner(a * b + a * c, a * u + c * c + c * a) == 100
    @test dot(2im * a + b, a) == -4im
    @test dot(a, 2im * a + b) == 4im
    @test inner(2im + 0 * a, 3 + 0 * a) === -6im
    @test inner(im * a + b, c + d + UKet(5)) === 0im
    @test inner(im * a + b + c, d + UKet(5)) === 0im
end

@testitem "VectorInterface: inner with operators that act on keyed kets" begin
    using LinearAlgebra, VectorInterface
    import NonCommutativeProducts: @nc, mul_effect, contraction_key, with_autosort, AddTerms, Swap

    # a†|n⟩ = |n+1⟩ and a|n⟩ = n|n-1⟩, which are each other's adjoints when ⟨n|n⟩ = n!
    struct RKet
        n::Int
        bra::Bool
    end
    RKet(n) = RKet(n, false)
    Base.adjoint(k::RKet) = RKet(k.n, !k.bra)
    struct Ladder
        up::Bool
    end
    Base.adjoint(a::Ladder) = Ladder(!a.up)
    @nc RKet Ladder
    mul_effect(a::RKet, b::RKet) = a.bra && !b.bra ? (a.n == b.n ? factorial(a.n) : 0) : nothing
    mul_effect(a::Ladder, k::RKet) = k.bra ? nothing : a.up ? RKet(k.n + 1) : k.n * RKet(k.n - 1)
    mul_effect(k::RKet, a::Ladder) = !k.bra ? nothing : a.up ? k.n * RKet(k.n - 1, true) : RKet(k.n + 1, true)
    mul_effect(a::Ladder, b::Ladder) = !a.up && b.up ? AddTerms((Swap(1), 1)) : nothing
    contraction_key(k::RKet) = k.bra ? nothing : k.n

    up, dn = Ladder(true), Ladder(false)
    k1, k2, k3 = RKet.(1:3)
    with_autosort(false) do
        @test inner(up * k1, k2) == inner(k2, up * k1) == 2
        @test inner(dn * k2, k1) == 2
        @test inner(up * k1, up * k1) == 2
        @test inner(im * up * k1 + k3, k2 + 2k3) == 12 - 2im # -2im + 2 * 3!
        @test dot(k2 + 2k3, im * up * k1 + k3) == 12 + 2im
    end
end

@testitem "mul_terms!! hook" begin
    import NonCommutativeProducts: @nc, mul_effect, with_autosort, add!!, terms, nterms, NCMul, NCAdd

    struct HKet
        n::Int
        bra::Bool
    end
    HKet(n) = HKet(n, false)
    Base.adjoint(k::HKet) = HKet(k.n, !k.bra)
    @nc HKet
    mul_effect(a::HKet, b::HKet) = a.bra && !b.bra ? Int(a.n == b.n) : nothing

    H = sum(HKet(n) * HKet(n + 1)' + 2 * HKet(n + 1) * HKet(n)' for n in 1:5)
    ψ = sum(n * HKet(n) for n in 1:6) + HKet(1) * HKet(2)
    ref = H * ψ
    unsorted = with_autosort(() -> H * ψ, false)

    const calls = Ref(0)
    # the generic path through products of terms; autosort is off in the hook, so it is enabled for the reference
    function NonCommutativeProducts.mul_terms!!(c, x::NCMul{<:Any,HKet}, b::NCAdd)
        calls[] += 1
        return add!!(c, with_autosort(() -> sum(x * y for y in terms(b))))
    end
    @test H * ψ == ref
    @test calls[] == nterms(H)
    calls[] = 0
    @test with_autosort(() -> H * ψ, false) == unsorted
    @test calls[] == 0
end

@testitem "VectorInterface: the type of inner doesn't depend on autosort" setup = [Fermions] begin
    using VectorInterface: inner
    f = Fermion(1)
    x, y = zero(1.0im * f + 0), zero(f + 0)
    results = map((true, false)) do sorted
        NonCommutativeProducts.with_autosort(sorted) do
            (inner(x, y), inner(y, x), inner(2 + 0 * f, 1.0 + 0 * f))
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


@testitem "Truncation with an error bound" begin
    using NonCommutativeProducts, LinearAlgebra
    import NonCommutativeProducts as NC
    # Jordan-Wigner fermions on two modes, a spin-1/2 component and an atom without a bound
    struct JWFermion
        mode::Int
        dagger::Bool
    end
    struct HalfSpin end
    struct Unbounded end
    NC.@nc JWFermion HalfSpin Unbounded
    NC.mul_effect(::Union{JWFermion,HalfSpin,Unbounded}, ::Union{JWFermion,HalfSpin,Unbounded}) = nothing
    Base.adjoint(f::JWFermion) = JWFermion(f.mode, !f.dagger)
    NC.opnorm_bound(::JWFermion) = 1
    NC.opnorm_bound(::HalfSpin) = 1 // 2

    # a representation on two fermion modes and a spin, by 8×8 matrices
    lower = [0 1; 0 0]
    Z = [1 0; 0 -1]
    rep(f::JWFermion) = (c = f.mode == 1 ? kron(lower, I(2), I(2)) : kron(Z, lower, I(2)); f.dagger ? c' : c)
    rep(::HalfSpin) = kron(I(2), I(2), [1 0; 0 -1] / 2)
    rep(::Unbounded) = 3.0 * I(8)
    rep(x) = sum((c * prod(rep, fs; init=Matrix(1.0I, 8, 8)) for (fs, c) in NC.termpairs(x)); init=zeros(8, 8))

    c1, c2 = JWFermion(1, false), JWFermion(2, false)
    S, u = HalfSpin(), Unbounded()
    bounded = 0.5 * c1' * c2 + 0.25 * c2' * c1 + 0.125 * c1 + 0.125 * S + 0.03125 * c1' * c1 * c2 + 1
    x = bounded + u * c1 / 64
    # weights: 0.03125 < 0.0625 (the spin) < 0.125 < 0.25 < 0.5 < 1 (the scalar), and no weight for u * c1
    @test NC.opnorm_bound(bounded) == 1.96875
    @test NC.opnorm_bound(x) === nothing
    @test opnorm(rep(bounded)) <= NC.opnorm_bound(bounded)

    y, err = NC.truncate_terms(x; atol=0.1)
    @test y == x - 0.125 * S - 0.03125 * c1' * c1 * c2
    @test err == 0.09375
    @test opnorm(rep(x) - rep(y)) <= err
    # rtol is relative to the bound of the terms that can be dropped
    @test NC.truncate_terms(x; rtol=0.05) == (; expr=y, err=err)
    # everything but the term without a bound, including the scalar
    y, err = NC.truncate_terms(x; atol=2)
    @test y == u * c1 / 64 && err == 1.96875
    @test opnorm(rep(x) - rep(y)) <= err
    @test NC.truncate_terms(x) == (; expr=x, err=0)
    @test NC.truncate_terms(x; atol=Inf).expr == u * c1 / 64
    # terms with the same weight are dropped together or not at all
    @test NC.truncate_terms(0.25 * c1 + 0.25 * c2 + 1; atol=0.3) == (; expr=0.25 * c1 + 0.25 * c2 + 1, err=0)
    @test NC.truncate_terms(0.25 * c1 + 0.25 * c2 + 1; atol=0.5) == (; expr=1, err=0.5)
    # also a product and its adjoint, whose bounds multiply in opposite orders
    struct Bounded
        bound::Float64
        dagger::Bool
    end
    NC.@nc Bounded
    NC.mul_effect(::Bounded, ::Bounded) = nothing
    Base.adjoint(a::Bounded) = Bounded(a.bound, !a.dagger)
    NC.opnorm_bound(a::Bounded) = a.bound
    a, b = Bounded(sqrt(2.0), false), Bounded(sqrt(3.0), false)
    H = 0.1 * a * b + 0.1 * b' * a'
    @test NC.truncate_terms(H; atol=0.3) == (; expr=H, err=0)
    @test_throws ArgumentError NC.truncate_terms(Bounded(-1.0, false); atol=1)
    # rtol = 1 drops everything
    @test NC.truncate_terms(0.1 * c1 + 0.2 * c2 + 0.3 * c1' + 0.6 * c2'; rtol=1).expr == 0
    # NaN and infinite coefficients are never dropped, and make the bound infinite
    @test NC.truncate_terms(NaN * c1 + Inf * c2 + 0.5 * c1' * c2; atol=1).err == 0.5
    @test NC.opnorm_bound(NaN * c1 + c2) == Inf

    # products, atoms and numbers
    @test NC.truncate_terms(0.25 * c1 * S; atol=0.2) == (; expr=0, err=0.125)
    @test NC.truncate_terms(0.25 * c1 * S; atol=0.1) == (; expr=0.25 * c1 * S, err=0)
    @test NC.truncate_terms(c1; atol=1) == (; expr=0, err=1)
    @test NC.truncate_terms(u; atol=Inf) == (; expr=u, err=0)
    @test NC.truncate_terms(-0.5; atol=1) == (; expr=0, err=0.5)
    @test NC.truncate_terms(0.5 + 1im; atol=1) == (; expr=0.5 + 1im, err=0)
    @test NC.truncate_terms(0 * c1 + 0) == (; expr=0, err=0)
    # factors in a Tuple or an abstractly typed Vector
    @test NC.opnorm_bound(NC.NCMul(1, (c1, c2, c1'))) == 1
    @test NC.opnorm_bound(NC.NCMul(1, Any[c1, c2, c1'])) == 1
    @test NC.opnorm_bound(NC.NCMul(1, Union{JWFermion,HalfSpin}[c1, c2, S])) == 0.5

    @test_throws ArgumentError NC.truncate_terms(x; atol=-1)
    @test_throws ArgumentError NC.truncate_terms(x; rtol=-0.1)
    @test_throws ArgumentError NC.truncate_terms(x; atol=NaN)
end

@testitem "Public API" begin
    using NonCommutativeProducts
    public_names = [:NCMul, :NCAdd, Symbol("@nc"), Symbol("@commutative"), :mul_effect, :Swap, :AddTerms, :autosort, :with_autosort,
        :ncmap, :ncmapreduce, :scalar, :add!!, :terms, :factors, :coefficient]
    for name in public_names
        @test Base.ispublic(NonCommutativeProducts, name)
        @test Base.Docs.hasdoc(NonCommutativeProducts, name)
    end
    for name in [:enable_autosort!, :disable_autosort!, :bubble_sort, :MulAdd, :_autosort, :grade, :mul!!, :anyadd, :isscalar,
        :atom_group, :foreach_segment, :split_by_group, :segment_ranges]
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
    # load the extensions so that they are checked too
    using KrylovKit: KrylovKit
    using Symbolics: Symbolics
    test_explicit_imports(NonCommutativeProducts;
        # Base internals without a public equivalent, KrylovKit operator hooks, and package internals used by the extensions
        all_qualified_accesses_are_public=(; ignore=(:_growat!, :inferencebarrier, :literal_pow, :power_by_squaring, :promote_op,
            :apply, :apply_normal, :apply_adjoint, :MulAdd, :map_coefficients, :_unwrap_symbolic, :_coefficient_bound,
            :SizeUnknown, :EltypeUnknown, :Experimental, :register_error_hint)),
        # qualified names in macro-generated code
        no_self_qualified_accesses=(; ignore=(:NCAdd, :NCMul, :_has_scalartype, :_is_atom, :_ordered_before, :add!!, :anyadd,
            :mul_effect, :ncmapreduce, :terms, :termpairs, :map_coefficients, :coefficient, :factors, :nterms)),
    )
end
