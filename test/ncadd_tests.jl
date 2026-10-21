@testitem "NCAdd: failed in-place updates leave the target unchanged" setup = [Fermions] begin
    using VectorInterface
    import NonCommutativeProducts as NC
    f1, f2 = Fermion(1), Fermion(2)

    # Int coefficients scaled by 0.5: the coefficient 2 * 0.5 fits in Int, the term f1 * 0.5 doesn't
    x = 2 + f1 + 2 * f2
    @test_throws InexactError VectorInterface.scale!(x, 0.5)
    @test x == 2 + f1 + 2 * f2
    @test_throws InexactError NC.scale!(x, 0.5)
    @test x == 2 + f1 + 2 * f2
    @test_throws InexactError NC.add!(x, 1, 1, 0.5)
    @test x == 2 + f1 + 2 * f2
    # coefficient that doesn't fit
    @test_throws InexactError NC.add!(x, 0.5)
    @test x == 2 + f1 + 2 * f2

    y = 2 + f1
    @test_throws InexactError VectorInterface.scale!(y, 2 + f1 + 2 * f2, 0.5)
    @test y == 2 + f1
    @test_throws InexactError VectorInterface.scale!(y, NC.NCMul(1, [f2]), 0.5)
    @test y == 2 + f1
    @test_throws InexactError VectorInterface.add!(y, 2 + f1 + 2 * f2, 0.5, 1)
    @test y == 2 + f1
    @test_throws InexactError VectorInterface.add!(y, 2 + 2 * f2, 1, 0.5)
    @test y == 2 + f1

    # and when they succeed, they mutate
    y = 2 + f1
    @test VectorInterface.scale!(y, 2 + f1 + 2 * f2, 2) === y
    @test y == 4 + 2 * f1 + 4 * f2
    @test VectorInterface.scale!(y, 3 * f2, 2) === y
    @test y == 6 * f2
    @test VectorInterface.add!(y, 2 + f1, 1, 2) === y
    @test y == 2 + f1 + 12 * f2
end

@testitem "NCAdd: set_coeff!!" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    f1 = Fermion(1)
    a = f1 + 1
    b = NC.set_coeff!!(a, 3)
    @test b === a
    @test a == f1 + 3
    c = NC.set_coeff!!(a, 0.5)
    @test c !== a
    @test c == f1 + 0.5
    @test NC.additive_coeff(c) isa Float64
    @test a == f1 + 3
end

@testitem "NCAdd: constructor without filtering" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    f1 = Fermion(1)
    k = NC.NCMul(1, [f1])
    x = NC.NCAdd(0, Dict(k => 0); filter_zeros=false, filter_scalars=false)
    @test length(x.dict) == 1
    @test only(values(x.dict)) == 0
    @test NC.NCAdd(0, Dict(k => 0)) == 0
    # the terms are always a Dict with the coefficient type as value type
    @test typeof(x.dict) == Dict{typeof(k),Int}
end

@testitem "NCAdd: zero terms are filtered on all paths" setup = [Fermions, Bosons] begin
    using VectorInterface
    import NonCommutativeProducts as NC
    import NonCommutativeProducts: add!!, add!, scale!!
    f1, f2 = Fermion(1), Fermion(2)
    b1 = Boson()

    function iscanonical(y, ref)
        y == ref && hash(y) == hash(ref) && !any(iszero, values(y.dict))
    end

    @test iscanonical(add!!(f1 + 1, -f1), 1)
    @test iscanonical(add!!(f1 + 1, f1, -1, 1), 1)
    @test iscanonical(add!!(f1 + 1, -f1 + 0), 1)
    @test iscanonical(add!!(f1 + 1, -1 * (f1 + 1)), 0)
    @test iscanonical(add!!(f1 + 1, f1, 1, 0), f1)
    @test iscanonical(add!!(f1 + f2, -1.0 * f1), f2) # widening coefficients
    # a coefficient that rounds to zero in the value type (promote_type(Float16, Rational{Int}) == Float16)
    @test iscanonical(add!!(Float16(1) * f1 + 0, (1 // 10^8) * f2 + 0), Float16(1) * f1 + 0)
    @test iscanonical(add!!(Float16(1) * f1 + 0, (1 // 10^8) * f2), Float16(1) * f1 + 0)
    @test iscanonical(add!!(1 + f2, NC.NCMul(0, [f1])), 1 + f2) # no new zero term
    @test iscanonical(add!!(f1 + 1, NC.NCMul(2, typeof(f1)[])), f1 + 3) # scalar products go to the coefficient
    @test iscanonical(add!(f1 + 1, 1, 1, 0), 1)
    @test iscanonical(NC.scale!(f1 + 1, 0), 0)
    @test iscanonical(scale!!(f1 + 1, 0), 0)
    @test iscanonical(scale!!(f1 + 1, 0.0), 0)
    @test iscanonical(VectorInterface.scale!!(f1 + 1, VectorInterface.Zero()), 0)
    @test iscanonical(VectorInterface.scale!(f1 + 1, f1 + f2, 0), 0)
    @test iscanonical(VectorInterface.add!(f1 + 1, f1 + 0, -1, 1), 1)
    @test iscanonical(VectorInterface.add!!(f1 + 1, f1 + 0, -1, 1), 1)
    # scaling by a nonzero number that underflows the terms to zero
    @test iscanonical(add!(1e-200 * f1 + 1.0, 0, 1, 1e-200), 1e-200)
    @test iscanonical(NC.scale!(1e-200 * f1 + 1.0, 1e-200), 1e-200)
    @test iscanonical(VectorInterface.scale!!(1e-200 * f1 + 1.0, 1e-200), 1e-200)
    @test iscanonical(f1 + 1 - f1, 1)
    @test iscanonical((f1 + 1) + (-f1 + 0), 1)
    @test iscanonical((f1 + 2 * f2) + (-1.0 * f1 + 0), 2.0 * f2 + 0) # widening coefficients
    @test iscanonical((true * f1 + 0) + (true * f1 + 0), 2 * f1 + 0) # Bool + Bool isa Int
    @test iscanonical((1 * f1 + b1) + (-b1 + 0), f1 + 0) # keys of different types
    # coefficients that round to zero in the promoted type (promote_type(Float16, Rational{Int}) == Float16)
    @test iscanonical(((1 // 10^8) * f1 + 0) + (Float16(1) * f2 + 0), Float16(1) * f2 + 0)
    @test iscanonical((Float16(1) * f2 + 0) + ((1 // 10^8) * f1 + 0), Float16(1) * f2 + 0)
    @test iscanonical(0 * (f1 + 1), 0)
    # mixed species: the keys of the terms have different types
    @test iscanonical(add!!(1 * f1 + b1, -b1), f1 + 0)
    @test iscanonical(add!!(f1 + 1, b1, 1, 0), b1 + 0)
    @test iscanonical(add!!(1 * f1 + b1 + 1, -f1 - b1), 1)
    @test iscanonical(1 * f1 + b1 - b1, f1 + 0)
end

@testitem "NCAdd: UniformScaling + NCAdd copies" setup = [Fermions] begin
    using LinearAlgebra
    import NonCommutativeProducts: add!!
    f1, f2 = Fermion(1), Fermion(2)
    for op in (x -> I + x, x -> x + I, x -> 2I + x, x -> x - I)
        x = f1 + 1
        y = op(x)
        @test y !== x
        @test y.dict !== x.dict
        add!!(y, f2)
        @test x == f1 + 1
    end
end

@testitem "NCAdd: show prints all terms in sorted order" setup = [Fermions] begin
    f = Fermion.(1:5)
    x = 1 + f[3] - 2 * f[1] + f[2] + f[4] + f[5]
    @test string(x) == "1I - 2*c[1] + c[2] + c[3] + c[4] + c[5]"
    @test string(x) == string(1 + f[5] + f[4] + f[3] + f[2] - 2 * f[1])
    @test !occursin("...", string(x))
    @test sprint(show, MIME"text/plain"(), x) == string(x)

    limited = IOContext(IOBuffer(), :limit => true)
    show(limited, MIME"text/plain"(), x)
    @test String(take!(limited.io)) == "Sum with 6 terms: \n1I-2*c[1]+c[2]+c[3] + ..."
    short = 1 + f[1]
    show(limited, MIME"text/plain"(), short)
    @test String(take!(limited.io)) == string(short)
end

@testitem "NCAdd: outer and inner constructors agree" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    f1, f2 = Fermion(1), Fermion(2)
    base = (f1 + 2 * f2).dict
    K = keytype(base)
    scalar_key = K(1, empty(first(keys(base)).factors))
    f1key = first(k for k in keys(base) if k == NC.NCMul(1, [f1]))
    f2key = first(k for k in keys(base) if k == NC.NCMul(1, [f2]))

    dicts = (
        Dict(f1key => 1, f2key => 2),
        Dict(f1key => 1, scalar_key => 3),
        Dict(f1key => 1, f2key => 0),
        Dict(f1key => 0, scalar_key => 0),
    )
    for d in dicts, c in (0, 5)
        x = NC.NCAdd(c, copy(d))
        T, K2 = typeof(x).parameters
        y = NC.NCAdd{T,K2}(c, copy(d))
        @test typeof(x) == typeof(y)
        @test x == y
        @test x.dict == y.dict
        @test all(!iszero, values(x.dict))
        @test all(!NC.isscalar, keys(x.dict))
    end

    # the scalar term is absorbed into the coefficient
    @test NC.NCAdd(5, Dict(f1key => 1, scalar_key => 3)).coeff == 8

    # the coefficient type follows arithmetic: false + false is an Int
    z = NC.NCAdd(false, Dict(f1key => false, scalar_key => false))
    @test z.coeff isa Int
    @test z == NC.NCAdd{Int,K}(false, Dict(f1key => false, scalar_key => false))
end

@testitem "NCAdd: add_term!! inserts, updates and deletes like a Dict" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    using Random: Xoshiro
    # the reference uses only the public Dict API
    function ref_add_term(dict, key, coeff)
        newval = get(dict, key, false) + coeff
        T = promote_type(valtype(dict), typeof(newval))
        newdict = Dict{keytype(dict),T}(dict)
        iszero(newval) ? delete!(newdict, key) : (newdict[key] = newval)
        return newdict
    end
    keys_ = [NC.NCMul(1, [Fermion(n) for n in 1:k]) for k in 1:40]
    rng = Xoshiro(1)
    for (V, coeffs) in ((Int, -2:2), (Float64, [-1.0, -0.5, 0.0, 0.5, 1.0]), (Bool, [true, false]))
        dict = Dict{eltype(keys_),V}()
        ref = copy(dict)
        for _ in 1:2000
            key, coeff = rand(rng, keys_), rand(rng, coeffs)
            dict = NC.add_term!!(dict, key, coeff)
            ref = ref_add_term(ref, key, coeff)
            @test dict == ref
        end
        @test !any(iszero, values(dict))
        @test all(k -> haskey(dict, k), keys(ref))
    end
    # values that don't fit the value type widen the Dict
    k = first(keys_)
    d = Dict(k => true)
    d2 = NC.add_term!!(d, k, true)
    @test d2 == Dict(k => 2) && valtype(d2) == Int
    d = Dict(k => 1)
    @test NC.add_term!!(d, last(keys_), 0.5) == Dict(k => 1.0, last(keys_) => 0.5)
    # values that fit the value type are converted, also when they round to zero
    d = NC.add_term!!(Dict(k => 1.0), k, 1)
    @test d == Dict(k => 2.0) && valtype(d) == Float64
    d = NC.add_term!!(Dict(k => Float16(1)), last(keys_), 1 // 10^8)
    @test d == Dict(k => 1) && valtype(d) == Float16
    # keys of another type are converted
    tuplekey = NC.NCMul(1, (Fermion(1), Fermion(2)))
    d = NC.add_term!!(Dict(k => 1), tuplekey, 2)
    @test d == Dict(k => 1, keys_[2] => 2) && keytype(d) == typeof(k)
    @test NC.add_term!!(d, tuplekey, -2) == Dict(k => 1)
end
